/**
 * FA Extended Cross-Validation Report: Carm vs R
 * Covers: PAF oblimin/quartimin, ML oblimin/varimax, CFA fit indices, Diagnostics (KMO, Bartlett)
 *
 * Reference: R psych::fa(), lavaan::cfa(), psych::KMO(), psych::cortest.bartlett()
 * Run: npx tsx validation/ts-harness/fa-extended-report.ts && open validation/reports/fa-extended-report.html
 */
import { readFileSync, writeFileSync, mkdirSync } from 'fs'
import { runEFA, runCFA, runFADiagnostics } from 'carm'

// ── Types ────────────────────────────────────────────────────────────────

interface CfaStdLoading {
  factor: string; item: string; est: number; se: number; pvalue: number
}

interface FAExtRef {
  seed: number
  n: number
  p: number
  nf: number
  variableNames: string[]
  data: Record<string, number[]>   // { x1: [], x2: [], ... } column arrays
  paf_oblimin: {
    loadings: number[][]; communalities: number[]; phi: number[][]
  } | null
  paf_quartimin: {
    loadings: number[][]; communalities: number[]
  } | null
  ml_oblimin: {
    loadings: number[][]
    fit: { chiSq: number; df: number; pValue: number; rmsea: number; cfi: number; tli: number; srmr: number }
  } | null
  ml_varimax: {
    loadings: number[][]
  } | null
  cfa: {
    model: string
    fitMeasures: { chisq: number; df: number; pvalue: number; rmsea: number; cfi: number; tli: number; srmr: number }
    standardizedLoadings: CfaStdLoading[]
  } | null
  diagnostics: {
    kmo: number
    bartlett: { chiSq: number; pValue: number }
  }
}

// ── Metric tracking ──────────────────────────────────────────────────────

interface MetricDef {
  label: string
  tolerance: number
  errors: number[]
  category: string
}

function makeDefs(): Record<string, MetricDef> {
  const m = (label: string, tolerance: number, category: string): MetricDef =>
    ({ label, tolerance, errors: [], category })
  return {
    // PAF + oblimin
    paf_obl_load:    m('PAF Oblimin Loadings MAE',    0.05,  'PAF Oblimin'),
    paf_obl_comm:    m('PAF Oblimin Communalities MAE',0.05, 'PAF Oblimin'),
    // PAF + quartimin
    paf_qrt_load:    m('PAF Quartimin Loadings MAE',  0.05,  'PAF Quartimin'),
    // ML + oblimin
    ml_obl_load:     m('ML Oblimin Loadings MAE',     0.05,  'ML Oblimin'),
    ml_obl_chisq:    m('ML Oblimin χ² |Δ|',           10,    'ML Oblimin Fit'),
    ml_obl_rmsea:    m('ML Oblimin RMSEA |Δ|',        0.02,  'ML Oblimin Fit'),
    ml_obl_cfi:      m('ML Oblimin CFI |Δ|',          0.02,  'ML Oblimin Fit'),
    ml_obl_tli:      m('ML Oblimin TLI |Δ|',          0.02,  'ML Oblimin Fit'),
    ml_obl_srmr:     m('ML Oblimin SRMR |Δ|',         0.01,  'ML Oblimin Fit'),
    // ML + varimax
    ml_var_load:     m('ML Varimax Loadings MAE',     0.05,  'ML Varimax'),
    // CFA
    cfa_chisq:       m('CFA χ² |Δ|',                  10,    'CFA Fit'),
    cfa_rmsea:       m('CFA RMSEA |Δ|',               0.05,  'CFA Fit'),
    cfa_cfi:         m('CFA CFI |Δ|',                 0.05,  'CFA Fit'),
    cfa_tli:         m('CFA TLI |Δ|',                 0.05,  'CFA Fit'),
    cfa_srmr:        m('CFA SRMR |Δ|',                0.02,  'CFA Fit'),
    cfa_load:        m('CFA Loadings MAE',            0.1,   'CFA Loadings'),
    // Diagnostics
    diag_kmo:        m('KMO |Δ|',                     1e-4,  'Diagnostics'),
    diag_bart_chi2:  m('Bartlett χ² |Δ|',             0.01,  'Diagnostics'),
    diag_bart_p:     m('Bartlett p |Δ|',              1e-3,  'Diagnostics'),
  }
}

// ── Helpers ──────────────────────────────────────────────────────────────

function push(m: MetricDef, error: number) { m.errors.push(error) }

// ── Factor matching with exhaustive permutation + sign search ────────────

function* permutations(n: number): Generator<number[]> {
  const a = Array.from({ length: n }, (_, i) => i)
  const c = new Array(n).fill(0) as number[]
  yield [...a]
  let i = 0
  while (i < n) {
    if (c[i]! < i) {
      if (i % 2 === 0) { const t = a[0]!; a[0] = a[i]!; a[i] = t }
      else { const t = a[c[i]!]!; a[c[i]!] = a[i]!; a[i] = t }
      yield [...a]
      c[i]!++
      i = 0
    } else { c[i] = 0; i++ }
  }
}

function matchFactors(rL: number[][], cL: number[][]): { aligned: number[][]; mae: number } {
  const p = rL.length, k = rL[0]!.length
  let bestPerm: number[] = Array.from({ length: k }, (_, i) => i)
  let bestSigns: number[] = new Array(k).fill(1) as number[]
  let bestMAE = Infinity

  // Precompute cross-correlation for sign decisions
  const corr: number[][] = []
  for (let rj = 0; rj < k; rj++) {
    corr[rj] = []
    for (let cj = 0; cj < k; cj++) {
      let dot = 0
      for (let i = 0; i < p; i++) dot += (rL[i]?.[rj] ?? 0) * (cL[i]?.[cj] ?? 0)
      corr[rj]![cj] = dot
    }
  }

  for (const perm of permutations(k)) {
    const signs: number[] = []
    let totalErr = 0
    for (let rj = 0; rj < k; rj++) {
      const cj = perm[rj]!
      const sign = corr[rj]![cj]! >= 0 ? 1 : -1
      signs.push(sign)
      for (let i = 0; i < p; i++) {
        totalErr += Math.abs((rL[i]?.[rj] ?? 0) - (cL[i]?.[cj] ?? 0) * sign)
      }
    }
    const mae = totalErr / (p * k)
    if (mae < bestMAE) { bestMAE = mae; bestPerm = [...perm]; bestSigns = [...signs] }
  }

  const aligned = rL.map((_, i) =>
    bestPerm.map((cj, rj) => (cL[i]?.[cj] ?? 0) * bestSigns[rj]!)
  )
  return { aligned, mae: bestMAE }
}

function mae(a: number[], b: number[]): number {
  const n = Math.min(a.length, b.length)
  if (n === 0) return NaN
  let sum = 0
  for (let i = 0; i < n; i++) sum += Math.abs(a[i]! - b[i]!)
  return sum / n
}

/** Convert column-dict data to row-oriented */
function toRows(data: Record<string, number[]>, varNames: string[]): number[][] {
  const n = data[varNames[0]!]!.length
  const rows: number[][] = []
  for (let i = 0; i < n; i++) {
    const row: number[] = []
    for (const v of varNames) {
      row.push(data[v]![i]!)
    }
    rows.push(row)
  }
  return rows
}

/** Parse CFA model string (lavaan syntax) to Carm model format */
function parseCFAModel(modelStr: string, varNames: string[]): Record<string, number[]> {
  const model: Record<string, number[]> = {}
  const lines = modelStr.split('\n').map(l => l.trim()).filter(l => l.length > 0)
  for (const line of lines) {
    const match = line.match(/^(\w+)\s*=~\s*(.+)$/)
    if (!match) continue
    const factorName = match[1]!
    const items = match[2]!.split('+').map(s => s.trim())
    const indices: number[] = []
    for (const item of items) {
      const idx = varNames.indexOf(item)
      if (idx >= 0) indices.push(idx)
    }
    if (indices.length > 0) model[factorName] = indices
  }
  return model
}

interface AggStat {
  label: string; mean: number; median: number; p95: number; max: number
  threshold: number; passRate: number; n: number; category: string
}

function aggStat(m: MetricDef): AggStat {
  const sorted = [...m.errors].filter(v => !isNaN(v)).sort((a, b) => a - b)
  const n = sorted.length
  if (n === 0) return { label: m.label, mean: NaN, median: NaN, p95: NaN, max: NaN,
    threshold: m.tolerance, passRate: 0, n: 0, category: m.category }
  return {
    label: m.label,
    mean: sorted.reduce((a, b) => a + b, 0) / n,
    median: n % 2 === 0 ? (sorted[n / 2 - 1]! + sorted[n / 2]!) / 2 : sorted[Math.floor(n / 2)]!,
    p95: sorted[Math.min(Math.floor(n * 0.95), n - 1)]!,
    max: sorted[n - 1]!,
    threshold: m.tolerance,
    passRate: sorted.filter(v => v <= m.tolerance).length / n,
    n,
    category: m.category,
  }
}

function fmt(v: number, digits = 6): string {
  if (isNaN(v)) return 'NaN'
  if (v === 0) return '0'
  if (Math.abs(v) < 1e-12) return v.toExponential(2)
  if (Math.abs(v) < 0.001) return v.toExponential(2)
  return v.toFixed(digits)
}

function cc(val: number, threshold: number): string {
  if (isNaN(val)) return 'nan'
  if (val <= threshold) return 'pass'
  if (val <= threshold * 3) return 'warn'
  return 'fail'
}

function pr(rate: number): string {
  if (rate >= 0.99) return 'pass'
  if (rate >= 0.90) return 'warn'
  return 'fail'
}

// ── Load reference data ──────────────────────────────────────────────────

console.log('Loading FA extended reference data...')
const ref: FAExtRef[] = JSON.parse(readFileSync('validation/data/fa-extended-ref.json', 'utf-8'))
console.log(`Loaded ${ref.length} datasets`)

const metrics = makeDefs()

// ── Per-dataset tracking ─────────────────────────────────────────────────

interface DatasetResult {
  seed: number; n: number; p: number; nf: number; status: 'ok' | 'fail'; error?: string
  pafOblLoadMAE: number; pafOblCommMAE: number
  pafQrtLoadMAE: number
  mlOblLoadMAE: number; mlOblChisqΔ: number; mlOblRmseaΔ: number; mlOblCfiΔ: number; mlOblTliΔ: number; mlOblSrmrΔ: number
  mlVarLoadMAE: number
  cfaChisqΔ: number; cfaRmseaΔ: number; cfaCfiΔ: number; cfaTliΔ: number; cfaSrmrΔ: number; cfaLoadMAE: number
  diagKmoΔ: number; diagBartChi2Δ: number; diagBartPΔ: number
}

const results: DatasetResult[] = []
let nOk = 0; let nFail = 0
let nPafObl = 0; let nPafQrt = 0; let nMlObl = 0; let nMlVar = 0; let nCfa = 0; let nDiag = 0

for (const r of ref) {
  const dr: Partial<DatasetResult> = { seed: r.seed, n: r.n, p: r.p, nf: r.nf }

  try {
    const rowData = toRows(r.data, r.variableNames)

    // ── PAF + oblimin ──
    if (r.paf_oblimin) {
      try {
        const cFA = runEFA(rowData, {
          nFactors: r.nf, extraction: 'paf', rotation: 'oblimin',
          variableNames: [...r.variableNames], seed: 42,
        })
        const { mae: loadMAE } = matchFactors(r.paf_oblimin.loadings, cFA.loadings as number[][])
        dr.pafOblLoadMAE = loadMAE
        push(metrics.paf_obl_load!, loadMAE)

        const commMAE = mae(r.paf_oblimin.communalities, [...cFA.communalities])
        dr.pafOblCommMAE = commMAE
        push(metrics.paf_obl_comm!, commMAE)
        nPafObl++
      } catch (err) {
        dr.pafOblLoadMAE = NaN; dr.pafOblCommMAE = NaN
      }
    } else {
      dr.pafOblLoadMAE = NaN; dr.pafOblCommMAE = NaN
    }

    // ── PAF + quartimin ──
    if (r.paf_quartimin) {
      try {
        const cFA = runEFA(rowData, {
          nFactors: r.nf, extraction: 'paf', rotation: 'quartimin',
          variableNames: [...r.variableNames], seed: 42,
        })
        const { mae: loadMAE } = matchFactors(r.paf_quartimin.loadings, cFA.loadings as number[][])
        dr.pafQrtLoadMAE = loadMAE
        push(metrics.paf_qrt_load!, loadMAE)
        nPafQrt++
      } catch {
        dr.pafQrtLoadMAE = NaN
      }
    } else {
      dr.pafQrtLoadMAE = NaN
    }

    // ── ML + oblimin ──
    if (r.ml_oblimin) {
      try {
        const cFA = runEFA(rowData, {
          nFactors: r.nf, extraction: 'ml', rotation: 'oblimin',
          variableNames: [...r.variableNames], seed: 42,
        })
        const { mae: loadMAE } = matchFactors(r.ml_oblimin.loadings, cFA.loadings as number[][])
        dr.mlOblLoadMAE = loadMAE
        push(metrics.ml_obl_load!, loadMAE)

        const rFit = r.ml_oblimin.fit
        dr.mlOblChisqΔ = Math.abs(cFA.fit.chiSq - rFit.chiSq)
        push(metrics.ml_obl_chisq!, dr.mlOblChisqΔ)
        dr.mlOblRmseaΔ = Math.abs(cFA.fit.rmsea - rFit.rmsea)
        push(metrics.ml_obl_rmsea!, dr.mlOblRmseaΔ)
        dr.mlOblCfiΔ = Math.abs(cFA.fit.cfi - rFit.cfi)
        push(metrics.ml_obl_cfi!, dr.mlOblCfiΔ)
        dr.mlOblTliΔ = Math.abs(cFA.fit.tli - rFit.tli)
        push(metrics.ml_obl_tli!, dr.mlOblTliΔ)
        dr.mlOblSrmrΔ = Math.abs(cFA.fit.srmr - rFit.srmr)
        push(metrics.ml_obl_srmr!, dr.mlOblSrmrΔ)
        nMlObl++
      } catch {
        dr.mlOblLoadMAE = NaN; dr.mlOblChisqΔ = NaN; dr.mlOblRmseaΔ = NaN
        dr.mlOblCfiΔ = NaN; dr.mlOblTliΔ = NaN; dr.mlOblSrmrΔ = NaN
      }
    } else {
      dr.mlOblLoadMAE = NaN; dr.mlOblChisqΔ = NaN; dr.mlOblRmseaΔ = NaN
      dr.mlOblCfiΔ = NaN; dr.mlOblTliΔ = NaN; dr.mlOblSrmrΔ = NaN
    }

    // ── ML + varimax ──
    if (r.ml_varimax) {
      try {
        const cFA = runEFA(rowData, {
          nFactors: r.nf, extraction: 'ml', rotation: 'varimax',
          variableNames: [...r.variableNames], seed: 42,
        })
        const { mae: loadMAE } = matchFactors(r.ml_varimax.loadings, cFA.loadings as number[][])
        dr.mlVarLoadMAE = loadMAE
        push(metrics.ml_var_load!, loadMAE)
        nMlVar++
      } catch {
        dr.mlVarLoadMAE = NaN
      }
    } else {
      dr.mlVarLoadMAE = NaN
    }

    // ── CFA ──
    if (r.cfa) {
      try {
        const cfaModel = parseCFAModel(r.cfa.model, r.variableNames)
        const cCFA = runCFA(rowData, cfaModel, {
          variableNames: [...r.variableNames],
        })

        const rFit = r.cfa.fitMeasures
        dr.cfaChisqΔ = Math.abs(cCFA.fit.chiSq - rFit.chisq)
        push(metrics.cfa_chisq!, dr.cfaChisqΔ)
        dr.cfaRmseaΔ = Math.abs(cCFA.fit.rmsea - rFit.rmsea)
        push(metrics.cfa_rmsea!, dr.cfaRmseaΔ)
        dr.cfaCfiΔ = Math.abs(cCFA.fit.cfi - rFit.cfi)
        push(metrics.cfa_cfi!, dr.cfaCfiΔ)
        dr.cfaTliΔ = Math.abs(cCFA.fit.tli - rFit.tli)
        push(metrics.cfa_tli!, dr.cfaTliΔ)
        dr.cfaSrmrΔ = Math.abs(cCFA.fit.srmr - rFit.srmr)
        push(metrics.cfa_srmr!, dr.cfaSrmrΔ)

        // Compare CFA standardized loadings
        // R ref has per-item standardized loadings, Carm has loadings matrix
        // Build flattened vector of R loadings (ordered by factor then item)
        const rLoadFlat = r.cfa.standardizedLoadings.map(l => l.est)
        // Carm CFA standardized loadings are in the same item order
        // CFA loadings are in model order: factor → items within that factor
        const cLoadFlat: number[] = []
        const factorNames = Object.keys(cfaModel)
        for (let fi = 0; fi < factorNames.length; fi++) {
          const itemIndices = cfaModel[factorNames[fi]!]!
          for (const idx of itemIndices) {
            // standardizedLoadings[item][factor]
            cLoadFlat.push(Math.abs(cCFA.standardizedLoadings[idx]?.[fi] ?? 0))
          }
        }
        // R loadings are absolute standardized estimates
        const rLoadAbs = rLoadFlat.map(v => Math.abs(v))
        dr.cfaLoadMAE = mae(rLoadAbs, cLoadFlat)
        push(metrics.cfa_load!, dr.cfaLoadMAE)
        nCfa++
      } catch {
        dr.cfaChisqΔ = NaN; dr.cfaRmseaΔ = NaN; dr.cfaCfiΔ = NaN
        dr.cfaTliΔ = NaN; dr.cfaSrmrΔ = NaN; dr.cfaLoadMAE = NaN
      }
    } else {
      dr.cfaChisqΔ = NaN; dr.cfaRmseaΔ = NaN; dr.cfaCfiΔ = NaN
      dr.cfaTliΔ = NaN; dr.cfaSrmrΔ = NaN; dr.cfaLoadMAE = NaN
    }

    // ── Diagnostics (KMO, Bartlett) ──
    try {
      const cDiag = runFADiagnostics(rowData, { seed: 42, parallelIterations: 100 })
      dr.diagKmoΔ = Math.abs(cDiag.kmo - r.diagnostics.kmo)
      push(metrics.diag_kmo!, dr.diagKmoΔ)
      dr.diagBartChi2Δ = Math.abs(cDiag.bartlett.chiSq - r.diagnostics.bartlett.chiSq)
      push(metrics.diag_bart_chi2!, dr.diagBartChi2Δ)
      dr.diagBartPΔ = Math.abs(cDiag.bartlett.pValue - r.diagnostics.bartlett.pValue)
      push(metrics.diag_bart_p!, dr.diagBartPΔ)
      nDiag++
    } catch {
      dr.diagKmoΔ = NaN; dr.diagBartChi2Δ = NaN; dr.diagBartPΔ = NaN
    }

    dr.status = 'ok'
    nOk++
  } catch (err) {
    dr.status = 'fail'
    dr.error = err instanceof Error ? err.message : String(err)
    dr.pafOblLoadMAE = NaN; dr.pafOblCommMAE = NaN
    dr.pafQrtLoadMAE = NaN
    dr.mlOblLoadMAE = NaN; dr.mlOblChisqΔ = NaN; dr.mlOblRmseaΔ = NaN
    dr.mlOblCfiΔ = NaN; dr.mlOblTliΔ = NaN; dr.mlOblSrmrΔ = NaN
    dr.mlVarLoadMAE = NaN
    dr.cfaChisqΔ = NaN; dr.cfaRmseaΔ = NaN; dr.cfaCfiΔ = NaN
    dr.cfaTliΔ = NaN; dr.cfaSrmrΔ = NaN; dr.cfaLoadMAE = NaN
    dr.diagKmoΔ = NaN; dr.diagBartChi2Δ = NaN; dr.diagBartPΔ = NaN
    nFail++
  }
  results.push(dr as DatasetResult)
}

// ── Aggregate ────────────────────────────────────────────────────────────

const allStats = Object.values(metrics).filter(m => m.errors.length > 0).map(aggStat)
const categories = [...new Set(allStats.map(s => s.category))]

// ── Console summary ──────────────────────────────────────────────────────

console.log('\n═══ FA Extended Cross-Validation Summary ═══\n')
console.log(`Datasets: ${ref.length} total, ${nOk} ok, ${nFail} failed`)
console.log(`Successful: PAF oblimin=${nPafObl}, PAF quartimin=${nPafQrt}, ML oblimin=${nMlObl}, ML varimax=${nMlVar}, CFA=${nCfa}, Diagnostics=${nDiag}\n`)

for (const cat of categories) {
  console.log(`── ${cat} ──`)
  const catStats = allStats.filter(s => s.category === cat)
  for (const s of catStats) {
    const status = s.passRate >= 0.99 ? '✓' : s.passRate >= 0.90 ? '~' : '✗'
    console.log(`  ${status} ${s.label.padEnd(32)} mean=${fmt(s.mean, 8).padStart(14)} max=${fmt(s.max, 8).padStart(14)} pass=${(s.passRate * 100).toFixed(1).padStart(5)}% (n=${s.n})`)
  }
}

const totalPass = allStats.filter(s => s.passRate >= 0.99).length
const totalMetrics = allStats.length
console.log(`\nOverall: ${totalPass}/${totalMetrics} metrics at ≥99% pass rate`)

// ── HTML Report ──────────────────────────────────────────────────────────

const timestamp = new Date().toISOString().slice(0, 19)

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>FA Extended Cross-Validation: Carm vs R</title>
<style>
  :root { --bg: #0d1117; --fg: #c9d1d9; --card: #161b22; --border: #30363d;
    --green: #3fb950; --yellow: #d29922; --red: #f85149; --blue: #58a6ff; --purple: #bc8cff; }
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: 'SF Mono', 'Cascadia Code', Consolas, monospace; font-size: 13px;
    background: var(--bg); color: var(--fg); padding: 24px 32px; line-height: 1.5; }
  h1 { font-size: 22px; font-weight: 700; margin-bottom: 4px; color: #fff; }
  h2 { font-size: 16px; font-weight: 600; margin: 32px 0 12px; color: var(--blue);
    border-bottom: 1px solid var(--border); padding-bottom: 6px; }
  h3 { font-size: 14px; font-weight: 600; margin: 20px 0 8px; color: var(--purple); }
  .subtitle { color: #8b949e; margin-bottom: 20px; }
  .summary-bar { display: flex; gap: 14px; margin-bottom: 24px; flex-wrap: wrap; }
  .summary-card { background: var(--card); border: 1px solid var(--border);
    border-radius: 8px; padding: 12px 18px; min-width: 130px; }
  .summary-card .num { font-size: 28px; font-weight: 700; }
  .summary-card .label { font-size: 11px; color: #8b949e; text-transform: uppercase; letter-spacing: 0.5px; }
  .pass { color: var(--green); }
  .warn { color: var(--yellow); }
  .fail { color: var(--red); }
  .nan { color: #6e7681; }
  table { width: 100%; border-collapse: collapse; margin-bottom: 16px; }
  th { text-align: left; padding: 8px 10px; background: var(--card); border: 1px solid var(--border);
    font-weight: 600; font-size: 11px; text-transform: uppercase; letter-spacing: 0.3px;
    color: #8b949e; position: sticky; top: 0; z-index: 1; }
  td { padding: 5px 10px; border: 1px solid var(--border); font-variant-numeric: tabular-nums; }
  tr:nth-child(even) td { background: rgba(22, 27, 34, 0.5); }
  tr:hover td { background: rgba(56, 139, 253, 0.08); }
  td.pass { background: rgba(63, 185, 80, 0.1); }
  td.warn { background: rgba(210, 153, 34, 0.1); }
  td.fail { background: rgba(248, 81, 73, 0.1); }
  .threshold { font-size: 10px; color: #6e7681; }
  .pct { font-weight: 700; }
  .section-note { color: #8b949e; font-size: 12px; margin-bottom: 12px; }
  .badge { display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 11px; font-weight: 600; margin-right: 6px; }
  .badge.ok { background: rgba(63, 185, 80, 0.15); color: var(--green); }
  .badge.err { background: rgba(248, 81, 73, 0.15); color: var(--red); }
</style>
</head>
<body>

<h1>FA Extended Cross-Validation</h1>
<p class="subtitle">Carm (TypeScript) vs R psych::fa() / lavaan::cfa() — ${ref.length} datasets — ${timestamp}</p>

<div class="summary-bar">
  <div class="summary-card"><div class="num">${ref.length}</div><div class="label">Datasets</div></div>
  <div class="summary-card"><div class="num">${totalMetrics}</div><div class="label">Metrics Tested</div></div>
  <div class="summary-card"><div class="num ${totalPass === totalMetrics ? 'pass' : totalPass >= totalMetrics * 0.9 ? 'warn' : 'fail'}">${totalPass}/${totalMetrics}</div><div class="label">≥99% Pass Rate</div></div>
  <div class="summary-card"><div class="num ${nOk === ref.length ? 'pass' : 'fail'}">${nOk}/${ref.length}</div><div class="label">Ran OK</div></div>
  <div class="summary-card"><div class="num">${nPafObl}</div><div class="label">PAF Oblimin</div></div>
  <div class="summary-card"><div class="num">${nMlObl}</div><div class="label">ML Oblimin</div></div>
  <div class="summary-card"><div class="num">${nCfa}</div><div class="label">CFA</div></div>
  <div class="summary-card"><div class="num">${nDiag}</div><div class="label">Diagnostics</div></div>
</div>

${categories.map(cat => {
  const catStats = allStats.filter(s => s.category === cat)
  return `
<h2>${cat}</h2>
<p class="section-note">Green = within threshold, Yellow = within 3x, Red = exceeds 3x.</p>
<table>
<thead><tr><th>Metric</th><th>Threshold</th><th>n</th><th>Mean |e|</th><th>Median |e|</th><th>P95 |e|</th><th>Max |e|</th><th>Pass Rate</th></tr></thead>
<tbody>
${catStats.map(s => `<tr>
  <td>${s.label}</td>
  <td class="threshold">${fmt(s.threshold, 4)}</td>
  <td>${s.n}</td>
  <td class="${cc(s.mean, s.threshold)}">${fmt(s.mean)}</td>
  <td class="${cc(s.median, s.threshold)}">${fmt(s.median)}</td>
  <td class="${cc(s.p95, s.threshold)}">${fmt(s.p95)}</td>
  <td class="${cc(s.max, s.threshold)}">${fmt(s.max)}</td>
  <td class="${pr(s.passRate)}"><span class="pct">${(s.passRate * 100).toFixed(1)}%</span></td>
</tr>`).join('\n')}
</tbody>
</table>`
}).join('\n')}

<h2>Per-Dataset: EFA Loadings (top 50 by max loading error)</h2>
<p class="section-note">All loading MAEs after permutation + sign alignment. Sorted by worst loading error.</p>
<table>
<thead><tr><th>Seed</th><th>n</th><th>p</th><th>nf</th><th>PAF Oblimin</th><th>PAF Quartimin</th><th>ML Oblimin</th><th>ML Varimax</th><th>KMO |Δ|</th></tr></thead>
<tbody>
${results
  .filter(r => r.status === 'ok')
  .sort((a, b) => {
    const maxA = Math.max(a.pafOblLoadMAE || 0, a.pafQrtLoadMAE || 0, a.mlOblLoadMAE || 0, a.mlVarLoadMAE || 0)
    const maxB = Math.max(b.pafOblLoadMAE || 0, b.pafQrtLoadMAE || 0, b.mlOblLoadMAE || 0, b.mlVarLoadMAE || 0)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => `<tr>
  <td>${d.seed}</td><td>${d.n}</td><td>${d.p}</td><td>${d.nf}</td>
  <td class="${cc(d.pafOblLoadMAE, 0.05)}">${fmt(d.pafOblLoadMAE, 4)}</td>
  <td class="${cc(d.pafQrtLoadMAE, 0.05)}">${fmt(d.pafQrtLoadMAE, 4)}</td>
  <td class="${cc(d.mlOblLoadMAE, 0.05)}">${fmt(d.mlOblLoadMAE, 4)}</td>
  <td class="${cc(d.mlVarLoadMAE, 0.05)}">${fmt(d.mlVarLoadMAE, 4)}</td>
  <td class="${cc(d.diagKmoΔ, 1e-4)}">${fmt(d.diagKmoΔ)}</td>
</tr>`).join('\n')}
</tbody>
</table>

<h2>Per-Dataset: ML Oblimin Fit Indices (top 50 by max fit error)</h2>
<p class="section-note">Sorted by worst fit index discrepancy.</p>
<table>
<thead><tr><th>Seed</th><th>n</th><th>p</th><th>nf</th><th>χ² |Δ|</th><th>RMSEA |Δ|</th><th>CFI |Δ|</th><th>TLI |Δ|</th><th>SRMR |Δ|</th></tr></thead>
<tbody>
${results
  .filter(r => r.status === 'ok' && !isNaN(r.mlOblChisqΔ))
  .sort((a, b) => {
    const maxA = Math.max(a.mlOblChisqΔ || 0, (a.mlOblRmseaΔ || 0) * 100, (a.mlOblCfiΔ || 0) * 100, (a.mlOblTliΔ || 0) * 100, (a.mlOblSrmrΔ || 0) * 100)
    const maxB = Math.max(b.mlOblChisqΔ || 0, (b.mlOblRmseaΔ || 0) * 100, (b.mlOblCfiΔ || 0) * 100, (b.mlOblTliΔ || 0) * 100, (b.mlOblSrmrΔ || 0) * 100)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => `<tr>
  <td>${d.seed}</td><td>${d.n}</td><td>${d.p}</td><td>${d.nf}</td>
  <td class="${cc(d.mlOblChisqΔ, 10)}">${fmt(d.mlOblChisqΔ, 2)}</td>
  <td class="${cc(d.mlOblRmseaΔ, 0.02)}">${fmt(d.mlOblRmseaΔ, 4)}</td>
  <td class="${cc(d.mlOblCfiΔ, 0.02)}">${fmt(d.mlOblCfiΔ, 4)}</td>
  <td class="${cc(d.mlOblTliΔ, 0.02)}">${fmt(d.mlOblTliΔ, 4)}</td>
  <td class="${cc(d.mlOblSrmrΔ, 0.01)}">${fmt(d.mlOblSrmrΔ, 4)}</td>
</tr>`).join('\n')}
</tbody>
</table>

<h2>Per-Dataset: CFA (top 50 by max error)</h2>
<p class="section-note">Sorted by worst CFA fit discrepancy.</p>
<table>
<thead><tr><th>Seed</th><th>n</th><th>p</th><th>nf</th><th>χ² |Δ|</th><th>RMSEA |Δ|</th><th>CFI |Δ|</th><th>TLI |Δ|</th><th>SRMR |Δ|</th><th>Load MAE</th></tr></thead>
<tbody>
${results
  .filter(r => r.status === 'ok' && !isNaN(r.cfaChisqΔ))
  .sort((a, b) => {
    const maxA = Math.max(a.cfaChisqΔ || 0, (a.cfaRmseaΔ || 0) * 100, (a.cfaCfiΔ || 0) * 100, (a.cfaTliΔ || 0) * 100)
    const maxB = Math.max(b.cfaChisqΔ || 0, (b.cfaRmseaΔ || 0) * 100, (b.cfaCfiΔ || 0) * 100, (b.cfaTliΔ || 0) * 100)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => `<tr>
  <td>${d.seed}</td><td>${d.n}</td><td>${d.p}</td><td>${d.nf}</td>
  <td class="${cc(d.cfaChisqΔ, 10)}">${fmt(d.cfaChisqΔ, 2)}</td>
  <td class="${cc(d.cfaRmseaΔ, 0.05)}">${fmt(d.cfaRmseaΔ, 4)}</td>
  <td class="${cc(d.cfaCfiΔ, 0.05)}">${fmt(d.cfaCfiΔ, 4)}</td>
  <td class="${cc(d.cfaTliΔ, 0.05)}">${fmt(d.cfaTliΔ, 4)}</td>
  <td class="${cc(d.cfaSrmrΔ, 0.02)}">${fmt(d.cfaSrmrΔ, 4)}</td>
  <td class="${cc(d.cfaLoadMAE, 0.1)}">${fmt(d.cfaLoadMAE, 4)}</td>
</tr>`).join('\n')}
</tbody>
</table>

<h2>Per-Dataset: Diagnostics (top 50 by max error)</h2>
<table>
<thead><tr><th>Seed</th><th>n</th><th>p</th><th>nf</th><th>KMO |Δ|</th><th>Bartlett χ² |Δ|</th><th>Bartlett p |Δ|</th></tr></thead>
<tbody>
${results
  .filter(r => r.status === 'ok' && !isNaN(r.diagKmoΔ))
  .sort((a, b) => {
    const maxA = Math.max(a.diagKmoΔ || 0, a.diagBartChi2Δ || 0)
    const maxB = Math.max(b.diagKmoΔ || 0, b.diagBartChi2Δ || 0)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => `<tr>
  <td>${d.seed}</td><td>${d.n}</td><td>${d.p}</td><td>${d.nf}</td>
  <td class="${cc(d.diagKmoΔ, 1e-4)}">${fmt(d.diagKmoΔ)}</td>
  <td class="${cc(d.diagBartChi2Δ, 0.01)}">${fmt(d.diagBartChi2Δ)}</td>
  <td class="${cc(d.diagBartPΔ, 1e-3)}">${fmt(d.diagBartPΔ)}</td>
</tr>`).join('\n')}
</tbody>
</table>

${nFail > 0 ? `
<h2>Failed Datasets</h2>
<table>
<thead><tr><th>Seed</th><th>n</th><th>p</th><th>nf</th><th>Error</th></tr></thead>
<tbody>
${results.filter(r => r.status === 'fail').map(r => `<tr>
  <td>${r.seed}</td><td>${r.n}</td><td>${r.p}</td><td>${r.nf}</td>
  <td style="color:var(--red);">${r.error ?? 'Unknown'}</td>
</tr>`).join('\n')}
</tbody>
</table>
` : ''}

<p style="margin-top:40px;color:#6e7681;font-size:11px;">
  Generated ${timestamp} · Carm vs R psych/lavaan · ${ref.length} datasets · ${totalMetrics} metrics ·
  PAF oblimin=${nPafObl}, PAF quartimin=${nPafQrt}, ML oblimin=${nMlObl}, ML varimax=${nMlVar}, CFA=${nCfa}
</p>

</body>
</html>`

mkdirSync('validation/reports', { recursive: true })
writeFileSync('validation/reports/fa-extended-report.html', html)
console.log(`\nReport saved to validation/reports/fa-extended-report.html`)
