/**
 * Full FA Cross-Validation Report: Carm vs R
 * Covers: Promax, Geomin, Diagnostics, Real dataset
 *
 * Run: npx tsx validation/ts-harness/fa-full-report.ts && open validation/reports/fa-crossval-report.html
 */
import { readFileSync, writeFileSync } from 'fs'
import { runEFA, runFADiagnostics } from 'carm'
type FAResult = import('carm').FAResult
type FADiagnostics = import('carm').FADiagnostics

// ── Types ────────────────────────────────────────────────────────────────

interface DatasetRef {
  id: number
  params: { n: number; p: number; k: number; loadingStrength: number; seed: number }
  data: number[][]
  variableNames: string[]
  rResults: {
    eigenvalues: number[]
    loadings: number[][]
    communalities: number[]
    uniqueness: number[]
    fit: {
      chiSq: number; df: number; pValue: number
      rmsea: number; rmseaLower: number; rmseaUpper: number
      cfi: number; tli: number; srmr: number; aic: number; bic: number
    }
    kmo: { overall: number; perItem: number[] }
    bartlett: { chiSq: number; df: number; pValue: number }
    nFactors: number
  }
}

interface GeominRef {
  id: number; n: number; p: number; k: number
  loadings: number[][]; phi: number[][]; converged: boolean
}

interface RealGeominRef {
  [k: string]: { k: number; loadings: number[][]; phi: number[][]; converged: boolean }
}

// ── Factor matching ──────────────────────────────────────────────────────

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

  const corr: number[][] = []
  for (let rj = 0; rj < k; rj++) {
    corr[rj] = []
    for (let cj = 0; cj < k; cj++) {
      let dot = 0
      for (let i = 0; i < p; i++) dot += rL[i]![rj]! * cL[i]![cj]!
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
        totalErr += Math.abs(rL[i]![rj]! - cL[i]![cj]! * sign)
      }
    }
    const mae = totalErr / (p * k)
    if (mae < bestMAE) { bestMAE = mae; bestPerm = [...perm]; bestSigns = [...signs] }
  }

  const aligned = rL.map((_, i) =>
    bestPerm.map((cj, rj) => cL[i]![cj]! * bestSigns[rj]!)
  )
  return { aligned, mae: bestMAE }
}

function mae(a: number[], b: number[]): number {
  let sum = 0
  for (let i = 0; i < a.length; i++) sum += Math.abs(a[i]! - b[i]!)
  return sum / a.length
}

function maxErr(a: number[], b: number[]): number {
  let mx = 0
  for (let i = 0; i < a.length; i++) mx = Math.max(mx, Math.abs(a[i]! - b[i]!))
  return mx
}

// ── Load data ────────────────────────────────────────────────────────────

const refData = JSON.parse(readFileSync('validation/data/fa-crossval-data.json', 'utf-8')) as { datasets: DatasetRef[] }
const datasets = refData.datasets
const geominRefs = JSON.parse(readFileSync('validation/data/fa-geomin-ref.json', 'utf-8')) as GeominRef[]
const geominRefMap = new Map(geominRefs.map(r => [r.id, r]))
const realGeominRef = JSON.parse(readFileSync('validation/data/fa-geomin-real-ref.json', 'utf-8')) as RealGeominRef

console.log(`Loaded ${datasets.length} synthetic datasets`)

// ── 1. Promax cross-validation ───────────────────────────────────────────

interface PromaxResult {
  id: number; n: number; p: number; k: number; lam: number
  status: 'ok' | 'fail'; error?: string
  loadMAE: number; loadMax: number; commMAE: number; uniqMAE: number
  chiSqDiff: number; dfDiff: number; rmseaDiff: number; cfiDiff: number; tliDiff: number; srmrDiff: number
  eigMAE: number; kmoΔ: number; bartΔ: number
}

const promaxResults: PromaxResult[] = []

for (const ds of datasets) {
  const r: Partial<PromaxResult> = {
    id: ds.id, n: ds.params.n, p: ds.params.p, k: ds.params.k, lam: ds.params.loadingStrength,
  }
  try {
    const diag: FADiagnostics = runFADiagnostics(ds.data, { seed: 42, parallelIterations: 100 })
    const fa: FAResult = runEFA(ds.data, {
      nFactors: ds.params.k, extraction: 'ml', rotation: 'promax',
      variableNames: [...ds.variableNames], seed: 42,
    })

    const { aligned, mae: lMAE } = matchFactors(ds.rResults.loadings, fa.loadings as number[][])
    r.loadMAE = lMAE
    r.loadMax = maxErr(ds.rResults.loadings.flat(), aligned.flat())
    r.commMAE = mae(ds.rResults.communalities, [...fa.communalities])
    r.uniqMAE = mae(ds.rResults.uniqueness, [...fa.uniqueness])
    r.chiSqDiff = Math.abs(ds.rResults.fit.chiSq - fa.fit.chiSq)
    r.dfDiff = Math.abs(ds.rResults.fit.df - fa.fit.df)
    r.rmseaDiff = Math.abs(ds.rResults.fit.rmsea - fa.fit.rmsea)
    r.cfiDiff = Math.abs(ds.rResults.fit.cfi - fa.fit.cfi)
    r.tliDiff = Math.abs(ds.rResults.fit.tli - fa.fit.tli)
    r.srmrDiff = Math.abs(ds.rResults.fit.srmr - fa.fit.srmr)

    const rEig = ds.rResults.eigenvalues
    const cEig = [...diag.parallelEigenvalues]
    r.eigMAE = mae(rEig.slice(0, Math.min(rEig.length, cEig.length)), cEig.slice(0, Math.min(rEig.length, cEig.length)))
    r.kmoΔ = Math.abs(ds.rResults.kmo.overall - diag.kmo)
    r.bartΔ = Math.abs(ds.rResults.bartlett.chiSq - diag.bartlett.chiSq)
    r.status = 'ok'
  } catch (err) {
    r.status = 'fail'
    r.error = err instanceof Error ? err.message : String(err)
    r.loadMAE = NaN; r.loadMax = NaN; r.commMAE = NaN; r.uniqMAE = NaN
    r.chiSqDiff = NaN; r.dfDiff = NaN; r.rmseaDiff = NaN; r.cfiDiff = NaN; r.tliDiff = NaN; r.srmrDiff = NaN
    r.eigMAE = NaN; r.kmoΔ = NaN; r.bartΔ = NaN
  }
  promaxResults.push(r as PromaxResult)
}

const promaxPass = promaxResults.filter(r => r.status === 'ok' && r.loadMAE <= 0.05).length
const promaxFail = promaxResults.length - promaxPass
console.log(`Promax: ${promaxPass}/${datasets.length} pass`)

// ── 2. Geomin cross-validation ───────────────────────────────────────────

interface GeominResult {
  id: number; n: number; p: number; k: number
  loadMAE: number; pass: boolean
}

const geominResults: GeominResult[] = []

for (const ds of datasets) {
  const ref = geominRefMap.get(ds.id)
  if (!ref) continue

  const fa: FAResult = runEFA(ds.data, {
    nFactors: ds.params.k, extraction: 'ml', rotation: 'geomin',
    geominDelta: 0.01, variableNames: [...ds.variableNames], seed: 42,
  })

  const { mae: lMAE } = matchFactors(ref.loadings, fa.loadings as number[][])
  const pass = lMAE <= 0.05
  geominResults.push({ id: ds.id, n: ref.n, p: ref.p, k: ref.k, loadMAE: lMAE, pass })
}

const geominPass = geominResults.filter(r => r.pass).length
console.log(`Geomin: ${geominPass}/${geominResults.length} pass`)

// ── 3. Real dataset ──────────────────────────────────────────────────────

interface RealResult {
  k: number; rotation: string; loadMAE: number; pass: boolean
  firstVars: { name: string; carm: number[]; r: number[] }[]
}

const realResults: RealResult[] = []

let realData: number[][] | null = null
let realHeader: string[] = []
try {
  const csv = readFileSync('/Users/mohammedsaqr/Downloads/rraw_dataaw_data.csv', 'utf-8')
  const lines = csv.trim().split('\n')
  realHeader = lines[0]!.split(',')
  realData = lines.slice(1).map(l => l.split(',').map(Number))
  console.log(`Real dataset: ${realData.length}×${realHeader.length}`)
} catch { console.log('Real dataset not found, skipping') }

if (realData) {
  // Load real promax reference
  let realPromaxRef: Record<string, { k: number; loadings: number[][] }> | null = null
  try {
    realPromaxRef = JSON.parse(readFileSync('validation/data/fa-real-crossval-ref.json', 'utf-8'))
  } catch { /* no promax ref */ }

  // Geomin real
  for (const kStr of ['3', '5']) {
    const rr = realGeominRef[kStr]
    if (!rr) continue
    const fa = runEFA(realData, {
      nFactors: rr.k, extraction: 'ml', rotation: 'geomin',
      geominDelta: 0.01, variableNames: realHeader,
    })
    const { aligned, mae: lMAE } = matchFactors(rr.loadings, fa.loadings as number[][])
    const firstVars = realHeader.slice(0, 5).map((name, i) => ({
      name,
      carm: aligned[i]!,
      r: rr.loadings[i]!,
    }))
    realResults.push({ k: rr.k, rotation: 'geomin', loadMAE: lMAE, pass: lMAE <= 0.05, firstVars })
  }

  // Promax real
  if (realPromaxRef) {
    for (const kStr of ['3', '4', '5', '6']) {
      const rr = realPromaxRef[kStr]
      if (!rr) continue
      const fa = runEFA(realData, {
        nFactors: rr.k, extraction: 'ml', rotation: 'promax',
        variableNames: realHeader,
      })
      const { aligned, mae: lMAE } = matchFactors(rr.loadings, fa.loadings as number[][])
      const firstVars = realHeader.slice(0, 5).map((name, i) => ({
        name,
        carm: aligned[i]!,
        r: rr.loadings[i]!,
      }))
      realResults.push({ k: rr.k, rotation: 'promax', loadMAE: lMAE, pass: lMAE <= 0.05, firstVars })
    }
  }
}

// ── Aggregate stats ──────────────────────────────────────────────────────

interface AggStat { label: string; mean: number; median: number; p95: number; max: number; threshold: number; passRate: number }

function aggStat(label: string, vals: number[], threshold: number): AggStat {
  const sorted = [...vals].filter(v => !isNaN(v)).sort((a, b) => a - b)
  const n = sorted.length
  if (n === 0) return { label, mean: NaN, median: NaN, p95: NaN, max: NaN, threshold, passRate: 0 }
  const mean = sorted.reduce((s, v) => s + v, 0) / n
  const median = n % 2 === 0 ? (sorted[n / 2 - 1]! + sorted[n / 2]!) / 2 : sorted[Math.floor(n / 2)]!
  const p95 = sorted[Math.min(Math.floor(n * 0.95), n - 1)]!
  const max = sorted[n - 1]!
  const passCount = sorted.filter(v => v <= threshold).length
  return { label, mean, median, p95, max, threshold, passRate: passCount / n }
}

const okPromax = promaxResults.filter(r => r.status === 'ok')
const promaxStats: AggStat[] = [
  aggStat('Eigenvalue MAE', okPromax.map(c => c.eigMAE), 1e-8),
  aggStat('KMO Overall |Δ|', okPromax.map(c => c.kmoΔ), 1e-8),
  aggStat('Bartlett χ² |Δ|', okPromax.map(c => c.bartΔ), 0.01),
  aggStat('Loading MAE', okPromax.map(c => c.loadMAE), 0.05),
  aggStat('Loading Max |Δ|', okPromax.map(c => c.loadMax), 0.20),
  aggStat('Communality MAE', okPromax.map(c => c.commMAE), 0.05),
  aggStat('Uniqueness MAE', okPromax.map(c => c.uniqMAE), 0.05),
  aggStat('χ² |Δ|', okPromax.map(c => c.chiSqDiff), 10.0),
  aggStat('RMSEA |Δ|', okPromax.map(c => c.rmseaDiff), 0.02),
  aggStat('CFI |Δ|', okPromax.map(c => c.cfiDiff), 0.02),
  aggStat('TLI |Δ|', okPromax.map(c => c.tliDiff), 0.02),
  aggStat('SRMR |Δ|', okPromax.map(c => c.srmrDiff), 0.01),
]

const geominStats: AggStat[] = [
  aggStat('Loading MAE', geominResults.map(c => c.loadMAE), 0.05),
]

// ── HTML ─────────────────────────────────────────────────────────────────

function fmt(v: number, d = 6): string {
  if (isNaN(v)) return '<span class="nan">NaN</span>'
  if (v === 0) return '0'
  if (Math.abs(v) < 1e-12) return v.toExponential(2)
  if (Math.abs(v) < 0.001) return v.toExponential(2)
  return v.toFixed(d)
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

const totalObs = datasets.map(d => d.params.n).reduce((a, b) => a + b, 0)
const timestamp = new Date().toISOString().slice(0, 19)

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>FA Cross-Validation: Carm vs R — Full Report</title>
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
  .real-cmp { font-size: 12px; margin-bottom: 6px; }
  .real-cmp .varname { color: var(--blue); font-weight: 600; min-width: 60px; display: inline-block; }
  .real-cmp .vals { color: #8b949e; }
  .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
  @media (max-width: 900px) { .grid-2 { grid-template-columns: 1fr; } }
</style>
</head>
<body>

<h1>Factor Analysis Cross-Validation</h1>
<p class="subtitle">Carm (TypeScript) vs R psych::fa() + GPArotation::geominQ() — ${datasets.length} synthetic + real dataset</p>

<div class="summary-bar">
  <div class="summary-card"><div class="num">${datasets.length}</div><div class="label">Synthetic Datasets</div></div>
  <div class="summary-card"><div class="num">${totalObs.toLocaleString()}</div><div class="label">Total Observations</div></div>
  <div class="summary-card"><div class="num ${promaxPass === datasets.length ? 'pass' : 'fail'}">${promaxPass}/${datasets.length}</div><div class="label">Promax Pass</div></div>
  <div class="summary-card"><div class="num ${geominPass === geominResults.length ? 'pass' : 'fail'}">${geominPass}/${geominResults.length}</div><div class="label">Geomin Pass</div></div>
  <div class="summary-card"><div class="num ${realResults.every(r => r.pass) ? 'pass' : 'fail'}">${realResults.filter(r => r.pass).length}/${realResults.length}</div><div class="label">Real Dataset</div></div>
</div>

<!-- ═══ SECTION 1: PROMAX ═══ -->
<h2>1. Promax Rotation — Aggregate (vs R psych::fa, ML + promax)</h2>
<p class="section-note">Green = within threshold, Yellow = within 3x, Red = exceeds 3x. Threshold = maximum acceptable |R - Carm|.</p>
<table>
<thead><tr><th>Metric</th><th>Threshold</th><th>Mean e</th><th>Median e</th><th>P95 e</th><th>Max e</th><th>Pass Rate</th></tr></thead>
<tbody>
${promaxStats.map(s => `<tr>
  <td>${s.label}</td>
  <td class="threshold">${fmt(s.threshold, 4)}</td>
  <td class="${cc(s.mean, s.threshold)}">${fmt(s.mean)}</td>
  <td class="${cc(s.median, s.threshold)}">${fmt(s.median)}</td>
  <td class="${cc(s.p95, s.threshold)}">${fmt(s.p95)}</td>
  <td class="${cc(s.max, s.threshold)}">${fmt(s.max)}</td>
  <td class="${pr(s.passRate)}"><span class="pct">${(s.passRate * 100).toFixed(1)}%</span></td>
</tr>`).join('\n')}
</tbody>
</table>

<h3>Per-Dataset: Promax Loadings + Fit Indices</h3>
<table>
<thead><tr><th>#</th><th>n</th><th>p</th><th>k</th><th>Load MAE</th><th>Load Max</th><th>Comm MAE</th><th>chi2 |D|</th><th>RMSEA |D|</th><th>CFI |D|</th><th>TLI |D|</th><th>SRMR |D|</th></tr></thead>
<tbody>
${okPromax.map(c => `<tr>
  <td>${c.id}</td><td>${c.n}</td><td>${c.p}</td><td>${c.k}</td>
  <td class="${cc(c.loadMAE, 0.05)}">${fmt(c.loadMAE, 4)}</td>
  <td class="${cc(c.loadMax, 0.20)}">${fmt(c.loadMax, 4)}</td>
  <td class="${cc(c.commMAE, 0.05)}">${fmt(c.commMAE, 4)}</td>
  <td class="${cc(c.chiSqDiff, 10)}">${fmt(c.chiSqDiff, 4)}</td>
  <td class="${cc(c.rmseaDiff, 0.02)}">${fmt(c.rmseaDiff, 4)}</td>
  <td class="${cc(c.cfiDiff, 0.02)}">${fmt(c.cfiDiff, 4)}</td>
  <td class="${cc(c.tliDiff, 0.02)}">${fmt(c.tliDiff, 4)}</td>
  <td class="${cc(c.srmrDiff, 0.01)}">${fmt(c.srmrDiff, 4)}</td>
</tr>`).join('\n')}
</tbody>
</table>

<!-- ═══ SECTION 2: GEOMIN ═══ -->
<h2>2. Geomin Rotation — Aggregate (vs R GPArotation::geominQ, delta=0.01)</h2>
<p class="section-note">Geomin is an oblique rotation. R reference: geominQ() applied to ML-extracted unrotated loadings.</p>
<table>
<thead><tr><th>Metric</th><th>Threshold</th><th>Mean e</th><th>Median e</th><th>P95 e</th><th>Max e</th><th>Pass Rate</th></tr></thead>
<tbody>
${geominStats.map(s => `<tr>
  <td>${s.label}</td>
  <td class="threshold">${fmt(s.threshold, 4)}</td>
  <td class="${cc(s.mean, s.threshold)}">${fmt(s.mean)}</td>
  <td class="${cc(s.median, s.threshold)}">${fmt(s.median)}</td>
  <td class="${cc(s.p95, s.threshold)}">${fmt(s.p95)}</td>
  <td class="${cc(s.max, s.threshold)}">${fmt(s.max)}</td>
  <td class="${pr(s.passRate)}"><span class="pct">${(s.passRate * 100).toFixed(1)}%</span></td>
</tr>`).join('\n')}
</tbody>
</table>

<h3>Per-Dataset: Geomin Loadings</h3>
<table>
<thead><tr><th>#</th><th>n</th><th>p</th><th>k</th><th>Load MAE</th><th>Status</th></tr></thead>
<tbody>
${geominResults.map(r => `<tr>
  <td>${r.id}</td><td>${r.n}</td><td>${r.p}</td><td>${r.k}</td>
  <td class="${cc(r.loadMAE, 0.05)}">${fmt(r.loadMAE, 4)}</td>
  <td><span class="badge ${r.pass ? 'ok' : 'err'}">${r.pass ? 'PASS' : 'FAIL'}</span></td>
</tr>`).join('\n')}
</tbody>
</table>

<!-- ═══ SECTION 3: REAL DATASET ═══ -->
<h2>3. Real Dataset (${realData ? `${realData.length}x${realHeader.length}` : 'N/A'})</h2>
<p class="section-note">${realData ? `Survey data: ${realHeader.slice(0, 5).join(', ')}... (${realHeader.length} variables). Both promax and geomin tested at multiple k.` : 'Real dataset file not found.'}</p>

${realResults.length > 0 ? `
<table>
<thead><tr><th>Rotation</th><th>k</th><th>Load MAE</th><th>Status</th></tr></thead>
<tbody>
${realResults.map(r => `<tr>
  <td>${r.rotation}</td><td>${r.k}</td>
  <td class="${cc(r.loadMAE, 0.05)}">${fmt(r.loadMAE, 4)}</td>
  <td><span class="badge ${r.pass ? 'ok' : 'err'}">${r.pass ? 'PASS' : 'FAIL'}</span></td>
</tr>`).join('\n')}
</tbody>
</table>

<h3>Sample Loadings Comparison (first 5 variables)</h3>
${realResults.map(r => `
<div style="margin-bottom:16px;">
<div style="font-weight:600;margin-bottom:4px;">${r.rotation} k=${r.k} (MAE=${r.loadMAE.toFixed(4)})</div>
${r.firstVars.map(v => `<div class="real-cmp">
  <span class="varname">${v.name}</span>
  <span class="vals">Carm=[${v.carm.map(x => x.toFixed(4)).join(', ')}]</span><br>
  <span class="varname">&nbsp;</span>
  <span class="vals">R   =[${v.r.map(x => x.toFixed(4)).join(', ')}]</span>
</div>`).join('\n')}
</div>
`).join('\n')}
` : '<p style="color:#8b949e;">No real dataset results available.</p>'}

<!-- ═══ SECTION 4: DIAGNOSTICS ═══ -->
<h2>4. Diagnostics (Eigenvalues, KMO, Bartlett)</h2>
<p class="section-note">Deterministic functions of the correlation matrix — should match to machine precision (~1e-12).</p>
<table>
<thead><tr><th>#</th><th>n</th><th>p</th><th>k</th><th>Eigenvalue MAE</th><th>KMO |D|</th><th>Bartlett chi2 |D|</th></tr></thead>
<tbody>
${okPromax.map(c => `<tr>
  <td>${c.id}</td><td>${c.n}</td><td>${c.p}</td><td>${c.k}</td>
  <td class="${cc(c.eigMAE, 1e-8)}">${fmt(c.eigMAE)}</td>
  <td class="${cc(c.kmoΔ, 1e-8)}">${fmt(c.kmoΔ)}</td>
  <td class="${cc(c.bartΔ, 0.01)}">${fmt(c.bartΔ)}</td>
</tr>`).join('\n')}
</tbody>
</table>

<p style="margin-top:40px;color:#6e7681;font-size:11px;">
  Generated ${timestamp} — Carm vs R psych/GPArotation — ${datasets.length} synthetic datasets + real (${realData ? `${realData.length}x${realHeader.length}` : 'N/A'})
</p>
</body>
</html>`

writeFileSync('validation/reports/fa-crossval-report.html', html)
console.log(`\nFull HTML report saved to validation/reports/fa-crossval-report.html`)
