/**
 * PCA Cross-Validation Report: Carm vs R
 * Covers: Eigenvalues, Variance Explained, Loadings, Scores, Varimax Rotation
 *
 * Reference: R prcomp() + varimax()
 * Run: npx tsx validation/ts-harness/pca-report.ts && open validation/reports/pca-report.html
 */
import { readFileSync, writeFileSync, mkdirSync } from 'fs'
import { runPCA, varimaxRotation } from 'carm'

// ── Types ────────────────────────────────────────────────────────────────

interface PCARef {
  seed: number
  n: number
  p: number
  data: { x1: number[]; x2: number[]; x3: number[]; x4: number[] }
  eigenvalues: number[]
  varianceExplained: number[]
  cumulativeVariance: number[]
  loadings: number[][]     // p × p
  scores: number[][]       // 5 × p (first 5 rows)
  varimax: { k: number; loadings: number[][] }  // p × k
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
    eigenvalues:     m('Eigenvalues MAE',          1e-8,  'Eigendecomposition'),
    var_explained:   m('Variance Explained MAE',   1e-8,  'Eigendecomposition'),
    cum_variance:    m('Cumulative Variance MAE',  1e-8,  'Eigendecomposition'),
    loading_c1:      m('Loading PC1 MAE',          1e-4,  'Loadings'),
    loading_c2:      m('Loading PC2 MAE',          1e-4,  'Loadings'),
    loading_c3:      m('Loading PC3 MAE',          1e-4,  'Loadings'),
    loading_c4:      m('Loading PC4 MAE',          1e-4,  'Loadings'),
    loading_all:     m('Loadings Overall MAE',     1e-4,  'Loadings'),
    score_c1:        m('Score PC1 MAE (5 rows)',   1e-4,  'Scores'),
    score_c2:        m('Score PC2 MAE (5 rows)',   1e-4,  'Scores'),
    score_c3:        m('Score PC3 MAE (5 rows)',   1e-4,  'Scores'),
    score_c4:        m('Score PC4 MAE (5 rows)',   1e-4,  'Scores'),
    score_all:       m('Scores Overall MAE',       1e-4,  'Scores'),
    varimax_all:     m('Varimax Loadings MAE',     1e-3,  'Varimax'),
  }
}

// ── Helpers ──────────────────────────────────────────────────────────────

function push(m: MetricDef, error: number) { m.errors.push(error) }

function mae(a: number[], b: number[]): number {
  const n = Math.min(a.length, b.length)
  if (n === 0) return NaN
  let sum = 0
  for (let i = 0; i < n; i++) sum += Math.abs(a[i]! - b[i]!)
  return sum / n
}

/** Align sign indeterminacy: if dot product of R column j and Carm column j is negative, flip Carm column. */
function alignSigns(rMatrix: number[][], cMatrix: number[][]): number[][] {
  const nRows = cMatrix.length
  const nCols = cMatrix[0]!.length
  const aligned = cMatrix.map(row => [...row])
  for (let j = 0; j < nCols; j++) {
    let dot = 0
    for (let i = 0; i < nRows; i++) dot += (rMatrix[i]?.[j] ?? 0) * (cMatrix[i]?.[j] ?? 0)
    if (dot < 0) {
      for (let i = 0; i < nRows; i++) aligned[i]![j] *= -1
    }
  }
  return aligned
}

/** Column-wise MAE after sign alignment */
function columnMAE(rMatrix: number[][], cMatrix: number[][], col: number): number {
  const n = Math.min(rMatrix.length, cMatrix.length)
  let sum = 0
  for (let i = 0; i < n; i++) sum += Math.abs((rMatrix[i]?.[col] ?? 0) - (cMatrix[i]?.[col] ?? 0))
  return sum / n
}

/** Overall matrix MAE */
function matrixMAE(rM: number[][], cM: number[][]): number {
  const rows = Math.min(rM.length, cM.length)
  if (rows === 0) return NaN
  const cols = Math.min(rM[0]!.length, cM[0]!.length)
  let sum = 0, count = 0
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      sum += Math.abs((rM[i]?.[j] ?? 0) - (cM[i]?.[j] ?? 0))
      count++
    }
  }
  return count > 0 ? sum / count : NaN
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

console.log('Loading PCA reference data...')
const ref: PCARef[] = JSON.parse(readFileSync('validation/data/pca-ref.json', 'utf-8'))
console.log(`Loaded ${ref.length} datasets`)

const metrics = makeDefs()

// ── Per-dataset tracking ─────────────────────────────────────────────────

interface DatasetResult {
  seed: number; n: number; status: 'ok' | 'fail'; error?: string
  eigMAE: number; varExpMAE: number; cumVarMAE: number
  loadMAE: number; scoreMAE: number; varimaxMAE: number
}

const results: DatasetResult[] = []
let nOk = 0
let nFail = 0

for (const r of ref) {
  const dr: Partial<DatasetResult> = { seed: r.seed, n: r.n }
  try {
    // Transpose column arrays to row-oriented data
    const n = r.data.x1.length
    const data: number[][] = []
    for (let i = 0; i < n; i++) {
      data.push([r.data.x1[i]!, r.data.x2[i]!, r.data.x3[i]!, r.data.x4[i]!])
    }

    const cPCA = runPCA(data)

    // ── Eigenvalues ──
    const eigErr = mae(r.eigenvalues, [...cPCA.eigenvalues])
    push(metrics.eigenvalues!, eigErr)
    dr.eigMAE = eigErr

    // ── Variance Explained ──
    // NOTE: R's summary(prcomp)$importance rounds to 5dp (round(vars, 5L)).
    // Recompute full-precision values from the stored eigenvalues.
    const rEigTotal = r.eigenvalues.reduce((a: number, b: number) => a + b, 0)
    const rVarExpFull = r.eigenvalues.map((v: number) => v / rEigTotal)
    const rCumVarFull = rVarExpFull.reduce((acc: number[], v: number, i: number) => {
      acc.push((acc[i - 1] ?? 0) + v)
      return acc
    }, [] as number[])

    const varExpErr = mae(rVarExpFull, [...cPCA.varianceExplained])
    push(metrics.var_explained!, varExpErr)
    dr.varExpMAE = varExpErr

    // ── Cumulative Variance ──
    const cumVarErr = mae(rCumVarFull, [...cPCA.cumulativeVariance])
    push(metrics.cum_variance!, cumVarErr)
    dr.cumVarMAE = cumVarErr

    // ── Loadings (sign-aligned) ──
    // Carm loadings: [nVars × nComponents] — same shape as R's p × p
    const cLoadings = cPCA.loadings.map(row => [...row])
    const alignedLoadings = alignSigns(r.loadings, cLoadings)

    const nComps = Math.min(r.loadings[0]!.length, alignedLoadings[0]!.length)
    for (let j = 0; j < Math.min(nComps, 4); j++) {
      const colErr = columnMAE(r.loadings, alignedLoadings, j)
      const key = `loading_c${j + 1}` as keyof typeof metrics
      if (metrics[key]) push(metrics[key]!, colErr)
    }
    const loadAllErr = matrixMAE(r.loadings, alignedLoadings)
    push(metrics.loading_all!, loadAllErr)
    dr.loadMAE = loadAllErr

    // ── Scores (first 5 rows, sign-aligned) ──
    // cPCA.scores is [nObs × nComponents], r.scores is [5 × p]
    const cScores5 = cPCA.scores.slice(0, 5).map(row => [...row])
    const alignedScores = alignSigns(r.scores, cScores5)

    const nScoreCols = Math.min(r.scores[0]!.length, alignedScores[0]!.length)
    for (let j = 0; j < Math.min(nScoreCols, 4); j++) {
      const colErr = columnMAE(r.scores, alignedScores, j)
      const key = `score_c${j + 1}` as keyof typeof metrics
      if (metrics[key]) push(metrics[key]!, colErr)
    }
    const scoreAllErr = matrixMAE(r.scores, alignedScores)
    push(metrics.score_all!, scoreAllErr)
    dr.scoreMAE = scoreAllErr

    // ── Varimax ──
    const varimaxK = r.varimax.k
    // Extract first varimaxK components from PCA loadings
    const loadingsForVarimax = cPCA.loadings.map(row => [...row].slice(0, varimaxK))
    const { rotatedLoadings } = varimaxRotation(loadingsForVarimax)
    const alignedVarimax = alignSigns(r.varimax.loadings, rotatedLoadings)
    const varimaxErr = matrixMAE(r.varimax.loadings, alignedVarimax)
    push(metrics.varimax_all!, varimaxErr)
    dr.varimaxMAE = varimaxErr

    dr.status = 'ok'
    nOk++
  } catch (err) {
    dr.status = 'fail'
    dr.error = err instanceof Error ? err.message : String(err)
    dr.eigMAE = NaN; dr.varExpMAE = NaN; dr.cumVarMAE = NaN
    dr.loadMAE = NaN; dr.scoreMAE = NaN; dr.varimaxMAE = NaN
    nFail++
  }
  results.push(dr as DatasetResult)
}

// ── Aggregate ────────────────────────────────────────────────────────────

const allStats = Object.values(metrics).filter(m => m.errors.length > 0).map(aggStat)
const categories = [...new Set(allStats.map(s => s.category))]

// ── Console summary ──────────────────────────────────────────────────────

console.log('\n═══ PCA Cross-Validation Summary ═══\n')
console.log(`Datasets: ${ref.length} total, ${nOk} ok, ${nFail} failed\n`)

for (const cat of categories) {
  console.log(`── ${cat} ──`)
  const catStats = allStats.filter(s => s.category === cat)
  for (const s of catStats) {
    const status = s.passRate >= 0.99 ? '✓' : s.passRate >= 0.90 ? '~' : '✗'
    console.log(`  ${status} ${s.label.padEnd(28)} mean=${fmt(s.mean, 8).padStart(14)} max=${fmt(s.max, 8).padStart(14)} pass=${(s.passRate * 100).toFixed(1).padStart(5)}% (n=${s.n})`)
  }
}

const totalPass = allStats.filter(s => s.passRate >= 0.99).length
const totalMetrics = allStats.length
console.log(`\nOverall: ${totalPass}/${totalMetrics} metrics at ≥99% pass rate`)

// ── HTML Report ──────────────────────────────────────────────────────────

const timestamp = new Date().toISOString().slice(0, 19)
const okResults = results.filter(r => r.status === 'ok')

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>PCA Cross-Validation: Carm vs R</title>
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

<h1>PCA Cross-Validation</h1>
<p class="subtitle">Carm (TypeScript) vs R prcomp() + varimax() — ${ref.length} datasets — ${timestamp}</p>

<div class="summary-bar">
  <div class="summary-card"><div class="num">${ref.length}</div><div class="label">Datasets</div></div>
  <div class="summary-card"><div class="num">${totalMetrics}</div><div class="label">Metrics Tested</div></div>
  <div class="summary-card"><div class="num ${totalPass === totalMetrics ? 'pass' : totalPass >= totalMetrics * 0.9 ? 'warn' : 'fail'}">${totalPass}/${totalMetrics}</div><div class="label">≥99% Pass Rate</div></div>
  <div class="summary-card"><div class="num ${nOk === ref.length ? 'pass' : 'fail'}">${nOk}/${ref.length}</div><div class="label">Ran Successfully</div></div>
  <div class="summary-card"><div class="num ${nFail === 0 ? 'pass' : 'fail'}">${nFail}</div><div class="label">Failures</div></div>
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
  <td class="threshold">${fmt(s.threshold, 8)}</td>
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

<h2>Per-Dataset Details (top 50 by max error)</h2>
<p class="section-note">Sorted by largest single error across all metrics.</p>
<table>
<thead><tr><th>Seed</th><th>n</th><th>Status</th><th>Eigenval MAE</th><th>VarExp MAE</th><th>CumVar MAE</th><th>Loadings MAE</th><th>Scores MAE</th><th>Varimax MAE</th></tr></thead>
<tbody>
${results
  .sort((a, b) => {
    const maxA = Math.max(a.eigMAE || 0, a.varExpMAE || 0, a.cumVarMAE || 0, a.loadMAE || 0, a.scoreMAE || 0, a.varimaxMAE || 0)
    const maxB = Math.max(b.eigMAE || 0, b.varExpMAE || 0, b.cumVarMAE || 0, b.loadMAE || 0, b.scoreMAE || 0, b.varimaxMAE || 0)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => `<tr>
  <td>${d.seed}</td>
  <td>${d.n}</td>
  <td><span class="badge ${d.status === 'ok' ? 'ok' : 'err'}">${d.status === 'ok' ? 'OK' : 'FAIL'}</span>${d.error ? ` <span style="color:#6e7681;font-size:11px;">${d.error.slice(0, 60)}</span>` : ''}</td>
  <td class="${cc(d.eigMAE, 1e-8)}">${fmt(d.eigMAE)}</td>
  <td class="${cc(d.varExpMAE, 1e-8)}">${fmt(d.varExpMAE)}</td>
  <td class="${cc(d.cumVarMAE, 1e-8)}">${fmt(d.cumVarMAE)}</td>
  <td class="${cc(d.loadMAE, 1e-4)}">${fmt(d.loadMAE)}</td>
  <td class="${cc(d.scoreMAE, 1e-4)}">${fmt(d.scoreMAE)}</td>
  <td class="${cc(d.varimaxMAE, 1e-3)}">${fmt(d.varimaxMAE)}</td>
</tr>`).join('\n')}
</tbody>
</table>

${nFail > 0 ? `
<h2>Failed Datasets</h2>
<table>
<thead><tr><th>Seed</th><th>n</th><th>Error</th></tr></thead>
<tbody>
${results.filter(r => r.status === 'fail').map(r => `<tr>
  <td>${r.seed}</td><td>${r.n}</td><td style="color:var(--red);">${r.error ?? 'Unknown'}</td>
</tr>`).join('\n')}
</tbody>
</table>
` : ''}

<p style="margin-top:40px;color:#6e7681;font-size:11px;">
  Generated ${timestamp} · Carm vs R prcomp/varimax · ${ref.length} datasets · ${totalMetrics} metrics
</p>

</body>
</html>`

mkdirSync('validation/reports', { recursive: true })
writeFileSync('validation/reports/pca-report.html', html)
console.log(`\nReport saved to validation/reports/pca-report.html`)
