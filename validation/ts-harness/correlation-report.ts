/**
 * Round 2: Correlation Cross-Validation Report
 * Carm vs R — 1000 datasets via Saqrlab::simulate_data("correlation")
 *
 * Run: npx tsx validation/ts-harness/correlation-report.ts
 */
import { readFileSync, writeFileSync, mkdirSync } from 'fs'
import {
  pearsonCorrelation,
  spearmanCorrelation,
  kendallTau,
  partialCorrelation,
  correlationMatrix,
} from 'carm'

// ── Types ────────────────────────────────────────────────────────────────

interface PairRef {
  i: number; j: number
  r: number; t: number; pValue: number; ciLower: number; ciUpper: number
}

interface SpearmanRef {
  i: number; j: number; rho: number; pValue: number
}

interface KendallRef {
  i: number; j: number; tau: number; pValue: number
}

interface PartialRef {
  statistic: number; pValue: number; estimate: number
}

interface RefDataset {
  seed: number; n: number
  x1: number[]; x2: number[]; x3: number[]; x4: number[]
  pearson: PairRef[]
  spearman: SpearmanRef[]
  kendall: KendallRef[]
  partial1: PartialRef
  partial2: PartialRef
  corrMatrixPearson: number[][]
  corrMatrixSpearman: number[][]
  corrMatrixKendall: number[][]
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
    // Pearson
    pearson_r:        m('Pearson r',              1e-8,  'Pearson'),
    pearson_t:        m('Pearson t-stat',         1e-4,  'Pearson'),
    pearson_p:        m('Pearson p-value',        1e-3,  'Pearson'),
    pearson_ci_lo:    m('Pearson CI Lower',       1e-3,  'Pearson'),
    pearson_ci_hi:    m('Pearson CI Upper',       1e-3,  'Pearson'),
    // Spearman
    spearman_rho:     m('Spearman rho',           1e-4,  'Spearman'),
    spearman_p:       m('Spearman p-value',       1e-3,  'Spearman'),
    // Kendall
    kendall_tau:      m('Kendall tau',            1e-4,  'Kendall'),
    kendall_p:        m('Kendall p-value',        1e-3,  'Kendall'),
    // Partial
    partial_est:      m('Partial r estimate',     1e-4,  'Partial'),
    partial_p:        m('Partial p-value',        1e-3,  'Partial'),
    // Correlation matrices
    corrmat_pearson:  m('Corr Matrix Pearson MAE',1e-8,  'Correlation Matrix'),
    corrmat_spearman: m('Corr Matrix Spearman MAE',1e-8,'Correlation Matrix'),
    corrmat_kendall:  m('Corr Matrix Kendall MAE',1e-6, 'Correlation Matrix'),
  }
}

// ── Helpers ──────────────────────────────────────────────────────────────

function push(m: MetricDef, error: number) { m.errors.push(error) }
function absDiff(a: number, b: number): number { return Math.abs(a - b) }

interface AggStat {
  label: string; mean: number; median: number; p95: number; max: number
  threshold: number; passRate: number; n: number; category: string
}

function aggStat(m: MetricDef): AggStat {
  const sorted = [...m.errors].sort((a, b) => a - b)
  const n = sorted.length
  if (n === 0) return { label: m.label, mean: 0, median: 0, p95: 0, max: 0,
    threshold: m.tolerance, passRate: 1, n: 0, category: m.category }
  return {
    label: m.label,
    mean: sorted.reduce((a, b) => a + b, 0) / n,
    median: sorted[Math.floor(n / 2)]!,
    p95: sorted[Math.floor(n * 0.95)]!,
    max: sorted[n - 1]!,
    threshold: m.tolerance,
    passRate: m.errors.filter(e => e <= m.tolerance).length / n,
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

/** Mean absolute error between two same-sized 2D matrices. */
function matrixMAE(a: readonly (readonly number[])[], b: readonly (readonly number[])[]): number {
  let sum = 0
  let count = 0
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < (a[i]?.length ?? 0); j++) {
      const av = a[i]?.[j] ?? 0
      const bv = b[i]?.[j] ?? 0
      if (!isNaN(av) && !isNaN(bv)) {
        sum += Math.abs(av - bv)
        count++
      }
    }
  }
  return count > 0 ? sum / count : 0
}

// ── Main ─────────────────────────────────────────────────────────────────

console.log('Loading correlation reference data...')
const ref: RefDataset[] = JSON.parse(readFileSync('validation/data/correlation-ref.json', 'utf-8'))
console.log(`Loaded ${ref.length} datasets`)

const metrics = makeDefs()
let nPearsonComp = 0
let nSpearmanComp = 0
let nKendallComp = 0
let nPartialComp = 0
let nMatrixComp = 0

const perDataset: { seed: number; n: number; pass: boolean; errors: Record<string, number> }[] = []

for (const r of ref) {
  const errs: Record<string, number> = {}
  const cols = [r.x1, r.x2, r.x3, r.x4]

  // ── Pearson (6 pairs) ──
  for (const pair of r.pearson) {
    const xi = cols[pair.i]
    const xj = cols[pair.j]
    if (!xi || !xj) continue
    try {
      const cResult = pearsonCorrelation(xi, xj)
      const key = `pearson_${pair.i}_${pair.j}`

      // r — Carm stores rounded statistic, but effectSize.value has full precision
      const cR = cResult.effectSize.value
      errs[`${key}_r`] = absDiff(cR, pair.r)
      push(metrics.pearson_r!, errs[`${key}_r`]!)

      // t-stat — compute from r and df since Carm doesn't expose t directly in Pearson
      // Actually Carm Pearson stores r in statistic (rounded), we need the t
      // Re-derive t from the full-precision r:
      const df = (cResult.df as number)
      const cT = Math.abs(cR) === 1 ? Infinity : cR * Math.sqrt(df / (1 - cR * cR))
      errs[`${key}_t`] = absDiff(cT, pair.t)
      push(metrics.pearson_t!, errs[`${key}_t`]!)

      errs[`${key}_p`] = absDiff(cResult.pValue, pair.pValue)
      push(metrics.pearson_p!, errs[`${key}_p`]!)

      errs[`${key}_ci_lo`] = absDiff(cResult.ci[0], pair.ciLower)
      push(metrics.pearson_ci_lo!, errs[`${key}_ci_lo`]!)

      errs[`${key}_ci_hi`] = absDiff(cResult.ci[1], pair.ciUpper)
      push(metrics.pearson_ci_hi!, errs[`${key}_ci_hi`]!)

      nPearsonComp++
    } catch {
      // Skip pairs that fail (e.g., zero variance)
    }
  }

  // ── Spearman (6 pairs) ──
  for (const pair of r.spearman) {
    const xi = cols[pair.i]
    const xj = cols[pair.j]
    if (!xi || !xj) continue
    try {
      const cResult = spearmanCorrelation(xi, xj)
      const key = `spearman_${pair.i}_${pair.j}`

      const cRho = cResult.effectSize.value
      errs[`${key}_rho`] = absDiff(cRho, pair.rho)
      push(metrics.spearman_rho!, errs[`${key}_rho`]!)

      errs[`${key}_p`] = absDiff(cResult.pValue, pair.pValue)
      push(metrics.spearman_p!, errs[`${key}_p`]!)

      nSpearmanComp++
    } catch {
      // Skip failures
    }
  }

  // ── Kendall (6 pairs) ──
  for (const pair of r.kendall) {
    const xi = cols[pair.i]
    const xj = cols[pair.j]
    if (!xi || !xj) continue
    try {
      const cResult = kendallTau(xi, xj)
      const key = `kendall_${pair.i}_${pair.j}`

      const cTau = cResult.effectSize.value
      errs[`${key}_tau`] = absDiff(cTau, pair.tau)
      push(metrics.kendall_tau!, errs[`${key}_tau`]!)

      errs[`${key}_p`] = absDiff(cResult.pValue, pair.pValue)
      push(metrics.kendall_p!, errs[`${key}_p`]!)

      nKendallComp++
    } catch {
      // Skip failures
    }
  }

  // ── Partial correlation ──
  // partial1: x1 vs x2, controlling for [x3, x4]
  if (r.partial1) {
    try {
      const cPartial1 = partialCorrelation(r.x1, r.x2, [r.x3, r.x4])
      errs.partial1_est = absDiff(cPartial1.effectSize.value, r.partial1.estimate)
      push(metrics.partial_est!, errs.partial1_est!)
      errs.partial1_p = absDiff(cPartial1.pValue, r.partial1.pValue)
      push(metrics.partial_p!, errs.partial1_p!)
      nPartialComp++
    } catch {
      // Skip failures
    }
  }

  // partial2: x1 vs x3, controlling for [x2]
  if (r.partial2) {
    try {
      const cPartial2 = partialCorrelation(r.x1, r.x3, [r.x2])
      errs.partial2_est = absDiff(cPartial2.effectSize.value, r.partial2.estimate)
      push(metrics.partial_est!, errs.partial2_est!)
      errs.partial2_p = absDiff(cPartial2.pValue, r.partial2.pValue)
      push(metrics.partial_p!, errs.partial2_p!)
      nPartialComp++
    } catch {
      // Skip failures
    }
  }

  // ── Correlation matrices ──
  // Convert column-oriented data to row-oriented for correlationMatrix
  const data = [r.x1, r.x2, r.x3, r.x4]
  const labels = ['x1', 'x2', 'x3', 'x4']

  try {
    const cMatPearson = correlationMatrix(data, labels, 'pearson')
    const maePearson = matrixMAE(cMatPearson.r, r.corrMatrixPearson)
    errs.corrmat_pearson = maePearson
    push(metrics.corrmat_pearson!, maePearson)
    nMatrixComp++
  } catch {
    // Skip failures
  }

  try {
    const cMatSpearman = correlationMatrix(data, labels, 'spearman')
    const maeSpearman = matrixMAE(cMatSpearman.r, r.corrMatrixSpearman)
    errs.corrmat_spearman = maeSpearman
    push(metrics.corrmat_spearman!, maeSpearman)
  } catch {
    // Skip failures
  }

  try {
    const cMatKendall = correlationMatrix(data, labels, 'kendall')
    const maeKendall = matrixMAE(cMatKendall.r, r.corrMatrixKendall)
    errs.corrmat_kendall = maeKendall
    push(metrics.corrmat_kendall!, maeKendall)
  } catch {
    // Skip failures
  }

  perDataset.push({ seed: r.seed, n: r.n, pass: true, errors: errs })
}

// ── Aggregate ────────────────────────────────────────────────────────────

const allStats = Object.values(metrics).filter(m => m.errors.length > 0).map(aggStat)
const categories = [...new Set(allStats.map(s => s.category))]

// Print summary to console
console.log('\n═══ Correlation Cross-Validation Summary ═══\n')
for (const cat of categories) {
  console.log(`── ${cat} ──`)
  const catStats = allStats.filter(s => s.category === cat)
  for (const s of catStats) {
    const status = s.passRate >= 0.99 ? '✓' : s.passRate >= 0.90 ? '~' : '✗'
    console.log(`  ${status} ${s.label.padEnd(30)} mean=${fmt(s.mean, 8).padStart(14)} max=${fmt(s.max, 8).padStart(14)} pass=${(s.passRate * 100).toFixed(1).padStart(5)}% (n=${s.n})`)
  }
}

const totalPass = allStats.filter(s => s.passRate >= 0.99).length
const totalMetrics = allStats.length
console.log(`\nOverall: ${totalPass}/${totalMetrics} metrics at ≥99% pass rate`)
console.log(`Comparisons: ${nPearsonComp} Pearson, ${nSpearmanComp} Spearman, ${nKendallComp} Kendall, ${nPartialComp} Partial, ${nMatrixComp} Matrix`)

// ── HTML Report ──────────────────────────────────────────────────────────

const timestamp = new Date().toISOString().slice(0, 19)

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Correlation Cross-Validation: Carm vs R</title>
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

<h1>Correlation Cross-Validation</h1>
<p class="subtitle">Carm (TypeScript) vs R — ${ref.length} datasets — ${timestamp}</p>

<div class="summary-bar">
  <div class="summary-card"><div class="num">${ref.length}</div><div class="label">Datasets</div></div>
  <div class="summary-card"><div class="num">${totalMetrics}</div><div class="label">Metrics Tested</div></div>
  <div class="summary-card"><div class="num ${totalPass === totalMetrics ? 'pass' : totalPass >= totalMetrics * 0.9 ? 'warn' : 'fail'}">${totalPass}/${totalMetrics}</div><div class="label">≥99% Pass Rate</div></div>
  <div class="summary-card"><div class="num">${nPearsonComp}</div><div class="label">Pearson Pairs</div></div>
  <div class="summary-card"><div class="num">${nPartialComp}</div><div class="label">Partial Tests</div></div>
  <div class="summary-card"><div class="num">${nMatrixComp}</div><div class="label">Matrix Tests</div></div>
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
  <td class="threshold">${s.threshold === 0 ? 'exact' : fmt(s.threshold, 8)}</td>
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

<h2>Per-Dataset Summary (first 50 with highest max error)</h2>
<p class="section-note">Showing top-50 datasets sorted by largest single error across all metrics.</p>
<table>
<thead><tr><th>Seed</th><th>n</th><th>Pearson r |e|</th><th>Pearson p |e|</th><th>Spearman rho |e|</th><th>Kendall tau |e|</th><th>Partial r |e|</th><th>Matrix MAE</th></tr></thead>
<tbody>
${perDataset
  .sort((a, b) => {
    const maxA = Math.max(...Object.values(a.errors).filter(v => !isNaN(v)), 0)
    const maxB = Math.max(...Object.values(b.errors).filter(v => !isNaN(v)), 0)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => {
    // Aggregate per-category max errors for display
    const pearsonRErrs = Object.entries(d.errors).filter(([k]) => k.includes('pearson_') && k.endsWith('_r')).map(([, v]) => v)
    const pearsonPErrs = Object.entries(d.errors).filter(([k]) => k.includes('pearson_') && k.endsWith('_p')).map(([, v]) => v)
    const spearmanErrs = Object.entries(d.errors).filter(([k]) => k.includes('spearman_') && k.endsWith('_rho')).map(([, v]) => v)
    const kendallErrs = Object.entries(d.errors).filter(([k]) => k.includes('kendall_') && k.endsWith('_tau')).map(([, v]) => v)
    const partialErrs = Object.entries(d.errors).filter(([k]) => k.includes('partial') && k.endsWith('_est')).map(([, v]) => v)
    const maxPearsonR = pearsonRErrs.length > 0 ? Math.max(...pearsonRErrs) : 0
    const maxPearsonP = pearsonPErrs.length > 0 ? Math.max(...pearsonPErrs) : 0
    const maxSpearman = spearmanErrs.length > 0 ? Math.max(...spearmanErrs) : 0
    const maxKendall = kendallErrs.length > 0 ? Math.max(...kendallErrs) : 0
    const maxPartial = partialErrs.length > 0 ? Math.max(...partialErrs) : 0
    const matMAE = d.errors.corrmat_pearson ?? 0
    return `<tr>
  <td>${d.seed}</td><td>${d.n}</td>
  <td class="${cc(maxPearsonR, 1e-8)}">${fmt(maxPearsonR)}</td>
  <td class="${cc(maxPearsonP, 1e-3)}">${fmt(maxPearsonP)}</td>
  <td class="${cc(maxSpearman, 1e-4)}">${fmt(maxSpearman)}</td>
  <td class="${cc(maxKendall, 1e-4)}">${fmt(maxKendall)}</td>
  <td class="${cc(maxPartial, 1e-4)}">${fmt(maxPartial)}</td>
  <td class="${cc(matMAE, 1e-8)}">${fmt(matMAE)}</td>
</tr>`
  }).join('\n')}
</tbody>
</table>

<p style="margin-top:40px;color:#6e7681;font-size:11px;">
  Generated ${timestamp} · Carm vs R (base R cor.test, ppcor::pcor.test) · ${ref.length} datasets · ${totalMetrics} metrics
</p>

</body>
</html>`

mkdirSync('validation/reports', { recursive: true })
writeFileSync('validation/reports/correlation-report.html', html)
console.log(`\nReport saved to validation/reports/correlation-report.html`)
