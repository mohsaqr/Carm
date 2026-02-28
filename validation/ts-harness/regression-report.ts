/**
 * Round 4: Regression Cross-Validation Report
 * Carm vs R — 500 datasets via Saqrlab::simulate_data("regression")
 *
 * Run: npx tsx validation/ts-harness/regression-report.ts
 */
import { readFileSync, writeFileSync, mkdirSync } from 'fs'
import {
  linearRegression,
  multipleRegression,
  polynomialRegression,
  logisticRegression,
  regressionDiagnostics,
  median,
} from 'carm'

// ── Types ────────────────────────────────────────────────────────────────

interface SimpleRef {
  intercept: number; slope: number
  interceptSE: number; slopeSE: number
  interceptT: number; slopeT: number
  interceptP: number; slopeP: number
  r2: number; adjR2: number
  fStatistic: number; fPvalue: number
  interceptCILower: number; interceptCIUpper: number
  slopeCILower: number; slopeCIUpper: number
}

interface MultipleRef {
  coefficients: number[]; se: number[]; t: number[]; p: number[]
  r2: number; adjR2: number; fStatistic: number; fPvalue: number
  aic: number; bic: number; coeffNames: string[]
}

interface PolyRef {
  coefficients: number[]; se: number[]; t: number[]; p: number[]
  r2: number; adjR2: number; coeffNames: string[]
}

interface LogisticRef {
  coefficients: number[]; se: number[]; z: number[]; p: number[]
  aic: number; converged: boolean; coeffNames: string[]
}

interface DiagnosticsRef {
  leverage: number[]; cooksD: number[]
  vif: number[]; vifNames: string[]
}

interface RefDataset {
  seed: number; n: number
  y: number[]; x1: number[]; x2: number[]; x3: number[]; x4: number[]
  simple: SimpleRef
  multiple: MultipleRef
  polynomial: PolyRef
  logistic: LogisticRef | null
  diagnostics: DiagnosticsRef
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
    // Simple regression
    simple_intercept:   m('Simple intercept',       1e-3,  'Simple Regression'),
    simple_slope:       m('Simple slope',           1e-3,  'Simple Regression'),
    simple_intercept_se:m('Simple intercept SE',    1e-3,  'Simple Regression'),
    simple_slope_se:    m('Simple slope SE',        1e-3,  'Simple Regression'),
    simple_intercept_t: m('Simple intercept t',     1e-4,  'Simple Regression'),
    simple_slope_t:     m('Simple slope t',         1e-4,  'Simple Regression'),
    simple_intercept_p: m('Simple intercept p',     1e-3,  'Simple Regression'),
    simple_slope_p:     m('Simple slope p',         1e-3,  'Simple Regression'),
    simple_r2:          m('Simple R2',              1e-4,  'Simple Regression'),
    simple_adjr2:       m('Simple adj.R2',          1e-4,  'Simple Regression'),
    simple_f:           m('Simple F-stat',          1e-3,  'Simple Regression'),
    simple_fp:          m('Simple F p-value',       1e-3,  'Simple Regression'),
    simple_ci_int_lo:   m('Simple int. CI Lo',      1e-3,  'Simple Regression'),
    simple_ci_int_hi:   m('Simple int. CI Hi',      1e-3,  'Simple Regression'),
    simple_ci_slp_lo:   m('Simple slope CI Lo',     1e-3,  'Simple Regression'),
    simple_ci_slp_hi:   m('Simple slope CI Hi',     1e-3,  'Simple Regression'),
    // Multiple regression
    multi_coef_mae:     m('Multiple coef MAE',      1e-3,  'Multiple Regression'),
    multi_se_mae:       m('Multiple SE MAE',        1e-3,  'Multiple Regression'),
    multi_t_mae:        m('Multiple t MAE',         1e-3,  'Multiple Regression'),
    multi_p_mae:        m('Multiple p MAE',         1e-3,  'Multiple Regression'),
    multi_r2:           m('Multiple R2',            1e-4,  'Multiple Regression'),
    multi_adjr2:        m('Multiple adj.R2',        1e-4,  'Multiple Regression'),
    multi_f:            m('Multiple F-stat',        1e-3,  'Multiple Regression'),
    multi_aic:          m('Multiple AIC',           1e-1,  'Multiple Regression'),
    multi_bic:          m('Multiple BIC',           1e-1,  'Multiple Regression'),
    // Polynomial regression
    poly_coef_mae:      m('Polynomial coef MAE',    1e-3,  'Polynomial Regression'),
    poly_r2:            m('Polynomial R2',          1e-4,  'Polynomial Regression'),
    poly_adjr2:         m('Polynomial adj.R2',      1e-4,  'Polynomial Regression'),
    // Logistic regression
    log_coef_mae:       m('Logistic coef MAE',      1e-2,  'Logistic Regression'),
    log_se_mae:         m('Logistic SE MAE',        1e-2,  'Logistic Regression'),
    log_z_mae:          m('Logistic z MAE',         1e-2,  'Logistic Regression'),
    log_p_mae:          m('Logistic p MAE',         1e-2,  'Logistic Regression'),
    log_aic:            m('Logistic AIC',           1e-1,  'Logistic Regression'),
    // Diagnostics
    diag_leverage_mae:  m('Leverage MAE',           1e-3,  'Diagnostics'),
    diag_cooks_mae:     m("Cook's D MAE",           1e-3,  'Diagnostics'),
    diag_vif_mae:       m('VIF MAE',                1e-2,  'Diagnostics'),
  }
}

// ── Helpers ──────────────────────────────────────────────────────────────

function push(m: MetricDef, error: number) { m.errors.push(error) }
function absDiff(a: number, b: number): number { return Math.abs(a - b) }

/** Mean absolute error between two number arrays. */
function mae(a: readonly number[], b: readonly number[]): number {
  const n = Math.min(a.length, b.length)
  if (n === 0) return 0
  let sum = 0
  for (let i = 0; i < n; i++) {
    const av = a[i] ?? 0
    const bv = b[i] ?? 0
    if (!isNaN(av) && !isNaN(bv)) {
      sum += Math.abs(av - bv)
    }
  }
  return sum / n
}

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

// ── Main ─────────────────────────────────────────────────────────────────

console.log('Loading regression reference data...')
const ref: RefDataset[] = JSON.parse(readFileSync('validation/data/regression-ref.json', 'utf-8'))
console.log(`Loaded ${ref.length} datasets`)

const metrics = makeDefs()
let nSimpleComp = 0
let nMultiComp = 0
let nPolyComp = 0
let nLogisticComp = 0
let nDiagComp = 0

const perDataset: { seed: number; n: number; pass: boolean; errors: Record<string, number> }[] = []

for (const r of ref) {
  const errs: Record<string, number> = {}

  const predictors = [
    { name: 'x1', values: r.x1 },
    { name: 'x2', values: r.x2 },
    { name: 'x3', values: r.x3 },
    { name: 'x4', values: r.x4 },
  ]

  // ── Simple linear regression ──
  try {
    const cSimple = linearRegression(r.x1, r.y)
    const cInt = cSimple.coefficients[0]!
    const cSlp = cSimple.coefficients[1]!

    errs.simple_intercept = absDiff(cInt.estimate, r.simple.intercept)
    push(metrics.simple_intercept!, errs.simple_intercept!)
    errs.simple_slope = absDiff(cSlp.estimate, r.simple.slope)
    push(metrics.simple_slope!, errs.simple_slope!)

    errs.simple_intercept_se = absDiff(cInt.se, r.simple.interceptSE)
    push(metrics.simple_intercept_se!, errs.simple_intercept_se!)
    errs.simple_slope_se = absDiff(cSlp.se, r.simple.slopeSE)
    push(metrics.simple_slope_se!, errs.simple_slope_se!)

    errs.simple_intercept_t = absDiff(cInt.tValue, r.simple.interceptT)
    push(metrics.simple_intercept_t!, errs.simple_intercept_t!)
    errs.simple_slope_t = absDiff(cSlp.tValue, r.simple.slopeT)
    push(metrics.simple_slope_t!, errs.simple_slope_t!)

    errs.simple_intercept_p = absDiff(cInt.pValue, r.simple.interceptP)
    push(metrics.simple_intercept_p!, errs.simple_intercept_p!)
    errs.simple_slope_p = absDiff(cSlp.pValue, r.simple.slopeP)
    push(metrics.simple_slope_p!, errs.simple_slope_p!)

    errs.simple_r2 = absDiff(cSimple.r2, r.simple.r2)
    push(metrics.simple_r2!, errs.simple_r2!)
    errs.simple_adjr2 = absDiff(cSimple.adjR2, r.simple.adjR2)
    push(metrics.simple_adjr2!, errs.simple_adjr2!)

    errs.simple_f = absDiff(cSimple.fStatistic, r.simple.fStatistic)
    push(metrics.simple_f!, errs.simple_f!)
    errs.simple_fp = absDiff(cSimple.fPValue, r.simple.fPvalue)
    push(metrics.simple_fp!, errs.simple_fp!)

    errs.simple_ci_int_lo = absDiff(cInt.ci[0], r.simple.interceptCILower)
    push(metrics.simple_ci_int_lo!, errs.simple_ci_int_lo!)
    errs.simple_ci_int_hi = absDiff(cInt.ci[1], r.simple.interceptCIUpper)
    push(metrics.simple_ci_int_hi!, errs.simple_ci_int_hi!)
    errs.simple_ci_slp_lo = absDiff(cSlp.ci[0], r.simple.slopeCILower)
    push(metrics.simple_ci_slp_lo!, errs.simple_ci_slp_lo!)
    errs.simple_ci_slp_hi = absDiff(cSlp.ci[1], r.simple.slopeCIUpper)
    push(metrics.simple_ci_slp_hi!, errs.simple_ci_slp_hi!)

    nSimpleComp++
  } catch {
    // Skip simple regression failures
  }

  // ── Multiple regression ──
  try {
    const cMulti = multipleRegression(r.y, predictors)
    const cCoefs = cMulti.coefficients.map(c => c.estimate)
    const cSEs = cMulti.coefficients.map(c => c.se)
    const cTs = cMulti.coefficients.map(c => c.tValue)
    const cPs = cMulti.coefficients.map(c => c.pValue)

    errs.multi_coef_mae = mae(cCoefs, r.multiple.coefficients)
    push(metrics.multi_coef_mae!, errs.multi_coef_mae!)
    errs.multi_se_mae = mae(cSEs, r.multiple.se)
    push(metrics.multi_se_mae!, errs.multi_se_mae!)
    errs.multi_t_mae = mae(cTs, r.multiple.t)
    push(metrics.multi_t_mae!, errs.multi_t_mae!)
    errs.multi_p_mae = mae(cPs, r.multiple.p)
    push(metrics.multi_p_mae!, errs.multi_p_mae!)

    errs.multi_r2 = absDiff(cMulti.r2, r.multiple.r2)
    push(metrics.multi_r2!, errs.multi_r2!)
    errs.multi_adjr2 = absDiff(cMulti.adjR2, r.multiple.adjR2)
    push(metrics.multi_adjr2!, errs.multi_adjr2!)
    errs.multi_f = absDiff(cMulti.fStatistic, r.multiple.fStatistic)
    push(metrics.multi_f!, errs.multi_f!)
    errs.multi_aic = absDiff(cMulti.aic, r.multiple.aic)
    push(metrics.multi_aic!, errs.multi_aic!)
    errs.multi_bic = absDiff(cMulti.bic, r.multiple.bic)
    push(metrics.multi_bic!, errs.multi_bic!)

    nMultiComp++

    // ── Diagnostics (uses multiple regression result) ──
    try {
      const cDiag = regressionDiagnostics(cMulti, predictors)

      // Reference leverage/cooksD only has first 10 values
      const refLev = r.diagnostics.leverage
      const refCook = r.diagnostics.cooksD
      const nDiag = refLev.length  // typically 10
      const cLev = cDiag.leverage.slice(0, nDiag)
      const cCook = (cDiag.cooksDistance as readonly number[]).slice(0, nDiag)

      errs.diag_leverage_mae = mae([...cLev], refLev)
      push(metrics.diag_leverage_mae!, errs.diag_leverage_mae!)
      errs.diag_cooks_mae = mae([...cCook], refCook)
      push(metrics.diag_cooks_mae!, errs.diag_cooks_mae!)

      // VIF
      errs.diag_vif_mae = mae([...cDiag.vif], r.diagnostics.vif)
      push(metrics.diag_vif_mae!, errs.diag_vif_mae!)

      nDiagComp++
    } catch {
      // Skip diagnostics failures
    }
  } catch {
    // Skip multiple regression failures
  }

  // ── Polynomial regression ──
  try {
    const cPoly = polynomialRegression(r.x1, r.y, 2)
    const cPolyCoefs = cPoly.coefficients.map(c => c.estimate)

    errs.poly_coef_mae = mae(cPolyCoefs, r.polynomial.coefficients)
    push(metrics.poly_coef_mae!, errs.poly_coef_mae!)
    errs.poly_r2 = absDiff(cPoly.r2, r.polynomial.r2)
    push(metrics.poly_r2!, errs.poly_r2!)
    errs.poly_adjr2 = absDiff(cPoly.adjR2, r.polynomial.adjR2)
    push(metrics.poly_adjr2!, errs.poly_adjr2!)

    nPolyComp++
  } catch {
    // Skip polynomial failures
  }

  // ── Logistic regression ──
  if (r.logistic && r.logistic.converged) {
    try {
      // Binarize y at the median
      const med = median(r.y)
      const yBin = r.y.map(v => v > med ? 1 : 0)

      const cLogistic = logisticRegression(yBin, predictors)
      const cLogCoefs = cLogistic.coefficients.map(c => c.estimate)
      const cLogSEs = cLogistic.coefficients.map(c => c.se)
      const cLogZs = cLogistic.coefficients.map(c => c.tValue) // z-values stored as tValue
      const cLogPs = cLogistic.coefficients.map(c => c.pValue)

      errs.log_coef_mae = mae(cLogCoefs, r.logistic.coefficients)
      push(metrics.log_coef_mae!, errs.log_coef_mae!)
      errs.log_se_mae = mae(cLogSEs, r.logistic.se)
      push(metrics.log_se_mae!, errs.log_se_mae!)
      errs.log_z_mae = mae(cLogZs, r.logistic.z)
      push(metrics.log_z_mae!, errs.log_z_mae!)
      errs.log_p_mae = mae(cLogPs, r.logistic.p)
      push(metrics.log_p_mae!, errs.log_p_mae!)
      errs.log_aic = absDiff(cLogistic.aic, r.logistic.aic)
      push(metrics.log_aic!, errs.log_aic!)

      nLogisticComp++
    } catch {
      // Skip logistic failures
    }
  }

  perDataset.push({ seed: r.seed, n: r.n, pass: true, errors: errs })
}

// ── Aggregate ────────────────────────────────────────────────────────────

const allStats = Object.values(metrics).filter(m => m.errors.length > 0).map(aggStat)
const categories = [...new Set(allStats.map(s => s.category))]

// Print summary to console
console.log('\n═══ Regression Cross-Validation Summary ═══\n')
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
console.log(`Comparisons: ${nSimpleComp} Simple, ${nMultiComp} Multiple, ${nPolyComp} Polynomial, ${nLogisticComp} Logistic, ${nDiagComp} Diagnostics`)

// ── HTML Report ──────────────────────────────────────────────────────────

const timestamp = new Date().toISOString().slice(0, 19)

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Regression Cross-Validation: Carm vs R</title>
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

<h1>Regression Cross-Validation</h1>
<p class="subtitle">Carm (TypeScript) vs R — ${ref.length} datasets — ${timestamp}</p>

<div class="summary-bar">
  <div class="summary-card"><div class="num">${ref.length}</div><div class="label">Datasets</div></div>
  <div class="summary-card"><div class="num">${totalMetrics}</div><div class="label">Metrics Tested</div></div>
  <div class="summary-card"><div class="num ${totalPass === totalMetrics ? 'pass' : totalPass >= totalMetrics * 0.9 ? 'warn' : 'fail'}">${totalPass}/${totalMetrics}</div><div class="label">≥99% Pass Rate</div></div>
  <div class="summary-card"><div class="num">${nSimpleComp}</div><div class="label">Simple Reg</div></div>
  <div class="summary-card"><div class="num">${nLogisticComp}</div><div class="label">Logistic Reg</div></div>
  <div class="summary-card"><div class="num">${nDiagComp}</div><div class="label">Diag. Tests</div></div>
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
<thead><tr><th>Seed</th><th>n</th><th>Simple R2 |e|</th><th>Multi coef MAE</th><th>Multi AIC |e|</th><th>Poly R2 |e|</th><th>Logistic AIC |e|</th><th>VIF MAE</th></tr></thead>
<tbody>
${perDataset
  .sort((a, b) => {
    const maxA = Math.max(...Object.values(a.errors).filter(v => !isNaN(v)), 0)
    const maxB = Math.max(...Object.values(b.errors).filter(v => !isNaN(v)), 0)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => `<tr>
  <td>${d.seed}</td><td>${d.n}</td>
  <td class="${cc(d.errors.simple_r2 ?? 0, 1e-4)}">${fmt(d.errors.simple_r2 ?? 0)}</td>
  <td class="${cc(d.errors.multi_coef_mae ?? 0, 1e-3)}">${fmt(d.errors.multi_coef_mae ?? 0)}</td>
  <td class="${cc(d.errors.multi_aic ?? 0, 1e-1)}">${fmt(d.errors.multi_aic ?? 0)}</td>
  <td class="${cc(d.errors.poly_r2 ?? 0, 1e-4)}">${fmt(d.errors.poly_r2 ?? 0)}</td>
  <td class="${cc(d.errors.log_aic ?? 0, 1e-1)}">${fmt(d.errors.log_aic ?? 0)}</td>
  <td class="${cc(d.errors.diag_vif_mae ?? 0, 1e-2)}">${fmt(d.errors.diag_vif_mae ?? 0)}</td>
</tr>`).join('\n')}
</tbody>
</table>

<p style="margin-top:40px;color:#6e7681;font-size:11px;">
  Generated ${timestamp} · Carm vs R (base R lm, glm, hatvalues, cooks.distance, car::vif) · ${ref.length} datasets · ${totalMetrics} metrics
</p>

</body>
</html>`

mkdirSync('validation/reports', { recursive: true })
writeFileSync('validation/reports/regression-report.html', html)
console.log(`\nReport saved to validation/reports/regression-report.html`)
