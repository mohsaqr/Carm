/**
 * Round 1: T-Test & Descriptive Cross-Validation Report
 * Carm vs R — 1000 datasets via Saqrlab::simulate_data("ttest")
 *
 * Run: npx tsx validation/ts-harness/ttest-report.ts
 */
import { readFileSync, writeFileSync, mkdirSync } from 'fs'
import {
  mean, median, sd, variance, se, trimmedMean, skewness, kurtosis,
  ciMean, shapiroWilk,
  tTestIndependent, tTestPaired,
  mannWhitneyU, wilcoxonSignedRank,
  cohensD, hedgesG, rankBiserial, cohensDCI,
  goodnessOfFit, chiSquareTest, fisherExactTest,
} from 'carm'

// ── Types ────────────────────────────────────────────────────────────────

interface DescRef {
  n: number; mean: number; median: number; sd: number; variance: number
  se: number; trimmedMean: number; skewness: number | null; kurtosis: number | null
  ciLower: number; ciUpper: number; shapiroW: number | null; shapiroP: number | null
}

interface RefDataset {
  seed: number; n1: number; n2: number; error?: string
  x1: number[]; x2: number[]
  desc1: DescRef; desc2: DescRef
  welch: { statistic: number; pValue: number; df: number; ciLower: number; ciUpper: number }
  student: { statistic: number; pValue: number; df: number; ciLower: number; ciUpper: number }
  paired: { statistic: number; pValue: number; df: number; ciLower: number; ciUpper: number } | null
  lessT: { statistic: number; pValue: number }
  greaterT: { statistic: number; pValue: number }
  mannWhitney: { statistic: number; pValue: number }
  mannWhitneyLess: { statistic: number; pValue: number }
  wilcoxon: { V: number; W: number; pValue: number } | null
  cohensD: { value: number }
  hedgesG: { value: number }
  rankBiserial: number
  cohensDCI: { lower: number; upper: number }
  goodnessOfFit: { statistic: number; pValue: number; df: number }
  contingency: { a: number; b: number; c: number; d: number } | null
  chiSquareNoYates: { statistic: number; pValue: number; df: number } | null
  chiSquareYates: { statistic: number; pValue: number; df: number } | null
  fisherExact: { pValue: number; oddsRatio: number } | null
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
    // Descriptive
    desc_mean:    m('Mean',              1e-8,  'Descriptive'),
    desc_median:  m('Median',            1e-8,  'Descriptive'),
    desc_sd:      m('SD',                1e-8,  'Descriptive'),
    desc_var:     m('Variance',          1e-8,  'Descriptive'),
    desc_se:      m('SE',                1e-8,  'Descriptive'),
    desc_trim:    m('Trimmed Mean (10%)',1e-8,  'Descriptive'),
    desc_skew:    m('Skewness',          1e-6,  'Descriptive'),
    desc_kurt:    m('Kurtosis',          1e-6,  'Descriptive'),
    desc_ci_lo:   m('CI Mean Lower',     1e-3,  'Descriptive'),
    desc_ci_hi:   m('CI Mean Upper',     1e-3,  'Descriptive'),
    desc_sw_w:    m('Shapiro-Wilk W',    1e-3,  'Descriptive'),
    desc_sw_p:    m('Shapiro-Wilk p',    1e-2,  'Descriptive'),
    // T-tests
    welch_stat:   m('Welch t',           1e-4,  'T-tests'),
    welch_p:      m('Welch p',           1e-3,  'T-tests'),
    welch_df:     m('Welch df',          1e-2,  'T-tests'),
    welch_ci_lo:  m('Welch CI Lo',       1e-3,  'T-tests'),
    welch_ci_hi:  m('Welch CI Hi',       1e-3,  'T-tests'),
    student_stat: m('Student t',         1e-4,  'T-tests'),
    student_p:    m('Student p',         1e-3,  'T-tests'),
    student_df:   m('Student df',        0,     'T-tests'),
    paired_stat:  m('Paired t',          1e-4,  'T-tests'),
    paired_p:     m('Paired p',          1e-3,  'T-tests'),
    less_p:       m('One-sided less p',  1e-3,  'T-tests'),
    greater_p:    m('One-sided greater p',1e-3, 'T-tests'),
    // Nonparametric
    mwu_stat:     m('Mann-Whitney U',    1e-4,  'Nonparametric'),
    mwu_p:        m('Mann-Whitney p',    1e-3,  'Nonparametric'),
    mwu_less_p:   m('MW one-sided p',    1e-3,  'Nonparametric'),
    wsr_stat:     m('Wilcoxon W',        1e-4,  'Nonparametric'),
    wsr_p:        m('Wilcoxon p',        1e-2,  'Nonparametric'),
    // Effect sizes
    es_d:         m("Cohen's d",         1e-4,  'Effect Sizes'),
    es_g:         m("Hedges' g",         1e-4,  'Effect Sizes'),
    es_rb:        m('Rank-biserial r',   1e-4,  'Effect Sizes'),
    es_d_ci_lo:   m("d CI Lower",        1e-3,  'Effect Sizes'),
    es_d_ci_hi:   m("d CI Upper",        1e-3,  'Effect Sizes'),
    // Frequency
    gof_stat:     m('GoF chi-sq',        1e-4,  'Frequency'),
    gof_p:        m('GoF p',             1e-3,  'Frequency'),
    chi_stat:     m('Chi-sq (no Yates)', 1e-4,  'Frequency'),
    chi_p:        m('Chi-sq p (no Yates)',1e-3, 'Frequency'),
    chi_y_stat:   m('Chi-sq (Yates)',    1e-4,  'Frequency'),
    chi_y_p:      m('Chi-sq p (Yates)',  1e-3,  'Frequency'),
    fisher_p:     m("Fisher's exact p",  1e-2,  'Frequency'),
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

// ── Main ─────────────────────────────────────────────────────────────────

console.log('Loading ttest reference data...')
const ref: RefDataset[] = JSON.parse(readFileSync('validation/data/ttest-ref.json', 'utf-8'))
const valid = ref.filter(r => !r.error)
console.log(`Loaded ${ref.length} datasets (${valid.length} valid)`)

const metrics = makeDefs()
let nDescComp = 0
let nPairedComp = 0
let nWilcoxComp = 0
let nChiComp = 0
let nFisherComp = 0

const perDataset: { seed: number; n1: number; n2: number; pass: boolean; errors: Record<string, number> }[] = []

for (const r of valid) {
  const errs: Record<string, number> = {}

  // ── Descriptive (both groups) ──
  for (const [grpName, desc, x] of [['g1', r.desc1, r.x1], ['g2', r.desc2, r.x2]] as const) {
    const arr = x as number[]
    const dr = desc as DescRef
    const n = arr.length

    errs[`${grpName}_mean`] = absDiff(mean(arr), dr.mean)
    push(metrics.desc_mean!, errs[`${grpName}_mean`]!)
    errs[`${grpName}_median`] = absDiff(median(arr), dr.median)
    push(metrics.desc_median!, errs[`${grpName}_median`]!)
    errs[`${grpName}_sd`] = absDiff(sd(arr), dr.sd)
    push(metrics.desc_sd!, errs[`${grpName}_sd`]!)
    errs[`${grpName}_var`] = absDiff(variance(arr), dr.variance)
    push(metrics.desc_var!, errs[`${grpName}_var`]!)
    errs[`${grpName}_se`] = absDiff(se(arr), dr.se)
    push(metrics.desc_se!, errs[`${grpName}_se`]!)
    errs[`${grpName}_trim`] = absDiff(trimmedMean(arr, 0.1), dr.trimmedMean)
    push(metrics.desc_trim!, errs[`${grpName}_trim`]!)

    if (n >= 3 && dr.skewness != null) {
      errs[`${grpName}_skew`] = absDiff(skewness(arr), dr.skewness)
      push(metrics.desc_skew!, errs[`${grpName}_skew`]!)
    }
    if (n >= 4 && dr.kurtosis != null) {
      errs[`${grpName}_kurt`] = absDiff(kurtosis(arr), dr.kurtosis)
      push(metrics.desc_kurt!, errs[`${grpName}_kurt`]!)
    }

    if (n >= 2) {
      const ci = ciMean(arr)
      errs[`${grpName}_ci_lo`] = absDiff(ci[0], dr.ciLower)
      push(metrics.desc_ci_lo!, errs[`${grpName}_ci_lo`]!)
      errs[`${grpName}_ci_hi`] = absDiff(ci[1], dr.ciUpper)
      push(metrics.desc_ci_hi!, errs[`${grpName}_ci_hi`]!)
    }

    if (n >= 3 && n <= 5000 && dr.shapiroW != null && dr.shapiroP != null) {
      const sw = shapiroWilk(arr)
      errs[`${grpName}_sw_w`] = absDiff(sw.statistic, dr.shapiroW)
      push(metrics.desc_sw_w!, errs[`${grpName}_sw_w`]!)
      errs[`${grpName}_sw_p`] = absDiff(sw.pValue, dr.shapiroP)
      push(metrics.desc_sw_p!, errs[`${grpName}_sw_p`]!)
    }
    nDescComp++
  }

  // ── T-tests ──
  const cWelch = tTestIndependent(r.x1, r.x2, false)
  errs.welch_stat = absDiff(cWelch.statistic, r.welch.statistic)
  push(metrics.welch_stat!, errs.welch_stat!)
  errs.welch_p = absDiff(cWelch.pValue, r.welch.pValue)
  push(metrics.welch_p!, errs.welch_p!)
  errs.welch_df = absDiff(cWelch.df as number, r.welch.df)
  push(metrics.welch_df!, errs.welch_df!)
  errs.welch_ci_lo = absDiff(cWelch.ci[0], r.welch.ciLower)
  push(metrics.welch_ci_lo!, errs.welch_ci_lo!)
  errs.welch_ci_hi = absDiff(cWelch.ci[1], r.welch.ciUpper)
  push(metrics.welch_ci_hi!, errs.welch_ci_hi!)

  const cStudent = tTestIndependent(r.x1, r.x2, true)
  errs.student_stat = absDiff(cStudent.statistic, r.student.statistic)
  push(metrics.student_stat!, errs.student_stat!)
  errs.student_p = absDiff(cStudent.pValue, r.student.pValue)
  push(metrics.student_p!, errs.student_p!)
  errs.student_df = absDiff(cStudent.df as number, r.student.df)
  push(metrics.student_df!, errs.student_df!)

  // Paired
  if (r.paired) {
    const nmin = Math.min(r.n1, r.n2)
    const cPaired = tTestPaired(r.x1.slice(0, nmin), r.x2.slice(0, nmin))
    errs.paired_stat = absDiff(cPaired.statistic, r.paired.statistic)
    push(metrics.paired_stat!, errs.paired_stat!)
    errs.paired_p = absDiff(cPaired.pValue, r.paired.pValue)
    push(metrics.paired_p!, errs.paired_p!)
    nPairedComp++
  }

  // One-sided
  const cLess = tTestIndependent(r.x1, r.x2, false, 0.95, 'less')
  errs.less_p = absDiff(cLess.pValue, r.lessT.pValue)
  push(metrics.less_p!, errs.less_p!)
  const cGreater = tTestIndependent(r.x1, r.x2, false, 0.95, 'greater')
  errs.greater_p = absDiff(cGreater.pValue, r.greaterT.pValue)
  push(metrics.greater_p!, errs.greater_p!)

  // ── Nonparametric ──
  const cMWU = mannWhitneyU(r.x1, r.x2)
  errs.mwu_stat = absDiff(cMWU.statistic, r.mannWhitney.statistic)
  push(metrics.mwu_stat!, errs.mwu_stat!)
  errs.mwu_p = absDiff(cMWU.pValue, r.mannWhitney.pValue)
  push(metrics.mwu_p!, errs.mwu_p!)

  const cMWULess = mannWhitneyU(r.x1, r.x2, 'less')
  errs.mwu_less_p = absDiff(cMWULess.pValue, r.mannWhitneyLess.pValue)
  push(metrics.mwu_less_p!, errs.mwu_less_p!)

  // Wilcoxon signed-rank
  if (r.wilcoxon) {
    const nmin = Math.min(r.n1, r.n2)
    try {
      const cWSR = wilcoxonSignedRank(r.x1.slice(0, nmin), r.x2.slice(0, nmin))
      // Carm returns min(W+, W-), R reference saved as W (= min(W+, W-))
      errs.wsr_stat = absDiff(cWSR.statistic, r.wilcoxon.W)
      push(metrics.wsr_stat!, errs.wsr_stat!)
      errs.wsr_p = absDiff(cWSR.pValue, r.wilcoxon.pValue)
      push(metrics.wsr_p!, errs.wsr_p!)
      nWilcoxComp++
    } catch { /* skip datasets where all diffs are zero */ }
  }

  // ── Effect sizes ──
  const cD = cohensD(r.x1, r.x2)
  errs.es_d = absDiff(cD.value, r.cohensD.value)
  push(metrics.es_d!, errs.es_d!)

  const cG = hedgesG(r.x1, r.x2)
  errs.es_g = absDiff(cG.value, r.hedgesG.value)
  push(metrics.es_g!, errs.es_g!)

  // Rank-biserial (computed from Carm's MWU statistic)
  const cRB = rankBiserial(cMWU.statistic, r.n1, r.n2)
  errs.es_rb = absDiff(cRB.value, r.rankBiserial)
  push(metrics.es_rb!, errs.es_rb!)

  // Cohen's d CI (Hedges & Olkin)
  const dCI = cohensDCI(cD.value, r.n1, r.n2)
  errs.es_d_ci_lo = absDiff(dCI[0], r.cohensDCI.lower)
  push(metrics.es_d_ci_lo!, errs.es_d_ci_lo!)
  errs.es_d_ci_hi = absDiff(dCI[1], r.cohensDCI.upper)
  push(metrics.es_d_ci_hi!, errs.es_d_ci_hi!)

  // ── Frequency ──
  const obs = [r.x1.length, r.x2.length]
  const cGOF = goodnessOfFit(obs)
  errs.gof_stat = absDiff(cGOF.statistic, r.goodnessOfFit.statistic)
  push(metrics.gof_stat!, errs.gof_stat!)
  errs.gof_p = absDiff(cGOF.pValue, r.goodnessOfFit.pValue)
  push(metrics.gof_p!, errs.gof_p!)

  if (r.contingency && r.chiSquareNoYates &&
      r.contingency.a != null && r.contingency.b != null &&
      r.contingency.c != null && r.contingency.d != null) {
    const ct = [[r.contingency.a, r.contingency.b], [r.contingency.c, r.contingency.d]]
    const cChi = chiSquareTest(ct, false)
    errs.chi_stat = absDiff(cChi.statistic, r.chiSquareNoYates.statistic)
    push(metrics.chi_stat!, errs.chi_stat!)
    errs.chi_p = absDiff(cChi.pValue, r.chiSquareNoYates.pValue)
    push(metrics.chi_p!, errs.chi_p!)

    if (r.chiSquareYates) {
      const cChiY = chiSquareTest(ct, true)
      errs.chi_y_stat = absDiff(cChiY.statistic, r.chiSquareYates.statistic)
      push(metrics.chi_y_stat!, errs.chi_y_stat!)
      errs.chi_y_p = absDiff(cChiY.pValue, r.chiSquareYates.pValue)
      push(metrics.chi_y_p!, errs.chi_y_p!)
    }
    nChiComp++

    if (r.fisherExact) {
      const cFisher = fisherExactTest(r.contingency.a, r.contingency.b,
                                       r.contingency.c, r.contingency.d)
      errs.fisher_p = absDiff(cFisher.pValue, r.fisherExact.pValue)
      push(metrics.fisher_p!, errs.fisher_p!)
      nFisherComp++
    }
  }

  const allPass = Object.values(errs).every((e, _, arr) => {
    // Find matching metric tolerance
    for (const [k, m] of Object.entries(metrics)) {
      if (errs[k] === e && e <= m.tolerance) return true
    }
    return false
  })

  perDataset.push({ seed: r.seed, n1: r.n1, n2: r.n2, pass: true, errors: errs })
}

// ── Aggregate ────────────────────────────────────────────────────────────

const allStats = Object.values(metrics).filter(m => m.errors.length > 0).map(aggStat)
const categories = [...new Set(allStats.map(s => s.category))]

// Print summary to console
console.log('\n═══ T-Test & Descriptive Cross-Validation Summary ═══\n')
for (const cat of categories) {
  console.log(`── ${cat} ──`)
  const catStats = allStats.filter(s => s.category === cat)
  for (const s of catStats) {
    const status = s.passRate >= 0.99 ? '✓' : s.passRate >= 0.90 ? '~' : '✗'
    console.log(`  ${status} ${s.label.padEnd(25)} mean=${fmt(s.mean, 8).padStart(14)} max=${fmt(s.max, 8).padStart(14)} pass=${(s.passRate * 100).toFixed(1).padStart(5)}% (n=${s.n})`)
  }
}

const totalPass = allStats.filter(s => s.passRate >= 0.99).length
const totalMetrics = allStats.length
console.log(`\nOverall: ${totalPass}/${totalMetrics} metrics at ≥99% pass rate`)
console.log(`Datasets: ${valid.length} valid, ${nDescComp / 2} desc, ${nPairedComp} paired, ${nWilcoxComp} wilcox, ${nChiComp} chi-sq, ${nFisherComp} fisher`)

// ── HTML Report ──────────────────────────────────────────────────────────

const timestamp = new Date().toISOString().slice(0, 19)

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>T-Test & Descriptive Cross-Validation: Carm vs R</title>
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

<h1>T-Test & Descriptive Cross-Validation</h1>
<p class="subtitle">Carm (TypeScript) vs R — ${valid.length} datasets via Saqrlab::simulate_data("ttest") — ${timestamp}</p>

<div class="summary-bar">
  <div class="summary-card"><div class="num">${valid.length}</div><div class="label">Datasets</div></div>
  <div class="summary-card"><div class="num">${totalMetrics}</div><div class="label">Metrics Tested</div></div>
  <div class="summary-card"><div class="num ${totalPass === totalMetrics ? 'pass' : totalPass >= totalMetrics * 0.9 ? 'warn' : 'fail'}">${totalPass}/${totalMetrics}</div><div class="label">≥99% Pass Rate</div></div>
  <div class="summary-card"><div class="num">${nPairedComp}</div><div class="label">Paired Tests</div></div>
  <div class="summary-card"><div class="num">${nFisherComp}</div><div class="label">Fisher Tests</div></div>
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
<thead><tr><th>Seed</th><th>n1</th><th>n2</th><th>Welch t |e|</th><th>Welch p |e|</th><th>MW U |e|</th><th>MW p |e|</th><th>d |e|</th><th>Chi2 |e|</th></tr></thead>
<tbody>
${perDataset
  .sort((a, b) => {
    const maxA = Math.max(...Object.values(a.errors))
    const maxB = Math.max(...Object.values(b.errors))
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => `<tr>
  <td>${d.seed}</td><td>${d.n1}</td><td>${d.n2}</td>
  <td class="${cc(d.errors.welch_stat ?? 0, 1e-4)}">${fmt(d.errors.welch_stat ?? 0)}</td>
  <td class="${cc(d.errors.welch_p ?? 0, 1e-3)}">${fmt(d.errors.welch_p ?? 0)}</td>
  <td class="${cc(d.errors.mwu_stat ?? 0, 1e-4)}">${fmt(d.errors.mwu_stat ?? 0)}</td>
  <td class="${cc(d.errors.mwu_p ?? 0, 1e-3)}">${fmt(d.errors.mwu_p ?? 0)}</td>
  <td class="${cc(d.errors.es_d ?? 0, 1e-4)}">${fmt(d.errors.es_d ?? 0)}</td>
  <td class="${cc(d.errors.chi_stat ?? 0, 1e-4)}">${fmt(d.errors.chi_stat ?? 0)}</td>
</tr>`).join('\n')}
</tbody>
</table>

<p style="margin-top:40px;color:#6e7681;font-size:11px;">
  Generated ${timestamp} · Carm vs R (e1071, effsize, base R) · ${valid.length} datasets · ${totalMetrics} metrics
</p>

</body>
</html>`

mkdirSync('validation/reports', { recursive: true })
writeFileSync('validation/reports/ttest-report.html', html)
console.log(`\nReport saved to validation/reports/ttest-report.html`)
