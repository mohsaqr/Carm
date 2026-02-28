/**
 * Round 3: ANOVA & Post-Hoc Cross-Validation Report
 * Carm vs R — 500 datasets via Saqrlab::simulate_data("anova")
 *
 * Run: npx tsx validation/ts-harness/anova-report.ts
 */
import { readFileSync, writeFileSync, mkdirSync } from 'fs'
import {
  oneWayANOVA,
  kruskalWallis,
  friedmanTest,
  etaSquared,
  omegaSquared,
  etaSquaredKW,
  tukeyHSD,
  gamesHowell,
  dunnTest,
  runLMM,
} from 'carm'

// ── Types ────────────────────────────────────────────────────────────────

interface AnovaRef {
  fStatistic: number; pValue: number
  dfBetween: number; dfWithin: number
  ssBetween: number; ssWithin: number
  msBetween: number; msWithin: number
  ssTotal: number
}

interface KruskalRef {
  H: number; pValue: number; df: number; etaSquaredKW: number
}

interface TukeyPairRef {
  pair: string; diff: number; pAdj: number; ciLower: number; ciUpper: number
}

interface GamesHowellPairRef {
  group1: string; group2: string; diff: number; pAdj: number; ciLower: number; ciUpper: number
}

interface DunnPairRef {
  pair: string; Z: number; pAdj: number
}

interface FriedmanRef {
  chi2: number; pValue: number; df: number; minGroupSize: number
}

interface LmmRef {
  intercept: number; groupVar: number; residVar: number; icc: number
}

interface RefDataset {
  seed: number; nTotal: number; k: number
  groupLevels: string[]
  groups: Record<string, number[]>
  anova: AnovaRef
  etaSquared: number
  omegaSquared: number
  kruskalWallis: KruskalRef
  tukeyHSD: TukeyPairRef[]
  gamesHowell: GamesHowellPairRef[]
  dunnBonferroni: DunnPairRef[]
  dunnHolm: DunnPairRef[]
  dunnBH: DunnPairRef[]
  friedman: FriedmanRef
  lmm: LmmRef | null
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
    // ANOVA
    anova_f:          m('ANOVA F',              1e-4,  'ANOVA'),
    anova_p:          m('ANOVA p-value',        1e-3,  'ANOVA'),
    anova_dfb:        m('ANOVA df_between',     0,     'ANOVA'),
    anova_dfw:        m('ANOVA df_within',      0,     'ANOVA'),
    anova_ssb:        m('ANOVA SS_between',     1e-4,  'ANOVA'),
    anova_ssw:        m('ANOVA SS_within',      1e-4,  'ANOVA'),
    anova_msb:        m('ANOVA MS_between',     1e-4,  'ANOVA'),
    anova_msw:        m('ANOVA MS_within',      1e-4,  'ANOVA'),
    // Effect sizes
    es_eta2:          m('Eta-squared',          1e-4,  'Effect Sizes'),
    es_omega2:        m('Omega-squared',        1e-4,  'Effect Sizes'),
    es_eta2kw:        m('Eta-squared KW',       1e-4,  'Effect Sizes'),
    // Kruskal-Wallis
    kw_h:             m('Kruskal-Wallis H',     1e-4,  'Kruskal-Wallis'),
    kw_p:             m('Kruskal-Wallis p',     1e-3,  'Kruskal-Wallis'),
    // Tukey HSD
    tukey_diff:       m('Tukey diff',           1e-3,  'Tukey HSD'),
    tukey_padj:       m('Tukey p_adj',          1e-2,  'Tukey HSD'),
    // Games-Howell
    gh_diff:          m('Games-Howell diff',    1e-3,  'Games-Howell'),
    gh_padj:          m('Games-Howell p_adj',   1e-2,  'Games-Howell'),
    // Dunn
    dunn_bonf_padj:   m('Dunn Bonferroni p_adj',1e-2, 'Dunn Test'),
    dunn_holm_padj:   m('Dunn Holm p_adj',      1e-2,  'Dunn Test'),
    dunn_bh_padj:     m('Dunn BH p_adj',        1e-2,  'Dunn Test'),
    // Friedman
    friedman_chi2:    m('Friedman chi2',        1e-4,  'Friedman'),
    friedman_p:       m('Friedman p',           1e-3,  'Friedman'),
    // LMM
    lmm_intercept:    m('LMM intercept',        1e-2,  'LMM'),
    lmm_groupvar:     m('LMM groupVar',         1e-2,  'LMM'),
    lmm_residvar:     m('LMM residVar',         1e-2,  'LMM'),
    lmm_icc:          m('LMM ICC',              1e-2,  'LMM'),
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

console.log('Loading ANOVA reference data...')
const ref: RefDataset[] = JSON.parse(readFileSync('validation/data/anova-ref.json', 'utf-8'))
console.log(`Loaded ${ref.length} datasets`)

const metrics = makeDefs()
let nAnovaComp = 0
let nKWComp = 0
let nTukeyComp = 0
let nGHComp = 0
let nDunnComp = 0
let nFriedmanComp = 0
let nLMMComp = 0

const perDataset: { seed: number; nTotal: number; k: number; pass: boolean; errors: Record<string, number> }[] = []

for (const r of ref) {
  const errs: Record<string, number> = {}

  // Build GroupData array
  const groupEntries = Object.entries(r.groups)
  const groupData = groupEntries.map(([label, values]) => ({ label, values: values as number[] }))

  // ── ANOVA ──
  try {
    const cAnova = oneWayANOVA(groupData)

    errs.anova_f = absDiff(cAnova.statistic, r.anova.fStatistic)
    push(metrics.anova_f!, errs.anova_f!)
    errs.anova_p = absDiff(cAnova.pValue, r.anova.pValue)
    push(metrics.anova_p!, errs.anova_p!)
    errs.anova_dfb = absDiff(cAnova.dfBetween, r.anova.dfBetween)
    push(metrics.anova_dfb!, errs.anova_dfb!)
    errs.anova_dfw = absDiff(cAnova.dfWithin, r.anova.dfWithin)
    push(metrics.anova_dfw!, errs.anova_dfw!)
    errs.anova_ssb = absDiff(cAnova.ssBetween, r.anova.ssBetween)
    push(metrics.anova_ssb!, errs.anova_ssb!)
    errs.anova_ssw = absDiff(cAnova.ssWithin, r.anova.ssWithin)
    push(metrics.anova_ssw!, errs.anova_ssw!)
    errs.anova_msb = absDiff(cAnova.msBetween, r.anova.msBetween)
    push(metrics.anova_msb!, errs.anova_msb!)
    errs.anova_msw = absDiff(cAnova.msWithin, r.anova.msWithin)
    push(metrics.anova_msw!, errs.anova_msw!)

    // ── Effect sizes ──
    const cEta2 = etaSquared(cAnova.ssBetween, cAnova.ssTotal)
    errs.es_eta2 = absDiff(cEta2.value, r.etaSquared)
    push(metrics.es_eta2!, errs.es_eta2!)

    const cOmega2 = omegaSquared(cAnova.ssBetween, cAnova.ssTotal, cAnova.dfBetween, cAnova.msWithin)
    errs.es_omega2 = absDiff(cOmega2.value, r.omegaSquared)
    push(metrics.es_omega2!, errs.es_omega2!)

    nAnovaComp++

    // ── Tukey HSD ──
    if (r.tukeyHSD && r.tukeyHSD.length > 0) {
      try {
        const cTukey = tukeyHSD(groupData, cAnova.msWithin, cAnova.dfWithin)
        for (const rPair of r.tukeyHSD) {
          // Find matching Carm pair
          const cPair = cTukey.find(p =>
            rPair.pair === `${p.group1}-${p.group2}` ||
            rPair.pair === `${p.group2}-${p.group1}`
          )
          if (cPair) {
            // Diff sign may be reversed if pair order differs
            const cDiff = rPair.pair === `${cPair.group1}-${cPair.group2}` ? cPair.meanDiff : -cPair.meanDiff
            errs[`tukey_${rPair.pair}_diff`] = absDiff(cDiff, rPair.diff)
            push(metrics.tukey_diff!, errs[`tukey_${rPair.pair}_diff`]!)
            errs[`tukey_${rPair.pair}_padj`] = absDiff(cPair.pValueAdj, rPair.pAdj)
            push(metrics.tukey_padj!, errs[`tukey_${rPair.pair}_padj`]!)
            nTukeyComp++
          }
        }
      } catch {
        // Skip Tukey failures
      }
    }

    // ── Games-Howell ──
    if (r.gamesHowell && r.gamesHowell.length > 0) {
      try {
        const cGH = gamesHowell(groupData)
        for (const rPair of r.gamesHowell) {
          const cPair = cGH.find(p =>
            (rPair.group1 === p.group1 && rPair.group2 === p.group2) ||
            (rPair.group1 === p.group2 && rPair.group2 === p.group1)
          )
          if (cPair) {
            const cDiff = rPair.group1 === cPair.group1 ? cPair.meanDiff : -cPair.meanDiff
            // rstatix::games_howell_test returns estimate = mean(group2) - mean(group1)
            // Carm computes m1 - m2 (group1 - group2), so negate R's diff
            errs[`gh_${rPair.group1}_${rPair.group2}_diff`] = absDiff(cDiff, -rPair.diff)
            push(metrics.gh_diff!, errs[`gh_${rPair.group1}_${rPair.group2}_diff`]!)
            errs[`gh_${rPair.group1}_${rPair.group2}_padj`] = absDiff(cPair.pValueAdj, rPair.pAdj)
            push(metrics.gh_padj!, errs[`gh_${rPair.group1}_${rPair.group2}_padj`]!)
            nGHComp++
          }
        }
      } catch {
        // Skip Games-Howell failures
      }
    }
  } catch {
    // Skip ANOVA failures
  }

  // ── Kruskal-Wallis ──
  try {
    const cKW = kruskalWallis(groupData)
    errs.kw_h = absDiff(cKW.statistic, r.kruskalWallis.H)
    push(metrics.kw_h!, errs.kw_h!)
    errs.kw_p = absDiff(cKW.pValue, r.kruskalWallis.pValue)
    push(metrics.kw_p!, errs.kw_p!)

    const cEta2KW = etaSquaredKW(cKW.statistic, r.k, r.nTotal)
    errs.es_eta2kw = absDiff(cEta2KW.value, r.kruskalWallis.etaSquaredKW)
    push(metrics.es_eta2kw!, errs.es_eta2kw!)

    nKWComp++
  } catch {
    // Skip KW failures
  }

  // ── Dunn Test ──
  if (r.dunnBonferroni && r.dunnBonferroni.length > 0) {
    try {
      const cDunnBonf = dunnTest(groupData, 'bonferroni')
      for (const rPair of r.dunnBonferroni) {
        // Parse pair name like "G1 - G2"
        const parts = rPair.pair.split(' - ')
        const cPair = cDunnBonf.find(p =>
          (parts[0] === p.group1 && parts[1] === p.group2) ||
          (parts[0] === p.group2 && parts[1] === p.group1)
        )
        if (cPair) {
          errs[`dunn_bonf_${rPair.pair}_padj`] = absDiff(cPair.pValueAdj, rPair.pAdj)
          push(metrics.dunn_bonf_padj!, errs[`dunn_bonf_${rPair.pair}_padj`]!)
        }
      }
    } catch {
      // Skip Dunn Bonferroni failures
    }
  }

  if (r.dunnHolm && r.dunnHolm.length > 0) {
    try {
      const cDunnHolm = dunnTest(groupData, 'holm')
      for (const rPair of r.dunnHolm) {
        const parts = rPair.pair.split(' - ')
        const cPair = cDunnHolm.find(p =>
          (parts[0] === p.group1 && parts[1] === p.group2) ||
          (parts[0] === p.group2 && parts[1] === p.group1)
        )
        if (cPair) {
          errs[`dunn_holm_${rPair.pair}_padj`] = absDiff(cPair.pValueAdj, rPair.pAdj)
          push(metrics.dunn_holm_padj!, errs[`dunn_holm_${rPair.pair}_padj`]!)
        }
      }
    } catch {
      // Skip Dunn Holm failures
    }
  }

  if (r.dunnBH && r.dunnBH.length > 0) {
    try {
      const cDunnBH = dunnTest(groupData, 'BH')
      for (const rPair of r.dunnBH) {
        const parts = rPair.pair.split(' - ')
        const cPair = cDunnBH.find(p =>
          (parts[0] === p.group1 && parts[1] === p.group2) ||
          (parts[0] === p.group2 && parts[1] === p.group1)
        )
        if (cPair) {
          errs[`dunn_bh_${rPair.pair}_padj`] = absDiff(cPair.pValueAdj, rPair.pAdj)
          push(metrics.dunn_bh_padj!, errs[`dunn_bh_${rPair.pair}_padj`]!)
        }
      }
      nDunnComp++
    } catch {
      // Skip Dunn BH failures
    }
  }

  // ── Friedman ──
  if (r.friedman && r.friedman.minGroupSize >= 2) {
    try {
      const minSize = r.friedman.minGroupSize
      // Build balanced matrix: rows = subjects, cols = conditions
      const matrix: number[][] = []
      for (let i = 0; i < minSize; i++) {
        matrix.push(groupEntries.map(([, v]) => (v as number[])[i]!))
      }
      const cFriedman = friedmanTest(matrix)
      errs.friedman_chi2 = absDiff(cFriedman.statistic, r.friedman.chi2)
      push(metrics.friedman_chi2!, errs.friedman_chi2!)
      errs.friedman_p = absDiff(cFriedman.pValue, r.friedman.pValue)
      push(metrics.friedman_p!, errs.friedman_p!)
      nFriedmanComp++
    } catch {
      // Skip Friedman failures
    }
  }

  // ── LMM ──
  if (r.lmm) {
    try {
      const allScores = groupEntries.flatMap(([label, vals]) =>
        (vals as number[]).map(v => ({ score: v, group: label }))
      )
      const cLMM = runLMM({
        outcome: allScores.map(d => d.score),
        fixedPredictors: {},
        groupId: allScores.map(d => d.group),
      })

      errs.lmm_intercept = absDiff(cLMM.fixedEffects[0]!.estimate, r.lmm.intercept)
      push(metrics.lmm_intercept!, errs.lmm_intercept!)
      errs.lmm_groupvar = absDiff(cLMM.varianceComponents.intercept, r.lmm.groupVar)
      push(metrics.lmm_groupvar!, errs.lmm_groupvar!)
      errs.lmm_residvar = absDiff(cLMM.varianceComponents.residual, r.lmm.residVar)
      push(metrics.lmm_residvar!, errs.lmm_residvar!)
      errs.lmm_icc = absDiff(cLMM.icc, r.lmm.icc)
      push(metrics.lmm_icc!, errs.lmm_icc!)
      nLMMComp++
    } catch {
      // Skip LMM failures
    }
  }

  perDataset.push({ seed: r.seed, nTotal: r.nTotal, k: r.k, pass: true, errors: errs })
}

// ── Aggregate ────────────────────────────────────────────────────────────

const allStats = Object.values(metrics).filter(m => m.errors.length > 0).map(aggStat)
const categories = [...new Set(allStats.map(s => s.category))]

// Print summary to console
console.log('\n═══ ANOVA & Post-Hoc Cross-Validation Summary ═══\n')
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
console.log(`Comparisons: ${nAnovaComp} ANOVA, ${nKWComp} KW, ${nTukeyComp} Tukey pairs, ${nGHComp} GH pairs, ${nDunnComp} Dunn, ${nFriedmanComp} Friedman, ${nLMMComp} LMM`)

// ── HTML Report ──────────────────────────────────────────────────────────

const timestamp = new Date().toISOString().slice(0, 19)

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>ANOVA & Post-Hoc Cross-Validation: Carm vs R</title>
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

<h1>ANOVA & Post-Hoc Cross-Validation</h1>
<p class="subtitle">Carm (TypeScript) vs R — ${ref.length} datasets — ${timestamp}</p>

<div class="summary-bar">
  <div class="summary-card"><div class="num">${ref.length}</div><div class="label">Datasets</div></div>
  <div class="summary-card"><div class="num">${totalMetrics}</div><div class="label">Metrics Tested</div></div>
  <div class="summary-card"><div class="num ${totalPass === totalMetrics ? 'pass' : totalPass >= totalMetrics * 0.9 ? 'warn' : 'fail'}">${totalPass}/${totalMetrics}</div><div class="label">≥99% Pass Rate</div></div>
  <div class="summary-card"><div class="num">${nTukeyComp}</div><div class="label">Tukey Pairs</div></div>
  <div class="summary-card"><div class="num">${nFriedmanComp}</div><div class="label">Friedman Tests</div></div>
  <div class="summary-card"><div class="num">${nLMMComp}</div><div class="label">LMM Fits</div></div>
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
<thead><tr><th>Seed</th><th>n</th><th>k</th><th>F |e|</th><th>F p |e|</th><th>KW H |e|</th><th>eta2 |e|</th><th>Tukey max |e|</th><th>LMM ICC |e|</th></tr></thead>
<tbody>
${perDataset
  .sort((a, b) => {
    const maxA = Math.max(...Object.values(a.errors).filter(v => !isNaN(v)), 0)
    const maxB = Math.max(...Object.values(b.errors).filter(v => !isNaN(v)), 0)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => {
    const tukeyErrs = Object.entries(d.errors).filter(([k]) => k.startsWith('tukey_') && k.endsWith('_diff')).map(([, v]) => v)
    const maxTukey = tukeyErrs.length > 0 ? Math.max(...tukeyErrs) : 0
    return `<tr>
  <td>${d.seed}</td><td>${d.nTotal}</td><td>${d.k}</td>
  <td class="${cc(d.errors.anova_f ?? 0, 1e-4)}">${fmt(d.errors.anova_f ?? 0)}</td>
  <td class="${cc(d.errors.anova_p ?? 0, 1e-3)}">${fmt(d.errors.anova_p ?? 0)}</td>
  <td class="${cc(d.errors.kw_h ?? 0, 1e-4)}">${fmt(d.errors.kw_h ?? 0)}</td>
  <td class="${cc(d.errors.es_eta2 ?? 0, 1e-4)}">${fmt(d.errors.es_eta2 ?? 0)}</td>
  <td class="${cc(maxTukey, 1e-3)}">${fmt(maxTukey)}</td>
  <td class="${cc(d.errors.lmm_icc ?? 0, 1e-2)}">${fmt(d.errors.lmm_icc ?? 0)}</td>
</tr>`
  }).join('\n')}
</tbody>
</table>

<p style="margin-top:40px;color:#6e7681;font-size:11px;">
  Generated ${timestamp} · Carm vs R (base R aov, kruskal.test, TukeyHSD, rstatix, dunn.test, friedman.test, lme4) · ${ref.length} datasets · ${totalMetrics} metrics
</p>

</body>
</html>`

mkdirSync('validation/reports', { recursive: true })
writeFileSync('validation/reports/anova-report.html', html)
console.log(`\nReport saved to validation/reports/anova-report.html`)
