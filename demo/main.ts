/**
 * Carm demo application.
 * Demonstrates all major analysis and visualization modules
 * with synthetic data, publication-quality plots, and styled tables.
 */

import {
  describe as describeStats,
  oneWayANOVA,
  tukeyHSD,
  pearsonCorrelation,
  correlationMatrix,
  linearRegression,
  multipleRegression,
  regressionDiagnostics,
  runLMM,
  computeBLUPs,
  runPCA,
  renderHistogram,
  renderQQPlot,
  renderViolinBox,
  renderRaincloud,
  renderScatterStats,
  renderCorrelogram,
  renderBarStats,
  renderCoefPlot,
  renderResidualPanel,
  renderMixedPlot,
  renderPCAPlot,
  renderDistribution,
  // Batch 1
  renderDensity,
  renderBoxplot,
  renderLollipop,
  renderDotPlot,
  renderGroupedBar,
  renderLineChart,
  renderBubbleChart,
  // Batch 2
  renderPareto,
  renderFunnel,
  renderPieChart,
  renderAreaChart,
  renderForestPlot,
  renderROCCurve,
  renderStripPlot,
  // Batch 3
  renderSwarmPlot,
  renderMosaicPlot,
  renderPairPlot,
  renderRadarChart,
  renderParallelCoords,
  renderTreemap,
  renderWaffleChart,
  renderSparkline,
  // Hierarchical & Proportional
  renderSunburst,
  renderMarimekko,
  // Networks & Flows
  renderChordDiagram,
  renderAlluvialPlot,
  renderArcDiagram,
  renderEdgeBundling,
} from 'carm'

// ─── Seeded RNG (LCG + Box-Muller) ────────────────────────────────────────

function randn(seed: number): () => number {
  let s = seed
  const lcg = () => { s = (s * 1664525 + 1013904223) & 0xffffffff; return (s >>> 0) / 0xffffffff }
  let spare: number | null = null
  return () => {
    if (spare !== null) { const v = spare; spare = null; return v }
    const u = lcg(), v = lcg()
    const r = Math.sqrt(-2 * Math.log(u + 1e-10))
    spare = r * Math.sin(2 * Math.PI * v)
    return r * Math.cos(2 * Math.PI * v)
  }
}

function randu(seed: number): () => number {
  let s = seed
  return () => { s = (s * 1664525 + 1013904223) & 0xffffffff; return (s >>> 0) / 0xffffffff }
}

// ─── Statistical helpers ───────────────────────────────────────────────────

/** Pooled-SD Cohen's d (unsigned). */
function cohensD(m1: number, s1: number, n1: number, m2: number, s2: number, n2: number): number {
  const sp = Math.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
  return sp > 0 ? Math.abs(m1 - m2) / sp : 0
}

/** Cohen's d verbal label (Cohen 1988). */
function interpretD(d: number): string {
  return d < 0.2 ? 'negligible' : d < 0.5 ? 'small' : d < 0.8 ? 'medium' : 'large'
}

/** Sample standard deviation. */
function sdOf(arr: readonly number[]): number {
  const m = arr.reduce((s, v) => s + v, 0) / arr.length
  return Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / (arr.length - 1))
}

/** Percentile via linear interpolation (p ∈ [0,1]). */
function pctile(arr: readonly number[], p: number): number {
  const sorted = [...arr].sort((a, b) => a - b)
  const pos = p * (sorted.length - 1)
  const lo = Math.floor(pos)
  const frac = pos - lo
  return lo + 1 < sorted.length ? sorted[lo]! + frac * (sorted[lo + 1]! - sorted[lo]!) : sorted[lo]!
}

// ─── Table utilities ───────────────────────────────────────────────────────

type Align = 'l' | 'r' | 'c'

function sig(p: number): string {
  const s = p < 0.001 ? '***' : p < 0.01 ? '**' : p < 0.05 ? '*' : p < 0.10 ? '†' : ''
  return s ? `<sup class="stat-sig">${s}</sup>` : ''
}

function fmtP(p: number): string {
  if (p < 0.001) return '&lt;&thinsp;.001'
  return p.toFixed(3).replace(/^0\./, '.')
}

function fmtCI(lo: number, hi: number, d = 2): string {
  return `[${lo.toFixed(d)},&thinsp;${hi.toFixed(d)}]`
}

function fmtR(r: number): string {
  const s = r.toFixed(3)
  return r < 0 ? s.replace('-0.', '&minus;.') : s.replace('0.', '.')
}

function fmtD(d: number): string {
  return d.toFixed(2).replace(/^0\./, '.')
}

function tbl(
  container: HTMLElement,
  opts: {
    title?: string
    subtitle?: string
    note?: string
    cols: Array<{ h: string; a?: Align }>
    rows: string[][]
  }
): void {
  const a = (i: number): Align => opts.cols[i]?.a ?? 'r'
  const hdr = opts.cols.map((c, i) => `<th class="${a(i)}">${c.h}</th>`).join('')
  const bdy = opts.rows.map(row =>
    `<tr>${row.map((cell, i) => `<td class="${a(i)}">${cell}</td>`).join('')}</tr>`
  ).join('')
  container.insertAdjacentHTML('beforeend', `
    <div class="stat-table-wrap">
      ${opts.title    ? `<div class="stat-table-title">${opts.title}</div>` : ''}
      ${opts.subtitle ? `<div class="stat-table-subtitle">${opts.subtitle}</div>` : ''}
      <table class="stat-table">
        <thead><tr>${hdr}</tr></thead>
        <tbody>${bdy}</tbody>
      </table>
      ${opts.note ? `<div class="stat-table-note"><em>Note.</em> ${opts.note}</div>` : ''}
    </div>`)
}


// ─── Synthetic datasets ────────────────────────────────────────────────────

const rng  = randn(42)
const rng2 = randn(99)
const rng3 = randn(7)
const rng4 = randn(13)
const rng5 = randn(55)

// Group comparison: 3 groups n=30 each
const groupA = Array.from({ length: 30 }, () => 5.0 + rng() * 1.5)
const groupB = Array.from({ length: 30 }, () => 7.0 + rng() * 1.8)
const groupC = Array.from({ length: 30 }, () => 6.0 + rng() * 1.2)

// Correlation / scatter: positively correlated pair
const x1 = Array.from({ length: 50 }, (_, i) => i * 0.1 + rng() * 0.5)
const y1 = x1.map(x => 2.5 + 1.8 * x + rng() * 0.8)

// Multiple regression: y = 1 + 2·x1 + 0.5·x2 + noise
const pred1 = Array.from({ length: 40 }, () => rng2() * 3)
const pred2 = Array.from({ length: 40 }, () => rng2() * 2 + 1)
const yReg  = pred1.map((x, i) => 1 + 2 * x + 0.5 * (pred2[i]!) + rng2() * 0.8)

// LMM: 3 groups × 10 obs
const lmmGroups = Array.from({ length: 30 }, (_, i) => Math.floor(i / 10) + 1)
const lmmX = lmmGroups.map(() => rng3() * 2 + 1)
const groupEffects: Record<number, number> = { 1: -2, 2: 0, 3: 3 }
const lmmY = lmmX.map((x, i) => (groupEffects[lmmGroups[i]!] ?? 0) + 1.2 * x + rng3() * 0.5)

// PCA: 5 variables, 60 obs
const pcaData = Array.from({ length: 60 }, () => {
  const f1 = rng4(), f2 = rng4()
  return [
    f1 + rng4() * 0.3,
    f1 * 0.8  + rng4() * 0.4,
    f1 * -0.6 + rng4() * 0.5,
    f2 + rng4() * 0.3,
    f2 * 0.7  + rng4() * 0.4,
  ]
})

// Time-series data: 24 months, 3 series
const months = Array.from({ length: 24 }, (_, i) => i)
const trendA = months.map(m => 10 + m * 0.8 + rng5() * 2.0)
const trendB = months.map(m => 15 + m * 0.3 + rng5() * 2.5)
const trendC = months.map(m => 8  + m * 1.2 + rng5() * 1.5)

// Ranking data
const rankLabels = ['Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Eta', 'Theta']
const rankValues = [42, 78, 31, 91, 55, 67, 23, 84]

// Composition data
const pieSlices = [
  { label: 'Category A', value: 35 },
  { label: 'Category B', value: 28 },
  { label: 'Category C', value: 18 },
  { label: 'Category D', value: 12 },
  { label: 'Category E', value: 7 },
]

// Bubble data
const bubblePoints = Array.from({ length: 20 }, (_, i) => ({
  x: rng5() * 10 + 5,
  y: rng5() * 8 + 4,
  r: Math.abs(rng5()) * 8 + 2,
  label: `P${i + 1}`,
  group: i < 7 ? 'Group A' : i < 14 ? 'Group B' : 'Group C',
}))

// Forest plot
const forestStudies = [
  { label: 'Smith et al. (2018)', estimate: 0.45, ciLow: 0.12, ciHigh: 0.78, weight: 3.2 },
  { label: 'Jones et al. (2019)', estimate: 0.62, ciLow: 0.31, ciHigh: 0.93, weight: 4.1 },
  { label: 'Brown et al. (2020)', estimate: 0.28, ciLow: -0.05, ciHigh: 0.61, weight: 2.8 },
  { label: 'Davis et al. (2021)', estimate: 0.71, ciLow: 0.45, ciHigh: 0.97, weight: 5.0 },
  { label: 'Wilson et al. (2022)', estimate: 0.53, ciLow: 0.22, ciHigh: 0.84, weight: 3.7 },
]
const pooledEstimate = { estimate: 0.54, ciLow: 0.38, ciHigh: 0.70 }

// ROC curve
const rocFpr: number[] = [0]
const rocTpr: number[] = [0]
const rocRng = randu(123)
for (let t = 1; t >= 0; t -= 0.02) {
  rocFpr.push(Math.min(1, Math.max(0, Math.pow(t, 2.5) * (0.9 + rocRng() * 0.2))))
  rocTpr.push(Math.min(1, Math.max(0, 1 - Math.pow(1 - t, 1.8) * (0.85 + rocRng() * 0.1))))
}
rocFpr.push(1); rocTpr.push(1)
const rocAuc = 0.82

// Radar
const radarAxes = ['Speed', 'Strength', 'Agility', 'Endurance', 'Technique']
const radarSeries = [
  { label: 'Athlete A', values: [85, 72, 90, 65, 88] },
  { label: 'Athlete B', values: [70, 91, 75, 82, 77] },
  { label: 'Athlete C', values: [92, 68, 83, 78, 71] },
]

// Parallel coordinates
const pcAxes = ['Sepal L', 'Sepal W', 'Petal L', 'Petal W']
const pcRng = randn(31)
const pcRows = Array.from({ length: 45 }, (_, i) => {
  const g = Math.floor(i / 15)
  const bases = [[5.0, 3.4, 1.5, 0.2], [5.9, 2.8, 4.3, 1.3], [6.5, 3.0, 5.5, 2.0]]
  return bases[g]!.map(b => Math.max(0.1, b + pcRng() * 0.6))
})
const pcGroups = Array.from({ length: 45 }, (_, i) => Math.floor(i / 15))

// Mosaic
const mosaicTable = [[30, 20, 10], [25, 40, 20], [10, 25, 35]]
const mosaicRowLabels = ['Low Ed.', 'Mid Ed.', 'High Ed.']
const mosaicColLabels = ['Low Inc.', 'Mid Inc.', 'High Inc.']

// Treemap
const treemapChildren = [
  { label: 'Electronics', value: 420, group: 'Tech' },
  { label: 'Software', value: 310, group: 'Tech' },
  { label: 'Hardware', value: 180, group: 'Tech' },
  { label: 'Clothing', value: 250, group: 'Retail' },
  { label: 'Footwear', value: 190, group: 'Retail' },
  { label: 'Accessories', value: 120, group: 'Retail' },
  { label: 'Food', value: 200, group: 'FMCG' },
  { label: 'Beverages', value: 150, group: 'FMCG' },
]

const waffleSlices = [
  { label: 'Product A', value: 35 },
  { label: 'Product B', value: 28 },
  { label: 'Product C', value: 20 },
  { label: 'Product D', value: 17 },
]

const funnelStages = [
  { label: 'Awareness',      value: 10000 },
  { label: 'Interest',       value: 6500  },
  { label: 'Consideration',  value: 3800  },
  { label: 'Intent',         value: 2100  },
  { label: 'Purchase',       value: 980   },
]

const pairData = pcaData.slice(0, 30).map(row => row.slice(0, 4))
const pairLabels = ['V1', 'V2', 'V3', 'V4']

const swarmRng = randn(77)
const swarmGroups = [
  { label: 'Control',   values: Array.from({ length: 25 }, () => 5.0 + swarmRng() * 1.2) },
  { label: 'Treatment', values: Array.from({ length: 25 }, () => 6.5 + swarmRng() * 1.4) },
  { label: 'Placebo',   values: Array.from({ length: 25 }, () => 5.3 + swarmRng() * 1.1) },
]

// ═══════════════════════════════════════════════════════════════════════════
// DESCRIPTIVE SECTION
// ═══════════════════════════════════════════════════════════════════════════

const allValues    = [...groupA, ...groupB, ...groupC]
const descriptives = describeStats(allValues)

renderHistogram(
  document.getElementById('histogram-plot')!,
  { values: allValues, descriptives },
  {
    title: 'Distribution of Group Values',
    xLabel: 'Value', showNormalCurve: true, width: 650, height: 380,
    caption: `N = ${descriptives.n}. Curve: fitted normal distribution. Shapiro–Wilk W = ${descriptives.shapiroWilk.statistic.toFixed(3)}, p = ${descriptives.shapiroWilk.pValue.toFixed(3)}.`,
  }
)
renderQQPlot(
  document.getElementById('qq-plot')!,
  allValues,
  {
    title: 'Normal QQ Plot',
    width: 450, height: 380,
    caption: 'Points along the diagonal indicate normality. Shaded band: 95% pointwise confidence envelope.',
  }
)
renderDensity(
  document.getElementById('density-plot')!,
  {
    series: [
      { label: 'Group A', values: groupA },
      { label: 'Group B', values: groupB },
      { label: 'Group C', values: groupC },
    ],
  },
  {
    title: 'KDE Density by Group',
    xLabel: 'Value', yLabel: 'Density', showRug: true, width: 650, height: 340,
    caption: 'Kernel density estimation (Silverman bandwidth). Rug marks: individual observations.',
  }
)
renderBoxplot(
  document.getElementById('boxplot-plot')!,
  { groups: [
    { label: 'Group A', values: groupA },
    { label: 'Group B', values: groupB },
    { label: 'Group C', values: groupC },
  ] },
  {
    title: 'Standalone Box Plot',
    xLabel: 'Group', yLabel: 'Value', width: 500, height: 380,
    caption: 'Box: IQR (Q1–Q3). Whiskers: 1.5 × IQR. Circles: outliers.',
  }
)

const dTables = document.getElementById('descriptive-tables')!
const d = descriptives
const normConclusion = d.shapiroWilk.pValue > 0.05
  ? '<span style="color:#198754;font-weight:600">Normal</span>'
  : '<span style="color:#dc3545;font-weight:600">Non-normal</span>'

tbl(dTables, {
  title: 'Table 1.  Descriptive Statistics',
  subtitle: `N = ${d.n} observations`,
  cols: [
    { h: 'Statistic',  a: 'l' },
    { h: 'Value',      a: 'r' },
    { h: 'Statistic',  a: 'l' },
    { h: 'Value',      a: 'r' },
  ],
  rows: [
    ['Mean (<em>M</em>)',                        d.mean.toFixed(3),
     'Trimmed Mean (5%)',                         d.trimmedMean.toFixed(3)],
    ['Standard Deviation (<em>SD</em>)',         d.sd.toFixed(3),
     'Variance (<em>s</em>&sup2;)',               d.variance.toFixed(3)],
    ['Standard Error (<em>SE</em>)',             d.se.toFixed(3),
     '95% CI for <em>M</em>',                    fmtCI(d.ci[0], d.ci[1])],
    ['Median',                                    d.median.toFixed(3),
     'IQR',                                       d.iqr.toFixed(3)],
    ['Q1 (25th pctile)',                         d.q1.toFixed(3),
     'Q3 (75th pctile)',                         d.q3.toFixed(3)],
    ['P5',                                        pctile(allValues, 0.05).toFixed(3),
     'P95',                                       pctile(allValues, 0.95).toFixed(3)],
    ['Min',                                       d.min.toFixed(3),
     'Max',                                       d.max.toFixed(3)],
    ['Skewness (<em>g</em><sub>1</sub>)',        d.skewness.toFixed(3),
     'Excess Kurtosis (<em>g</em><sub>2</sub>)', d.kurtosis.toFixed(3)],
  ],
  note: 'Trimmed mean excludes top and bottom 5%. Skewness and kurtosis are moment-based estimators (G1/G2).',
})

tbl(dTables, {
  title: 'Table 2.  Normality Assessment',
  subtitle: 'Shapiro–Wilk test of normality',
  cols: [
    { h: 'Statistic', a: 'l' },
    { h: 'Value',     a: 'r' },
  ],
  rows: [
    ['<em>W</em>',                d.shapiroWilk.statistic.toFixed(4)],
    ['<em>p</em>-value',          fmtP(d.shapiroWilk.pValue) + sig(d.shapiroWilk.pValue)],
    ['<em>n</em>',                `${d.n}`],
    ['Conclusion',                normConclusion],
  ],
  note: 'H₀: data are drawn from a normal distribution. '
      + 'W close to 1 indicates normality. '
      + '*&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
})

// ═══════════════════════════════════════════════════════════════════════════
// COMPARISON SECTION
// ═══════════════════════════════════════════════════════════════════════════

const groups = [
  { label: 'Group A', values: groupA },
  { label: 'Group B', values: groupB },
  { label: 'Group C', values: groupC },
]
const anovaResult  = oneWayANOVA(groups)
const tukeyResults = tukeyHSD(groups, anovaResult.msWithin, anovaResult.dfWithin)
const groupDescs   = groups.map(g => describeStats([...g.values]))
const groupDescMap = new Map(groups.map((g, i) => [g.label, groupDescs[i]!]))

const nTotal = groups.reduce((s, g) => s + g.values.length, 0)
const anovaSubtitle = `${anovaResult.formatted}  •  N = ${nTotal}`

renderViolinBox(
  document.getElementById('violin-plot')!,
  { groups, testResult: anovaResult, pairwise: tukeyResults },
  {
    title: 'Group Comparison',
    xLabel: 'Group', yLabel: 'Value', showBrackets: true, width: 650, height: 500,
    caption: `N = ${nTotal} (${groups.map(g => g.values.length).join(' / ')} per group). One-way ANOVA with Tukey HSD post-hoc. Brackets: p < .05.`,
  }
)
renderRaincloud(
  document.getElementById('raincloud-plot')!,
  { groups, testResult: anovaResult },
  {
    title: 'Raincloud Plot',
    xLabel: 'Group', yLabel: 'Value', width: 650, height: 480,
    caption: 'Half-violin (KDE density) + box (IQR, whiskers) + jittered individual points.',
  }
)
renderGroupedBar(
  document.getElementById('grouped-bar-plot')!,
  {
    categories: ['Q1', 'Q2', 'Q3', 'Q4'],
    series: [
      { label: 'Product A', values: [42, 58, 71, 65] },
      { label: 'Product B', values: [35, 49, 62, 78] },
      { label: 'Product C', values: [28, 41, 55, 69] },
    ],
  },
  { type: 'grouped', title: 'Quarterly Sales by Product (Grouped)', xLabel: 'Quarter', yLabel: 'Units', width: 650, height: 400 }
)
renderGroupedBar(
  document.getElementById('stacked-bar-plot')!,
  {
    categories: ['Q1', 'Q2', 'Q3', 'Q4'],
    series: [
      { label: 'Product A', values: [42, 58, 71, 65] },
      { label: 'Product B', values: [35, 49, 62, 78] },
      { label: 'Product C', values: [28, 41, 55, 69] },
    ],
  },
  { type: 'stacked', title: 'Quarterly Sales by Product (Stacked)', xLabel: 'Quarter', yLabel: 'Cumulative Units', width: 650, height: 400 }
)
renderBubbleChart(
  document.getElementById('bubble-plot')!,
  { points: bubblePoints },
  { title: 'Bubble Chart', xLabel: 'X Variable', yLabel: 'Y Variable', width: 650, height: 480 }
)

const cTables = document.getElementById('comparison-tables')!

tbl(cTables, {
  title: 'Table 3.  Descriptive Statistics by Group',
  subtitle: `N = ${nTotal} total`,
  cols: [
    { h: 'Group',                            a: 'l' },
    { h: '<em>n</em>',                       a: 'c' },
    { h: '<em>M</em>',                       a: 'r' },
    { h: '<em>SD</em>',                      a: 'r' },
    { h: '<em>SE</em>',                      a: 'r' },
    { h: '95% CI',                           a: 'c' },
    { h: 'Mdn',                              a: 'r' },
    { h: 'Min &ndash; Max',                  a: 'c' },
    { h: 'Skew',                             a: 'r' },
  ],
  rows: groups.map((g, i) => {
    const gd = groupDescs[i]!
    return [
      `<strong>${g.label}</strong>`,
      `${gd.n}`,
      gd.mean.toFixed(3),
      gd.sd.toFixed(3),
      gd.se.toFixed(3),
      fmtCI(gd.ci[0], gd.ci[1]),
      gd.median.toFixed(3),
      `${gd.min.toFixed(2)} &ndash; ${gd.max.toFixed(2)}`,
      gd.skewness.toFixed(2),
    ]
  }),
})

tbl(cTables, {
  title: 'Table 4.  One-Way ANOVA Summary',
  cols: [
    { h: 'Source',                        a: 'l' },
    { h: '<em>SS</em>',                   a: 'r' },
    { h: '<em>df</em>',                   a: 'c' },
    { h: '<em>MS</em>',                   a: 'r' },
    { h: '<em>F</em>',                    a: 'r' },
    { h: '<em>p</em>',                    a: 'r' },
    { h: '&eta;&sup2;',                   a: 'r' },
    { h: '&omega;&sup2;',                 a: 'r' },
  ],
  rows: [
    [
      'Between groups',
      anovaResult.ssBetween.toFixed(3),
      `${anovaResult.dfBetween}`,
      anovaResult.msBetween.toFixed(3),
      anovaResult.statistic.toFixed(3) + sig(anovaResult.pValue),
      fmtP(anovaResult.pValue),
      // η² = SS_between / SS_total
      (anovaResult.ssBetween / anovaResult.ssTotal).toFixed(3),
      anovaResult.effectSize.value.toFixed(3),
    ],
    ['Within groups', anovaResult.ssWithin.toFixed(3), `${anovaResult.dfWithin}`, anovaResult.msWithin.toFixed(3), '', '', '', ''],
    ['<strong>Total</strong>', anovaResult.ssTotal.toFixed(3), `${anovaResult.dfBetween + anovaResult.dfWithin}`, '', '', '', '', ''],
  ],
  note: `η² = SS_between / SS_total (unadjusted). `
      + `ω² = (SS_B − df_B·MS_W) / (SS_T + MS_W) (bias-corrected). `
      + `Effect size: ω² = ${anovaResult.effectSize.value.toFixed(3)} (${anovaResult.effectSize.interpretation}). `
      + '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
})

tbl(cTables, {
  title: 'Table 5.  Post-hoc Pairwise Comparisons (Tukey HSD)',
  subtitle: 'Family-wise error rate controlled at α = .05',
  cols: [
    { h: 'Comparison',                          a: 'l' },
    { h: '<em>MD</em>',                         a: 'r' },
    { h: '<em>SE</em>',                         a: 'r' },
    { h: '<em>q</em>',                          a: 'r' },
    { h: '<em>p</em><sub>adj</sub>',            a: 'r' },
    { h: '95% CI',                              a: 'c' },
    { h: "<em>d</em>",                          a: 'r' },
    { h: 'Size',                                a: 'l' },
    { h: '',                                    a: 'c' },
  ],
  rows: tukeyResults.map(pw => {
    const d1 = groupDescMap.get(pw.group1)!
    const d2 = groupDescMap.get(pw.group2)!
    const d_val = cohensD(d1.mean, d1.sd, d1.n, d2.mean, d2.sd, d2.n)
    return [
      `${pw.group1} vs ${pw.group2}`,
      pw.meanDiff.toFixed(3),
      pw.se.toFixed(3),
      pw.statistic.toFixed(3),
      fmtP(pw.pValueAdj),
      fmtCI(pw.ci[0], pw.ci[1]),
      fmtD(d_val),
      `<span style="color:#6c757d;font-size:10px">${interpretD(d_val)}</span>`,
      sig(pw.pValueAdj) || '<span class="stat-dim">ns</span>',
    ]
  }),
  note: 'MD = mean difference. d = Cohen\'s d (pooled SD). '
      + 'Size benchmarks (Cohen 1988): d &lt; .20 negligible, .20 small, .50 medium, .80 large. '
      + 'ns = not significant.',
})

// ═══════════════════════════════════════════════════════════════════════════
// CORRELATION SECTION
// ═══════════════════════════════════════════════════════════════════════════

const corrResult    = pearsonCorrelation(x1, y1)
const regForScatter = linearRegression(x1, y1)
const r2_scatter    = corrResult.statistic ** 2

renderScatterStats(
  document.getElementById('scatter-plot')!,
  { x: x1, y: y1, correlationResult: corrResult, regressionResult: regForScatter },
  {
    title: 'Scatter with Regression Line',
    xLabel: 'X', yLabel: 'Y', width: 650, height: 480,
    caption: `n = ${corrResult.n}. Pearson r = ${corrResult.statistic.toFixed(3)}, p ${corrResult.pValue < .001 ? '< .001' : '= ' + corrResult.pValue.toFixed(3)}. R² = ${r2_scatter.toFixed(3)} (${(r2_scatter * 100).toFixed(1)}% variance explained). Shaded band: 95% CI.`,
  }
)

const cm4 = correlationMatrix(
  Array.from({ length: 4 }, (_, vi) => pcaData.map(row => row[vi] ?? 0)),
  ['V1', 'V2', 'V3', 'V4']
)
renderCorrelogram(
  document.getElementById('correlogram-plot')!,
  cm4,
  {
    title: 'Correlation Matrix (V1–V4)',
    width: 500, height: 500,
    caption: 'Color scale: red = positive correlation, blue = negative. * p < .05  ** p < .01  *** p < .001',
  }
)

const corrTables = document.getElementById('correlation-tables')!
const cr   = corrResult
const crDf = typeof cr.df === 'number' ? cr.df : cr.df[0]
const r2   = cr.statistic ** 2

tbl(corrTables, {
  title: 'Table 6.  Pearson Correlation',
  subtitle: 'X and Y (n = 50)',
  cols: [
    { h: '<em>r</em>',         a: 'r' },
    { h: '<em>r</em>&sup2;',  a: 'r' },
    { h: '<em>df</em>',        a: 'c' },
    { h: '<em>t</em>',         a: 'r' },
    { h: '<em>p</em>',         a: 'r' },
    { h: '95% CI',             a: 'c' },
    { h: '% Var. Explained',   a: 'r' },
    { h: 'Interpretation',     a: 'l' },
  ],
  rows: [[
    fmtR(cr.statistic) + sig(cr.pValue),
    fmtR(r2),
    `${crDf}`,
    (cr.statistic * Math.sqrt(cr.n - 2) / Math.sqrt(1 - cr.statistic ** 2)).toFixed(3),
    fmtP(cr.pValue),
    fmtCI(cr.ci[0], cr.ci[1], 3),
    `${(r2 * 100).toFixed(1)}%`,
    cr.effectSize.interpretation.charAt(0).toUpperCase() + cr.effectSize.interpretation.slice(1),
  ]],
  note: 'r² = proportion of variance in Y explained by X. '
      + 'Benchmarks (Cohen 1988): |r| &ge; .10 small, .30 medium, .50 large. '
      + '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
})

const cmLabels = cm4.labels
tbl(corrTables, {
  title: 'Table 7.  Correlation Matrix',
  subtitle: 'Variables V1–V4 (n = 60). Lower triangle: Pearson r. Upper triangle: r².',
  cols: [
    { h: 'Variable', a: 'l' },
    ...cmLabels.map(l => ({ h: `<strong>${l}</strong>`, a: 'c' as Align })),
  ],
  rows: cmLabels.map((label, i) => [
    `<strong>${label}</strong>`,
    ...cmLabels.map((_, j): string => {
      if (i === j) return '<span class="stat-dim">—</span>'
      const r = cm4.r[i]?.[j] ?? 0
      const p = cm4.pValues[i]?.[j] ?? 1
      if (j < i) return fmtR(r) + sig(p)          // lower: r
      return `<span class="stat-dim">${fmtR(r * r)}</span>`   // upper: r²
    }),
  ]),
  note: 'Lower triangle: Pearson r. Upper triangle (dim): r² (coefficient of determination). '
      + '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
})

// ═══════════════════════════════════════════════════════════════════════════
// REGRESSION SECTION
// ═══════════════════════════════════════════════════════════════════════════

const multiReg = multipleRegression(yReg, [
  { name: 'Predictor 1', values: pred1 },
  { name: 'Predictor 2', values: pred2 },
])

renderCoefPlot(
  document.getElementById('coef-plot')!,
  multiReg.coefficients,
  {
    title: 'Regression Coefficients (95% CI)',
    xLabel: 'Unstandardised Estimate', width: 550, height: 280,
    caption: `R² = ${multiReg.r2.toFixed(3)}, adj. R² = ${multiReg.adjR2.toFixed(3)}, F(${multiReg.fDf[0]},${multiReg.fDf[1]}) = ${multiReg.fStatistic.toFixed(2)}, p ${multiReg.fPValue < .001 ? '< .001' : '= ' + multiReg.fPValue.toFixed(3)}.`,
  }
)

const diag = regressionDiagnostics(multiReg, [
  { name: 'Predictor 1', values: pred1 },
  { name: 'Predictor 2', values: pred2 },
])
renderResidualPanel(
  document.getElementById('residual-plot')!,
  multiReg,
  diag.leverage,
  {
    title: 'Regression Diagnostics',
    width: 700, height: 580,
    caption: 'Top-left: residuals vs. fitted. Top-right: normal QQ of residuals. Bottom-left: scale-location. Bottom-right: leverage (Cook\'s distance).',
  }
)

// Compute standardised betas: β_std = β * (SD_x / SD_y)
const sdYReg  = sdOf(yReg)
const sdPreds: Record<string, number> = {
  'Predictor 1': sdOf(pred1),
  'Predictor 2': sdOf(pred2),
}

const regTables = document.getElementById('regression-tables')!
tbl(regTables, {
  title: 'Table 8.  Model Summary',
  subtitle: `N = ${multiReg.n}, predictors = 2`,
  cols: [{ h: 'Index', a: 'l' }, { h: 'Value', a: 'r' }, { h: 'Index', a: 'l' }, { h: 'Value', a: 'r' }],
  rows: [
    ['<em>R</em>',                     Math.sqrt(multiReg.r2).toFixed(4),  '<em>F</em>',  multiReg.fStatistic.toFixed(3) + sig(multiReg.fPValue)],
    ['<em>R</em>&sup2;',               multiReg.r2.toFixed(4),              'df<sub>1</sub>, df<sub>2</sub>', `${multiReg.fDf[0]},&thinsp;${multiReg.fDf[1]}`],
    ['Adj. <em>R</em>&sup2;',          multiReg.adjR2.toFixed(4),           '<em>p</em>',  fmtP(multiReg.fPValue)],
    ['<em>SE</em><sub>est</sub>',      Math.sqrt(multiReg.residualVariance ?? (multiReg.ssResidual / (multiReg.n - 3))).toFixed(4), 'AIC', multiReg.aic.toFixed(2)],
    ['<em>n</em>',                     `${multiReg.n}`,                    'BIC',  multiReg.bic.toFixed(2)],
  ],
  note: `Dependent variable: Y. Predictors: Predictor 1, Predictor 2. Method: OLS. `
      + `R² = ${multiReg.r2.toFixed(4)} indicates ${(multiReg.r2 * 100).toFixed(1)}% of variance explained.`,
})

tbl(regTables, {
  title: 'Table 9.  Regression Coefficients',
  subtitle: 'Unstandardised and standardised estimates with 95% confidence intervals',
  cols: [
    { h: 'Predictor',          a: 'l' },
    { h: '<em>B</em>',         a: 'r' },
    { h: '<em>&beta;</em>*',   a: 'r' },
    { h: '<em>SE</em>',        a: 'r' },
    { h: '<em>t</em>',         a: 'r' },
    { h: '<em>p</em>',         a: 'r' },
    { h: '95% CI',             a: 'c' },
    { h: '',                   a: 'c' },
  ],
  rows: multiReg.coefficients.map(c => {
    const sdX = sdPreds[c.name] ?? null
    const stdBeta = sdX !== null && sdYReg > 0 ? (c.estimate * sdX / sdYReg).toFixed(3) : '—'
    return [
      c.name === '(Intercept)' ? '<em>Intercept</em>' : `<strong>${c.name}</strong>`,
      c.estimate.toFixed(3),
      stdBeta,
      c.se.toFixed(3),
      c.tValue.toFixed(3),
      fmtP(c.pValue),
      fmtCI(c.ci[0], c.ci[1]),
      sig(c.pValue) || '<span class="stat-dim">ns</span>',
    ]
  }),
  note: '<em>B</em> = unstandardised coefficient. <em>&beta;*</em> = standardised coefficient (β × SD_x / SD_y). '
      + 'OLS regression. ns = not significant. '
      + '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
})

// ═══════════════════════════════════════════════════════════════════════════
// LMM SECTION
// ═══════════════════════════════════════════════════════════════════════════

const lmmResult = runLMM({ outcome: lmmY, fixedPredictors: { x: lmmX }, groupId: lmmGroups })
const blups = computeBLUPs({ outcome: lmmY, fixedPredictors: { x: lmmX }, groupId: lmmGroups }, lmmResult)

renderMixedPlot(
  document.getElementById('lmm-plot')!,
  lmmResult, blups,
  {
    title: 'Random Effects — BLUPs Caterpillar Plot',
    width: 550, height: 280,
    caption: `k = ${lmmResult.nGroups} groups, N = ${lmmResult.nObs}. BLUPs = Best Linear Unbiased Predictors. Bars: 95% CI.`,
  }
)

const lmmTables = document.getElementById('lmm-tables')!
tbl(lmmTables, {
  title: 'Table 10.  Fixed Effects',
  subtitle: `Linear Mixed Model (REML). N = ${lmmResult.nObs} observations, k = ${lmmResult.nGroups} groups`,
  cols: [
    { h: 'Effect',          a: 'l' },
    { h: '<em>&beta;</em>', a: 'r' },
    { h: '<em>SE</em>',    a: 'r' },
    { h: '<em>t</em>',     a: 'r' },
    { h: '<em>p</em>',     a: 'r' },
    { h: '95% CI',         a: 'c' },
    { h: '',               a: 'c' },
  ],
  rows: lmmResult.fixedEffects.map(fe => [
    fe.name === '(Intercept)' ? '<em>Intercept</em>' : `<strong>${fe.name}</strong>`,
    fe.estimate.toFixed(3),
    fe.se.toFixed(3),
    fe.tValue.toFixed(3),
    fmtP(fe.pValue),
    fmtCI(fe.ci[0], fe.ci[1]),
    sig(fe.pValue) || '<span class="stat-dim">ns</span>',
  ]),
  note: 'Estimated via REML. Degrees of freedom: approximate (Satterthwaite). ns = not significant.',
})

const vc = lmmResult.varianceComponents
const totalVar = vc.intercept + vc.residual
renderBarStats(
  document.getElementById('lmm-variance-plot')!,
  {
    rows: [
      { value: 'Random intercept (σ²_b)', count: vc.intercept, relative: vc.intercept / totalVar, cumulative: vc.intercept / totalVar },
      { value: 'Residual (σ²_e)',          count: vc.residual,  relative: vc.residual  / totalVar, cumulative: 1 },
    ],
  },
  { title: 'Variance Components', yLabel: 'σ²', showPercentages: true, width: 460, height: 320 }
)

const lmmGroupData = [1, 2, 3].map(g => ({
  label: `Group ${g}`,
  values: lmmY.filter((_, i) => lmmGroups[i] === g),
}))
renderViolinBox(
  document.getElementById('lmm-groups-plot')!,
  { groups: lmmGroupData },
  { title: 'Observed Values by Group', xLabel: 'Group', yLabel: 'Y', width: 500, height: 320 }
)

tbl(lmmTables, {
  title: 'Table 11.  Variance Components and Model Fit',
  cols: [
    { h: 'Component / Index',  a: 'l' },
    { h: '&sigma;&sup2;',      a: 'r' },
    { h: '&sigma;',            a: 'r' },
    { h: '% Total',            a: 'r' },
  ],
  rows: [
    ['Between-group (<em>u</em><sub>0</sub>)', vc.intercept.toFixed(4), Math.sqrt(vc.intercept).toFixed(4), `${(100 * vc.intercept / totalVar).toFixed(1)}%`],
    ['Within-group (residual &epsilon;)',       vc.residual.toFixed(4),  Math.sqrt(vc.residual).toFixed(4),  `${(100 * vc.residual / totalVar).toFixed(1)}%`],
    ['<strong>Total</strong>',                  totalVar.toFixed(4),     Math.sqrt(totalVar).toFixed(4),     '100.0%'],
    ['ICC',                                     lmmResult.icc.toFixed(4), '—', '—'],
    ['log<em>L</em>',                           lmmResult.logLik.toFixed(3), '—', '—'],
    ['AIC',                                     lmmResult.aic.toFixed(2), '—', '—'],
    ['BIC',                                     lmmResult.bic.toFixed(2), '—', '—'],
  ],
  note: `ICC = σ²_b / (σ²_b + σ²_e) = ${lmmResult.icc.toFixed(4)}: ${(lmmResult.icc * 100).toFixed(1)}% of total variance is between groups. `
      + 'AIC/BIC: lower is better. Estimated via REML.',
})

// ═══════════════════════════════════════════════════════════════════════════
// PCA SECTION
// ═══════════════════════════════════════════════════════════════════════════

const pca = runPCA(pcaData, 3)
const varLabels = ['V1', 'V2', 'V3', 'V4', 'V5']

renderPCAPlot(
  document.getElementById('scree-plot')!,
  pca,
  {
    title: 'Scree Plot — Eigenvalues',
    type: 'scree', width: 500, height: 350, variableLabels: varLabels,
    caption: 'Kaiser criterion: retain components with λ ≥ 1 (above dashed line). Elbow rule: look for the bend.',
  }
)
renderPCAPlot(
  document.getElementById('biplot')!,
  pca,
  {
    title: 'PCA Biplot',
    type: 'biplot', width: 550, height: 500, variableLabels: varLabels,
    caption: `n = 60, p = 5 variables. Arrows: variable loadings on PC1/PC2. Scores: individual observations.`,
  }
)

const pcaTables = document.getElementById('pca-tables')!
const nTotalVars = pcaData[0]!.length
tbl(pcaTables, {
  title: 'Table 12.  Eigenvalues and Variance Explained',
  subtitle: `PCA via SVD (n = 60, variables = ${nTotalVars})`,
  cols: [
    { h: 'Component',   a: 'l' },
    { h: '&lambda;',    a: 'r' },
    { h: 'Var%',        a: 'r' },
    { h: 'Cum%',        a: 'r' },
    { h: '&lambda; &ge; 1?', a: 'c' },
  ],
  rows: Array.from({ length: pca.nComponents }, (_, i) => {
    const ev    = pca.eigenvalues[i] ?? 0
    const varPc = ((pca.varianceExplained[i] ?? 0) * 100).toFixed(1)
    const cumPc = ((pca.cumulativeVariance[i] ?? 0) * 100).toFixed(1)
    const keep  = ev >= 1
      ? '<span style="color:#198754;font-weight:600">&#10003;</span>'
      : '<span style="color:#dc3545">&#10005;</span>'
    return [`<strong>PC${i + 1}</strong>`, ev.toFixed(4), `${varPc}%`, `${cumPc}%`, keep]
  }),
  note: 'Kaiser criterion: retain components with λ ≥ 1. Variance extracted via SVD on the standardised data matrix.',
})

tbl(pcaTables, {
  title: 'Table 13.  Component Loadings',
  subtitle: 'Unrotated. Values with |loading| ≥ .30 in bold.',
  cols: [
    { h: 'Variable', a: 'l' },
    ...Array.from({ length: pca.nComponents }, (_, i) => ({ h: `<strong>PC${i + 1}</strong>`, a: 'r' as Align })),
    { h: 'Communality', a: 'r' },
  ],
  rows: varLabels.map((label, vi) => {
    const loadings = Array.from({ length: pca.nComponents }, (_, ci) => pca.loadings[ci]?.[vi] ?? 0)
    const communality = loadings.reduce((s, l) => s + l ** 2, 0)
    return [
      `<strong>${label}</strong>`,
      ...loadings.map(loading => {
        const fmtd = loading.toFixed(3).replace(/^-0\./, '&minus;.').replace(/^0\./, '.')
        return Math.abs(loading) >= 0.30 ? `<strong>${fmtd}</strong>` : `<span class="stat-dim">${fmtd}</span>`
      }),
      communality.toFixed(3),
    ]
  }),
  note: 'Communality = sum of squared loadings across retained components (h²). Higher values indicate more variance captured by the solution.',
})

// ═══════════════════════════════════════════════════════════════════════════
// DISTRIBUTIONS SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderDistribution(
  document.getElementById('dist-normal')!,
  { distribution: 'normal', params: { mean: 0, sd: 1 }, highlightX: 1.96 },
  { title: 'Standard Normal — N(0,1)', width: 550, height: 320,
    caption: 'Shaded area: P(Z > 1.96) = .025. Dashed line at z = 1.96 (two-tailed α = .05 critical value).' }
)
renderDistribution(
  document.getElementById('dist-t')!,
  { distribution: 't', params: { df: 10 }, highlightX: 2.228 },
  { title: 't-Distribution (df = 10)', width: 550, height: 320,
    caption: 'Critical value t(10) = 2.228 at α = .05 (two-tailed). Heavy tails vs. normal: uncertainty from estimating σ.' }
)

// ═══════════════════════════════════════════════════════════════════════════
// RANKING SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderLollipop(
  document.getElementById('lollipop-plot')!,
  { labels: rankLabels, values: rankValues },
  { title: 'Lollipop Chart — Category Scores', xLabel: 'Score', sorted: true, width: 600, height: 400 }
)

const rng6b = randn(44)
const dotGroup1 = rankValues.map(v => v + rng6b() * 5)
renderDotPlot(
  document.getElementById('dot-plot-plot')!,
  { labels: rankLabels, values: rankValues, group2: dotGroup1, group1Label: 'Before', group2Label: 'After' },
  { title: 'Cleveland Dot Plot — Before vs. After', xLabel: 'Score', width: 600, height: 400 }
)

renderPareto(
  document.getElementById('pareto-plot')!,
  { labels: rankLabels, values: rankValues },
  {
    title: 'Pareto Chart', xLabel: 'Category', yLabel: 'Count', width: 650, height: 420,
    caption: 'Bars: counts (descending). Line: cumulative %. Horizontal rule: 80% threshold (Pareto principle).',
  }
)

renderFunnel(
  document.getElementById('funnel-plot')!,
  { stages: funnelStages },
  {
    title: 'Sales Funnel', width: 600, height: 440,
    caption: 'Width proportional to value. Labels: stage value + % of previous stage.',
  }
)

// ═══════════════════════════════════════════════════════════════════════════
// COMPOSITION SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderPieChart(
  document.getElementById('pie-plot')!,
  { slices: pieSlices },
  { title: 'Market Share — Pie Chart', showPercentages: true, donut: false, width: 500, height: 420 }
)
renderPieChart(
  document.getElementById('donut-plot')!,
  { slices: pieSlices },
  { title: 'Market Share — Donut Chart', showPercentages: true, donut: true, width: 500, height: 420 }
)

renderAreaChart(
  document.getElementById('area-plot')!,
  { series: [
    { label: 'Series A', x: months, y: trendA },
    { label: 'Series B', x: months, y: trendB },
    { label: 'Series C', x: months, y: trendC },
  ] },
  { title: 'Area Chart — Monthly Trends', xLabel: 'Month', yLabel: 'Value', stacked: false, width: 700, height: 400 }
)

renderTreemap(
  document.getElementById('treemap-plot')!,
  { children: treemapChildren },
  { title: 'Treemap — Revenue by Category', width: 700, height: 420,
    caption: 'Cell area proportional to value. Colour by product group.' }
)

renderWaffleChart(
  document.getElementById('waffle-plot')!,
  { slices: waffleSlices },
  { title: 'Waffle Chart — Market Share', width: 400, height: 380,
    caption: '10 × 10 grid = 100 cells. Each cell = 1% of total. Largest-remainder rounding.' }
)

// ═══════════════════════════════════════════════════════════════════════════
// TIME SERIES SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderLineChart(
  document.getElementById('linechart-plot')!,
  { series: [
    { label: 'Series A', x: months, y: trendA },
    { label: 'Series B', x: months, y: trendB },
    { label: 'Series C', x: months, y: trendC },
  ] },
  { title: 'Multi-line Time Series', xLabel: 'Month', yLabel: 'Value', showArea: false, width: 700, height: 400 }
)

const spkRng = randu(91)
const spkA = Array.from({ length: 20 }, (_, i) => 50 + i * 2 + (spkRng() - 0.5) * 10)
const spkB = Array.from({ length: 20 }, (_, i) => 80 - i * 1.5 + (spkRng() - 0.5) * 8)
const spkC = Array.from({ length: 20 }, () => 60 + (spkRng() - 0.5) * 20)

renderSparkline(document.getElementById('sparkline-a')!, { values: spkA }, { title: 'Revenue', showArea: true, width: 300, height: 70 })
renderSparkline(document.getElementById('sparkline-b')!, { values: spkB }, { title: 'Costs',   showArea: true, width: 300, height: 70 })
renderSparkline(document.getElementById('sparkline-c')!, { values: spkC }, { title: 'Margin',  showArea: false, width: 300, height: 70 })

// ═══════════════════════════════════════════════════════════════════════════
// STATISTICAL SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderForestPlot(
  document.getElementById('forest-plot')!,
  { studies: forestStudies, pooled: pooledEstimate },
  {
    title: 'Forest Plot — Meta-Analysis',
    xLabel: 'Effect Size (Hedges g)', width: 740, height: 420,
    caption: `k = ${forestStudies.length} studies. Square size proportional to weight. Diamond: pooled estimate ${pooledEstimate.estimate.toFixed(2)} [${pooledEstimate.ciLow.toFixed(2)}, ${pooledEstimate.ciHigh.toFixed(2)}]. Vertical line: null effect (g = 0).`,
  }
)

renderROCCurve(
  document.getElementById('roc-plot')!,
  { fpr: rocFpr, tpr: rocTpr, auc: rocAuc },
  {
    title: 'ROC Curve — Classifier Performance',
    width: 500, height: 500,
    caption: `AUC = ${rocAuc} (area under the ROC curve). Diagonal dashed line: random classifier (AUC = 0.50). Higher AUC = better discrimination.`,
  }
)

// ═══════════════════════════════════════════════════════════════════════════
// MULTIVARIATE SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderRadarChart(
  document.getElementById('radar-plot')!,
  { series: radarSeries, axes: radarAxes },
  {
    title: 'Radar Chart — Athlete Profiles',
    width: 600, height: 500,
    caption: 'Each polygon = one athlete. Values normalised to [0,1] per axis (min–max across athletes).',
  }
)

renderParallelCoords(
  document.getElementById('parallel-plot')!,
  { rows: pcRows, axes: pcAxes, groups: pcGroups },
  {
    title: 'Parallel Coordinates — Iris-like Data',
    width: 700, height: 420,
    caption: 'Each polyline = one observation. Colour by species group. Axes independently min–max scaled.',
  }
)

renderPairPlot(
  document.getElementById('pair-plot')!,
  { data: pairData, labels: pairLabels },
  {
    title: 'Scatter Matrix (Pair Plot)',
    width: 700, height: 700,
    caption: 'Off-diagonal: scatter plots with Pearson r. Diagonal: histogram of each variable.',
  }
)

// ═══════════════════════════════════════════════════════════════════════════
// CATEGORICAL SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderStripPlot(
  document.getElementById('strip-plot')!,
  { groups: swarmGroups },
  {
    title: 'Strip Plot — Jittered Data Points',
    xLabel: 'Group', yLabel: 'Value', width: 600, height: 420,
    caption: 'Jittered points (seeded pseudo-random). Horizontal line: group mean. n = 25 per group.',
  }
)

renderSwarmPlot(
  document.getElementById('swarm-plot')!,
  { groups: swarmGroups },
  {
    title: 'Beeswarm — Collision-free Layout',
    xLabel: 'Group', yLabel: 'Value', width: 600, height: 440,
    caption: 'Greedy collision-avoidance: points displaced left/right to prevent overlap while preserving y = exact value.',
  }
)

renderMosaicPlot(
  document.getElementById('mosaic-plot')!,
  { table: mosaicTable, rowLabels: mosaicRowLabels, colLabels: mosaicColLabels },
  {
    title: 'Mosaic Plot — Education × Income',
    width: 620, height: 460,
    caption: 'Column width ∝ column marginals. Row height ∝ conditional proportion within column. Colour: Pearson residual (blue = positive, red = negative association).',
  }
)

// ═══════════════════════════════════════════════════════════════════════════
// HIERARCHICAL & PROPORTIONAL SECTION
// ═══════════════════════════════════════════════════════════════════════════

// ── Sunburst: academic publication hierarchy ────────────────────────────────
renderSunburst(
  document.getElementById('sunburst-plot')!,
  {
    root: {
      name: 'All Publications',
      children: [
        {
          name: 'Quantitative',
          children: [
            { name: 'Experiments', value: 142 },
            { name: 'Surveys',     value: 98  },
            { name: 'Regression',  value: 87  },
            { name: 'SEM',         value: 54  },
          ],
        },
        {
          name: 'Qualitative',
          children: [
            { name: 'Interviews',   value: 76  },
            { name: 'Case Studies', value: 58  },
            { name: 'Ethnography',  value: 34  },
          ],
        },
        {
          name: 'Mixed Methods',
          children: [
            { name: 'Sequential',   value: 61  },
            { name: 'Concurrent',   value: 49  },
            { name: 'Transformative', value: 23 },
          ],
        },
        {
          name: 'Reviews',
          children: [
            { name: 'Systematic',   value: 88  },
            { name: 'Meta-Analysis',value: 67  },
            { name: 'Scoping',      value: 41  },
          ],
        },
      ],
    },
  },
  {
    title: 'Publication Methodology Breakdown',
    caption: 'N = 878 publications across 4 methodology families and 12 sub-categories. Ring 1: family, Ring 2: sub-category.',
    width: 520, height: 520,
    innerRadius: 0.22,
  }
)

// ── Marimekko: market share by region and product segment ─────────────────
renderMarimekko(
  document.getElementById('marimekko-plot')!,
  {
    categories: ['North America', 'Europe', 'Asia-Pacific', 'LatAm', 'MEA'],
    series: [
      { label: 'Premium',    values: [320, 280, 190,  55,  40] },
      { label: 'Mid-range',  values: [410, 350, 520, 120,  90] },
      { label: 'Budget',     values: [180, 140, 640, 180, 130] },
      { label: 'Enterprise', values: [290, 230, 150,  45,  60] },
    ],
  },
  {
    title: 'Market Share by Region & Segment',
    xLabel: 'Region (column width ∝ total market size)',
    yLabel: 'Segment share (%)',
    caption: 'Column width proportional to regional market size. Cell height = segment share within region. N = 5,173 M units.',
    width: 680, height: 500,
    showPercentLabels: true,
  }
)

// ═══════════════════════════════════════════════════════════════════════════
// NETWORKS & FLOWS SECTION
// ═══════════════════════════════════════════════════════════════════════════

// ── Chord diagram: inter-disciplinary collaboration matrix ─────────────────
// Rows/cols = Academic School. Cell [i][j] = joint publications between schools.
renderChordDiagram(
  document.getElementById('chord-plot')!,
  {
    labels: ['Education', 'Psychology', 'Comp Sci', 'Statistics', 'Linguistics'],
    matrix: [
      [  0, 42, 18, 12,  8],
      [ 42,  0, 29, 31, 54],
      [ 18, 29,  0, 47, 11],
      [ 12, 31, 47,  0, 22],
      [  8, 54, 11, 22,  0],
    ],
  },
  {
    title: 'Inter-Disciplinary Collaboration',
    caption: 'Chord width ∝ number of joint publications between departments. N = 5 academic schools, 2019–2024.',
    width: 620, height: 540,
  }
)

// ── Alluvial: student academic pathway across 3 checkpoints ───────────────
renderAlluvialPlot(
  document.getElementById('alluvial-plot')!,
  {
    nodes: [
      // Enrolment stage
      { id: 'stem-0',     label: 'STEM',        stage: 0 },
      { id: 'hum-0',      label: 'Humanities',  stage: 0 },
      { id: 'bus-0',      label: 'Business',    stage: 0 },
      // Year 2 stage
      { id: 'pass-1',     label: 'Passing',     stage: 1 },
      { id: 'cond-1',     label: 'Conditional', stage: 1 },
      { id: 'fail-1',     label: 'Failed',      stage: 1 },
      // Graduation stage
      { id: 'hon-2',      label: 'Honours',     stage: 2 },
      { id: 'ord-2',      label: 'Ordinary',    stage: 2 },
      { id: 'drop-2',     label: 'Dropout',     stage: 2 },
    ],
    flows: [
      // Enrolment → Year 2
      { source: 'stem-0', target: 'pass-1', value: 180 },
      { source: 'stem-0', target: 'cond-1', value:  55 },
      { source: 'stem-0', target: 'fail-1', value:  25 },
      { source: 'hum-0',  target: 'pass-1', value: 120 },
      { source: 'hum-0',  target: 'cond-1', value:  60 },
      { source: 'hum-0',  target: 'fail-1', value:  40 },
      { source: 'bus-0',  target: 'pass-1', value: 130 },
      { source: 'bus-0',  target: 'cond-1', value:  50 },
      { source: 'bus-0',  target: 'fail-1', value:  20 },
      // Year 2 → Graduation
      { source: 'pass-1', target: 'hon-2',  value: 290 },
      { source: 'pass-1', target: 'ord-2',  value: 100 },
      { source: 'cond-1', target: 'ord-2',  value: 110 },
      { source: 'cond-1', target: 'drop-2', value:  55 },
      { source: 'fail-1', target: 'drop-2', value:  85 },
    ],
    stageLabels: ['Enrolment', 'Year 2', 'Graduation'],
  },
  {
    title: 'Student Academic Pathway',
    caption: 'N = 680 students tracked across 3 checkpoints. Flow width ∝ number of students following each pathway.',
    width: 700, height: 460,
  }
)

// ── Arc diagram: co-authorship network ─────────────────────────────────────
renderArcDiagram(
  document.getElementById('arc-plot')!,
  {
    nodes: [
      { id: 'Adams',    group: 'Quant'  },
      { id: 'Brown',    group: 'Qual'   },
      { id: 'Chen',     group: 'Quant'  },
      { id: 'Davies',   group: 'Mixed'  },
      { id: 'Evans',    group: 'Quant'  },
      { id: 'Flores',   group: 'Qual'   },
      { id: 'Garcia',   group: 'Mixed'  },
      { id: 'Harris',   group: 'Quant'  },
      { id: 'Ibrahim',  group: 'Qual'   },
      { id: 'Jensen',   group: 'Mixed'  },
    ],
    edges: [
      { source: 'Adams',   target: 'Chen',    value: 3 },
      { source: 'Adams',   target: 'Evans',   value: 2 },
      { source: 'Adams',   target: 'Harris',  value: 1 },
      { source: 'Brown',   target: 'Flores',  value: 4 },
      { source: 'Brown',   target: 'Ibrahim', value: 2 },
      { source: 'Chen',    target: 'Davies',  value: 1 },
      { source: 'Chen',    target: 'Evans',   value: 3 },
      { source: 'Davies',  target: 'Garcia',  value: 2 },
      { source: 'Davies',  target: 'Jensen',  value: 1 },
      { source: 'Evans',   target: 'Harris',  value: 4 },
      { source: 'Flores',  target: 'Garcia',  value: 2 },
      { source: 'Garcia',  target: 'Jensen',  value: 3 },
      { source: 'Harris',  target: 'Jensen',  value: 1 },
      { source: 'Ibrahim', target: 'Jensen',  value: 2 },
    ],
  },
  {
    title: 'Co-authorship Network',
    sortByGroup: true,
    caption: 'Nodes grouped by methodology (Quantitative / Qualitative / Mixed). Arc height ∝ node distance. Arc weight ∝ joint publications.',
    width: 700, height: 420,
  }
)

// ── Edge bundling: software module dependencies ────────────────────────────
renderEdgeBundling(
  document.getElementById('edge-bundling-plot')!,
  {
    nodes: [
      // Root
      { id: 'root',     label: 'Carm',      parent: '' },
      // Inner
      { id: 'core',     label: 'core',      parent: 'root' },
      { id: 'stats',    label: 'stats',     parent: 'root' },
      { id: 'viz',      label: 'viz',       parent: 'root' },
      // core leaves
      { id: 'types',    label: 'types',     parent: 'core',  group: 'core'  },
      { id: 'matrix',   label: 'matrix',    parent: 'core',  group: 'core'  },
      { id: 'math',     label: 'math',      parent: 'core',  group: 'core'  },
      { id: 'apa',      label: 'apa',       parent: 'core',  group: 'core'  },
      // stats leaves
      { id: 'desc',     label: 'descriptive', parent: 'stats', group: 'stats' },
      { id: 'comp',     label: 'comparison',  parent: 'stats', group: 'stats' },
      { id: 'corr',     label: 'correlation', parent: 'stats', group: 'stats' },
      { id: 'reg',      label: 'regression',  parent: 'stats', group: 'stats' },
      { id: 'mixed',    label: 'mixed',       parent: 'stats', group: 'stats' },
      { id: 'pca',      label: 'pca',         parent: 'stats', group: 'stats' },
      // viz leaves
      { id: 'violin',   label: 'violin-box',  parent: 'viz',   group: 'viz'   },
      { id: 'scatter',  label: 'scatter',     parent: 'viz',   group: 'viz'   },
      { id: 'hist',     label: 'histogram',   parent: 'viz',   group: 'viz'   },
      { id: 'corrgram', label: 'correlogram', parent: 'viz',   group: 'viz'   },
    ],
    edges: [
      // stats → core
      { source: 'desc',     target: 'types'  },
      { source: 'desc',     target: 'math'   },
      { source: 'comp',     target: 'types'  },
      { source: 'comp',     target: 'math'   },
      { source: 'corr',     target: 'types'  },
      { source: 'corr',     target: 'math'   },
      { source: 'reg',      target: 'matrix' },
      { source: 'reg',      target: 'math'   },
      { source: 'mixed',    target: 'matrix' },
      { source: 'mixed',    target: 'math'   },
      { source: 'pca',      target: 'matrix' },
      // viz → stats
      { source: 'violin',   target: 'comp'   },
      { source: 'scatter',  target: 'corr'   },
      { source: 'scatter',  target: 'reg'    },
      { source: 'hist',     target: 'desc'   },
      { source: 'corrgram', target: 'corr'   },
      // viz → core
      { source: 'violin',   target: 'apa'    },
      { source: 'scatter',  target: 'apa'    },
      { source: 'hist',     target: 'types'  },
    ],
  },
  {
    title: 'Carm Module Dependency Graph',
    caption: 'Hierarchical edge bundling (β = 0.85). Edges routed through lowest common ancestor — shared dependencies cluster visually.',
    width: 620, height: 620,
    beta: 0.85,
  }
)
