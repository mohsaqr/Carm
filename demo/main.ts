/**
 * Carm demo application.
 * Demonstrates all major analysis and visualization modules
 * with synthetic data, publication-quality plots, and styled tables.
 */

import {
  describe as describeStats,
  tTestIndependent,
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
const rng6 = randu(17)

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

// Time-series data: monthly values over 2 years, 3 series
const months = Array.from({ length: 24 }, (_, i) => i)
const trendA = months.map(m => 10 + m * 0.8 + rng5() * 2.0)
const trendB = months.map(m => 15 + m * 0.3 + rng5() * 2.5)
const trendC = months.map(m => 8  + m * 1.2 + rng5() * 1.5)

// Ranking data: 8 categories
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

// Forest plot: 5 studies + pooled
const forestStudies = [
  { label: 'Smith et al. (2018)', estimate: 0.45, ciLow: 0.12, ciHigh: 0.78, weight: 3.2 },
  { label: 'Jones et al. (2019)', estimate: 0.62, ciLow: 0.31, ciHigh: 0.93, weight: 4.1 },
  { label: 'Brown et al. (2020)', estimate: 0.28, ciLow: -0.05, ciHigh: 0.61, weight: 2.8 },
  { label: 'Davis et al. (2021)', estimate: 0.71, ciLow: 0.45, ciHigh: 0.97, weight: 5.0 },
  { label: 'Wilson et al. (2022)', estimate: 0.53, ciLow: 0.22, ciHigh: 0.84, weight: 3.7 },
]
const pooledEstimate = { estimate: 0.54, ciLow: 0.38, ciHigh: 0.70 }

// ROC curve: synthetic FPR/TPR from a decent classifier
const rocFpr: number[] = [0]
const rocTpr: number[] = [0]
const rocRng = randu(123)
for (let t = 1; t >= 0; t -= 0.02) {
  const fpr = Math.pow(t, 2.5) * (0.9 + rocRng() * 0.2)
  const tpr = 1 - Math.pow(1 - t, 1.8) * (0.85 + rocRng() * 0.1)
  rocFpr.push(Math.min(1, Math.max(0, fpr)))
  rocTpr.push(Math.min(1, Math.max(0, tpr)))
}
rocFpr.push(1); rocTpr.push(1)
const rocAuc = 0.82

// Radar data: 3 athletes, 5 skills
const radarAxes = ['Speed', 'Strength', 'Agility', 'Endurance', 'Technique']
const radarSeries = [
  { label: 'Athlete A', values: [85, 72, 90, 65, 88] },
  { label: 'Athlete B', values: [70, 91, 75, 82, 77] },
  { label: 'Athlete C', values: [92, 68, 83, 78, 71] },
]

// Parallel coordinates: iris-like data (4 vars, 3 groups)
const pcAxes = ['Sepal L', 'Sepal W', 'Petal L', 'Petal W']
const pcRng = randn(31)
const pcRows = Array.from({ length: 45 }, (_, i) => {
  const g = Math.floor(i / 15)
  const bases = [[5.0, 3.4, 1.5, 0.2], [5.9, 2.8, 4.3, 1.3], [6.5, 3.0, 5.5, 2.0]]
  return bases[g]!.map(b => Math.max(0.1, b + pcRng() * 0.6))
})
const pcGroups = Array.from({ length: 45 }, (_, i) => Math.floor(i / 15))

// Mosaic: 3×3 contingency table (education × income)
const mosaicTable = [
  [30, 20, 10],
  [25, 40, 20],
  [10, 25, 35],
]
const mosaicRowLabels = ['Low Ed.', 'Mid Ed.', 'High Ed.']
const mosaicColLabels = ['Low Inc.', 'Mid Inc.', 'High Inc.']

// Treemap: product categories
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

// Waffle: market share
const waffleSlices = [
  { label: 'Product A', value: 35 },
  { label: 'Product B', value: 28 },
  { label: 'Product C', value: 20 },
  { label: 'Product D', value: 17 },
]

// Funnel: sales funnel
const funnelStages = [
  { label: 'Awareness', value: 10000 },
  { label: 'Interest', value: 6500 },
  { label: 'Consideration', value: 3800 },
  { label: 'Intent', value: 2100 },
  { label: 'Purchase', value: 980 },
]

// Pair plot variables (use pcaData first 4 cols, 30 rows)
const pairData = pcaData.slice(0, 30).map(row => row.slice(0, 4))
const pairLabels = ['V1', 'V2', 'V3', 'V4']

// Swarm + strip: 3 groups, n=25 each
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
  { title: 'Distribution of Group Values', xLabel: 'Value', showNormalCurve: true, width: 650, height: 380 }
)
renderQQPlot(
  document.getElementById('qq-plot')!,
  allValues,
  { title: 'Normal QQ Plot', width: 450, height: 380 }
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
  { title: 'KDE Density by Group', xLabel: 'Value', yLabel: 'Density', showRug: true, width: 650, height: 340 }
)
renderBoxplot(
  document.getElementById('boxplot-plot')!,
  { groups: [
    { label: 'Group A', values: groupA },
    { label: 'Group B', values: groupB },
    { label: 'Group C', values: groupC },
  ] },
  { title: 'Standalone Box Plot', xLabel: 'Group', yLabel: 'Value', width: 500, height: 380 }
)

const dTables = document.getElementById('descriptive-tables')!
const d = descriptives
const normConclusion = d.shapiroWilk.pValue > 0.05
  ? '<span style="color:#198754;font-weight:600">Normal</span>'
  : '<span style="color:#dc3545;font-weight:600">Non-normal</span>'

tbl(dTables, {
  title: 'Summary Statistics',
  subtitle: `N = ${d.n} observations`,
  cols: [
    { h: 'Statistic', a: 'l' },
    { h: 'Value',     a: 'r' },
  ],
  rows: [
    ['<em>n</em>',                    `${d.n}`],
    ['Mean (<em>M</em>)',              d.mean.toFixed(3)],
    ['Standard Deviation (<em>SD</em>)', d.sd.toFixed(3)],
    ['Standard Error (<em>SE</em>)',  d.se.toFixed(3)],
    ['Median',                         d.median.toFixed(3)],
    ['IQR',                            d.iqr.toFixed(3)],
    ['Min &ndash; Max',                `${d.min.toFixed(2)} &ndash; ${d.max.toFixed(2)}`],
    ['Skewness',                       d.skewness.toFixed(3)],
    ['Excess Kurtosis',                d.kurtosis.toFixed(3)],
    [`95% CI for <em>M</em>`,         fmtCI(d.ci[0], d.ci[1])],
  ],
  note: 'Skewness and kurtosis are moment-based estimators.',
})

tbl(dTables, {
  title: 'Normality Test',
  subtitle: 'Shapiro–Wilk',
  cols: [
    { h: 'Statistic', a: 'l' },
    { h: 'Value',     a: 'r' },
  ],
  rows: [
    ['<em>W</em>',      d.shapiroWilk.statistic.toFixed(4)],
    ['<em>p</em>-value', fmtP(d.shapiroWilk.pValue) + sig(d.shapiroWilk.pValue)],
    ['Conclusion',       normConclusion],
  ],
  note: 'Null hypothesis: data are normally distributed.',
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

renderViolinBox(
  document.getElementById('violin-plot')!,
  { groups, testResult: anovaResult, pairwise: tukeyResults },
  { title: 'Group Comparison (ANOVA)', xLabel: 'Group', yLabel: 'Value', showBrackets: true, width: 650, height: 480 }
)
renderRaincloud(
  document.getElementById('raincloud-plot')!,
  { groups, testResult: anovaResult },
  { title: 'Raincloud Plot', xLabel: 'Group', yLabel: 'Value', width: 650, height: 480 }
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
  { type: 'grouped', title: 'Quarterly Sales by Product (Grouped)', xLabel: 'Quarter', yLabel: 'Sales', width: 650, height: 400 }
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
  { type: 'stacked', title: 'Quarterly Sales by Product (Stacked)', xLabel: 'Quarter', yLabel: 'Sales', width: 650, height: 400 }
)
renderBubbleChart(
  document.getElementById('bubble-plot')!,
  { points: bubblePoints },
  { title: 'Bubble Chart', xLabel: 'X Variable', yLabel: 'Y Variable', width: 650, height: 480 }
)

const cTables = document.getElementById('comparison-tables')!
const groupDescs = groups.map(g => describeStats([...g.values]))
tbl(cTables, {
  title: 'Table 1.  Descriptive Statistics by Group',
  cols: [
    { h: 'Group',   a: 'l' },
    { h: '<em>n</em>', a: 'c' },
    { h: '<em>M</em>', a: 'r' },
    { h: '<em>SD</em>', a: 'r' },
    { h: '<em>SE</em>', a: 'r' },
    { h: 'Min', a: 'r' },
    { h: 'Max', a: 'r' },
    { h: '95% CI', a: 'c' },
  ],
  rows: groups.map((g, i) => {
    const gd = groupDescs[i]!
    return [
      `<strong>${g.label}</strong>`,
      `${gd.n}`,
      gd.mean.toFixed(2),
      gd.sd.toFixed(2),
      gd.se.toFixed(2),
      gd.min.toFixed(2),
      gd.max.toFixed(2),
      fmtCI(gd.ci[0], gd.ci[1]),
    ]
  }),
})

const ar = anovaResult
tbl(cTables, {
  title: 'Table 2.  One-Way ANOVA',
  cols: [
    { h: 'Source',  a: 'l' },
    { h: '<em>SS</em>', a: 'r' },
    { h: '<em>df</em>', a: 'c' },
    { h: '<em>MS</em>', a: 'r' },
    { h: '<em>F</em>',  a: 'r' },
    { h: '<em>p</em>',  a: 'r' },
    { h: '&omega;&sup2;', a: 'r' },
  ],
  rows: [
    ['Between groups', ar.ssBetween.toFixed(3), `${ar.dfBetween}`, ar.msBetween.toFixed(3), ar.statistic.toFixed(3) + sig(ar.pValue), fmtP(ar.pValue), ar.effectSize.value.toFixed(3)],
    ['Within groups',  ar.ssWithin.toFixed(3),  `${ar.dfWithin}`,  ar.msWithin.toFixed(3),  '', '', ''],
    ['<strong>Total</strong>', ar.ssTotal.toFixed(3), `${ar.dfBetween + ar.dfWithin}`, '', '', '', ''],
  ],
  note: `&omega;&sup2; = ${ar.effectSize.value.toFixed(3)} (${ar.effectSize.interpretation} effect). `
      + '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
})

tbl(cTables, {
  title: 'Table 3.  Post-hoc Pairwise Comparisons (Tukey HSD)',
  cols: [
    { h: 'Comparison',  a: 'l' },
    { h: '<em>MD</em>', a: 'r' },
    { h: '<em>SE</em>', a: 'r' },
    { h: '<em>q</em>',  a: 'r' },
    { h: '<em>p</em><sub>adj</sub>', a: 'r' },
    { h: '95% CI',      a: 'c' },
    { h: '',            a: 'c' },
  ],
  rows: tukeyResults.map(pw => [
    `${pw.group1} vs ${pw.group2}`,
    pw.meanDiff.toFixed(3),
    pw.se.toFixed(3),
    pw.statistic.toFixed(3),
    fmtP(pw.pValueAdj),
    fmtCI(pw.ci[0], pw.ci[1]),
    sig(pw.pValueAdj) || '<span class="stat-dim">ns</span>',
  ]),
  note: 'MD = mean difference. ns = not significant.',
})

// ═══════════════════════════════════════════════════════════════════════════
// CORRELATION SECTION
// ═══════════════════════════════════════════════════════════════════════════

const corrResult = pearsonCorrelation(x1, y1)
const regForScatter = linearRegression(x1, y1)

renderScatterStats(
  document.getElementById('scatter-plot')!,
  { x: x1, y: y1, correlationResult: corrResult, regressionResult: regForScatter },
  { title: 'Scatter with Regression Line', xLabel: 'X', yLabel: 'Y', width: 650, height: 480 }
)

const cm4 = correlationMatrix(
  Array.from({ length: 4 }, (_, vi) => pcaData.map(row => row[vi] ?? 0)),
  ['V1', 'V2', 'V3', 'V4']
)
renderCorrelogram(
  document.getElementById('correlogram-plot')!,
  cm4,
  { title: 'Correlation Matrix (V1–V4)', width: 500, height: 500 }
)

const corrTables = document.getElementById('correlation-tables')!
const cr = corrResult
const crDf = typeof cr.df === 'number' ? cr.df : cr.df[0]
tbl(corrTables, {
  title: 'Table 4.  Pearson Correlation',
  subtitle: 'X and Y (n = 50)',
  cols: [
    { h: '<em>r</em>',  a: 'r' },
    { h: '<em>df</em>', a: 'c' },
    { h: '<em>t</em>',  a: 'r' },
    { h: '<em>p</em>',  a: 'r' },
    { h: '95% CI',      a: 'c' },
    { h: 'Interpretation', a: 'l' },
  ],
  rows: [[
    fmtR(cr.statistic) + sig(cr.pValue),
    `${crDf}`,
    (cr.statistic * Math.sqrt(cr.n - 2) / Math.sqrt(1 - cr.statistic ** 2)).toFixed(3),
    fmtP(cr.pValue),
    fmtCI(cr.ci[0], cr.ci[1], 3),
    cr.effectSize.interpretation.charAt(0).toUpperCase() + cr.effectSize.interpretation.slice(1),
  ]],
  note: '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
})

const cmLabels = cm4.labels
tbl(corrTables, {
  title: 'Table 5.  Correlation Matrix',
  subtitle: 'Variables V1–V4 (n = 60). Lower triangle: Pearson r.',
  cols: [
    { h: 'Variable', a: 'l' },
    ...cmLabels.map(l => ({ h: `<strong>${l}</strong>`, a: 'c' as Align })),
  ],
  rows: cmLabels.map((label, i) => [
    `<strong>${label}</strong>`,
    ...cmLabels.map((_, j): string => {
      if (i === j) return '<span class="stat-dim">—</span>'
      if (j > i)   return ''
      const r = cm4.r[i]?.[j] ?? 0
      const p = cm4.pValues[i]?.[j] ?? 1
      return fmtR(r) + sig(p)
    }),
  ]),
  note: '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
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
  { title: 'Regression Coefficients', xLabel: 'Estimate (95% CI)', width: 550, height: 280 }
)

const diag = regressionDiagnostics(multiReg, [
  { name: 'Predictor 1', values: pred1 },
  { name: 'Predictor 2', values: pred2 },
])
renderResidualPanel(
  document.getElementById('residual-plot')!,
  multiReg,
  diag.leverage,
  { title: 'Regression Diagnostics', width: 700, height: 580 }
)

const regTables = document.getElementById('regression-tables')!
tbl(regTables, {
  title: 'Table 6.  Model Summary',
  subtitle: `N = ${multiReg.n}, predictors = 2`,
  cols: [{ h: 'Statistic', a: 'l' }, { h: 'Value', a: 'r' }],
  rows: [
    ['<em>R</em>&sup2;',       multiReg.r2.toFixed(4)],
    ['Adj. <em>R</em>&sup2;',  multiReg.adjR2.toFixed(4)],
    ['<em>F</em>',             multiReg.fStatistic.toFixed(3) + sig(multiReg.fPValue)],
    ['df<sub>1</sub>, df<sub>2</sub>', `${multiReg.fDf[0]},&thinsp;${multiReg.fDf[1]}`],
    ['<em>p</em>',             fmtP(multiReg.fPValue)],
    ['AIC',                    multiReg.aic.toFixed(2)],
    ['BIC',                    multiReg.bic.toFixed(2)],
  ],
})

tbl(regTables, {
  title: 'Table 7.  Regression Coefficients with Confidence Intervals',
  cols: [
    { h: 'Predictor',  a: 'l' },
    { h: '<em>&beta;</em>', a: 'r' },
    { h: '<em>SE</em>',    a: 'r' },
    { h: '<em>t</em>',     a: 'r' },
    { h: '<em>p</em>',     a: 'r' },
    { h: '95% CI',         a: 'c' },
    { h: '',               a: 'c' },
  ],
  rows: multiReg.coefficients.map(c => [
    c.name === '(Intercept)' ? '<em>Intercept</em>' : `<strong>${c.name}</strong>`,
    c.estimate.toFixed(3),
    c.se.toFixed(3),
    c.tValue.toFixed(3),
    fmtP(c.pValue),
    fmtCI(c.ci[0], c.ci[1]),
    sig(c.pValue) || '<span class="stat-dim">ns</span>',
  ]),
  note: 'OLS regression. ns = not significant.',
})

// ═══════════════════════════════════════════════════════════════════════════
// LMM SECTION
// ═══════════════════════════════════════════════════════════════════════════

const lmmResult = runLMM({ outcome: lmmY, fixedPredictors: { x: lmmX }, groupId: lmmGroups })
const blups = computeBLUPs({ outcome: lmmY, fixedPredictors: { x: lmmX }, groupId: lmmGroups }, lmmResult)

renderMixedPlot(
  document.getElementById('lmm-plot')!,
  lmmResult, blups,
  { title: 'Random Effects (BLUPs)', width: 550, height: 280 }
)

const lmmTables = document.getElementById('lmm-tables')!
tbl(lmmTables, {
  title: 'Table 8.  Fixed Effects',
  subtitle: `Linear Mixed Model (REML). N = ${lmmResult.nObs}, groups = ${lmmResult.nGroups}`,
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
  note: 'ns = not significant.',
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
  title: 'Table 9.  Variance Components',
  cols: [
    { h: 'Component',     a: 'l' },
    { h: '&sigma;&sup2;', a: 'r' },
    { h: '&sigma;',       a: 'r' },
    { h: '% Total',       a: 'r' },
  ],
  rows: [
    ['Random intercept', vc.intercept.toFixed(4), Math.sqrt(vc.intercept).toFixed(4), `${(100 * vc.intercept / totalVar).toFixed(1)}%`],
    ['Residual',         vc.residual.toFixed(4),  Math.sqrt(vc.residual).toFixed(4),  `${(100 * vc.residual / totalVar).toFixed(1)}%`],
    ['<strong>Total</strong>', totalVar.toFixed(4), Math.sqrt(totalVar).toFixed(4), '100%'],
  ],
  note: `ICC = ${lmmResult.icc.toFixed(4)}.`,
})

// ═══════════════════════════════════════════════════════════════════════════
// PCA SECTION
// ═══════════════════════════════════════════════════════════════════════════

const pca = runPCA(pcaData, 3)
const varLabels = ['V1', 'V2', 'V3', 'V4', 'V5']

renderPCAPlot(
  document.getElementById('scree-plot')!,
  pca,
  { title: 'Scree Plot', type: 'scree', width: 500, height: 350, variableLabels: varLabels }
)
renderPCAPlot(
  document.getElementById('biplot')!,
  pca,
  { title: 'PCA Biplot', type: 'biplot', width: 550, height: 500, variableLabels: varLabels }
)

const pcaTables = document.getElementById('pca-tables')!
const nTotalVars = pcaData[0]!.length
tbl(pcaTables, {
  title: 'Table 10.  Eigenvalues and Variance Explained',
  subtitle: `PCA via SVD (n = 60, variables = ${nTotalVars}). Showing retained components.`,
  cols: [
    { h: 'Component', a: 'l' },
    { h: '&lambda;',  a: 'r' },
    { h: 'Var%',      a: 'r' },
    { h: 'Cum%',      a: 'r' },
    { h: 'Retain?',   a: 'c' },
  ],
  rows: Array.from({ length: pca.nComponents }, (_, i) => {
    const ev    = pca.eigenvalues[i] ?? 0
    const varPc = ((pca.varianceExplained[i] ?? 0) * 100).toFixed(1)
    const cumPc = ((pca.cumulativeVariance[i] ?? 0) * 100).toFixed(1)
    const keep  = ev >= 1
      ? '<span style="color:#198754;font-weight:600">&#10003; Yes</span>'
      : '<span style="color:#dc3545">&#10005; No</span>'
    return [`<strong>PC${i + 1}</strong>`, ev.toFixed(4), `${varPc}%`, `${cumPc}%`, keep]
  }),
  note: 'Kaiser criterion: retain components with &lambda; &ge; 1.',
})

tbl(pcaTables, {
  title: 'Table 11.  Component Loadings',
  subtitle: 'Values &gt; |.30| in bold.',
  cols: [
    { h: 'Variable', a: 'l' },
    ...Array.from({ length: pca.nComponents }, (_, i) => ({ h: `<strong>PC${i + 1}</strong>`, a: 'r' as Align })),
  ],
  rows: varLabels.map((label, vi) => [
    `<strong>${label}</strong>`,
    ...Array.from({ length: pca.nComponents }, (_, ci) => {
      const loading = pca.loadings[ci]?.[vi] ?? 0
      const fmtd = loading.toFixed(3).replace(/^-0\./, '&minus;.').replace(/^0\./, '.')
      return Math.abs(loading) >= 0.30 ? `<strong>${fmtd}</strong>` : `<span class="stat-dim">${fmtd}</span>`
    }),
  ]),
  note: 'Loadings are unrotated.',
})

// ═══════════════════════════════════════════════════════════════════════════
// DISTRIBUTIONS SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderDistribution(
  document.getElementById('dist-normal')!,
  { distribution: 'normal', params: { mean: 0, sd: 1 }, highlightX: 1.96 },
  { title: 'Standard Normal Distribution', width: 550, height: 320 }
)
renderDistribution(
  document.getElementById('dist-t')!,
  { distribution: 't', params: { df: 10 }, highlightX: 2.228 },
  { title: 't-Distribution (df = 10)', width: 550, height: 320 }
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
  {
    labels: rankLabels,
    values: rankValues,
    group2: dotGroup1,
    group1Label: 'Before',
    group2Label: 'After',
  },
  { title: 'Cleveland Dot Plot — Before vs After', xLabel: 'Score', width: 600, height: 400 }
)

renderPareto(
  document.getElementById('pareto-plot')!,
  { labels: rankLabels, values: rankValues },
  { title: 'Pareto Chart — Category Scores', xLabel: 'Category', yLabel: 'Count', width: 650, height: 420 }
)

renderFunnel(
  document.getElementById('funnel-plot')!,
  { stages: funnelStages },
  { title: 'Sales Funnel', width: 600, height: 440 }
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
  {
    series: [
      { label: 'Series A', x: months, y: trendA },
      { label: 'Series B', x: months, y: trendB },
      { label: 'Series C', x: months, y: trendC },
    ],
  },
  { title: 'Area Chart — Monthly Trends', xLabel: 'Month', yLabel: 'Value', stacked: false, width: 700, height: 400 }
)

renderTreemap(
  document.getElementById('treemap-plot')!,
  { children: treemapChildren },
  { title: 'Treemap — Revenue by Category', width: 700, height: 420 }
)

renderWaffleChart(
  document.getElementById('waffle-plot')!,
  { slices: waffleSlices },
  { title: 'Waffle Chart — Market Share', width: 400, height: 380 }
)

// ═══════════════════════════════════════════════════════════════════════════
// TIME SERIES SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderLineChart(
  document.getElementById('linechart-plot')!,
  {
    series: [
      { label: 'Series A', x: months, y: trendA },
      { label: 'Series B', x: months, y: trendB },
      { label: 'Series C', x: months, y: trendC },
    ],
  },
  { title: 'Multi-line Time Series', xLabel: 'Month', yLabel: 'Value', showArea: false, width: 700, height: 400 }
)

const spkRng = randu(91)
const spkA = Array.from({ length: 20 }, (_, i) => 50 + i * 2 + (spkRng() - 0.5) * 10)
const spkB = Array.from({ length: 20 }, (_, i) => 80 - i * 1.5 + (spkRng() - 0.5) * 8)
const spkC = Array.from({ length: 20 }, () => 60 + (spkRng() - 0.5) * 20)

renderSparkline(
  document.getElementById('sparkline-a')!,
  { values: spkA },
  { title: 'Revenue', showArea: true, width: 300, height: 70 }
)
renderSparkline(
  document.getElementById('sparkline-b')!,
  { values: spkB },
  { title: 'Costs', showArea: true, width: 300, height: 70 }
)
renderSparkline(
  document.getElementById('sparkline-c')!,
  { values: spkC },
  { title: 'Margin', showArea: false, width: 300, height: 70 }
)

// ═══════════════════════════════════════════════════════════════════════════
// STATISTICAL SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderForestPlot(
  document.getElementById('forest-plot')!,
  { studies: forestStudies, pooled: pooledEstimate },
  { title: 'Forest Plot — Meta-Analysis', xLabel: 'Effect Size (Hedges g)', width: 700, height: 420 }
)

renderROCCurve(
  document.getElementById('roc-plot')!,
  { fpr: rocFpr, tpr: rocTpr, auc: rocAuc },
  { title: 'ROC Curve — Classifier Performance', width: 500, height: 500 }
)

// ═══════════════════════════════════════════════════════════════════════════
// MULTIVARIATE SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderRadarChart(
  document.getElementById('radar-plot')!,
  { series: radarSeries, axes: radarAxes },
  { title: 'Radar Chart — Athlete Profiles', width: 600, height: 500 }
)

renderParallelCoords(
  document.getElementById('parallel-plot')!,
  { rows: pcRows, axes: pcAxes, groups: pcGroups },
  { title: 'Parallel Coordinates — Iris-like Data', width: 700, height: 420 }
)

renderPairPlot(
  document.getElementById('pair-plot')!,
  { data: pairData, labels: pairLabels },
  { title: 'Scatter Matrix (Pair Plot)', width: 700, height: 700 }
)

// ═══════════════════════════════════════════════════════════════════════════
// CATEGORICAL SECTION
// ═══════════════════════════════════════════════════════════════════════════

renderStripPlot(
  document.getElementById('strip-plot')!,
  { groups: swarmGroups },
  { title: 'Strip Plot — Jittered Data Points', xLabel: 'Group', yLabel: 'Value', width: 600, height: 400 }
)

renderSwarmPlot(
  document.getElementById('swarm-plot')!,
  { groups: swarmGroups },
  { title: 'Beeswarm Plot — Collision-free Layout', xLabel: 'Group', yLabel: 'Value', width: 600, height: 420 }
)

renderMosaicPlot(
  document.getElementById('mosaic-plot')!,
  { table: mosaicTable, rowLabels: mosaicRowLabels, colLabels: mosaicColLabels },
  { title: 'Mosaic Plot — Education × Income', width: 600, height: 440 }
)
