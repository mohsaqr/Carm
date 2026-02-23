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

// ─── Table utilities ───────────────────────────────────────────────────────

type Align = 'l' | 'r' | 'c'

/** APA-style significance stars as superscript HTML. */
function sig(p: number): string {
  const s = p < 0.001 ? '***' : p < 0.01 ? '**' : p < 0.05 ? '*' : p < 0.10 ? '†' : ''
  return s ? `<sup class="stat-sig">${s}</sup>` : ''
}

/** APA p-value: no leading zero, "< .001" for tiny values. */
function fmtP(p: number): string {
  if (p < 0.001) return '&lt;&thinsp;.001'
  return p.toFixed(3).replace(/^0\./, '.')
}

/** Format a confidence interval. */
function fmtCI(lo: number, hi: number, d = 2): string {
  return `[${lo.toFixed(d)},&thinsp;${hi.toFixed(d)}]`
}

/** Format correlation r — no leading zero, handles negatives cleanly. */
function fmtR(r: number): string {
  const s = r.toFixed(3)
  return r < 0 ? s.replace('-0.', '&minus;.') : s.replace('0.', '.')
}

/**
 * Render a single publication-quality table into a container.
 * Appends HTML — call multiple times for multiple tables in one container.
 */
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

// ═══════════════════════════════════════════════════════════════════════════
// DESCRIPTIVE SECTION
// ═══════════════════════════════════════════════════════════════════════════

const allValues   = [...groupA, ...groupB, ...groupC]
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

const cTables = document.getElementById('comparison-tables')!

// Table 1: Group descriptives
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

// Table 2: ANOVA summary
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
    [
      'Between groups',
      ar.ssBetween.toFixed(3),
      `${ar.dfBetween}`,
      ar.msBetween.toFixed(3),
      ar.statistic.toFixed(3) + sig(ar.pValue),
      fmtP(ar.pValue),
      ar.effectSize.value.toFixed(3),
    ],
    [
      'Within groups',
      ar.ssWithin.toFixed(3),
      `${ar.dfWithin}`,
      ar.msWithin.toFixed(3),
      '', '', '',
    ],
    [
      '<strong>Total</strong>',
      ar.ssTotal.toFixed(3),
      `${ar.dfBetween + ar.dfWithin}`,
      '', '', '', '',
    ],
  ],
  note: `&omega;&sup2; = ${ar.effectSize.value.toFixed(3)} (${ar.effectSize.interpretation} effect). `
      + '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
})

// Table 3: Tukey post-hoc
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
  note: 'MD = mean difference. '
      + '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001. '
      + 'ns = not significant.',
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

// Build a correct 4-variable correlation matrix (4 arrays of 60 values)
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

// Table 1: single Pearson result
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
    // t = r * sqrt(n-2) / sqrt(1-r²)
    (cr.statistic * Math.sqrt(cr.n - 2) / Math.sqrt(1 - cr.statistic ** 2)).toFixed(3),
    fmtP(cr.pValue),
    fmtCI(cr.ci[0], cr.ci[1], 3),
    cr.effectSize.interpretation.charAt(0).toUpperCase() + cr.effectSize.interpretation.slice(1),
  ]],
  note: '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001.',
})

// Table 2: correlation matrix (lower triangle)
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
    cols: [
      { h: 'Statistic', a: 'l' },
      { h: 'Value',     a: 'r' },
    ],
    rows: [
      ['<em>R</em>&sup2;',       multiReg.r2.toFixed(4)],
      ['Adj. <em>R</em>&sup2;',  multiReg.adjR2.toFixed(4)],
      ['<em>F</em>',             multiReg.fStatistic.toFixed(3) + sig(multiReg.fPValue)],
      ['df<sub>1</sub>, df<sub>2</sub>',
        `${multiReg.fDf[0]},&thinsp;${multiReg.fDf[1]}`],
      ['<em>p</em>',             fmtP(multiReg.fPValue)],
      ['AIC',                    multiReg.aic.toFixed(2)],
      ['BIC',                    multiReg.bic.toFixed(2)],
    ],
  })

  tbl(regTables, {
    title: 'Coefficient Summary',
    cols: [
      { h: 'Predictor', a: 'l' },
      { h: '<em>&beta;</em>', a: 'r' },
      { h: '<em>SE</em>',    a: 'r' },
      { h: '<em>t</em>',     a: 'r' },
      { h: '<em>p</em>',     a: 'r' },
    ],
    rows: multiReg.coefficients.map(c => [
      c.name === '(Intercept)' ? '<em>Intercept</em>' : c.name,
      c.estimate.toFixed(3),
      c.se.toFixed(3),
      c.tValue.toFixed(3) + sig(c.pValue),
      fmtP(c.pValue),
    ]),
  })

// Full coefficient table with CIs
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
  note: 'OLS regression. '
      + '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001. '
      + 'ns = not significant.',
})

// ═══════════════════════════════════════════════════════════════════════════
// LMM SECTION
// ═══════════════════════════════════════════════════════════════════════════

const lmmResult = runLMM({
  outcome: lmmY,
  fixedPredictors: { x: lmmX },
  groupId: lmmGroups,
})
const blups = computeBLUPs(
  { outcome: lmmY, fixedPredictors: { x: lmmX }, groupId: lmmGroups },
  lmmResult
)

renderMixedPlot(
  document.getElementById('lmm-plot')!,
  lmmResult,
  blups,
  { title: 'Random Effects (BLUPs)', width: 550, height: 280 }
)

const lmmTables = document.getElementById('lmm-tables')!

// Table 1: fixed effects
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
  note: '&dagger;&thinsp;<em>p</em> &lt; .10, *&thinsp;<em>p</em> &lt; .05, **&thinsp;<em>p</em> &lt; .01, ***&thinsp;<em>p</em> &lt; .001. ns = not significant.',
})

const vc = lmmResult.varianceComponents
const totalVar = vc.intercept + vc.residual

// Variance components bar chart
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

// Group-level observed distributions
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
      [
        'Random intercept',
        vc.intercept.toFixed(4),
        Math.sqrt(vc.intercept).toFixed(4),
        `${(100 * vc.intercept / totalVar).toFixed(1)}%`,
      ],
      [
        'Residual',
        vc.residual.toFixed(4),
        Math.sqrt(vc.residual).toFixed(4),
        `${(100 * vc.residual / totalVar).toFixed(1)}%`,
      ],
      [
        '<strong>Total</strong>',
        totalVar.toFixed(4),
        Math.sqrt(totalVar).toFixed(4),
        '100%',
      ],
    ],
    note: `ICC = ${lmmResult.icc.toFixed(4)}.`,
  })

  tbl(lmmTables, {
    title: 'Model Fit',
    cols: [
      { h: 'Index',   a: 'l' },
      { h: 'Value',   a: 'r' },
    ],
    rows: [
      ['log<em>L</em>',  lmmResult.logLik.toFixed(4)],
      ['AIC',             lmmResult.aic.toFixed(2)],
      ['BIC',             lmmResult.bic.toFixed(2)],
      ['ICC',             lmmResult.icc.toFixed(4)],
      ['<em>N</em> obs',  `${lmmResult.nObs}`],
      ['<em>k</em> groups', `${lmmResult.nGroups}`],
    ],
    note: 'Estimated via REML.',
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
const nTotalVars = pcaData[0]!.length  // 5

// Reconstruct all 5 eigenvalues (variance explained ratios * total variance)
// pca only stores nComponents, so we show what we have
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
    return [
      `<strong>PC${i + 1}</strong>`,
      ev.toFixed(4),
      `${varPc}%`,
      `${cumPc}%`,
      keep,
    ]
  }),
  note: 'Kaiser criterion: retain components with &lambda; &ge; 1.',
})

// Component loadings table
tbl(pcaTables, {
  title: 'Table 11.  Component Loadings',
  subtitle: 'Values &gt; |.30| in bold.',
  cols: [
    { h: 'Variable', a: 'l' },
    ...Array.from({ length: pca.nComponents }, (_, i) => ({
      h: `<strong>PC${i + 1}</strong>`,
      a: 'r' as Align,
    })),
  ],
  rows: varLabels.map((label, vi) => [
    `<strong>${label}</strong>`,
    ...Array.from({ length: pca.nComponents }, (_, ci) => {
      const loading = pca.loadings[ci]?.[vi] ?? 0
      const fmtd = loading.toFixed(3).replace(/^-0\./, '&minus;.').replace(/^0\./, '.')
      return Math.abs(loading) >= 0.30 ? `<strong>${fmtd}</strong>` : `<span class="stat-dim">${fmtd}</span>`
    }),
  ]),
  note: 'Loadings are unrotated. Large |loadings| (&ge; .30) indicate strong association with that component.',
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
