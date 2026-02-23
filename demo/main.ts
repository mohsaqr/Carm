/**
 * Carm demo application.
 * Demonstrates all major analysis and visualization modules
 * with synthetic data designed to produce interesting, readable plots.
 */

import {
  describe as describeStats,
  tTestIndependent,
  oneWayANOVA,
  tukeyHSD,
  pearsonCorrelation,
  linearRegression,
  multipleRegression,
  regressionDiagnostics,
  runLMM,
  computeBLUPs,
  runPCA,
  correlationMatrix,
  renderHistogram,
  renderQQPlot,
  renderViolinBox,
  renderRaincloud,
  renderScatterStats,
  renderCorrelogram,
  renderCoefPlot,
  renderResidualPanel,
  renderMixedPlot,
  renderPCAPlot,
  renderDistribution,
} from 'carm'

// ─── Synthetic datasets ───────────────────────────────────────────────────

function randn(seed: number): () => number {
  // Box-Muller with LCG seeding
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

const rng = randn(42)

// Group comparison: 3 groups, n=30 each
const groupA = Array.from({ length: 30 }, () => 5 + rng() * 1.5)
const groupB = Array.from({ length: 30 }, () => 7 + rng() * 1.8)
const groupC = Array.from({ length: 30 }, () => 6 + rng() * 1.2)

// Correlation: positively correlated pair
const x1 = Array.from({ length: 50 }, (_, i) => i * 0.1 + rng() * 0.5)
const y1 = x1.map(x => 2.5 + 1.8 * x + rng() * 0.8)

// Multiple regression: y = 1 + 2*x1 + 0.5*x2 + noise
const rng2 = randn(99)
const pred1 = Array.from({ length: 40 }, () => rng2() * 3)
const pred2 = Array.from({ length: 40 }, () => rng2() * 2 + 1)
const yReg = pred1.map((x, i) => 1 + 2 * x + 0.5 * (pred2[i]!) + rng2() * 0.8)

// LMM data: 3 groups × 10 obs
const rng3 = randn(7)
const lmmGroups = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
const lmmX = lmmGroups.map(() => rng3() * 2 + 1)
const groupEffects: Record<number, number> = { 1: -2, 2: 0, 3: 3 }
const lmmY = lmmX.map((x, i) => {
  const g = lmmGroups[i]!
  return (groupEffects[g] ?? 0) + 1.2 * x + rng3() * 0.5
})

// PCA data: 5 variables, 60 obs
const rng4 = randn(13)
const pcaData = Array.from({ length: 60 }, () => {
  const f1 = rng4(), f2 = rng4()
  return [f1 + rng4() * 0.3, f1 * 0.8 + rng4() * 0.4, f1 * -0.6 + rng4() * 0.5, f2 + rng4() * 0.3, f2 * 0.7 + rng4() * 0.4]
})

// ─── Descriptive section ──────────────────────────────────────────────────

const allGroupValues = [...groupA, ...groupB, ...groupC]
const descriptives = describeStats(allGroupValues)

renderHistogram(
  document.getElementById('histogram-plot')!,
  { values: allGroupValues, descriptives },
  { title: 'Distribution of Group Values', xLabel: 'Value', showNormalCurve: true, width: 650, height: 380 }
)

renderQQPlot(
  document.getElementById('qq-plot')!,
  allGroupValues,
  { title: 'Normal QQ Plot', width: 450, height: 380 }
)

document.getElementById('descriptive-stats')!.textContent = [
  `n = ${descriptives.n}`,
  `Mean = ${descriptives.mean.toFixed(3)}, Median = ${descriptives.median.toFixed(3)}`,
  `SD = ${descriptives.sd.toFixed(3)}, SE = ${descriptives.se.toFixed(3)}`,
  `Skewness = ${descriptives.skewness.toFixed(3)}, Kurtosis = ${descriptives.kurtosis.toFixed(3)}`,
  `95% CI: [${descriptives.ci[0].toFixed(3)}, ${descriptives.ci[1].toFixed(3)}]`,
  `Shapiro-Wilk: W = ${descriptives.shapiroWilk.statistic.toFixed(3)}, p = ${descriptives.shapiroWilk.pValue.toFixed(3)}`,
].join('\n')

// ─── Comparison section ───────────────────────────────────────────────────

const groups = [
  { label: 'Group A', values: groupA },
  { label: 'Group B', values: groupB },
  { label: 'Group C', values: groupC },
]

const anovaResult = oneWayANOVA(groups)
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

document.getElementById('comparison-stats')!.textContent = anovaResult.formatted

// ─── Correlation section ──────────────────────────────────────────────────

const corrResult = pearsonCorrelation(x1, y1)
const regForScatter = linearRegression(x1, y1)

renderScatterStats(
  document.getElementById('scatter-plot')!,
  { x: x1, y: y1, correlationResult: corrResult, regressionResult: regForScatter },
  { title: 'Scatter with Regression', xLabel: 'X', yLabel: 'Y', width: 650, height: 480 }
)

const corrMatrix = correlationMatrix(pcaData.map(row => row.slice(0, 4)).reduce<number[][]>(
  (acc, _, i) => { acc.push(pcaData.map(row => row[i] ?? 0)); return acc }, []
), ['V1', 'V2', 'V3', 'V4'])

renderCorrelogram(
  document.getElementById('correlogram-plot')!,
  corrMatrix,
  { title: 'Correlation Matrix', width: 500, height: 500 }
)

// ─── Regression section ───────────────────────────────────────────────────

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

// ─── LMM section ─────────────────────────────────────────────────────────

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

document.getElementById('lmm-stats')!.textContent = [
  lmmResult.formatted,
  '',
  'Fixed Effects:',
  ...lmmResult.fixedEffects.map(fe =>
    `  ${fe.name}: β = ${fe.estimate.toFixed(3)}, SE = ${fe.se.toFixed(3)}, t = ${fe.tValue.toFixed(2)}, p = ${fe.pValue.toFixed(3)}`
  ),
].join('\n')

// ─── PCA section ─────────────────────────────────────────────────────────

const pca = runPCA(pcaData, 3)

renderPCAPlot(
  document.getElementById('scree-plot')!,
  pca,
  { title: 'Scree Plot', type: 'scree', width: 500, height: 350,
    variableLabels: ['V1', 'V2', 'V3', 'V4', 'V5'] }
)

renderPCAPlot(
  document.getElementById('biplot')!,
  pca,
  { title: 'PCA Biplot', type: 'biplot', width: 550, height: 500,
    variableLabels: ['V1', 'V2', 'V3', 'V4', 'V5'] }
)

// ─── Distribution section ─────────────────────────────────────────────────

renderDistribution(
  document.getElementById('dist-normal')!,
  { distribution: 'normal', params: { mean: 0, sd: 1 }, highlightX: 1.96 },
  { title: 'Standard Normal Distribution', width: 550, height: 320 }
)

renderDistribution(
  document.getElementById('dist-t')!,
  { distribution: 't', params: { df: 10 }, highlightX: 2.228 },
  { title: 't-Distribution (df=10)', width: 550, height: 320 }
)
