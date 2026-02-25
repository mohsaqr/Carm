/**
 * Detailed side-by-side report: Carm vs mclust on engagement data
 * Run: npx tsx tmp/engagement_report.ts
 */
import { readFileSync } from 'fs'
import { fitGMM } from '../src/stats/clustering.js'

const ref = JSON.parse(readFileSync('tmp/engagement_mclust_ref.json', 'utf-8'))
const data: number[][] = ref.data
const varNames: string[] = ref.varNames

// ── Fit Carm (multi-seed) ──
let best: ReturnType<typeof fitGMM> | null = null
for (const seed of [42, 1, 7, 13, 99, 123, 256, 500, 1000, 2024, 3141, 9999]) {
  try {
    const res = fitGMM(data, { k: 3, model: 'VVI', seed, tol: 1e-8, maxIter: 1000 })
    if (!best || res.diagnostics.logLikelihood > best.diagnostics.logLikelihood) best = res
  } catch { /* skip */ }
}
const carm = best!

// ── Sort both by first mean (ascending) for alignment ──
const rMeansRaw: number[][] = ref.means
const rWeightsRaw: number[] = ref.weights
const rClassRaw: number[] = ref.classification
const rPosteriorsRaw: number[][] = ref.posteriors
const rAvePPRaw: number[] = ref.avepp
const rCaseEntropyRaw: number[] = ref.caseEntropy

// R cluster indices (1-indexed) → compute sizes
const rSizes = [0, 0, 0]
for (const c of rClassRaw) rSizes[c - 1]++

// Sort R clusters by first mean
const rOrder = [0, 1, 2].sort((a, b) => rMeansRaw[a]![0]! - rMeansRaw[b]![0]!)
const rMeans = rOrder.map(i => rMeansRaw[i]!)
const rWeights = rOrder.map(i => rWeightsRaw[i]!)
const rSizesSorted = rOrder.map(i => rSizes[i]!)

// Sort Carm clusters by first mean
const cOrder = [0, 1, 2].sort((a, b) => carm.means[a]![0]! - carm.means[b]![0]!)
const cMeans = cOrder.map(i => carm.means[i]!)
const cWeights = cOrder.map(i => carm.weights[i]!)
const cSizes = [0, 0, 0]
for (const l of carm.labels) cSizes[l]++
const cSizesSorted = cOrder.map(i => cSizes[i]!)

// R AvePP recomputed per sorted cluster
const rAvePPSorted = rOrder.map(i => rAvePPRaw[i]!)
// Carm AvePP per sorted cluster
const cAvePPSorted = cOrder.map(i => carm.diagnostics.avepp[i]!)

// ── Compute per-cluster entropy for Carm ──
function perClusterEntropy(posteriors: number[][], labels: number[], k: number, order: number[]) {
  const result: { mean: number; sd: number; min: number; max: number }[] = []
  for (const ci of order) {
    const eis: number[] = []
    for (let i = 0; i < posteriors.length; i++) {
      if (labels[i] === ci) {
        let rowSum = 0
        for (const p of posteriors[i]!) {
          if (p > 1e-300) rowSum += p * Math.log(p)
        }
        eis.push(1 + rowSum / Math.log(k))
      }
    }
    const mean = eis.reduce((a, b) => a + b, 0) / eis.length
    const sd = Math.sqrt(eis.reduce((s, e) => s + (e - mean) ** 2, 0) / (eis.length - 1))
    result.push({ mean, sd, min: Math.min(...eis), max: Math.max(...eis) })
  }
  return result
}

// R per-cluster entropy
function rPerClusterEntropy() {
  const result: { mean: number; sd: number; min: number; max: number }[] = []
  for (const ci of rOrder) {
    const eis: number[] = []
    for (let i = 0; i < rClassRaw.length; i++) {
      if (rClassRaw[i]! - 1 === ci) {
        eis.push(rCaseEntropyRaw[i]!)
      }
    }
    const mean = eis.reduce((a, b) => a + b, 0) / eis.length
    const sd = Math.sqrt(eis.reduce((s, e) => s + (e - mean) ** 2, 0) / (eis.length - 1))
    result.push({ mean, sd, min: Math.min(...eis), max: Math.max(...eis) })
  }
  return result
}

const cClusterEntropy = perClusterEntropy(carm.posteriors as number[][], carm.labels as number[], 3, cOrder)
const rClusterEntropy = rPerClusterEntropy()

// ── Compute per-cluster variances (diagonal for VVI) ──
function perClusterVar(data: number[][], labels: number[], means: number[][], order: number[]) {
  const d = data[0]!.length
  return order.map(ci => {
    const rows = data.filter((_, i) => labels[i] === ci)
    const n = rows.length
    const vars = new Array(d).fill(0) as number[]
    for (let j = 0; j < d; j++) {
      const m = means[ci]![j]!
      for (const row of rows) vars[j] += (row[j]! - m) ** 2
      vars[j] /= (n - 1)
    }
    return vars
  })
}
const cVars = perClusterVar(data, carm.labels as number[], carm.means as number[][], cOrder)

// ── Helpers ──
const f = (v: number, dp = 4) => v.toFixed(dp)
const pad = (s: string, w: number) => s.padStart(w)

// ═══════════════════════════════════════════════════════════════════════
// REPORT
// ═══════════════════════════════════════════════════════════════════════

console.log('╔══════════════════════════════════════════════════════════════════════════╗')
console.log('║   GMM Cross-Validation Report: Carm vs R mclust                        ║')
console.log('║   Dataset: School Engagement (N=717, D=3)                               ║')
console.log('║   Model: VVI (diagonal, variable volume/shape), K=3                     ║')
console.log('║   R: mclust 6.1.2 with priorControl()                                  ║')
console.log('║   Carm: pure MLE, K-Means++ init, multi-seed (12 seeds)                ║')
console.log('╚══════════════════════════════════════════════════════════════════════════╝')
console.log()

// ── 1. Global Fit Statistics ──
console.log('┌─────────────────────────────────────────────────────────────┐')
console.log('│ 1. GLOBAL FIT STATISTICS                                   │')
console.log('├──────────────────┬─────────────┬─────────────┬─────────────┤')
console.log('│ Metric           │    R mclust │        Carm │          |Δ|│')
console.log('├──────────────────┼─────────────┼─────────────┼─────────────┤')

const rows = [
  ['Log-Likelihood', ref.logLikelihood, carm.diagnostics.logLikelihood],
  ['BIC', ref.bic, carm.diagnostics.bic],
  ['AIC', 2 * ref.df - 2 * ref.logLikelihood, carm.diagnostics.aic],
  ['ICL', ref.icl, carm.diagnostics.icl],
  ['Entropy (norm)', ref.entropy, carm.diagnostics.entropy],
  ['DF', ref.df, carm.diagnostics.df],
  ['Converged', 1, carm.diagnostics.converged ? 1 : 0],
  ['Iterations', NaN, carm.diagnostics.iterations],
]

for (const [label, rVal, cVal] of rows) {
  const rv = typeof rVal === 'number' && !isNaN(rVal as number) ? f(rVal as number, label === 'DF' || label === 'Converged' || label === 'Iterations' ? 0 : 4) : '—'
  const cv = typeof cVal === 'number' ? f(cVal as number, label === 'DF' || label === 'Converged' || label === 'Iterations' ? 0 : 4) : String(cVal)
  const delta = typeof rVal === 'number' && typeof cVal === 'number' && !isNaN(rVal as number) ? f(Math.abs((rVal as number) - (cVal as number)), 4) : '—'
  console.log(`│ ${(label as string).padEnd(16)} │ ${pad(rv, 11)} │ ${pad(cv, 11)} │ ${pad(delta, 11)} │`)
}
console.log('└──────────────────┴─────────────┴─────────────┴─────────────┘')
console.log()

// ── 2. Cluster Means ──
console.log('┌──────────────────────────────────────────────────────────────────────────┐')
console.log('│ 2. CLUSTER MEANS (sorted by Emotional Engagement)                       │')
console.log('├─────────┬──────────────────────┬──────────────────────┬──────────────────┤')
console.log('│ Cluster │      R mclust        │         Carm         │      |Δ|         │')
console.log('├─────────┼──────────────────────┼──────────────────────┼──────────────────┤')

const clusterLabels = ['Low', 'Mid', 'High']
for (let k = 0; k < 3; k++) {
  const rStr = rMeans[k]!.map(v => f(v)).join(', ')
  const cStr = cMeans[k]!.map(v => f(v)).join(', ')
  const dStr = rMeans[k]!.map((v, j) => f(Math.abs(v - cMeans[k]![j]!))).join(', ')
  console.log(`│ ${(k + 1 + ' (' + clusterLabels[k] + ')').padEnd(7)} │ ${pad(rStr, 20)} │ ${pad(cStr, 20)} │ ${pad(dStr, 16)} │`)
}
console.log('├─────────┴──────────────────────┴──────────────────────┴──────────────────┤')
console.log('│ Variables: ZPRE_ENG_EMOC, ZPRE_ENG_COGN, ZPRE_ENG_COND                  │')
console.log('└──────────────────────────────────────────────────────────────────────────┘')
console.log()

// ── 3. Cluster Weights & Sizes ──
console.log('┌──────────────────────────────────────────────────────────────────────┐')
console.log('│ 3. CLUSTER WEIGHTS & SIZES                                         │')
console.log('├─────────┬────────────┬────────────┬──────────┬──────────┬───────────┤')
console.log('│ Cluster │  R weight  │ Carm weight│   R size │Carm size │   Δ size  │')
console.log('├─────────┼────────────┼────────────┼──────────┼──────────┼───────────┤')
for (let k = 0; k < 3; k++) {
  console.log(`│ ${(k + 1 + ' (' + clusterLabels[k] + ')').padEnd(7)} │ ${pad(f(rWeights[k]!), 10)} │ ${pad(f(cWeights[k]!), 10)} │ ${pad(String(rSizesSorted[k]), 8)} │ ${pad(String(cSizesSorted[k]), 8)} │ ${pad(String(Math.abs(rSizesSorted[k]! - cSizesSorted[k]!)), 9)} │`)
}
console.log(`│ Total   │ ${pad(f(rWeights.reduce((a, b) => a + b, 0)), 10)} │ ${pad(f(cWeights.reduce((a, b) => a + b, 0)), 10)} │ ${pad(String(rSizesSorted.reduce((a, b) => a + b, 0)), 8)} │ ${pad(String(cSizesSorted.reduce((a, b) => a + b, 0)), 8)} │           │`)
console.log('└─────────┴────────────┴────────────┴──────────┴──────────┴───────────┘')
console.log()

// ── 4. Entropy Detail ──
console.log('┌──────────────────────────────────────────────────────────────────────┐')
console.log('│ 4. ENTROPY DIAGNOSTICS                                             │')
console.log('├──────────────────┬─────────────┬─────────────┬─────────────────────┤')
console.log('│ Metric           │    R mclust │        Carm │ Note                │')
console.log('├──────────────────┼─────────────┼─────────────┼─────────────────────┤')
console.log(`│ Entropy (norm)   │ ${pad(f(ref.entropy, 6), 11)} │ ${pad(f(carm.diagnostics.entropy, 6), 11)} │ 1=perfect, >0.8 OK  │`)
console.log(`│ |Δ Entropy|      │             │             │ ${pad(f(Math.abs(ref.entropy - carm.diagnostics.entropy), 6), 19)} │`)
console.log('├──────────────────┴─────────────┴─────────────┴─────────────────────┤')
console.log('│ Per-Cluster Entropy (case-specific Ei averaged per cluster)        │')
console.log('├─────────┬────────────────────────┬────────────────────────┐        │')

console.log('│ Cluster │ R: mean (sd) [min,max] │ Carm: mean (sd) [min,max]      │')
console.log('├─────────┼────────────────────────┼────────────────────────────────┤')
for (let k = 0; k < 3; k++) {
  const re = rClusterEntropy[k]!
  const ce = cClusterEntropy[k]!
  const rStr = `${f(re.mean, 3)} (${f(re.sd, 3)}) [${f(re.min, 3)}, ${f(re.max, 3)}]`
  const cStr = `${f(ce.mean, 3)} (${f(ce.sd, 3)}) [${f(ce.min, 3)}, ${f(ce.max, 3)}]`
  console.log(`│ ${(k + 1 + ' (' + clusterLabels[k] + ')').padEnd(7)} │ ${rStr.padEnd(22)} │ ${cStr.padEnd(30)} │`)
}
console.log('└─────────┴────────────────────────┴────────────────────────────────┘')
console.log()

// ── 5. Average Posterior Probabilities ──
console.log('┌──────────────────────────────────────────────────────────────┐')
console.log('│ 5. AVERAGE POSTERIOR PROBABILITIES (AvePP)                  │')
console.log('│ Cutoff: > 0.7 = acceptable assignment certainty            │')
console.log('├─────────┬──────────┬──────────┬──────────┬─────────────────┤')
console.log('│ Cluster │  R AvePP │Carm AvePP│     |Δ|  │ Status          │')
console.log('├─────────┼──────────┼──────────┼──────────┼─────────────────┤')
for (let k = 0; k < 3; k++) {
  const ra = rAvePPSorted[k]!
  const ca = cAvePPSorted[k]!
  const status = ca > 0.8 ? 'Good (>0.8)' : ca > 0.7 ? 'Acceptable' : 'POOR (<0.7)'
  console.log(`│ ${(k + 1 + ' (' + clusterLabels[k] + ')').padEnd(7)} │ ${pad(f(ra), 8)} │ ${pad(f(ca), 8)} │ ${pad(f(Math.abs(ra - ca)), 8)} │ ${status.padEnd(15)} │`)
}
console.log('└─────────┴──────────┴──────────┴──────────┴─────────────────┘')
console.log()

// ── 6. Per-Variable Means Table ──
console.log('┌──────────────────────────────────────────────────────────────────────────┐')
console.log('│ 6. PER-VARIABLE CLUSTER MEANS (detailed)                                │')
console.log('├─────────────────────┬──────────┬──────────┬──────────┬──────────┬────────┤')
console.log('│ Variable            │ Cluster  │  R mean  │Carm mean │   |Δ|    │  Carm σ│')
console.log('├─────────────────────┼──────────┼──────────┼──────────┼──────────┼────────┤')
for (let j = 0; j < varNames.length; j++) {
  for (let k = 0; k < 3; k++) {
    const vLabel = k === 0 ? varNames[j]! : ''
    const rm = rMeans[k]![j]!
    const cm = cMeans[k]![j]!
    const cv = Math.sqrt(cVars[k]![j]!)
    console.log(`│ ${vLabel.padEnd(19)} │ ${(k + 1 + ' (' + clusterLabels[k] + ')').padEnd(8)} │ ${pad(f(rm), 8)} │ ${pad(f(cm), 8)} │ ${pad(f(Math.abs(rm - cm)), 8)} │ ${pad(f(cv, 3), 6)} │`)
  }
  if (j < varNames.length - 1) {
    console.log('│                     │          │          │          │          │        │')
  }
}
console.log('└─────────────────────┴──────────┴──────────┴──────────┴──────────┴────────┘')
console.log()

// ── 7. Summary ──
const maxMeanDiff = Math.max(...[0, 1, 2].flatMap(k => [0, 1, 2].map(j => Math.abs(rMeans[k]![j]! - cMeans[k]![j]!))))
const maxAvePPDiff = Math.max(...[0, 1, 2].map(k => Math.abs(rAvePPSorted[k]! - cAvePPSorted[k]!)))
const deltaLL = Math.abs(ref.logLikelihood - carm.diagnostics.logLikelihood)
const deltaE = Math.abs(ref.entropy - carm.diagnostics.entropy)

console.log('┌──────────────────────────────────────────────────────────────┐')
console.log('│ 7. SUMMARY                                                  │')
console.log('├──────────────────────────────────────────────────────────────┤')
console.log(`│ Max |Δ mean|:     ${pad(f(maxMeanDiff, 6), 10)}                            │`)
console.log(`│ Max |Δ AvePP|:    ${pad(f(maxAvePPDiff, 6), 10)}                            │`)
console.log(`│ |Δ LL|:           ${pad(f(deltaLL, 6), 10)}                            │`)
console.log(`│ |Δ Entropy|:      ${pad(f(deltaE, 6), 10)}                            │`)
console.log(`│ |Δ BIC|:          ${pad(f(Math.abs(ref.bic - carm.diagnostics.bic), 6), 10)}                            │`)
console.log(`│ |Δ ICL|:          ${pad(f(Math.abs(ref.icl - carm.diagnostics.icl), 6), 10)}                            │`)
console.log('├──────────────────────────────────────────────────────────────┤')
const pass = deltaLL < 5 && deltaE < 0.02 && maxMeanDiff < 0.1 && maxAvePPDiff < 0.03
console.log(`│ Verdict: ${pass ? '✅ PASS — numerical equivalence confirmed' : '❌ FAIL — see details above'}        │`)
console.log('│ Note: R uses conjugate prior (priorControl), Carm uses MLE  │')
console.log('└──────────────────────────────────────────────────────────────┘')
