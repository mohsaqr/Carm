/**
 * Detailed side-by-side report: Carm vs mclust (NO prior) on engagement data.
 * Run: npx tsx tmp/engagement_report2.ts
 */
import { readFileSync } from 'fs'
import { fitGMM } from '../src/stats/clustering.js'

const ref = JSON.parse(readFileSync('tmp/engagement_noprior_ref.json', 'utf-8'))
const data: number[][] = ref.data
const varNames: string[] = ref.varNames
const K = 3
const D = 3

// ── Fit Carm ──
let best: ReturnType<typeof fitGMM> | null = null
for (const seed of [42, 1, 7, 13, 99, 123, 256, 500, 1000, 2024, 3141, 9999]) {
  try {
    const res = fitGMM(data, { k: K, model: 'VVI', seed, tol: 1e-8, maxIter: 1000 })
    if (!best || res.diagnostics.logLikelihood > best.diagnostics.logLikelihood) best = res
  } catch { /* skip */ }
}
const carm = best!

// ── Sort clusters by first mean (ascending) → Low, Mid, High ──
const rOrder = [0, 1, 2].sort((a, b) => (ref.means as number[][])[a]![0]! - (ref.means as number[][])[b]![0]!)
const cOrder = [0, 1, 2].sort((a, b) => carm.means[a]![0]! - carm.means[b]![0]!)
const labels = ['Low', 'Mid', 'High']

const rM = rOrder.map(i => (ref.means as number[][])[i]!)
const cM = cOrder.map(i => carm.means[i]!)
const rW = rOrder.map(i => (ref.weights as number[])[i]!)
const cW = cOrder.map(i => carm.weights[i]!)
const rPC = rOrder.map(i => (ref.perCluster as any[])[i])
const rAP = rOrder.map(i => (ref.avepp as number[])[i]!)
const cAP = cOrder.map(i => carm.diagnostics.avepp[i]!)

// Carm sizes
const cSizesAll = new Array(K).fill(0) as number[]
for (const l of carm.labels) cSizesAll[l]++
const cSizes = cOrder.map(i => cSizesAll[i]!)
const rSizes = rOrder.map(i => rPC[rOrder.indexOf(i)].n) as number[]

// Carm per-cluster stats
function carmClusterStats(clusterIdx: number) {
  const idx: number[] = []
  for (let i = 0; i < data.length; i++) if (carm.labels[i] === clusterIdx) idx.push(i)
  const n = idx.length
  const means = new Array(D).fill(0) as number[]
  const sds = new Array(D).fill(0) as number[]
  for (let j = 0; j < D; j++) {
    let s = 0, s2 = 0
    for (const i of idx) { const v = data[i]![j]!; s += v; s2 += v * v }
    means[j] = s / n
    sds[j] = Math.sqrt((s2 - s * s / n) / (n - 1))
  }
  // Case entropy
  const eis = idx.map(i => {
    let rowSum = 0
    for (const p of carm.posteriors[i]!) { if (p > 1e-300) rowSum += p * Math.log(p) }
    return 1 + rowSum / Math.log(K)
  })
  const eMean = eis.reduce((a, b) => a + b, 0) / eis.length
  const eSd = Math.sqrt(eis.reduce((s, e) => s + (e - eMean) ** 2, 0) / (eis.length - 1))
  return { n, means, sds, eMean, eSd, eMin: Math.min(...eis), eMax: Math.max(...eis) }
}
const cStats = cOrder.map(i => carmClusterStats(i))

// ── Helpers ──
const f = (v: number, dp = 4) => v.toFixed(dp)
const p = (s: string, w: number) => s.padStart(w)
const pl = (s: string, w: number) => s.padEnd(w)

// ═══════════════════════════════════════════════════════════════════════
console.log()
console.log('══════════════════════════════════════════════════════════════════════════')
console.log('  GMM Cross-Validation: Carm vs R mclust (NO prior)')
console.log('  Data: School Engagement — N = 717, D = 3 (z-scored)')
console.log('  Model: VVI, K = 3')
console.log('══════════════════════════════════════════════════════════════════════════')
console.log()

// ── TABLE 1: Global Fit ──
console.log('TABLE 1. Global Fit Statistics')
console.log('─────────────────────────────────────────────────────────')
console.log(`  ${'Metric'.padEnd(18)} ${'R mclust'.padStart(12)} ${'Carm'.padStart(12)} ${'|Δ|'.padStart(12)}`)
console.log('─────────────────────────────────────────────────────────')
const fitRows: [string, number, number][] = [
  ['Log-Likelihood', ref.logLikelihood, carm.diagnostics.logLikelihood],
  ['BIC', ref.bic, carm.diagnostics.bic],
  ['AIC', ref.aic, carm.diagnostics.aic],
  ['ICL', ref.icl, carm.diagnostics.icl],
  ['DF', ref.df, carm.diagnostics.df],
  ['Entropy', ref.entropy, carm.diagnostics.entropy],
]
for (const [label, rv, cv] of fitRows) {
  const dp = label === 'DF' ? 0 : label === 'Entropy' ? 6 : 3
  console.log(`  ${pl(label, 18)} ${p(f(rv, dp), 12)} ${p(f(cv, dp), 12)} ${p(f(Math.abs(rv - cv), dp === 0 ? 0 : dp), 12)}`)
}
console.log('─────────────────────────────────────────────────────────')
console.log()

// ── TABLE 2: Cluster Sizes & Weights ──
console.log('TABLE 2. Cluster Sizes and Mixing Weights')
console.log('──────────────────────────────────────────────────────────────────')
console.log(`  ${'Cluster'.padEnd(10)} ${'R n'.padStart(6)} ${'Carm n'.padStart(8)} ${'Δn'.padStart(5)}  ${'R π'.padStart(8)} ${'Carm π'.padStart(8)} ${'|Δπ|'.padStart(8)}`)
console.log('──────────────────────────────────────────────────────────────────')
for (let k = 0; k < K; k++) {
  const rn = rPC[k].n as number
  const cn = cStats[k]!.n
  console.log(`  ${pl(`${k+1} (${labels[k]})`, 10)} ${p(String(rn), 6)} ${p(String(cn), 8)} ${p(String(cn - rn), 5)}  ${p(f(rW[k]!), 8)} ${p(f(cW[k]!), 8)} ${p(f(Math.abs(rW[k]! - cW[k]!)), 8)}`)
}
console.log(`  ${pl('Total', 10)} ${p('717', 6)} ${p('717', 8)}`)
console.log('──────────────────────────────────────────────────────────────────')
console.log()

// ── TABLE 3: Cluster Means ──
console.log('TABLE 3. Cluster Means — M (SD)')
console.log('─────────────────────────────────────────────────────────────────────────────────────────')
console.log(`  ${'Variable'.padEnd(18)} ${''.padEnd(8)} ${'Cluster 1 (Low)'.padStart(18)} ${'Cluster 2 (Mid)'.padStart(18)} ${'Cluster 3 (High)'.padStart(18)}`)
console.log('─────────────────────────────────────────────────────────────────────────────────────────')
for (let j = 0; j < D; j++) {
  // R row
  const rStrs = [0, 1, 2].map(k => {
    const m = rM[k]![j]!
    const sd = (rPC[k].sds as number[])[j]!
    return `${f(m, 3)} (${f(sd, 3)})`
  })
  console.log(`  ${pl(varNames[j]!, 18)} ${'R'.padEnd(8)} ${rStrs.map(s => p(s, 18)).join(' ')}`)

  // Carm row
  const cStrs = [0, 1, 2].map(k => {
    const m = cM[k]![j]!
    const sd = cStats[k]!.sds[j]!
    return `${f(m, 3)} (${f(sd, 3)})`
  })
  console.log(`  ${''.padEnd(18)} ${'Carm'.padEnd(8)} ${cStrs.map(s => p(s, 18)).join(' ')}`)

  // Delta row
  const dStrs = [0, 1, 2].map(k => f(Math.abs(rM[k]![j]! - cM[k]![j]!), 4))
  console.log(`  ${''.padEnd(18)} ${'|Δ|'.padEnd(8)} ${dStrs.map(s => p(s, 18)).join(' ')}`)

  if (j < D - 1) console.log()
}
console.log('─────────────────────────────────────────────────────────────────────────────────────────')
console.log()

// ── TABLE 4: Model-Estimated Variances (VVI diagonal) ──
console.log('TABLE 4. Model-Estimated Variances (VVI diagonal σ²)')
console.log('──────────────────────────────────────────────────────────────────────────')
console.log(`  ${'Variable'.padEnd(18)} ${''.padEnd(6)} ${'Cluster 1'.padStart(12)} ${'Cluster 2'.padStart(12)} ${'Cluster 3'.padStart(12)}`)
console.log('──────────────────────────────────────────────────────────────────────────')

// Carm covariances (VVI = diagonal)
const cVars = cOrder.map(ci => {
  const cov = carm.covariances[ci]!
  return Array.from({ length: D }, (_, j) => cov.get(j, j))
})

for (let j = 0; j < D; j++) {
  const rStrs = [0, 1, 2].map(k => f((rPC[k].vars as number[])[j]!))
  console.log(`  ${pl(varNames[j]!, 18)} ${'R'.padEnd(6)} ${rStrs.map(s => p(s, 12)).join(' ')}`)

  const cStrs = [0, 1, 2].map(k => f(cVars[k]![j]!))
  console.log(`  ${''.padEnd(18)} ${'Carm'.padEnd(6)} ${cStrs.map(s => p(s, 12)).join(' ')}`)

  const dStrs = [0, 1, 2].map(k => f(Math.abs((rPC[k].vars as number[])[j]! - cVars[k]![j]!)))
  console.log(`  ${''.padEnd(18)} ${'|Δ|'.padEnd(6)} ${dStrs.map(s => p(s, 12)).join(' ')}`)
  if (j < D - 1) console.log()
}
console.log('──────────────────────────────────────────────────────────────────────────')
console.log()

// ── TABLE 5: Entropy & AvePP per Cluster ──
console.log('TABLE 5. Per-Cluster Entropy and Average Posterior Probability')
console.log('─────────────────────────────────────────────────────────────────────────────────')
console.log(`  ${'Cluster'.padEnd(12)} ${''.padEnd(6)} ${'Entropy M'.padStart(11)} ${'(SD)'.padStart(9)} ${'[min'.padStart(9)} ${'max]'.padStart(8)} ${'AvePP'.padStart(9)}`)
console.log('─────────────────────────────────────────────────────────────────────────────────')
for (let k = 0; k < K; k++) {
  // R
  const re = rPC[k]
  console.log(`  ${pl(`${k+1} (${labels[k]})`, 12)} ${'R'.padEnd(6)} ${p(f(re.entropy_mean, 4), 11)} ${p(`(${f(re.entropy_sd, 4)})`, 9)} ${p(`[${f(re.entropy_min, 3)}`, 9)} ${p(`${f(re.entropy_max, 3)}]`, 8)} ${p(f(rAP[k]!, 4), 9)}`)
  // Carm
  const ce = cStats[k]!
  console.log(`  ${''.padEnd(12)} ${'Carm'.padEnd(6)} ${p(f(ce.eMean, 4), 11)} ${p(`(${f(ce.eSd, 4)})`, 9)} ${p(`[${f(ce.eMin, 3)}`, 9)} ${p(`${f(ce.eMax, 3)}]`, 8)} ${p(f(cAP[k]!, 4), 9)}`)
  // Delta
  console.log(`  ${''.padEnd(12)} ${'|Δ|'.padEnd(6)} ${p(f(Math.abs(re.entropy_mean - ce.eMean), 4), 11)} ${''.padStart(9)} ${''.padStart(9)} ${''.padStart(8)} ${p(f(Math.abs(rAP[k]! - cAP[k]!), 4), 9)}`)
  if (k < K - 1) console.log()
}
console.log('─────────────────────────────────────────────────────────────────────────────────')
console.log('  Note: Entropy scale [0, 1]; 1 = perfect separation. AvePP > 0.7 = acceptable.')
console.log()

// ── TABLE 6: Summary ──
const maxMeanDiff = Math.max(...[0, 1, 2].flatMap(k => [0, 1, 2].map(j => Math.abs(rM[k]![j]! - cM[k]![j]!))))
const maxVarDiff = Math.max(...[0, 1, 2].flatMap(k => [0, 1, 2].map(j => Math.abs((rPC[k].vars as number[])[j]! - cVars[k]![j]!))))
const maxAPDiff = Math.max(...[0, 1, 2].map(k => Math.abs(rAP[k]! - cAP[k]!)))
const dLL = Math.abs(ref.logLikelihood - carm.diagnostics.logLikelihood)
const dE = Math.abs(ref.entropy - carm.diagnostics.entropy)
const dBIC = Math.abs(ref.bic - carm.diagnostics.bic)
const maxSizeDiff = Math.max(...[0, 1, 2].map(k => Math.abs(rPC[k].n - cStats[k]!.n)))

console.log('TABLE 6. Summary of Differences')
console.log('────────────────────────────────────────────────')
console.log(`  ${'Metric'.padEnd(24)} ${'Value'.padStart(12)}`)
console.log('────────────────────────────────────────────────')
console.log(`  ${pl('|Δ Log-Likelihood|', 24)} ${p(f(dLL, 4), 12)}`)
console.log(`  ${pl('|Δ BIC|', 24)} ${p(f(dBIC, 4), 12)}`)
console.log(`  ${pl('|Δ Entropy|', 24)} ${p(f(dE, 6), 12)}`)
console.log(`  ${pl('Max |Δ mean|', 24)} ${p(f(maxMeanDiff, 6), 12)}`)
console.log(`  ${pl('Max |Δ variance|', 24)} ${p(f(maxVarDiff, 6), 12)}`)
console.log(`  ${pl('Max |Δ AvePP|', 24)} ${p(f(maxAPDiff, 6), 12)}`)
console.log(`  ${pl('Max |Δ cluster size|', 24)} ${p(String(maxSizeDiff), 12)}`)
console.log('────────────────────────────────────────────────')
console.log()

const pass = dLL < 2 && dE < 0.01 && maxMeanDiff < 0.05 && maxAPDiff < 0.02
console.log(pass
  ? '  VERDICT: ✅ Numerical equivalence confirmed (no prior, pure MLE both sides)'
  : '  VERDICT: ❌ Differences exceed tolerance — investigate')
console.log()
