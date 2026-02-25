/**
 * Cross-validation: Carm GMM VVI K=3 vs mclust on engagement data.
 *
 * Run: npx tsx tmp/engagement_xval.ts
 */
import { readFileSync } from 'fs'
import { fitGMM } from '../src/stats/clustering.js'

// Load R reference
const ref = JSON.parse(readFileSync('tmp/engagement_mclust_ref.json', 'utf-8'))
const data: number[][] = ref.data
const n = ref.n
const d = ref.d

console.log(`Data: ${n} observations × ${d} variables`)
console.log(`Variables: ${ref.varNames.join(', ')}`)
console.log()

// ── Run Carm: multi-seed search for best LL ──
let best: ReturnType<typeof fitGMM> | null = null
const seeds = [42, 1, 7, 13, 99, 123, 256, 500, 1000, 2024, 3141, 9999]
for (const seed of seeds) {
  try {
    const res = fitGMM(data, { k: 3, model: 'VVI', seed, tol: 1e-8, maxIter: 1000 })
    if (!best || res.diagnostics.logLikelihood > best.diagnostics.logLikelihood) {
      best = res
    }
  } catch {
    // skip failed seeds
  }
}

if (!best) {
  console.error('All seeds failed!')
  process.exit(1)
}

const res = best

// ── Compare ──
console.log('=== R mclust Reference ===')
console.log(`  LL: ${ref.logLikelihood}`)
console.log(`  BIC: ${ref.bic}`)
console.log(`  Entropy: ${ref.entropy}`)
console.log(`  AvePP: ${ref.avepp.map((a: number) => a.toFixed(4)).join(', ')}`)
console.log(`  ICL: ${ref.icl}`)
console.log(`  Means:`)
for (let k = 0; k < 3; k++) {
  console.log(`    Cluster ${k + 1}: ${ref.means[k].map((m: number) => m.toFixed(4)).join(', ')}`)
}
console.log(`  Weights: ${ref.weights.map((w: number) => w.toFixed(4)).join(', ')}`)
console.log(`  Sizes: ${[124, 171, 422].join(', ')}`)

console.log()
console.log('=== Carm Result ===')
console.log(`  LL: ${res.diagnostics.logLikelihood.toFixed(4)}`)
console.log(`  BIC: ${res.diagnostics.bic.toFixed(4)}`)
console.log(`  Entropy: ${res.diagnostics.entropy.toFixed(6)}`)
console.log(`  AvePP: ${res.diagnostics.avepp.map(a => a.toFixed(4)).join(', ')}`)
console.log(`  ICL: ${res.diagnostics.icl.toFixed(4)}`)

// Sort means by first column for comparison
const carmMeans = [...res.means].sort((a, b) => a[0]! - b[0]!)
const rMeans = [...ref.means as number[][]].sort((a, b) => a[0]! - b[0]!)
console.log(`  Means (sorted by col1):`)
for (let k = 0; k < 3; k++) {
  console.log(`    Cluster ${k + 1}: ${carmMeans[k]!.map(m => m.toFixed(4)).join(', ')}`)
}
console.log(`  Weights: ${res.weights.map(w => w.toFixed(4)).join(', ')}`)

// Cluster sizes
const sizes = new Array(3).fill(0)
for (const l of res.labels) sizes[l]++
console.log(`  Sizes: ${sizes.join(', ')}`)

console.log()
console.log('=== Differences ===')
const deltaLL = Math.abs(res.diagnostics.logLikelihood - ref.logLikelihood)
const deltaBIC = Math.abs(res.diagnostics.bic - ref.bic)
const deltaEntropy = Math.abs(res.diagnostics.entropy - ref.entropy)
const deltaICL = Math.abs(res.diagnostics.icl - ref.icl)
console.log(`  |ΔLL|: ${deltaLL.toFixed(4)}`)
console.log(`  |ΔBIC|: ${deltaBIC.toFixed(4)}`)
console.log(`  |ΔEntropy|: ${deltaEntropy.toFixed(6)}`)
console.log(`  |ΔICL|: ${deltaICL.toFixed(4)}`)

// Means comparison
console.log(`  Means (max abs diff):`)
let maxMeanDiff = 0
for (let k = 0; k < 3; k++) {
  for (let j = 0; j < d; j++) {
    const diff = Math.abs(carmMeans[k]![j]! - rMeans[k]![j]!)
    if (diff > maxMeanDiff) maxMeanDiff = diff
  }
}
console.log(`    ${maxMeanDiff.toFixed(6)}`)

// AvePP comparison
const rAvePP = [...ref.avepp as number[]].sort()
const carmAvePP = [...res.diagnostics.avepp].sort()
console.log(`  AvePP (sorted, max abs diff): ${Math.max(...carmAvePP.map((a, i) => Math.abs(a - rAvePP[i]!))).toFixed(6)}`)

// Pass/fail
console.log()
const PASS = deltaLL < 5 && deltaEntropy < 0.05 && maxMeanDiff < 0.3
console.log(PASS ? '✅ PASS — Carm matches mclust within tolerance' : '❌ FAIL — differences too large')
if (!PASS) {
  console.log('  Note: mclust uses priorControl() (conjugate prior) which Carm does not yet support.')
  console.log('  Some differences are expected with real-world data.')
}
