/**
 * Cross-validation: Carm GMM (VVI, K=3) vs R mclust on real-world engagement data.
 *
 * Dataset: School Engagement (N=717, 3 z-scored variables: Emotional, Cognitive, Behavioral)
 * Source: https://github.com/lamethods/data/raw/main/3_engSRLach/Manuscript_School%20Engagment.csv
 *
 * R reference:
 * > library(mclust)
 * > data <- read.csv2("engagement.csv")
 * > x <- data[, c("ZPRE_ENG_EMOC", "ZPRE_ENG_COGN", "ZPRE_ENG_COND")]
 * > mod <- Mclust(x, modelNames = "VVI", G = 3, prior = priorControl())
 * > probs <- mod$z
 * > E <- 1 + sum(probs * log(probs)) / (mod$n * log(mod$G))
 * > # E = 0.6939932
 * > # AvePP = 0.8524, 0.8657, 0.8564
 * > # LL = -2782.981, BIC = 5697.463
 *
 * Note: mclust uses priorControl() (conjugate prior regularization) which
 * Carm does not implement. Small differences in LL/BIC are expected.
 */
import { describe, it, expect } from 'vitest'
import { readFileSync } from 'fs'
import { join } from 'path'
import { fitGMM } from '../../src/stats/clustering.js'

interface EngagementRef {
  data: number[][]
  n: number
  d: number
  varNames: string[]
  model: string
  k: number
  weights: number[]
  means: number[][]
  logLikelihood: number
  bic: number
  df: number
  classification: number[]
  entropy: number
  caseEntropy: number[]
  avepp: number[]
  icl: number
  posteriors: number[][]
}

const refPath = join(__dirname, '..', 'engagement_mclust_ref.json')
const ref: EngagementRef = JSON.parse(readFileSync(refPath, 'utf-8'))

// Multi-seed search for best LL
function fitBest() {
  let best: ReturnType<typeof fitGMM> | null = null
  for (const seed of [42, 1, 7, 13, 99, 123, 256, 500, 1000, 2024, 3141, 9999]) {
    try {
      const res = fitGMM(ref.data, { k: 3, model: 'VVI', seed, tol: 1e-8, maxIter: 1000 })
      if (!best || res.diagnostics.logLikelihood > best.diagnostics.logLikelihood) {
        best = res
      }
    } catch {
      // skip failed seeds
    }
  }
  return best!
}

// Order-invariant comparison: sort by first mean dimension
function sortByMean0<T extends { [0]: number }>(arr: readonly T[]): T[] {
  return [...arr].sort((a, b) => a[0] - b[0])
}

describe('Engagement data: Carm GMM VVI K=3 vs R mclust', () => {
  const res = fitBest()

  it('data dimensions are correct', () => {
    expect(ref.data.length).toBe(717)
    expect(ref.data[0].length).toBe(3)
  })

  it('log-likelihood within 5 of R', () => {
    // mclust uses priorControl(); LL may differ slightly
    expect(Math.abs(res.diagnostics.logLikelihood - ref.logLikelihood)).toBeLessThan(5)
  })

  it('BIC within 10 of R', () => {
    expect(Math.abs(res.diagnostics.bic - ref.bic)).toBeLessThan(10)
  })

  it('means match R within 0.1 (order-invariant)', () => {
    const carmSorted = sortByMean0(res.means)
    const rSorted = sortByMean0(ref.means)
    for (let k = 0; k < 3; k++) {
      for (let j = 0; j < ref.d; j++) {
        expect(carmSorted[k]![j]).toBeCloseTo(rSorted[k]![j]!, 1)
      }
    }
  })

  it('weights match R within 0.02 (sorted)', () => {
    const carmW = [...res.weights].sort()
    const rW = [...ref.weights].sort()
    for (let k = 0; k < 3; k++) {
      expect(Math.abs(carmW[k]! - rW[k]!)).toBeLessThan(0.02)
    }
  })

  it('normalized entropy matches R within 0.02', () => {
    expect(Math.abs(res.diagnostics.entropy - ref.entropy)).toBeLessThan(0.02)
  })

  it('entropy is in [0, 1]', () => {
    expect(res.diagnostics.entropy).toBeGreaterThanOrEqual(0)
    expect(res.diagnostics.entropy).toBeLessThanOrEqual(1)
  })

  it('entropy is moderate (real-world data, not perfectly separated)', () => {
    // This dataset has overlapping clusters; entropy should be ~0.65–0.75
    expect(res.diagnostics.entropy).toBeGreaterThan(0.5)
    expect(res.diagnostics.entropy).toBeLessThan(0.85)
  })

  it('case-specific entropies average to overall entropy', () => {
    const Ei = res.posteriors.map(row => {
      let rowSum = 0
      for (const p of row) {
        if (p > 1e-300) rowSum += p * Math.log(p)
      }
      return 1 + rowSum / Math.log(3)
    })
    const meanEi = Ei.reduce((a, b) => a + b, 0) / Ei.length
    expect(meanEi).toBeCloseTo(res.diagnostics.entropy, 10)
  })

  it('AvePP all above 0.7 (acceptable assignment certainty)', () => {
    for (const a of res.diagnostics.avepp) {
      expect(a).toBeGreaterThan(0.7)
    }
  })

  it('AvePP matches R within 0.03 (sorted)', () => {
    const carmPP = [...res.diagnostics.avepp].sort()
    const rPP = [...ref.avepp].sort()
    for (let k = 0; k < 3; k++) {
      expect(Math.abs(carmPP[k]! - rPP[k]!)).toBeLessThan(0.03)
    }
  })

  it('ICL matches R within 15', () => {
    expect(Math.abs(res.diagnostics.icl - ref.icl)).toBeLessThan(15)
  })

  it('ICL = BIC + 2 * rawEntropy (internal consistency)', () => {
    let rawE = 0
    for (const row of res.posteriors) {
      for (const p of row) {
        if (p > 1e-300) rawE -= p * Math.log(p)
      }
    }
    expect(res.diagnostics.icl).toBeCloseTo(res.diagnostics.bic + 2 * rawE, 6)
  })

  it('cluster sizes match R approximately', () => {
    const sizes = new Array(3).fill(0) as number[]
    for (const l of res.labels) sizes[l]++
    const carmSorted = sizes.sort((a, b) => a - b)
    const rSizes = [124, 171, 422].sort((a, b) => a - b)
    // Allow ±20 per cluster
    for (let k = 0; k < 3; k++) {
      expect(Math.abs(carmSorted[k]! - rSizes[k]!)).toBeLessThan(20)
    }
  })

  it('posteriors sum to 1 per observation', () => {
    for (const row of res.posteriors) {
      expect(row.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 6)
    }
  })

  it('converges', () => {
    expect(res.diagnostics.converged).toBe(true)
  })
})
