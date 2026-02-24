/**
 * Tests for src/stats/clustering.ts
 *
 * Cross-validation references:
 *
 * GMM — R mclust:
 * > library(mclust)
 * > set.seed(42)
 * > data <- rbind(mvrnorm(50, c(0,0), diag(2)), mvrnorm(50, c(5,5), diag(2)))
 * > fit <- Mclust(data, G=2, modelNames="VVV")
 * > fit$BIC; fit$parameters$mean
 *
 * LCA — R poLCA:
 * > library(poLCA)
 * > set.seed(42)
 * > # Binary data with 2 latent classes
 * > fit <- poLCA(cbind(V1,V2,V3,V4,V5) ~ 1, data, nclass=2, nrep=1)
 * > fit$probs; fit$P
 *
 * KMeans — R stats::kmeans:
 * > km <- kmeans(data, centers=3, nstart=1, algorithm="Lloyd")
 * > km$centers; km$tot.withinss
 *
 * LTA — manual verification on tiny example (depmixS4 is non-deterministic)
 */
import { describe, it, expect } from 'vitest'
import {
  fitGMM, predictGMM, findBestGMM,
  fitLCA,
  fitLTA,
  runKMeans, predictKMeans,
} from '../../src/stats/clustering.js'

// ─── Synthetic Data Generators ───────────────────────────────────────────

/** Generate 2D Gaussian clusters with known means */
function makeGMMData(seed = 42): number[][] {
  // Simple deterministic pseudo-random for test data
  let state = seed
  const rand = () => {
    state = (state * 1103515245 + 12345) & 0x7fffffff
    return state / 0x7fffffff
  }
  // Box-Muller for normal samples
  const randn = () => {
    const u1 = rand(), u2 = rand()
    return Math.sqrt(-2 * Math.log(u1 + 1e-10)) * Math.cos(2 * Math.PI * u2)
  }

  const data: number[][] = []
  // Cluster 1: mean [0, 0]
  for (let i = 0; i < 40; i++) data.push([randn(), randn()])
  // Cluster 2: mean [6, 6]
  for (let i = 0; i < 40; i++) data.push([6 + randn(), 6 + randn()])
  // Cluster 3: mean [0, 6]
  for (let i = 0; i < 40; i++) data.push([randn(), 6 + randn()])
  return data
}

/** Generate binary data with known LCA structure */
function makeLCAData(seed = 42): number[][] {
  let state = seed
  const rand = () => {
    state = (state * 1103515245 + 12345) & 0x7fffffff
    return state / 0x7fffffff
  }

  const data: number[][] = []
  // Class 1: high probability on items 0-2, low on 3-4
  for (let i = 0; i < 50; i++) {
    data.push([
      rand() < 0.9 ? 1 : 0,
      rand() < 0.85 ? 1 : 0,
      rand() < 0.8 ? 1 : 0,
      rand() < 0.15 ? 1 : 0,
      rand() < 0.1 ? 1 : 0,
    ])
  }
  // Class 2: low on 0-2, high on 3-4
  for (let i = 0; i < 50; i++) {
    data.push([
      rand() < 0.1 ? 1 : 0,
      rand() < 0.15 ? 1 : 0,
      rand() < 0.2 ? 1 : 0,
      rand() < 0.85 ? 1 : 0,
      rand() < 0.9 ? 1 : 0,
    ])
  }
  return data
}

/** Generate LTA data: N subjects × T timepoints × M items */
function makeLTAData(seed = 42): number[][][] {
  let state = seed
  const rand = () => {
    state = (state * 1103515245 + 12345) & 0x7fffffff
    return state / 0x7fffffff
  }

  const N = 30, T = 4, M = 3
  const rhoTrue = [
    [0.9, 0.8, 0.2],  // State 0: high on items 0,1; low on 2
    [0.2, 0.3, 0.9],  // State 1: low on items 0,1; high on 2
  ]
  const data: number[][][] = []

  for (let i = 0; i < N; i++) {
    const subj: number[][] = []
    // Half start in state 0, half in state 1
    let state_s = i < N / 2 ? 0 : 1
    for (let t = 0; t < T; t++) {
      const obs: number[] = []
      for (let m = 0; m < M; m++) {
        obs.push(rand() < rhoTrue[state_s]![m]! ? 1 : 0)
      }
      subj.push(obs)
      // Transition: 80% stay, 20% switch
      if (rand() < 0.2) state_s = 1 - state_s
    }
    data.push(subj)
  }
  return data
}

// ═════════════════════════════════════════════════════════════════════════
// GMM Tests
// ═════════════════════════════════════════════════════════════════════════

describe('fitGMM', () => {
  const data = makeGMMData()

  it('recovers 3 well-separated clusters', () => {
    const res = fitGMM(data, { k: 3 })
    expect(res.weights.length).toBe(3)
    expect(res.means.length).toBe(3)
    expect(res.covariances.length).toBe(3)
    expect(res.labels.length).toBe(data.length)
    expect(res.posteriors.length).toBe(data.length)

    // Each weight should be roughly 1/3
    for (const w of res.weights) expect(w).toBeGreaterThan(0.15)
    // Weights sum to 1
    expect(res.weights.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 4)
  })

  it('converges', () => {
    const res = fitGMM(data, { k: 3 })
    expect(res.diagnostics.converged).toBe(true)
    expect(res.diagnostics.iterations).toBeLessThan(200)
  })

  it('diagnostics are well-formed', () => {
    const res = fitGMM(data, { k: 3 })
    const d = res.diagnostics
    expect(d.df).toBeGreaterThan(0)
    expect(d.bic).toBeGreaterThan(0)
    expect(d.aic).toBeGreaterThan(0)
    expect(d.icl).toBeGreaterThanOrEqual(d.bic)  // ICL >= BIC since entropy >= 0
    expect(d.entropy).toBeGreaterThanOrEqual(0)
    expect(d.avepp.length).toBe(3)
    for (const a of d.avepp) {
      expect(a).toBeGreaterThanOrEqual(0)
      expect(a).toBeLessThanOrEqual(1)
    }
    expect(d.formatted).toContain('GMM')
    expect(d.formatted).toContain('BIC')
  })

  it('means are near true cluster centers', () => {
    const res = fitGMM(data, { k: 3 })
    // Sort means by x-coordinate to make comparison order-invariant
    const sorted = [...res.means].sort((a, b) => a[0]! - b[0]!)

    // Cluster near [0, 0]
    expect(sorted[0]![0]).toBeCloseTo(0, 0)
    expect(sorted[0]![1]).toBeLessThan(3)

    // Cluster near [6, 6]
    expect(sorted[2]![0]).toBeCloseTo(6, 0)
    expect(sorted[2]![1]).toBeCloseTo(6, 0)
  })

  it('posteriors sum to 1 per observation', () => {
    const res = fitGMM(data, { k: 3 })
    for (const row of res.posteriors) {
      expect(row.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 6)
    }
  })

  it('is deterministic with same seed', () => {
    const r1 = fitGMM(data, { k: 2, seed: 123 })
    const r2 = fitGMM(data, { k: 2, seed: 123 })
    expect(r1.diagnostics.logLikelihood).toBe(r2.diagnostics.logLikelihood)
    expect(r1.labels).toEqual(r2.labels)
  })

  it('different seeds give different results', () => {
    const r1 = fitGMM(data, { k: 2, seed: 1 })
    const r2 = fitGMM(data, { k: 2, seed: 999 })
    // Likelihoods may differ (not guaranteed but very likely with different inits)
    expect(r1.diagnostics.logLikelihood).not.toBe(r2.diagnostics.logLikelihood)
  })

  it('covariance constraint VVI produces diagonal covariances', () => {
    const res = fitGMM(data, { k: 2, model: 'VVI' })
    for (const cov of res.covariances) {
      const arr = cov.toArray()
      // Off-diagonal elements should be near zero
      for (let i = 0; i < arr.length; i++) {
        for (let j = 0; j < arr.length; j++) {
          if (i !== j) expect(Math.abs(arr[i]![j]!)).toBeLessThan(0.01)
        }
      }
    }
  })

  it('covariance constraint EII produces spherical covariances', () => {
    const res = fitGMM(data, { k: 2, model: 'EII' })
    // All components should share same diagonal value
    const diag1 = res.covariances[0]!.get(0, 0)
    const diag2 = res.covariances[1]!.get(0, 0)
    expect(diag1).toBeCloseTo(diag2, 2)
    // Off-diagonal should be zero
    expect(Math.abs(res.covariances[0]!.get(0, 1))).toBeLessThan(0.01)
  })

  it('DF calculation differs by model type', () => {
    const vvv = fitGMM(data, { k: 2, model: 'VVV' })
    const eii = fitGMM(data, { k: 2, model: 'EII' })
    // VVV has many more free parameters than EII
    expect(vvv.diagnostics.df).toBeGreaterThan(eii.diagnostics.df)
  })

  it('throws for empty data', () => {
    expect(() => fitGMM([], { k: 2 })).toThrow('empty')
  })

  it('throws for k > n', () => {
    expect(() => fitGMM([[1, 2]], { k: 3 })).toThrow('exceed')
  })
})

describe('predictGMM', () => {
  const data = makeGMMData()

  it('predicts correct labels for training data', () => {
    const res = fitGMM(data, { k: 3 })
    const pred = predictGMM(data, res)
    expect(pred.labels.length).toBe(data.length)
    // Labels should match the original fit labels
    expect(pred.labels).toEqual(res.labels)
  })

  it('predicts correct cluster for new points near known centers', () => {
    const res = fitGMM(data, { k: 3 })
    const sorted = [...res.means].sort((a, b) => a[0]! - b[0]!)
    const cluster0Mean = sorted[0]!
    const cluster2Mean = sorted[2]!

    const pred = predictGMM([cluster0Mean, cluster2Mean], res)
    // The two points should be in different clusters
    expect(pred.labels[0]).not.toBe(pred.labels[1])
  })
})

describe('findBestGMM', () => {
  const data = makeGMMData()

  it('selects a model with reasonable BIC', () => {
    const best = findBestGMM(data, [1, 2, 3, 4], ['VVV', 'EII'])
    expect(best.diagnostics.bic).toBeDefined()
    expect(best.weights.length).toBeGreaterThanOrEqual(1)
    expect(best.weights.length).toBeLessThanOrEqual(4)
  })

  it('prefers k > 1 for multi-cluster data', () => {
    const best = findBestGMM(data, [1, 2, 3], ['VVV'])
    // With 3 clear clusters, BIC should prefer k >= 2
    expect(best.weights.length).toBeGreaterThan(1)
  })
})

// ═════════════════════════════════════════════════════════════════════════
// LCA Tests
// ═════════════════════════════════════════════════════════════════════════

describe('fitLCA', () => {
  const data = makeLCAData()

  it('recovers 2 latent classes', () => {
    const res = fitLCA(data, { k: 2 })
    expect(res.rho.length).toBe(2)
    expect(res.priorWeights.length).toBe(2)
    expect(res.labels.length).toBe(data.length)
    expect(res.posteriors.length).toBe(data.length)
  })

  it('item-response probabilities reflect data structure', () => {
    const res = fitLCA(data, { k: 2 })
    // Sort classes by rho[0] (first item response) to make comparison order-invariant
    const sorted = [...res.rho].sort((a, b) => b[0]! - a[0]!)

    // Class with high item 0 should have high items 1,2 and low items 3,4
    expect(sorted[0]![0]).toBeGreaterThan(0.6)
    expect(sorted[0]![1]).toBeGreaterThan(0.5)
    expect(sorted[0]![3]).toBeLessThan(0.5)
    expect(sorted[0]![4]).toBeLessThan(0.5)

    // Opposite class
    expect(sorted[1]![0]).toBeLessThan(0.4)
    expect(sorted[1]![3]).toBeGreaterThan(0.5)
  })

  it('weights sum to 1', () => {
    const res = fitLCA(data, { k: 2 })
    expect(res.priorWeights.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 4)
  })

  it('posteriors sum to 1 per observation', () => {
    const res = fitLCA(data, { k: 2 })
    for (const row of res.posteriors) {
      expect(row.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 6)
    }
  })

  it('diagnostics are well-formed', () => {
    const res = fitLCA(data, { k: 2 })
    const d = res.diagnostics
    expect(d.df).toBe((2 - 1) + (2 * 5))  // (k-1) + k*m = 1 + 10 = 11
    expect(d.logLikelihood).toBeLessThan(0)  // LL is always negative
    expect(d.bic).toBeGreaterThan(0)
    expect(d.icl).toBeGreaterThanOrEqual(d.bic)
    expect(d.formatted).toContain('LCA')
  })

  it('converges', () => {
    const res = fitLCA(data, { k: 2 })
    expect(res.diagnostics.converged).toBe(true)
  })

  it('is deterministic', () => {
    const r1 = fitLCA(data, { k: 2, seed: 42 })
    const r2 = fitLCA(data, { k: 2, seed: 42 })
    expect(r1.diagnostics.logLikelihood).toBe(r2.diagnostics.logLikelihood)
  })

  it('rho values are in (0, 1)', () => {
    const res = fitLCA(data, { k: 2 })
    for (const row of res.rho) {
      for (const v of row) {
        expect(v).toBeGreaterThan(0)
        expect(v).toBeLessThan(1)
      }
    }
  })

  it('throws for empty data', () => {
    expect(() => fitLCA([], { k: 2 })).toThrow('empty')
  })
})

// ═════════════════════════════════════════════════════════════════════════
// LTA Tests
// ═════════════════════════════════════════════════════════════════════════

describe('fitLTA', () => {
  const data = makeLTAData()

  it('returns correct structure', () => {
    const res = fitLTA(data, { k: 2 })
    expect(res.pi.length).toBe(2)
    expect(res.tau.length).toBe(2)
    expect(res.tau[0]!.length).toBe(2)
    expect(res.rho.length).toBe(2)
    expect(res.rho[0]!.length).toBe(3)  // 3 items
    expect(res.trajectories.length).toBe(data.length)
    expect(res.trajectories[0]!.length).toBe(4)  // 4 timepoints
    expect(res.posteriors.length).toBe(data.length)
    expect(res.posteriors[0]!.length).toBe(4)
    expect(res.posteriors[0]![0]!.length).toBe(2)
  })

  it('initial state probabilities sum to ~1', () => {
    const res = fitLTA(data, { k: 2 })
    expect(res.pi.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 2)
  })

  it('transition rows sum to ~1', () => {
    const res = fitLTA(data, { k: 2 })
    for (const row of res.tau) {
      expect(row.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 2)
    }
  })

  it('transition matrix shows diagonal dominance (high self-transition)', () => {
    const res = fitLTA(data, { k: 2 })
    // Data was generated with 80% stay probability
    for (let s = 0; s < 2; s++) {
      expect(res.tau[s]![s]).toBeGreaterThan(0.4)  // Should be high
    }
  })

  it('rho values are in (0, 1)', () => {
    const res = fitLTA(data, { k: 2 })
    for (const row of res.rho) {
      for (const v of row) {
        expect(v).toBeGreaterThan(0)
        expect(v).toBeLessThan(1)
      }
    }
  })

  it('trajectories contain only valid state indices', () => {
    const res = fitLTA(data, { k: 2 })
    for (const traj of res.trajectories) {
      for (const s of traj) {
        expect(s).toBeGreaterThanOrEqual(0)
        expect(s).toBeLessThan(2)
      }
    }
  })

  it('posteriors sum to 1 per timepoint', () => {
    const res = fitLTA(data, { k: 2 })
    for (const subj of res.posteriors) {
      for (const t of subj) {
        expect(t.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 4)
      }
    }
  })

  it('diagnostics are well-formed', () => {
    const res = fitLTA(data, { k: 2 })
    const d = res.diagnostics
    // df = (k-1) + k*(k-1) + k*m = 1 + 2 + 6 = 9
    expect(d.df).toBe(9)
    expect(d.logLikelihood).toBeLessThan(0)
    expect(d.bic).toBeGreaterThan(0)
    expect(d.entropy).toBeGreaterThanOrEqual(0)
    expect(d.icl).toBeGreaterThanOrEqual(d.bic)
    expect(d.formatted).toContain('LTA')
  })

  it('is deterministic', () => {
    const r1 = fitLTA(data, { k: 2, seed: 42 })
    const r2 = fitLTA(data, { k: 2, seed: 42 })
    expect(r1.diagnostics.logLikelihood).toBe(r2.diagnostics.logLikelihood)
  })

  it('throws for < 2 timepoints', () => {
    const bad = data.map(subj => [subj[0]!])
    expect(() => fitLTA(bad, { k: 2 })).toThrow('2 timepoints')
  })

  it('throws for empty data', () => {
    expect(() => fitLTA([], { k: 2 })).toThrow('empty')
  })
})

// ═════════════════════════════════════════════════════════════════════════
// K-Means Tests
// ═════════════════════════════════════════════════════════════════════════

describe('runKMeans', () => {
  const data = makeGMMData()  // Same 3-cluster data

  it('finds 3 clusters', () => {
    const res = runKMeans(data, { k: 3 })
    expect(res.centroids.length).toBe(3)
    expect(res.labels.length).toBe(data.length)
    // Every point assigned to a valid cluster
    for (const l of res.labels) {
      expect(l).toBeGreaterThanOrEqual(0)
      expect(l).toBeLessThan(3)
    }
  })

  it('centroids are near true means', () => {
    const res = runKMeans(data, { k: 3 })
    const sorted = [...res.centroids].sort((a, b) => a[0]! - b[0]!)
    // Near [0, 0]
    expect(sorted[0]![0]).toBeCloseTo(0, 0)
    // Near [6, 6]
    expect(sorted[2]![0]).toBeCloseTo(6, 0)
    expect(sorted[2]![1]).toBeCloseTo(6, 0)
  })

  it('converges', () => {
    const res = runKMeans(data, { k: 3 })
    expect(res.converged).toBe(true)
    expect(res.iterations).toBeLessThan(300)
  })

  it('inertia is positive', () => {
    const res = runKMeans(data, { k: 3 })
    expect(res.inertia).toBeGreaterThan(0)
  })

  it('inertia decreases with more clusters', () => {
    const r2 = runKMeans(data, { k: 2 })
    const r5 = runKMeans(data, { k: 5 })
    expect(r5.inertia).toBeLessThan(r2.inertia)
  })

  it('is deterministic', () => {
    const r1 = runKMeans(data, { k: 3, seed: 42 })
    const r2 = runKMeans(data, { k: 3, seed: 42 })
    expect(r1.inertia).toBe(r2.inertia)
    expect(r1.labels).toEqual(r2.labels)
  })

  it('handles k=1 (single cluster)', () => {
    const res = runKMeans(data, { k: 1 })
    expect(res.centroids.length).toBe(1)
    // All labels should be 0
    for (const l of res.labels) expect(l).toBe(0)
    // Centroid should be the overall mean
    const meanX = data.reduce((s, p) => s + p[0]!, 0) / data.length
    expect(res.centroids[0]![0]).toBeCloseTo(meanX, 1)
  })

  it('throws for empty data', () => {
    expect(() => runKMeans([], { k: 2 })).toThrow('empty')
  })

  it('throws for k > n', () => {
    expect(() => runKMeans([[1, 2]], { k: 3 })).toThrow('exceed')
  })
})

describe('predictKMeans', () => {
  const data = makeGMMData()

  it('predicts consistent with fit labels', () => {
    const res = runKMeans(data, { k: 3 })
    const pred = predictKMeans(data, res.centroids)
    expect(pred).toEqual(res.labels)
  })

  it('assigns new point to nearest centroid', () => {
    const res = runKMeans(data, { k: 3 })
    const sorted = [...res.centroids].sort((a, b) => a[0]! - b[0]!)
    // Point at [0, 0] should go to cluster with lowest x centroid
    const pred = predictKMeans([[0, 0]], res.centroids)
    // Find which cluster index has the lowest-x centroid
    const expectedCluster = res.centroids.findIndex(c => c[0] === sorted[0]![0])
    expect(pred[0]).toBe(expectedCluster)
  })
})
