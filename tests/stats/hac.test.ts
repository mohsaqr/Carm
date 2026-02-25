/**
 * Tests for Hierarchical Agglomerative Clustering (HAC).
 *
 * Cross-validated with R:
 * > set.seed(42)
 * > data <- rbind(matrix(rnorm(20,0,0.5),ncol=2), matrix(rnorm(20,5,0.5),ncol=2), matrix(rnorm(20,10,0.5),ncol=2))
 * > hc <- hclust(dist(data), method="ward.D2")
 * > hc$height; hc$order; cutree(hc, k=3)
 * > cor(cophenetic(hc), dist(data))
 *
 * R merge convention: negative = original obs (1-indexed), positive = merged cluster (1-indexed)
 * Carm merge convention: a,b are cluster IDs (0..n-1 = obs, n+ = intermediate)
 */

import { describe, it, expect } from 'vitest'
import { runHierarchical, cutTree, cutTreeHeight, silhouetteScores } from '../../src/stats/clustering.js'

// Same synthetic data as DBSCAN tests: 3 clusters at (0,0), (5,5), (10,10)
const DATA: readonly (readonly number[])[] = [
  [0.68547922357, 0.65243482711],
  [-0.2823490857, 1.1433226964],
  [0.18156420567, -0.69443035056],
  [0.31643130248, -0.13939438341],
  [0.20213416157, -0.066660668197],
  [-0.053062258046, 0.31797519904],
  [0.75576099872, -0.14212646071],
  [-0.047329519207, -1.3282277105],
  [1.0092118569, -1.2202334643],
  [-0.031357049526, 0.66005667287],
  [4.846680703, 5.2277250616],
  [4.109345783, 5.3524186686],
  [4.9140413221, 5.517551761],
  [5.6073373496, 4.6955368123],
  [5.9475967306, 5.2524775616],
  [4.7847654342, 4.1414956605],
  [4.8713653086, 4.6077704958],
  [4.1184184574, 4.5745462029],
  [5.2300486774, 3.792896175],
  [4.680002562, 5.0180613034],
  [10.1029993001, 10.1609626326],
  [9.8194713507, 9.6080805296],
  [10.3790816178, 10.7878637599],
  [9.6366475865, 10.3214496529],
  [9.3158594778, 10.0448803233],
  [10.2164090129, 10.1382753736],
  [9.5943034119, 10.339644408],
  [10.7220506309, 10.0449164433],
  [9.7842768987, 8.5034549584],
  [10.3278239417, 10.1424414768],
]

// R reference merge heights (ward.D2) — from tests/r_hac_reference.json
const R_WARD_HEIGHTS = [
  0.04608772325, 0.11149279244, 0.1354770451, 0.1967236742,
  0.26784416021, 0.3427693843, 0.47424862712, 0.47801705387,
  0.48780568536, 0.57505871701, 0.63262982661, 0.65265568616,
  0.67386306524, 0.77792537323, 0.80480580139, 0.81650344172,
  0.8439843547, 0.85361382696, 0.99153416982, 1.114260972,
  1.173946936, 1.62943858, 1.6409056958, 1.7624021909,
  2.2123971257, 2.2654588279, 2.9423429285, 21.3336516599,
  38.8346507139,
]

// R reference merge heights (single linkage) — from tests/r_hac_reference.json
const R_SINGLE_HEIGHTS = [
  0.04608772325, 0.11149279244, 0.11565671051, 0.1354770451,
  0.26784416021, 0.29755162355, 0.3427693843, 0.40548352629,
  0.40611059405, 0.4393381912, 0.45272314676, 0.4615950205,
  0.47424862712, 0.49319367847, 0.54455766607, 0.56550753128,
  0.57118653576, 0.62134267352, 0.64745445626, 0.65265568616,
  0.66139549974, 0.67386306524, 0.71559929497, 0.71687679204,
  0.7411867318, 0.98054553646, 1.1051860938, 5.0288138508,
  5.2122959842,
]

// R reference merge heights (complete linkage) — from tests/r_hac_reference.json
const R_COMPLETE_HEIGHTS = [
  0.04608772325, 0.11149279244, 0.1354770451, 0.22558624223,
  0.26784416021, 0.3427693843, 0.42355118311, 0.47424862712,
  0.55160208339, 0.55874659789, 0.62983431803, 0.65265568616,
  0.67386306524, 0.76543210604, 0.77792537323, 0.81074457315,
  0.81828996031, 0.89032236731, 1.062046344, 1.0852016109,
  1.2338051813, 1.4061911535, 1.4324072586, 1.6264233867,
  1.9507649993, 2.3605754015, 2.6934229881, 8.8204998879,
  15.984671461,
]

// R reference merge heights (average linkage) — from tests/r_hac_reference.json
const R_AVERAGE_HEIGHTS = [
  0.04608772325, 0.11149279244, 0.1354770451, 0.17062147637,
  0.26784416021, 0.3427693843, 0.4145173547, 0.42457685347,
  0.47424862712, 0.49904239454, 0.51671098246, 0.65265568616,
  0.66533327299, 0.70058098709, 0.7050981864, 0.72283449586,
  0.72791494929, 0.74355475312, 0.87094099204, 0.87923657992,
  0.92639949304, 1.062046344, 1.0869911624, 1.146541102,
  1.2852835384, 1.6410383525, 1.7396302781, 6.8050261704,
  10.6652765774,
]

describe('runHierarchical — ward linkage', () => {
  const hac = runHierarchical(DATA, { linkage: 'ward' })

  it('produces n-1 merges', () => {
    expect(hac.merges.length).toBe(29)
    expect(hac.heights.length).toBe(29)
  })

  it('heights are monotonically non-decreasing', () => {
    for (let i = 1; i < hac.heights.length; i++) {
      expect(hac.heights[i]!).toBeGreaterThanOrEqual(hac.heights[i - 1]! - 1e-10)
    }
  })

  it('merge heights match R hclust(method="ward.D2") to 1e-4', () => {
    for (let i = 0; i < R_WARD_HEIGHTS.length; i++) {
      expect(hac.heights[i]).toBeCloseTo(R_WARD_HEIGHTS[i]!, 4)
    }
  })

  it('first merge height matches R exactly', () => {
    // R: 0.04608772325
    expect(hac.heights[0]).toBeCloseTo(0.04608772325, 8)
  })

  it('last merge height matches R (final join)', () => {
    // R: 38.8346507139
    expect(hac.heights[28]).toBeCloseTo(38.8346507139, 3)
  })

  it('cophenetic correlation matches R', () => {
    // R: cor(cophenetic(hc), dist(data)) = 0.8656
    expect(hac.copheneticCorrelation).toBeCloseTo(0.86556235091, 3)
  })

  it('leaf order has all n observations', () => {
    expect(hac.order.length).toBe(30)
    const sorted = [...hac.order].sort((a, b) => a - b)
    for (let i = 0; i < 30; i++) expect(sorted[i]).toBe(i)
  })

  it('formatted string includes method and cophenetic', () => {
    expect(hac.formatted).toContain('ward')
    expect(hac.formatted).toContain('30')
  })
})

describe('runHierarchical — single linkage', () => {
  const hac = runHierarchical(DATA, { linkage: 'single' })

  it('heights match R hclust(method="single") to 1e-4', () => {
    for (let i = 0; i < R_SINGLE_HEIGHTS.length; i++) {
      expect(hac.heights[i]).toBeCloseTo(R_SINGLE_HEIGHTS[i]!, 4)
    }
  })

  it('first merge is same as ward (smallest distance)', () => {
    expect(hac.heights[0]).toBeCloseTo(R_WARD_HEIGHTS[0]!, 8)
  })
})

describe('runHierarchical — complete linkage', () => {
  const hac = runHierarchical(DATA, { linkage: 'complete' })

  it('heights match R hclust(method="complete") to 1e-4', () => {
    for (let i = 0; i < R_COMPLETE_HEIGHTS.length; i++) {
      expect(hac.heights[i]).toBeCloseTo(R_COMPLETE_HEIGHTS[i]!, 4)
    }
  })
})

describe('runHierarchical — average linkage', () => {
  const hac = runHierarchical(DATA, { linkage: 'average' })

  it('heights match R hclust(method="average") to 1e-4', () => {
    for (let i = 0; i < R_AVERAGE_HEIGHTS.length; i++) {
      expect(hac.heights[i]).toBeCloseTo(R_AVERAGE_HEIGHTS[i]!, 4)
    }
  })
})

describe('cutTree', () => {
  const hac = runHierarchical(DATA, { linkage: 'ward' })

  it('k=3 produces correct cluster partition (matches R)', () => {
    // R: cutree(hc, k=3) → [1,1,...,2,2,...,3,3,...] (1-indexed)
    // Our labels are 0-indexed but should partition the same way
    const labels = cutTree(hac, 3)
    expect(labels.length).toBe(30)

    // Verify same partitioning: points 0-9 same cluster, 10-19 same, 20-29 same
    const c1 = labels[0]
    const c2 = labels[10]
    const c3 = labels[20]
    expect(new Set([c1, c2, c3]).size).toBe(3) // 3 distinct clusters

    for (let i = 0; i < 10; i++) expect(labels[i]).toBe(c1)
    for (let i = 10; i < 20; i++) expect(labels[i]).toBe(c2)
    for (let i = 20; i < 30; i++) expect(labels[i]).toBe(c3)
  })

  it('k=1 puts everything in one cluster', () => {
    const labels = cutTree(hac, 1)
    expect(new Set(labels).size).toBe(1)
  })

  it('k=30 gives each point its own cluster', () => {
    const labels = cutTree(hac, 30)
    expect(new Set(labels).size).toBe(30)
  })

  it('k=2 merges two of the three clusters', () => {
    const labels = cutTree(hac, 2)
    expect(new Set(labels).size).toBe(2)
  })

  it('throws on invalid k', () => {
    expect(() => cutTree(hac, 0)).toThrow()
    expect(() => cutTree(hac, 31)).toThrow()
  })

  it('silhouette of k=3 cut matches R', () => {
    const labels = cutTree(hac, 3)
    const sil = silhouetteScores(DATA, labels as number[])
    // R: mean silhouette = 0.8466298
    expect(sil.mean).toBeCloseTo(0.8466298, 4)
  })
})

describe('cutTreeHeight', () => {
  const hac = runHierarchical(DATA, { linkage: 'ward' })

  it('very large height → 1 cluster', () => {
    const labels = cutTreeHeight(hac, 1000)
    expect(new Set(labels).size).toBe(1)
  })

  it('very small height → n clusters', () => {
    const labels = cutTreeHeight(hac, 0)
    expect(new Set(labels).size).toBe(30)
  })

  it('height just above last merge → 1 cluster', () => {
    const labels = cutTreeHeight(hac, hac.heights[28]! + 0.001)
    expect(new Set(labels).size).toBe(1)
  })

  it('height between 2nd and 3rd merges → ~28 clusters', () => {
    // Between merge 1 (height ~0.111) and merge 2 (height ~0.135)
    const labels = cutTreeHeight(hac, 0.12)
    const nClusters = new Set(labels).size
    expect(nClusters).toBe(28) // 2 merges done → 30 - 2 = 28 clusters
  })

  it('height between cluster-level merges → 3 clusters', () => {
    // Ward heights: merge 26 = 2.94, merge 27 = 21.33, merge 28 = 38.83
    // Cutting at h=30.0 should give 2 clusters (merges 0-27 done, 28 not)
    const labels = cutTreeHeight(hac, 30.0)
    expect(new Set(labels).size).toBe(2)

    // Cutting at h=10.0 should give 3 clusters (merges 0-26 done, 27-28 not)
    const labels3 = cutTreeHeight(hac, 10.0)
    expect(new Set(labels3).size).toBe(3)
  })
})

describe('edge cases', () => {
  it('2 points — produces 1 merge', () => {
    const data = [[0, 0], [1, 1]]
    const hac = runHierarchical(data)
    expect(hac.merges.length).toBe(1)
    expect(hac.heights[0]).toBeCloseTo(Math.sqrt(2), 10)
  })

  it('3 points — correct number of merges', () => {
    const data = [[0, 0], [1, 0], [10, 0]]
    const hac = runHierarchical(data, { linkage: 'single' })
    expect(hac.merges.length).toBe(2)
    expect(hac.heights[0]).toBeCloseTo(1, 10) // closest pair: 0 and 1
    expect(hac.heights[1]).toBeCloseTo(9, 10) // remaining: (0,1) and 2, single linkage uses min dist
  })

  it('identical points', () => {
    const data = [[1, 1], [1, 1], [1, 1]]
    const hac = runHierarchical(data, { linkage: 'ward' })
    expect(hac.merges.length).toBe(2)
    expect(hac.heights[0]).toBeCloseTo(0, 10)
    expect(hac.heights[1]).toBeCloseTo(0, 10)
  })

  it('throws on single point', () => {
    expect(() => runHierarchical([[1, 1]])).toThrow()
  })

  it('default linkage is ward', () => {
    const hac = runHierarchical(DATA)
    expect(hac.formatted).toContain('ward')
  })

  it('all four linkages produce same k=3 partition for well-separated data', () => {
    // R confirms: all four linkages give identical k=3 labels for this data
    for (const linkage of ['single', 'complete', 'average', 'ward'] as const) {
      const hac = runHierarchical(DATA, { linkage })
      const labels = cutTree(hac, 3)
      // Same partition: 0-9 together, 10-19 together, 20-29 together
      const c1 = labels[0], c2 = labels[10], c3 = labels[20]
      expect(new Set([c1, c2, c3]).size).toBe(3)
      for (let i = 0; i < 10; i++) expect(labels[i]).toBe(c1)
      for (let i = 10; i < 20; i++) expect(labels[i]).toBe(c2)
      for (let i = 20; i < 30; i++) expect(labels[i]).toBe(c3)
    }
  })
})
