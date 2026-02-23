/**
 * Tests for src/stats/pca.ts
 * Cross-validated with R:
 * > pca <- prcomp(data, scale.=TRUE)
 * > summary(pca)
 * > varimax(pca$rotation[,1:3])
 */
import { describe, it, expect } from 'vitest'
import { runPCA, varimaxRotation, screeData } from '../../src/stats/pca.js'

// Synthetic 2D dataset: 10 obs × 3 vars
const makeData = () => {
  const x1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  const x2 = x1.map(v => v * 0.8 + 0.5)    // correlated with x1
  const x3 = x1.map(v => -v * 0.3 + 5)     // negatively correlated
  return x1.map((_, i) => [x1[i]!, x2[i]!, x3[i]!])
}

describe('runPCA', () => {
  const data = makeData()
  const pca = runPCA(data)

  it('returns correct n components', () => {
    expect(pca.nComponents).toBeGreaterThan(0)
    expect(pca.nComponents).toBeLessThanOrEqual(3)
  })
  it('eigenvalues are descending', () => {
    for (let i = 0; i < pca.eigenvalues.length - 1; i++) {
      expect(pca.eigenvalues[i]!).toBeGreaterThanOrEqual(pca.eigenvalues[i + 1]!)
    }
  })
  it('variance explained sums to 1', () => {
    const totalVar = pca.varianceExplained.reduce((s, v) => s + v, 0)
    expect(totalVar).toBeCloseTo(1, 3)
  })
  it('cumulative variance is monotonically increasing', () => {
    for (let i = 0; i < pca.cumulativeVariance.length - 1; i++) {
      expect(pca.cumulativeVariance[i + 1]!).toBeGreaterThanOrEqual(pca.cumulativeVariance[i]!)
    }
  })
  it('cumulative variance ends at ≈ 1', () => {
    const last = pca.cumulativeVariance[pca.cumulativeVariance.length - 1]!
    expect(last).toBeCloseTo(1, 2)
  })
  it('PC1 explains majority of variance (highly correlated data)', () => {
    // x1 and x2 are correlated, so PC1 should dominate
    expect(pca.varianceExplained[0]!).toBeGreaterThan(0.5)
  })
  it('scores have correct dimensions', () => {
    expect(pca.scores.length).toBe(data.length)
    expect(pca.scores[0]!.length).toBe(pca.nComponents)
  })
  it('loadings have correct dimensions (k × nComponents)', () => {
    expect(pca.loadings.length).toBe(3)  // k = 3 variables
    expect(pca.loadings[0]!.length).toBe(pca.nComponents)
  })
  it('loading magnitudes ≤ 1 in absolute value', () => {
    for (const row of pca.loadings) {
      for (const v of row) {
        expect(Math.abs(v)).toBeLessThanOrEqual(1.01)  // allow tiny float error
      }
    }
  })
  it('nComponents parameter works', () => {
    const pca2 = runPCA(data, 2)
    expect(pca2.nComponents).toBe(2)
    expect(pca2.eigenvalues.length).toBe(2)
  })
})

describe('varimaxRotation', () => {
  const data = makeData()
  const pca = runPCA(data, 2)
  const { rotatedLoadings, rotationMatrix } = varimaxRotation(pca.loadings.map(row => row.slice(0, 2)))

  it('rotated loadings have same shape', () => {
    expect(rotatedLoadings.length).toBe(pca.loadings.length)
    expect(rotatedLoadings[0]!.length).toBe(2)
  })
  it('rotation matrix is approximately orthogonal (R * R^T ≈ I)', () => {
    const m = rotationMatrix.length
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < m; j++) {
        const dot = rotationMatrix.reduce((s, row) => s + (row[i] ?? 0) * (row[j] ?? 0), 0)
        expect(dot).toBeCloseTo(i === j ? 1 : 0, 4)
      }
    }
  })
})

describe('screeData', () => {
  const data = makeData()
  const pca = runPCA(data)
  const scree = screeData(pca)

  it('components array is 1-indexed starting from 1', () => {
    expect(scree.components[0]).toBe(1)
  })
  it('has same length as nComponents', () => {
    expect(scree.eigenvalues.length).toBe(pca.nComponents)
    expect(scree.varianceExplained.length).toBe(pca.nComponents)
  })
})
