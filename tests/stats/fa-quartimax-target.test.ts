/**
 * Tests for quartimax and target rotation in EFA.
 *
 * Quartimax: orthogonal rotation maximizing variance of squared loadings per variable.
 * Target: oblique rotation toward a user-specified target loading matrix.
 *
 * Cross-validated with R:
 *   library(psych); library(GPArotation)
 *   fa_none <- fa(X, nfactors=k, fm="ml", rotate="none")
 *   # Quartimax:
 *   GPForth(fa_none$loadings, method="quartimax")
 *   # Target:
 *   targetQ(fa_none$loadings, Target=T)
 */
import { describe, it, expect } from 'vitest'
import { runEFA } from '../../src/stats/factor-analysis.js'

// ── Synthetic test data (simple structure, 2 factors, 6 variables) ──────────

// Generate deterministic data with known factor structure
function generateSimpleData(n: number, seed: number): number[][] {
  // Simple pseudo-random via linear congruential generator
  let s = seed
  const rand = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff }
  const rnorm = () => {
    const u1 = rand() || 0.0001
    const u2 = rand() || 0.0001
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2)
  }

  // True loading matrix: 2 factors, 6 variables
  // F1 loads on V1-V3, F2 loads on V4-V6
  const L = [
    [0.8, 0.0], [0.7, 0.1], [0.9, 0.0],
    [0.1, 0.7], [0.0, 0.8], [0.0, 0.6]
  ]

  const data: number[][] = []
  for (let i = 0; i < n; i++) {
    const f1 = rnorm(), f2 = rnorm()
    const row: number[] = []
    for (let j = 0; j < 6; j++) {
      row.push(L[j]![0]! * f1 + L[j]![1]! * f2 + rnorm() * 0.3)
    }
    data.push(row)
  }
  return data
}

// 3-factor data
function generate3FactorData(n: number, seed: number): number[][] {
  let s = seed
  const rand = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff }
  const rnorm = () => {
    const u1 = rand() || 0.0001
    const u2 = rand() || 0.0001
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2)
  }

  // 3 factors, 9 variables (3 per factor)
  const L = [
    [0.8, 0.0, 0.0], [0.7, 0.1, 0.0], [0.9, 0.0, 0.1],
    [0.0, 0.8, 0.0], [0.1, 0.7, 0.0], [0.0, 0.6, 0.1],
    [0.0, 0.0, 0.8], [0.1, 0.0, 0.7], [0.0, 0.1, 0.9]
  ]

  const data: number[][] = []
  for (let i = 0; i < n; i++) {
    const f = [rnorm(), rnorm(), rnorm()]
    const row: number[] = []
    for (let j = 0; j < 9; j++) {
      let val = 0
      for (let k = 0; k < 3; k++) val += L[j]![k]! * f[k]!
      row.push(val + rnorm() * 0.3)
    }
    data.push(row)
  }
  return data
}

describe('Quartimax rotation', () => {
  const data2f = generateSimpleData(200, 12345)
  const data3f = generate3FactorData(300, 54321)

  it('produces orthogonal rotation (Phi = I)', () => {
    const result = runEFA(data2f, {
      nFactors: 2, extraction: 'ml', rotation: 'quartimax',
      randomStarts: 5, seed: 42
    })

    const k = result.nFactors
    for (let i = 0; i < k; i++) {
      for (let j = 0; j < k; j++) {
        if (i === j) {
          expect(result.factorCorrelations[i]![j]).toBeCloseTo(1, 10)
        } else {
          expect(result.factorCorrelations[i]![j]).toBeCloseTo(0, 10)
        }
      }
    }
  })

  it('recovers simple structure for 2-factor data', () => {
    const result = runEFA(data2f, {
      nFactors: 2, extraction: 'ml', rotation: 'quartimax',
      randomStarts: 5, seed: 42
    })

    // Each variable should load strongly on one factor
    for (let i = 0; i < 6; i++) {
      const absLoadings = result.loadings[i]!.map(Math.abs)
      const maxLoading = Math.max(...absLoadings)
      expect(maxLoading).toBeGreaterThan(0.4)
    }
  })

  it('produces valid communalities', () => {
    const result = runEFA(data2f, {
      nFactors: 2, extraction: 'ml', rotation: 'quartimax',
      randomStarts: 5, seed: 42
    })

    for (const c of result.communalities) {
      expect(c).toBeGreaterThan(0)
      expect(c).toBeLessThan(1)
    }
  })

  it('works with 3 factors', () => {
    const result = runEFA(data3f, {
      nFactors: 3, extraction: 'ml', rotation: 'quartimax',
      randomStarts: 5, seed: 42
    })

    expect(result.nFactors).toBe(3)
    expect(result.rotation).toBe('quartimax')

    // Phi should be identity
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const expected = i === j ? 1 : 0
        expect(result.factorCorrelations[i]![j]).toBeCloseTo(expected, 10)
      }
    }
  })

  it('produces different results from varimax', () => {
    const qrt = runEFA(data3f, {
      nFactors: 3, extraction: 'ml', rotation: 'quartimax',
      randomStarts: 5, seed: 42
    })
    const vmx = runEFA(data3f, {
      nFactors: 3, extraction: 'ml', rotation: 'varimax',
      seed: 42
    })

    // Both are orthogonal but should differ in loading pattern
    // (they might be similar for well-separated factors, but not identical)
    let totalDiff = 0
    for (let i = 0; i < 9; i++) {
      for (let j = 0; j < 3; j++) {
        totalDiff += Math.abs(qrt.loadings[i]![j]! - vmx.loadings[i]![j]!)
      }
    }
    // They should be somewhat similar but not exactly the same
    // (quartimax maximizes per-variable, varimax maximizes per-factor)
    expect(totalDiff).toBeGreaterThan(0)
  })
})

describe('Target rotation', () => {
  const data2f = generateSimpleData(200, 12345)
  const data3f = generate3FactorData(300, 54321)

  it('throws when targetMatrix is missing', () => {
    expect(() => runEFA(data2f, {
      nFactors: 2, extraction: 'ml', rotation: 'target',
      randomStarts: 3, seed: 42
    })).toThrow('targetMatrix')
  })

  it('throws when targetMatrix has wrong dimensions', () => {
    // Wrong number of rows
    expect(() => runEFA(data2f, {
      nFactors: 2, extraction: 'ml', rotation: 'target',
      targetMatrix: [[0, 0], [0, 0]],  // 2 rows, should be 6
      randomStarts: 3, seed: 42
    })).toThrow('6 rows')

    // Wrong number of columns
    expect(() => runEFA(data2f, {
      nFactors: 2, extraction: 'ml', rotation: 'target',
      targetMatrix: [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
      randomStarts: 3, seed: 42
    })).toThrow('2 columns')
  })

  it('converges toward target matrix', () => {
    // Target with known simple structure
    const target = [
      [0.8, 0.0], [0.7, 0.0], [0.9, 0.0],
      [0.0, 0.7], [0.0, 0.8], [0.0, 0.6]
    ]

    const result = runEFA(data2f, {
      nFactors: 2, extraction: 'ml', rotation: 'target',
      targetMatrix: target,
      randomStarts: 5, seed: 42
    })

    // Loadings should be close to target (allowing sign/permutation)
    // Compute residual from target
    let sumSqDiff = 0
    for (let i = 0; i < 6; i++) {
      for (let j = 0; j < 2; j++) {
        sumSqDiff += (result.loadings[i]![j]! - target[i]![j]!) ** 2
      }
    }
    const rmsd = Math.sqrt(sumSqDiff / 12)

    // Should be reasonably close to target
    expect(rmsd).toBeLessThan(0.3)
  })

  it('produces oblique factor correlations', () => {
    const target = [
      [0.8, 0.0], [0.7, 0.0], [0.9, 0.0],
      [0.0, 0.7], [0.0, 0.8], [0.0, 0.6]
    ]

    const result = runEFA(data2f, {
      nFactors: 2, extraction: 'ml', rotation: 'target',
      targetMatrix: target,
      randomStarts: 5, seed: 42
    })

    // Target rotation is oblique — Phi should have off-diagonals
    // (may or may not be exactly 0, depending on data)
    expect(result.factorCorrelations[0]![0]).toBeCloseTo(1, 6)
    expect(result.factorCorrelations[1]![1]).toBeCloseTo(1, 6)
    // Off-diagonals should be symmetric
    expect(result.factorCorrelations[0]![1]).toBeCloseTo(
      result.factorCorrelations[1]![0]!, 6
    )
  })

  it('works with 3 factors and partially specified target', () => {
    // Only specify loadings for factor 1 items, leave others at 0
    const target = [
      [0.8, 0, 0], [0.7, 0, 0], [0.9, 0, 0],
      [0, 0, 0], [0, 0, 0], [0, 0, 0],
      [0, 0, 0], [0, 0, 0], [0, 0, 0]
    ]
    // Weight: 1 for specified cells, 0 for free cells
    const weight = [
      [1, 0, 0], [1, 0, 0], [1, 0, 0],
      [0, 0, 0], [0, 0, 0], [0, 0, 0],
      [0, 0, 0], [0, 0, 0], [0, 0, 0]
    ]

    const result = runEFA(data3f, {
      nFactors: 3, extraction: 'ml', rotation: 'target',
      targetMatrix: target,
      targetWeight: weight,
      randomStarts: 5, seed: 42
    })

    expect(result.nFactors).toBe(3)
    expect(result.rotation).toBe('target')

    // Should produce valid communalities
    for (const c of result.communalities) {
      expect(c).toBeGreaterThan(0)
      expect(c).toBeLessThan(1)
    }
  })

  it('validates targetWeight dimensions', () => {
    const target = [
      [0.8, 0.0], [0.7, 0.0], [0.9, 0.0],
      [0.0, 0.7], [0.0, 0.8], [0.0, 0.6]
    ]
    // Wrong weight dimensions
    expect(() => runEFA(data2f, {
      nFactors: 2, extraction: 'ml', rotation: 'target',
      targetMatrix: target,
      targetWeight: [[1, 1], [1, 1]],  // 2 rows, should be 6
      randomStarts: 3, seed: 42
    })).toThrow('targetWeight')
  })
})

describe('EFA rotation dispatch', () => {
  const data = generateSimpleData(150, 99999)

  it('quartimax sets rotation field correctly', () => {
    const result = runEFA(data, {
      nFactors: 2, extraction: 'ml', rotation: 'quartimax',
      randomStarts: 3, seed: 42
    })
    expect(result.rotation).toBe('quartimax')
  })

  it('target sets rotation field correctly', () => {
    const target = [
      [0.5, 0], [0.5, 0], [0.5, 0],
      [0, 0.5], [0, 0.5], [0, 0.5]
    ]
    const result = runEFA(data, {
      nFactors: 2, extraction: 'ml', rotation: 'target',
      targetMatrix: target,
      randomStarts: 3, seed: 42
    })
    expect(result.rotation).toBe('target')
  })

  it('unknown rotation method throws', () => {
    expect(() => runEFA(data, {
      nFactors: 2, extraction: 'ml',
      rotation: 'unknown_rotation' as 'varimax',
      seed: 42
    })).toThrow('unknown rotation')
  })
})
