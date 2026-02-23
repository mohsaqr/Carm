/**
 * Tests for src/core/math.ts
 * All expected values cross-validated with R:
 * > pnorm(1.96)  # 0.975002
 * > pt(-3, df=10) * 2  # 0.01380
 * > p.adjust(c(0.01, 0.04, 0.20, 0.01), method='bonferroni')
 */
import { describe, it, expect } from 'vitest'
import {
  normalCDF, normalQuantile,
  tDistCDF, tDistPValue, tDistQuantile,
  chiSqCDF, chiSqPValue,
  fDistCDF, fDistPValue,
  adjustPValues,
  mean, variance, sd, se, median, quantile, rank, cov,
  nelderMead,
  logGamma, gamma, incompleteBeta,
} from '../../src/core/math.js'
import gt from '../fixtures/ground_truth.json'

describe('Normal distribution', () => {
  it('normalCDF(0) = 0.5', () => {
    // A&S 7.1.26 polynomial has ~1e-9 error at z=0; 6 dp (tol=5e-7) is sufficient
    expect(normalCDF(0)).toBeCloseTo(0.5, 6)
  })
  it('normalCDF(1.96) ≈ 0.975 (R ground truth)', () => {
    expect(normalCDF(1.96)).toBeCloseTo(0.975002, 5)
  })
  it('normalCDF(-1.96) ≈ 0.025 (R ground truth)', () => {
    expect(normalCDF(-1.96)).toBeCloseTo(0.024998, 5)
  })
  it('normalCDF(3) ≈ 0.99865 (R ground truth)', () => {
    expect(normalCDF(3)).toBeCloseTo(0.99865, 4)
  })
  it('normalQuantile(0.5) = 0', () => {
    expect(normalQuantile(0.5)).toBeCloseTo(0, 5)
  })
  it('normalQuantile(0.975) ≈ 1.96', () => {
    expect(normalQuantile(0.975)).toBeCloseTo(1.96, 2)
  })
  it('normalQuantile is inverse of normalCDF', () => {
    // Acklam normalQuantile has max error ~1.15e-9; round-trip is limited by
    // the erf approximation in normalCDF (~1e-7), so 6 dp is the right bound.
    for (const p of [0.1, 0.25, 0.5, 0.75, 0.9, 0.99]) {
      expect(normalCDF(normalQuantile(p))).toBeCloseTo(p, 6)
    }
  })
  it('normalQuantile — high precision at tail regions (Acklam vs A&S improvement)', () => {
    // R: qnorm(c(0.001, 0.01, 0.025, 0.975, 0.99, 0.999))
    const cases: [number, number][] = [
      [0.001, -3.090232], [0.01, -2.326348], [0.025, -1.959964],
      [0.975,  1.959964], [0.99,  2.326348], [0.999,  3.090232],
    ]
    for (const [p, expected] of cases) {
      expect(normalQuantile(p)).toBeCloseTo(expected, 5)
    }
  })
})

describe('t-distribution', () => {
  it('tDistCDF(0, 10) = 0.5', () => {
    expect(tDistCDF(0, 10)).toBeCloseTo(0.5, 8)
  })
  it('tDistPValue(2.31, 48) ≈ 0.025 (R ground truth)', () => {
    expect(tDistPValue(2.31, 48)).toBeCloseTo(gt.math.t_pvalue[0]!.p_two_sided, 3)
  })
  it('tDistPValue(-3, 10) ≈ 0.013 (R ground truth)', () => {
    expect(tDistPValue(-3, 10)).toBeCloseTo(gt.math.t_pvalue[1]!.p_two_sided, 2)
  })
  it('tDistQuantile(0.975, 48) ≈ 2.011', () => {
    expect(tDistQuantile(0.975, 48)).toBeCloseTo(2.011, 2)
  })
  it('t critical for df=Inf approaches normal', () => {
    expect(tDistQuantile(0.975, 1000)).toBeCloseTo(1.96, 1)
  })
})

describe('Chi-square distribution', () => {
  it('chiSqCDF(0) = 0', () => {
    expect(chiSqCDF(0, 5)).toBeCloseTo(0, 8)
  })
  it('chiSqPValue(3.841, 1) ≈ 0.05 (R ground truth)', () => {
    expect(chiSqPValue(3.841, 1)).toBeCloseTo(0.05, 2)
  })
  it('chiSqPValue(9.488, 4) ≈ 0.05 (R ground truth)', () => {
    expect(chiSqPValue(9.488, 4)).toBeCloseTo(0.05, 2)
  })
})

describe('F-distribution', () => {
  it('fDistCDF(0, 2, 10) = 0', () => {
    expect(fDistCDF(0, 2, 10)).toBe(0)
  })
  it('fDistPValue(F=3.33, df1=2, df2=27) ≈ 0.05', () => {
    // R: qf(0.95, 2, 27) = 3.354
    expect(fDistPValue(3.354, 2, 27)).toBeCloseTo(0.05, 2)
  })
})

describe('p-value adjustment', () => {
  const raw = gt.math.p_adjust_bonferroni.raw as number[]
  const n = raw.length

  it('Bonferroni: multiply by n, cap at 1 (R ground truth)', () => {
    const adj = adjustPValues(raw, 'bonferroni')
    const expected = gt.math.p_adjust_bonferroni.adjusted as number[]
    adj.forEach((p, i) => expect(p).toBeCloseTo(expected[i]!, 5))
  })

  it('none: returns unchanged', () => {
    const adj = adjustPValues(raw, 'none')
    raw.forEach((p, i) => expect(adj[i]).toBe(p))
  })

  it('BH: Benjamini-Hochberg (R ground truth)', () => {
    const adj = adjustPValues(raw, 'BH')
    const expected = gt.math.p_adjust_BH.adjusted as number[]
    adj.forEach((p, i) => expect(p).toBeCloseTo(expected[i]!, 4))
  })

  it('holm: returns adjusted values ≤ 1', () => {
    const adj = adjustPValues(raw, 'holm')
    adj.forEach(p => {
      expect(p).toBeGreaterThanOrEqual(0)
      expect(p).toBeLessThanOrEqual(1)
    })
  })
})

describe('Descriptive utilities', () => {
  const x = [2, 4, 4, 4, 5, 5, 7, 9]

  it('mean([2,4,4,4,5,5,7,9]) = 5 (R ground truth)', () => {
    expect(mean(x)).toBeCloseTo(gt.descriptive.mean, 8)
  })
  it('variance([2,4,4,4,5,5,7,9]) = 4 (R ground truth)', () => {
    expect(variance(x)).toBeCloseTo(gt.descriptive.variance, 6)
  })
  it('sd([2,4,4,4,5,5,7,9]) = 2 (R ground truth)', () => {
    expect(sd(x)).toBeCloseTo(gt.descriptive.sd, 6)
  })
  it('median([2,4,4,4,5,5,7,9]) = 4.5 (R ground truth)', () => {
    expect(median(x)).toBeCloseTo(gt.descriptive.median, 6)
  })
  it('quantile with linear interpolation (R type=7)', () => {
    expect(quantile(x, 0.25)).toBeCloseTo(gt.descriptive.q1, 4)
    expect(quantile(x, 0.75)).toBeCloseTo(gt.descriptive.q3, 4)
  })

  it('throws on empty array', () => {
    expect(() => mean([])).toThrow()
    expect(() => variance([1])).toThrow()
    expect(() => median([])).toThrow()
  })

  it('mean — coercion: Number() guards against runtime string concatenation', () => {
    // TypeScript prevents this at compile time; guard ensures correct behaviour
    // if TypeScript safety is bypassed (e.g. untyped CSV data).
    // Without Number(): [1,2,'3',4] would concatenate → '334' / 4 = 83.5
    const mixed = [1, 2, '3', 4] as unknown as number[]
    expect(mean(mixed)).toBeCloseTo(2.5, 8)
  })

  it('sd — Infinity input returns NaN explicitly (not opaque arithmetic NaN)', () => {
    // mean([1, Infinity]) = Infinity; (Infinity - Infinity)^2 was silently NaN before fix
    expect(isNaN(sd([1, Infinity]))).toBe(true)
    expect(isNaN(variance([1, Infinity]))).toBe(true)
    // mean itself still propagates correctly
    expect(mean([1, Infinity])).toBe(Infinity)
  })
})

describe('rank', () => {
  it('ranks [1,2,3] correctly', () => {
    expect(rank([1, 2, 3])).toEqual([1, 2, 3])
  })
  it('average ties', () => {
    const r = rank([1, 2, 2, 3])
    expect(r[1]).toBeCloseTo(2.5, 8)
    expect(r[2]).toBeCloseTo(2.5, 8)
  })
})

describe('Nelder-Mead', () => {
  it('minimizes x² + y² to (0, 0)', () => {
    const result = nelderMead(([x, y]: readonly number[]) => (x ?? 0) ** 2 + (y ?? 0) ** 2, [3, 4])
    expect(result.x[0]).toBeCloseTo(0, 4)
    expect(result.x[1]).toBeCloseTo(0, 4)
    expect(result.fval).toBeCloseTo(0, 6)
    expect(result.converged).toBe(true)
  })

  it('minimizes Rosenbrock function', () => {
    // f(x,y) = (1-x)² + 100(y-x²)²; minimum at (1,1)
    const rosen = ([x, y]: readonly number[]) =>
      (1 - (x ?? 0)) ** 2 + 100 * ((y ?? 0) - (x ?? 0) ** 2) ** 2
    const result = nelderMead(rosen, [0, 0], { maxIter: 10000, tol: 1e-10 })
    expect(result.x[0]).toBeCloseTo(1, 2)
    expect(result.x[1]).toBeCloseTo(1, 2)
  })
})

describe('logGamma', () => {
  it('logGamma(1) = 0', () => expect(logGamma(1)).toBeCloseTo(0, 10))
  it('logGamma(2) = 0', () => expect(logGamma(2)).toBeCloseTo(0, 10))
  it('logGamma(5) = log(24)', () => expect(logGamma(5)).toBeCloseTo(Math.log(24), 8))
  it('gamma(0.5) = sqrt(π)', () => expect(gamma(0.5)).toBeCloseTo(Math.sqrt(Math.PI), 6))
})

describe('incompleteBeta', () => {
  it('I_0(a,b) = 0', () => expect(incompleteBeta(0, 2, 3)).toBe(0))
  it('I_1(a,b) = 1', () => expect(incompleteBeta(1, 2, 3)).toBe(1))
  it('I_0.5(2,2) = 0.5 (symmetric)', () => expect(incompleteBeta(0.5, 2, 2)).toBeCloseTo(0.5, 6))
})
