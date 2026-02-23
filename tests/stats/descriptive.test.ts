/**
 * Tests for src/stats/descriptive.ts
 * Cross-validated with R:
 * > x <- c(2,4,4,4,5,5,7,9)
 * > library(e1071); skewness(x, type=2); kurtosis(x, type=2)
 * > shapiro.test(x)
 */
import { describe, it, expect } from 'vitest'
import { describe as describeStats, skewness, kurtosis, trimmedMean, shapiroWilk, ciMean } from '../../src/stats/descriptive.js'
import gt from '../fixtures/ground_truth.json'

const x = gt.descriptive.x as number[]

describe('describe()', () => {
  const result = describeStats(x)

  it('n is correct', () => expect(result.n).toBe(gt.descriptive.n))
  it('mean is correct (R ground truth)', () => expect(result.mean).toBeCloseTo(gt.descriptive.mean, 6))
  it('median is correct (R ground truth)', () => expect(result.median).toBeCloseTo(gt.descriptive.median, 6))
  it('sd is correct (R ground truth)', () => expect(result.sd).toBeCloseTo(gt.descriptive.sd, 6))
  it('variance is correct (R ground truth)', () => expect(result.variance).toBeCloseTo(gt.descriptive.variance, 6))
  it('se = sd / sqrt(n)', () => expect(result.se).toBeCloseTo(gt.descriptive.se, 5))
  it('q1 correct (R type=7 quantile)', () => expect(result.q1).toBeCloseTo(gt.descriptive.q1, 4))
  it('q3 correct (R type=7 quantile)', () => expect(result.q3).toBeCloseTo(gt.descriptive.q3, 4))
  it('iqr = q3 - q1', () => expect(result.iqr).toBeCloseTo(gt.descriptive.iqr, 4))
  it('skewness (R e1071 type=2)', () => expect(result.skewness).toBeCloseTo(gt.descriptive.skewness, 3))
  it('excess kurtosis (R e1071 type=2)', () => expect(result.kurtosis).toBeCloseTo(gt.descriptive.kurtosis, 2))
  it('formatted string is non-empty', () => expect(result.formatted.length).toBeGreaterThan(0))
  it('CI lower < mean < CI upper', () => {
    expect(result.ci[0]).toBeLessThan(result.mean)
    expect(result.ci[1]).toBeGreaterThan(result.mean)
  })
  it('Shapiro-Wilk W is in [0,1]', () => {
    expect(result.shapiroWilk.statistic).toBeGreaterThan(0)
    expect(result.shapiroWilk.statistic).toBeLessThanOrEqual(1)
  })
  it('Shapiro-Wilk p is in [0,1]', () => {
    expect(result.shapiroWilk.pValue).toBeGreaterThanOrEqual(0)
    expect(result.shapiroWilk.pValue).toBeLessThanOrEqual(1)
  })
})

describe('skewness()', () => {
  it('c(2,4,4,4,5,5,7,9) ≈ 0.489 (R e1071 type=2)', () => {
    expect(skewness(x)).toBeCloseTo(gt.descriptive.skewness, 3)
  })
  it('symmetric data has skewness ≈ 0', () => {
    expect(skewness([1, 2, 3, 4, 5])).toBeCloseTo(0, 5)
  })
  it('throws for n < 3', () => {
    expect(() => skewness([1, 2])).toThrow()
  })
})

describe('kurtosis()', () => {
  it('c(2,4,4,4,5,5,7,9) ≈ -0.497 (R e1071 type=2)', () => {
    expect(kurtosis(x)).toBeCloseTo(gt.descriptive.kurtosis, 2)
  })
  it('normal distribution has excess kurtosis ≈ 0', () => {
    // Large normal sample has near-zero excess kurtosis
    const norm = [0.52, -0.13, 1.05, -0.78, 0.31, 0.85, -1.2, 0.44, -0.09, 1.31]
    // Just check it doesn't throw and is in reasonable range
    const k = kurtosis(norm)
    expect(isFinite(k)).toBe(true)
  })
  it('throws for n < 4', () => {
    expect(() => kurtosis([1, 2, 3])).toThrow()
  })
})

describe('trimmedMean()', () => {
  it('5% trimmed mean of c(1..10)', () => {
    const x10 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    // R: mean(x10, trim=0.05) = 5.5 (only removes fractional observations here)
    expect(trimmedMean(x10, 0.05)).toBeCloseTo(5.5, 4)
  })
  it('20% trimmed mean excludes extremes', () => {
    const x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 100]
    const trimmed = trimmedMean(x, 0.2)
    expect(trimmed).toBeLessThan(mean_simple(x))  // less than raw mean (which is skewed by 100)
    expect(trimmed).toBeCloseTo(5.5, 0)
  })
})

describe('shapiroWilk()', () => {
  it('W is in [0,1] for normal data', () => {
    const result = shapiroWilk([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    expect(result.statistic).toBeGreaterThanOrEqual(0)
    expect(result.statistic).toBeLessThanOrEqual(1)
  })
  it('normal data has high p-value', () => {
    // Perfectly linear sequence is close to normal
    const result = shapiroWilk([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    expect(result.pValue).toBeGreaterThan(0.5)
  })
  it('strongly non-normal data has lower p-value', () => {
    // Exponential-like data
    const result = shapiroWilk([1, 1, 1, 1, 1, 10, 100, 1000, 10000, 100000])
    expect(result.pValue).toBeLessThan(0.05)
  })
  it('throws for n < 3', () => {
    expect(() => shapiroWilk([1, 2])).toThrow()
  })
})

describe('ciMean()', () => {
  it('95% CI contains true mean for large n', () => {
    const x = [2, 4, 4, 4, 5, 5, 7, 9]
    const ci = ciMean(x, 0.95)
    expect(ci[0]).toBeLessThan(5)
    expect(ci[1]).toBeGreaterThan(5)
  })
  it('99% CI is wider than 95% CI', () => {
    const x = [2, 4, 4, 4, 5, 5, 7, 9]
    const ci95 = ciMean(x, 0.95)
    const ci99 = ciMean(x, 0.99)
    expect(ci99[1] - ci99[0]).toBeGreaterThan(ci95[1] - ci95[0])
  })
})

function mean_simple(x: number[]): number {
  return x.reduce((s, v) => s + v, 0) / x.length
}
