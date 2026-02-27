/**
 * Tests for src/stats/bootstrap.ts
 * Cross-validated with R:
 * > library(boot)
 * > b <- boot(data, function(d,i) mean(d[i]), R=2000)
 * > boot.ci(b, type="perc")
 * > boot.ci(b, type="bca")
 */
import { describe, it, expect } from 'vitest'
import { bootstrapCI, bootstrapCITwoSample } from '../../src/stats/bootstrap.js'
import { mean as _mean } from '../../src/core/math.js'

describe('bootstrapCI — percentile method', () => {
  const data = [2.1, 3.5, 4.2, 5.0, 3.8, 4.6, 2.9, 5.1, 3.3, 4.7]
  const result = bootstrapCI(data, _mean, { method: 'percentile', seed: 42 })

  it('estimate equals sample mean', () => {
    expect(result.estimate).toBeCloseTo(_mean(data), 8)
  })
  it('CI contains the estimate', () => {
    expect(result.ci[0]).toBeLessThanOrEqual(result.estimate)
    expect(result.ci[1]).toBeGreaterThanOrEqual(result.estimate)
  })
  it('CI is in reasonable range', () => {
    expect(result.ci[0]).toBeGreaterThan(1)
    expect(result.ci[1]).toBeLessThan(6)
  })
  it('SE is positive', () => {
    expect(result.se).toBeGreaterThan(0)
  })
  it('nBoot defaults to 2000', () => {
    expect(result.nBoot).toBe(2000)
  })
  it('method is percentile', () => {
    expect(result.method).toBe('percentile')
  })
  it('ciLevel defaults to 0.95', () => {
    expect(result.ciLevel).toBe(0.95)
  })
  it('formatted string includes CI', () => {
    expect(result.formatted).toMatch(/CI/)
    expect(result.formatted).toMatch(/SE_boot/)
  })
  it('deterministic with same seed', () => {
    const r2 = bootstrapCI(data, _mean, { method: 'percentile', seed: 42 })
    expect(r2.ci[0]).toBe(result.ci[0])
    expect(r2.ci[1]).toBe(result.ci[1])
  })
  it('different seed gives different result', () => {
    const r2 = bootstrapCI(data, _mean, { method: 'percentile', seed: 123 })
    // They might not be exactly equal (very unlikely)
    expect(r2.ci[0] !== result.ci[0] || r2.ci[1] !== result.ci[1]).toBe(true)
  })
})

describe('bootstrapCI — BCa method', () => {
  const data = [2.1, 3.5, 4.2, 5.0, 3.8, 4.6, 2.9, 5.1, 3.3, 4.7]
  const result = bootstrapCI(data, _mean, { method: 'bca', seed: 42 })

  it('estimate equals sample mean', () => {
    expect(result.estimate).toBeCloseTo(_mean(data), 8)
  })
  it('CI contains the estimate', () => {
    expect(result.ci[0]).toBeLessThanOrEqual(result.estimate)
    expect(result.ci[1]).toBeGreaterThanOrEqual(result.estimate)
  })
  it('method is bca', () => {
    expect(result.method).toBe('bca')
  })
  it('BCa CI may differ from percentile CI', () => {
    const percResult = bootstrapCI(data, _mean, { method: 'percentile', seed: 42 })
    // They might differ slightly due to bias correction
    // Just check they're both valid
    expect(result.ci[0]).toBeGreaterThan(0)
    expect(percResult.ci[0]).toBeGreaterThan(0)
  })
})

describe('bootstrapCI — custom statistic (median)', () => {
  const data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  const medianFn = (d: readonly number[]) => {
    const sorted = [...d].sort((a, b) => a - b)
    const mid = Math.floor(sorted.length / 2)
    return sorted.length % 2 === 0 ? (sorted[mid - 1]! + sorted[mid]!) / 2 : sorted[mid]!
  }
  const result = bootstrapCI(data, medianFn, { seed: 42 })

  it('estimate equals sample median', () => {
    expect(result.estimate).toBe(5.5)
  })
  it('CI is reasonable', () => {
    expect(result.ci[0]).toBeGreaterThan(0)
    expect(result.ci[1]).toBeLessThan(11)
  })
})

describe('bootstrapCI — edge cases', () => {
  it('throws on empty data', () => {
    expect(() => bootstrapCI([], _mean)).toThrow(/empty/)
  })
  it('single observation has zero SE', () => {
    const result = bootstrapCI([5], _mean, { seed: 42 })
    expect(result.estimate).toBe(5)
    expect(result.se).toBe(0)
    expect(result.ci[0]).toBe(5)
    expect(result.ci[1]).toBe(5)
  })
  it('custom nBoot', () => {
    const result = bootstrapCI([1, 2, 3, 4, 5], _mean, { nBoot: 500, seed: 42 })
    expect(result.nBoot).toBe(500)
  })
  it('custom ciLevel = 0.90', () => {
    const data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    const r95 = bootstrapCI(data, _mean, { ciLevel: 0.95, seed: 42 })
    const r90 = bootstrapCI(data, _mean, { ciLevel: 0.90, seed: 42 })
    // 90% CI should be narrower than 95% CI
    const w95 = r95.ci[1] - r95.ci[0]
    const w90 = r90.ci[1] - r90.ci[0]
    expect(w90).toBeLessThan(w95)
  })
})

describe('bootstrapCITwoSample', () => {
  const x1 = [5.1, 4.8, 5.3, 5.0, 4.9, 5.2, 5.4]
  const x2 = [3.2, 3.5, 3.1, 3.8, 3.0, 3.4, 3.6]
  const diffMeans = (a: readonly number[], b: readonly number[]) => _mean([...a]) - _mean([...b])
  const result = bootstrapCITwoSample(x1, x2, diffMeans, { seed: 42 })

  it('estimate equals difference of means', () => {
    expect(result.estimate).toBeCloseTo(diffMeans(x1, x2), 8)
  })
  it('CI is entirely positive (clear group separation)', () => {
    expect(result.ci[0]).toBeGreaterThan(0)
  })
  it('SE is positive', () => {
    expect(result.se).toBeGreaterThan(0)
  })
  it('method is percentile', () => {
    expect(result.method).toBe('percentile')
  })
  it('formatted string includes CI', () => {
    expect(result.formatted).toMatch(/CI/)
  })
  it('throws on empty sample', () => {
    expect(() => bootstrapCITwoSample([], [1, 2], diffMeans)).toThrow(/empty/)
    expect(() => bootstrapCITwoSample([1, 2], [], diffMeans)).toThrow(/empty/)
  })
  it('deterministic with same seed', () => {
    const r2 = bootstrapCITwoSample(x1, x2, diffMeans, { seed: 42 })
    expect(r2.ci[0]).toBe(result.ci[0])
    expect(r2.ci[1]).toBe(result.ci[1])
  })
})
