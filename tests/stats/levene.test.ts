/**
 * Tests for leveneTest (Levene's / Brown-Forsythe test for homogeneity of variances).
 *
 * Cross-validated with R:
 * > car::leveneTest(y ~ group, center = median)   # Brown-Forsythe
 * > car::leveneTest(y ~ group, center = mean)     # Original Levene
 * > # Manual equivalent:
 * > medians <- tapply(y, group, median)
 * > abs_dev <- abs(y - medians[group])
 * > anova(lm(abs_dev ~ group))
 */
import { describe, it, expect } from 'vitest'
import { leveneTest } from '../../src/stats/comparison.js'
import type { GroupData } from '../../src/core/types.js'
import ref from '../../tmp/levene-r-ref.json'

// ── Test 1: Equal variances (not significant) ──────────────────────────
// R: car::leveneTest → F(2,12) = 1.031, p = 0.386

const groups1: GroupData[] = (ref.test1_equal_var.groups as number[][]).map((v, i) => ({
  label: `G${i + 1}`,
  values: v,
}))

describe('leveneTest — equal variances (3 groups)', () => {
  const median = leveneTest(groups1, 'median')
  const mean = leveneTest(groups1, 'mean')

  it('Brown-Forsythe F matches R (tol=1e-6)', () => {
    expect(median.statistic).toBeCloseTo(ref.test1_equal_var.median.F, 6)
  })

  it('Brown-Forsythe p matches R (tol=1e-6)', () => {
    expect(median.pValue).toBeCloseTo(ref.test1_equal_var.median.p, 6)
  })

  it('Brown-Forsythe df correct', () => {
    expect(median.df).toEqual([2, 12])
  })

  it('original Levene (mean) F matches R (tol=1e-6)', () => {
    expect(mean.statistic).toBeCloseTo(ref.test1_equal_var.mean.F, 6)
  })

  it('original Levene (mean) p matches R (tol=1e-6)', () => {
    expect(mean.pValue).toBeCloseTo(ref.test1_equal_var.mean.p, 6)
  })

  it('not significant (p > .05) — variances are equal', () => {
    expect(median.pValue).toBeGreaterThan(0.05)
  })

  it('testName contains Brown-Forsythe for median center', () => {
    expect(median.testName).toMatch(/Brown-Forsythe/)
  })

  it('testName contains Levene for mean center', () => {
    expect(mean.testName).toMatch(/Levene/)
  })
})

// ── Test 2: Unequal variances (significant) ────────────────────────────
// R: F(2,12) = 7.268, p = 0.00855

const groups2: GroupData[] = (ref.test2_unequal_var.groups as number[][]).map((v, i) => ({
  label: `G${i + 1}`,
  values: v,
}))

describe('leveneTest — unequal variances (3 groups)', () => {
  const result = leveneTest(groups2)

  it('F matches R (tol=1e-6)', () => {
    expect(result.statistic).toBeCloseTo(ref.test2_unequal_var.median.F, 6)
  })

  it('p matches R (tol=1e-6)', () => {
    expect(result.pValue).toBeCloseTo(ref.test2_unequal_var.median.p, 6)
  })

  it('significant (p < .05) — variances differ', () => {
    expect(result.pValue).toBeLessThan(0.05)
  })

  it('df correct', () => {
    expect(result.df).toEqual([2, 12])
  })
})

// ── Test 3: Two groups ─────────────────────────────────────────────────
// R: F(1,8) = 8.249, p = 0.0208

const groups3: GroupData[] = (ref.test3_two_groups.groups as number[][]).map((v, i) => ({
  label: `G${i + 1}`,
  values: v,
}))

describe('leveneTest — 2 groups', () => {
  const result = leveneTest(groups3)

  it('F matches R (tol=1e-6)', () => {
    expect(result.statistic).toBeCloseTo(ref.test3_two_groups.median.F, 6)
  })

  it('p matches R (tol=1e-6)', () => {
    expect(result.pValue).toBeCloseTo(ref.test3_two_groups.median.p, 6)
  })

  it('df = [1, 8]', () => {
    expect(result.df).toEqual([1, 8])
  })
})

// ── Edge cases ─────────────────────────────────────────────────────────

describe('leveneTest — edge cases', () => {
  it('throws for fewer than 2 groups', () => {
    expect(() => leveneTest([{ label: 'A', values: [1, 2, 3] }])).toThrow('at least 2 groups')
  })

  it('zero variance groups → F = 0, p = 1', () => {
    const result = leveneTest([
      { label: 'A', values: [5, 5, 5] },
      { label: 'B', values: [10, 10, 10] },
    ])
    expect(result.statistic).toBe(0)
    expect(result.pValue).toBeCloseTo(1, 4)
  })

  it('default center is median', () => {
    const explicit = leveneTest(groups1, 'median')
    const implicit = leveneTest(groups1)
    expect(implicit.statistic).toBe(explicit.statistic)
    expect(implicit.pValue).toBe(explicit.pValue)
  })

  it('formatted string is readable', () => {
    const result = leveneTest(groups2)
    expect(result.formatted).toMatch(/Brown-Forsythe/)
    expect(result.formatted).toMatch(/F\(2, 12\)/)
    expect(result.formatted).toMatch(/p/)
  })
})
