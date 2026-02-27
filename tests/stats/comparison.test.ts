/**
 * Tests for src/stats/comparison.ts
 * Cross-validated with R: t.test(), oneway.test(), wilcox.test(), kruskal.test()
 */
import { describe, it, expect } from 'vitest'
import {
  tTestIndependent,
  tTestPaired,
  oneWayANOVA,
  mannWhitneyU,
  wilcoxonSignedRank,
  kruskalWallis,
  friedmanTest,
} from '../../src/stats/comparison.js'
import gt from '../fixtures/ground_truth.json'

describe('tTestIndependent (Welch)', () => {
  const { x1, x2 } = gt.t_test_independent
  const result = tTestIndependent(x1 as number[], x2 as number[])

  it('t statistic (R ground truth, tol=0.01)', () => {
    expect(result.statistic).toBeCloseTo(gt.t_test_independent.t, 2)
  })
  it('df (Welch-Satterthwaite, tol=0.1)', () => {
    expect(result.df as number).toBeCloseTo(gt.t_test_independent.df, 1)
  })
  it('p-value (R ground truth, tol=0.001)', () => {
    expect(result.pValue).toBeCloseTo(gt.t_test_independent.p, 3)
  })
  it('CI contains mean difference', () => {
    const diff = (x1 as number[]).reduce((s, v) => s + v, 0) / x1.length -
      (x2 as number[]).reduce((s, v) => s + v, 0) / x2.length
    expect(result.ci[0]).toBeLessThan(diff)
    expect(result.ci[1]).toBeGreaterThan(diff)
  })
  it('formatted string contains t and p', () => {
    expect(result.formatted).toMatch(/t\(/)
    expect(result.formatted).toMatch(/p/)
  })
})

describe('tTestIndependent (Student, equal var)', () => {
  const { x1, x2 } = gt.t_test_independent
  const result = tTestIndependent(x1 as number[], x2 as number[], true)

  it('df is integer n1+n2-2', () => {
    expect(result.df as number).toBeCloseTo((x1 as number[]).length + (x2 as number[]).length - 2, 8)
  })
})

describe('tTestPaired', () => {
  const { x1, x2 } = gt.t_test_paired
  const result = tTestPaired(x1 as number[], x2 as number[])

  it('t statistic (R ground truth, tol=0.01)', () => {
    expect(result.statistic).toBeCloseTo(gt.t_test_paired.t, 2)
  })
  it('df = n-1', () => {
    expect(result.df as number).toBeCloseTo(gt.t_test_paired.df, 8)
  })
  it('p-value (R ground truth, tol=0.001)', () => {
    expect(result.pValue).toBeCloseTo(gt.t_test_paired.p, 3)
  })
})

describe('oneWayANOVA', () => {
  const groups = Object.entries(gt.anova.groups).map(([label, values]) => ({
    label,
    values: values as number[],
  }))
  const result = oneWayANOVA(groups)

  it('F statistic (R ground truth, tol=0.05)', () => {
    expect(result.statistic).toBeCloseTo(gt.anova.F, 1)
  })
  it('p-value (R ground truth, tol=0.01)', () => {
    expect(result.pValue).toBeCloseTo(gt.anova.p, 2)
  })
  it('omega-squared in (0, 1)', () => {
    expect(result.effectSize.value).toBeGreaterThan(0)
    expect(result.effectSize.value).toBeLessThan(1)
  })
  it('dfBetween = k-1', () => {
    expect(result.dfBetween).toBe(groups.length - 1)
  })
  it('group stats populated', () => {
    expect(result.groups.length).toBe(groups.length)
    for (const g of result.groups) {
      expect(g.mean).toBeGreaterThan(0)
      expect(g.sd).toBeGreaterThanOrEqual(0)
    }
  })
  it('levene field present with F, p, df, homogeneous', () => {
    expect(result.levene).toBeDefined()
    expect(result.levene.statistic).toBeGreaterThanOrEqual(0)
    expect(result.levene.pValue).toBeGreaterThanOrEqual(0)
    expect(result.levene.pValue).toBeLessThanOrEqual(1)
    expect(result.levene.df[0]).toBe(groups.length - 1)
    expect(typeof result.levene.homogeneous).toBe('boolean')
  })
})

describe('mannWhitneyU', () => {
  const { x1, x2 } = gt.mann_whitney
  const result = mannWhitneyU(x1 as number[], x2 as number[])

  it('W statistic (R ground truth)', () => {
    expect(result.statistic).toBeCloseTo(gt.mann_whitney.W, 0)
  })
  it('p-value (R ground truth, tol=0.01)', () => {
    expect(result.pValue).toBeCloseTo(gt.mann_whitney.p, 2)
  })
  it('rank-biserial in (-1, 1)', () => {
    expect(result.effectSize.value).toBeGreaterThan(-1)
    expect(result.effectSize.value).toBeLessThan(1)
  })
})

describe('wilcoxonSignedRank', () => {
  // > wilcox.test(c(1,2,3,4,5), c(2,4,6,8,10), paired=TRUE)
  // V = 0, p = 0.0625
  const result = wilcoxonSignedRank([1, 2, 3, 4, 5], [2, 4, 6, 8, 10])

  it('W = 0 (all differences negative)', () => {
    expect(result.statistic).toBeCloseTo(0, 0)
  })
  it('p-value ≈ 0.0625 (R ground truth)', () => {
    expect(result.pValue).toBeCloseTo(0.0625, 2)
  })
})

describe('kruskalWallis', () => {
  const groups = (gt.kruskal_wallis.groups as number[][]).map((values, i) => ({
    label: `Group${i + 1}`,
    values,
  }))
  const result = kruskalWallis(groups)

  it('H statistic (R ground truth, tol=0.01)', () => {
    expect(result.statistic).toBeCloseTo(gt.kruskal_wallis.H, 1)
  })
  it('p-value (R ground truth, tol=0.01)', () => {
    expect(result.pValue).toBeCloseTo(gt.kruskal_wallis.p, 2)
  })
  it('df = k-1', () => {
    expect(result.df as number).toBe(groups.length - 1)
  })
})

describe('friedmanTest', () => {
  // > friedman.test(matrix(c(1,2,3, 4,5,6, 7,8,9), nrow=3))
  // chi-sq = 6, df = 2, p = 0.0498
  const data = [[1, 4, 7], [2, 5, 8], [3, 6, 9]]
  const result = friedmanTest(data)

  it('chi-sq ≈ 6 (R ground truth)', () => {
    expect(result.statistic).toBeCloseTo(6, 1)
  })
  it('df = k-1 = 2', () => {
    expect(result.df as number).toBe(2)
  })
  it('p ≈ 0.05 (R ground truth)', () => {
    expect(result.pValue).toBeCloseTo(0.0498, 2)
  })
})
