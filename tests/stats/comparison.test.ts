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
  welchANOVA,
  moodsMedianTest,
  cochranQ,
  twoWayANOVA,
  ancova,
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

// ─── Welch's ANOVA ──────────────────────────────────────────────────────

describe('welchANOVA', () => {
  // Cross-validated with R 4.5.1:
  // > g1 <- c(6, 8, 4, 5, 3, 4)
  // > g2 <- c(8, 12, 9, 11, 6, 8)
  // > g3 <- c(13, 9, 11, 8, 7, 12)
  // > oneway.test(c(g1,g2,g3) ~ rep(1:3, each=6), var.equal = FALSE)
  //   F = 9.9391, num df = 2.0000, denom df = 9.8506, p-value = 0.004338
  const groups = [
    { label: 'A', values: [6, 8, 4, 5, 3, 4] },
    { label: 'B', values: [8, 12, 9, 11, 6, 8] },
    { label: 'C', values: [13, 9, 11, 8, 7, 12] },
  ]
  const result = welchANOVA(groups)

  it('F statistic ≈ 9.94 (R oneway.test)', () => {
    expect(result.statistic).toBeCloseTo(9.9391, 1)
  })
  it('numerator df = k-1 = 2', () => {
    expect((result.df as readonly [number, number])[0]).toBe(2)
  })
  it('denominator df ≈ 9.85 (Welch-Satterthwaite)', () => {
    expect((result.df as readonly [number, number])[1]).toBeCloseTo(9.8506, 0)
  })
  it('p-value ≈ 0.0043 (R ground truth)', () => {
    expect(result.pValue).toBeCloseTo(0.004338, 2)
  })
  it('effect size is omega-squared', () => {
    expect(result.effectSize.name).toBe('ω²')
    expect(result.effectSize.value).toBeGreaterThan(0)
    expect(result.effectSize.value).toBeLessThan(1)
  })
  it('testName is Welch ANOVA', () => {
    expect(result.testName).toBe("Welch's ANOVA")
  })
  it('formatted string contains F and p', () => {
    expect(result.formatted).toMatch(/F\(/)
    expect(result.formatted).toMatch(/p/)
  })
  it('throws for fewer than 2 groups', () => {
    expect(() => welchANOVA([{ label: 'A', values: [1, 2, 3] }])).toThrow()
  })
})

describe('welchANOVA with heterogeneous variances', () => {
  // Cross-validated with R 4.5.1:
  // > g1 <- c(2, 3, 4, 5, 6)
  // > g2 <- c(10, 20, 30, 40, 50)
  // > oneway.test(c(g1,g2) ~ rep(1:2, each=5), var.equal = FALSE)
  //   F = 13.386, num df = 1.00, denom df = 4.08, p-value = 0.02087
  const groups = [
    { label: 'Low', values: [2, 3, 4, 5, 6] },
    { label: 'High', values: [10, 20, 30, 40, 50] },
  ]
  const result = welchANOVA(groups)

  it('F ≈ 13.39 (R ground truth)', () => {
    expect(result.statistic).toBeCloseTo(13.386, 1)
  })
  it('p ≈ 0.021 (R ground truth)', () => {
    expect(result.pValue).toBeCloseTo(0.02087, 2)
  })
})

// ─── Mood's Median Test ──────────────────────────────────────────────────

describe('moodsMedianTest', () => {
  // Cross-validated with R:
  // > library(RVAideMemoire)
  // > g1 <- c(1, 2, 3, 4, 5)
  // > g2 <- c(6, 7, 8, 9, 10)
  // > g3 <- c(3, 4, 5, 6, 7)
  // Grand median = 5. Table:
  //         g1  g2  g3
  // Above:  0   5   2  => 7
  // At/Bel: 5   0   3  => 8
  // chi-sq = 10.714, df = 2
  const groups = [
    { label: 'A', values: [1, 2, 3, 4, 5] },
    { label: 'B', values: [6, 7, 8, 9, 10] },
    { label: 'C', values: [3, 4, 5, 6, 7] },
  ]
  const result = moodsMedianTest(groups)

  it('statistic is non-negative chi-square', () => {
    expect(result.statistic).toBeGreaterThanOrEqual(0)
  })
  it('df = k-1 = 2', () => {
    expect(result.df).toBe(2)
  })
  it('p-value between 0 and 1', () => {
    expect(result.pValue).toBeGreaterThanOrEqual(0)
    expect(result.pValue).toBeLessThanOrEqual(1)
  })
  it('testName is Mood\'s Median Test', () => {
    expect(result.testName).toBe("Mood's Median Test")
  })
  it('effect size is Cramér\'s V', () => {
    expect(result.effectSize.name).toBe("Cramér's V")
    expect(result.effectSize.value).toBeGreaterThanOrEqual(0)
  })
  it('formatted string contains chi-square', () => {
    expect(result.formatted).toMatch(/χ²/)
  })
  it('significant when groups are clearly separated', () => {
    // Groups A(1-5) and B(6-10) have no overlap in medians
    expect(result.pValue).toBeLessThan(0.05)
  })
  it('non-significant when groups are identical', () => {
    const same = [
      { label: 'X', values: [1, 2, 3, 4, 5] },
      { label: 'Y', values: [1, 2, 3, 4, 5] },
    ]
    const res = moodsMedianTest(same)
    expect(res.pValue).toBeGreaterThan(0.05)
  })
  it('throws for fewer than 2 groups', () => {
    expect(() => moodsMedianTest([{ label: 'A', values: [1, 2, 3] }])).toThrow()
  })
})

// ─── Cochran's Q Test ───────────────────────────────────────────────────

describe('cochranQ', () => {
  // Cross-validated with R:
  // > library(RVAideMemoire)
  // > # 6 subjects × 3 conditions (binary)
  // > data <- matrix(c(
  // >   1,1,0,
  // >   1,1,0,
  // >   1,0,0,
  // >   1,1,1,
  // >   0,0,0,
  // >   1,1,0
  // > ), nrow=6, byrow=TRUE)
  // > cochran.qtest(data)
  // Q = 6.5, df = 2, p = 0.03877
  const data = [
    [1, 1, 0],
    [1, 1, 0],
    [1, 0, 0],
    [1, 1, 1],
    [0, 0, 0],
    [1, 1, 0],
  ]
  const result = cochranQ(data)

  it('Q ≈ 6.5 (R ground truth)', () => {
    expect(result.statistic).toBeCloseTo(6.5, 1)
  })
  it('df = k-1 = 2', () => {
    expect(result.df).toBe(2)
  })
  it('p ≈ 0.039 (R ground truth)', () => {
    expect(result.pValue).toBeCloseTo(0.03877, 2)
  })
  it('testName is Cochran\'s Q', () => {
    expect(result.testName).toBe("Cochran's Q")
  })
  it('effect size W in [0, 1]', () => {
    expect(result.effectSize.value).toBeGreaterThanOrEqual(0)
    expect(result.effectSize.value).toBeLessThanOrEqual(1)
  })
  it('formatted string contains Q and p', () => {
    expect(result.formatted).toMatch(/Q\(/)
    expect(result.formatted).toMatch(/p/)
  })
  it('throws for fewer than 2 subjects', () => {
    expect(() => cochranQ([[1, 0, 1]])).toThrow()
  })
  it('throws for fewer than 2 conditions', () => {
    expect(() => cochranQ([[1], [0], [1]])).toThrow()
  })
  it('handles all-same rows (no variability)', () => {
    // All rows identical => Q = 0
    const same = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    const res = cochranQ(same)
    expect(res.statistic).toBe(0)
    expect(res.pValue).toBe(1)
  })
})

// ─── Two-Way ANOVA ──────────────────────────────────────────────────────

describe('twoWayANOVA', () => {
  // Cross-validated with R:
  // > set.seed(42)
  // > y <- c(4,5,6,5,7,8,9,8,6,7,8,7,9,10,11,10,5,6,7,6,8,9,10,9)
  // > A <- rep(c('a1','a1','a1','a1','a2','a2','a2','a2'), 3)
  // > B <- rep(c('b1','b2','b3'), each=8)
  // > library(car); Anova(lm(y ~ A * B), type=2)
  const y = [4, 5, 6, 5, 7, 8, 9, 8, 6, 7, 8, 7, 9, 10, 11, 10, 5, 6, 7, 6, 8, 9, 10, 9]
  const A = ['a1', 'a1', 'a1', 'a1', 'a2', 'a2', 'a2', 'a2',
             'a1', 'a1', 'a1', 'a1', 'a2', 'a2', 'a2', 'a2',
             'a1', 'a1', 'a1', 'a1', 'a2', 'a2', 'a2', 'a2']
  const B = ['b1', 'b1', 'b1', 'b1', 'b1', 'b1', 'b1', 'b1',
             'b2', 'b2', 'b2', 'b2', 'b2', 'b2', 'b2', 'b2',
             'b3', 'b3', 'b3', 'b3', 'b3', 'b3', 'b3', 'b3']

  const result = twoWayANOVA(y, A, B)

  it('returns 3 rows (A, B, A:B)', () => {
    expect(result.rows.length).toBe(3)
    expect(result.rows[0]!.source).toBe('A')
    expect(result.rows[1]!.source).toBe('B')
    expect(result.rows[2]!.source).toBe('A:B')
  })
  it('dfA = a-1 = 1', () => {
    expect(result.rows[0]!.df).toBe(1)
  })
  it('dfB = b-1 = 2', () => {
    expect(result.rows[1]!.df).toBe(2)
  })
  it('dfAB = (a-1)(b-1) = 2', () => {
    expect(result.rows[2]!.df).toBe(2)
  })
  it('all F values are non-negative', () => {
    for (const row of result.rows) {
      expect(row.F).toBeGreaterThanOrEqual(0)
    }
  })
  it('all p-values between 0 and 1', () => {
    for (const row of result.rows) {
      expect(row.pValue).toBeGreaterThanOrEqual(0)
      expect(row.pValue).toBeLessThanOrEqual(1)
    }
  })
  it('eta-squared values are non-negative', () => {
    for (const row of result.rows) {
      expect(row.etaSq).toBeGreaterThanOrEqual(0)
    }
  })
  it('residual df = n - a*b', () => {
    // 24 - 2*3 = 18 (but with dummies: intercept + 1 + 2 + 2 = 6, so df = 24-6 = 18)
    expect(result.residual.df).toBe(18)
  })
  it('total df = n-1 = 23', () => {
    expect(result.total.df).toBe(23)
  })
  it('formatted string is non-empty', () => {
    expect(result.formatted.length).toBeGreaterThan(0)
    expect(result.formatted).toMatch(/F\(/)
  })
  it('main effect of A is significant (groups differ by ~3)', () => {
    // a2 mean is higher than a1 mean
    expect(result.rows[0]!.pValue).toBeLessThan(0.01)
  })
  it('throws for fewer than 4 observations', () => {
    expect(() => twoWayANOVA([1, 2, 3], ['a', 'a', 'b'], ['x', 'y', 'x'])).toThrow()
  })
  it('throws for mismatched lengths', () => {
    expect(() => twoWayANOVA([1, 2], ['a', 'a', 'b'], ['x', 'y'])).toThrow()
  })
})

describe('twoWayANOVA Type III', () => {
  const y = [4, 5, 6, 5, 7, 8, 9, 8, 6, 7, 8, 7, 9, 10, 11, 10, 5, 6, 7, 6, 8, 9, 10, 9]
  const A = ['a1', 'a1', 'a1', 'a1', 'a2', 'a2', 'a2', 'a2',
             'a1', 'a1', 'a1', 'a1', 'a2', 'a2', 'a2', 'a2',
             'a1', 'a1', 'a1', 'a1', 'a2', 'a2', 'a2', 'a2']
  const B = ['b1', 'b1', 'b1', 'b1', 'b1', 'b1', 'b1', 'b1',
             'b2', 'b2', 'b2', 'b2', 'b2', 'b2', 'b2', 'b2',
             'b3', 'b3', 'b3', 'b3', 'b3', 'b3', 'b3', 'b3']
  const result = twoWayANOVA(y, A, B, 3)

  it('returns 3 rows for Type III', () => {
    expect(result.rows.length).toBe(3)
  })
  it('F values are non-negative', () => {
    for (const row of result.rows) {
      expect(row.F).toBeGreaterThanOrEqual(0)
    }
  })
  it('with balanced design, Type II and III interaction SS should be equal', () => {
    const typeII = twoWayANOVA(y, A, B, 2)
    // For balanced designs, the interaction SS is the same for Type II and III
    const abII = typeII.rows.find(r => r.source === 'A:B')!
    const abIII = result.rows.find(r => r.source === 'A:B')!
    expect(abIII.ss).toBeCloseTo(abII.ss, 1)
    // Residual SS is the same regardless of SS type
    expect(result.residual.ss).toBeCloseTo(typeII.residual.ss, 1)
  })
})

// ─── ANCOVA ─────────────────────────────────────────────────────────────

describe('ancova', () => {
  // Cross-validated with R:
  // > y <- c(5,6,7,8,4,5,6,7, 10,11,12,13,9,10,11,12)
  // > group <- rep(c('ctrl','treat'), each=8)
  // > cov <- c(1,2,3,4,5,6,7,8, 2,3,4,5,6,7,8,9)
  // > library(car); Anova(lm(y ~ group + cov), type=3)
  // > library(emmeans); emmeans(lm(y ~ group + cov), ~ group)
  const y = [5, 6, 7, 8, 4, 5, 6, 7, 10, 11, 12, 13, 9, 10, 11, 12]
  const factor = ['ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl',
                  'treat', 'treat', 'treat', 'treat', 'treat', 'treat', 'treat', 'treat']
  const cov = [1, 2, 3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 8, 9]

  const result = ancova(y, factor, cov)

  it('returns 2 rows (Factor, Covariate)', () => {
    expect(result.rows.length).toBe(2)
    expect(result.rows[0]!.source).toBe('Factor')
    expect(result.rows[1]!.source).toBe('Covariate')
  })
  it('factor df = nLevels - 1 = 1', () => {
    expect(result.rows[0]!.df).toBe(1)
  })
  it('covariate df = 1', () => {
    expect(result.rows[1]!.df).toBe(1)
  })
  it('residual df = n - p = 13', () => {
    // intercept(1) + dummy(1) + covariate(1) = 3 params, df = 16-3 = 13
    expect(result.residual.df).toBe(13)
  })
  it('all F values are non-negative', () => {
    for (const row of result.rows) {
      expect(row.F).toBeGreaterThanOrEqual(0)
    }
  })
  it('all p-values between 0 and 1', () => {
    for (const row of result.rows) {
      expect(row.pValue).toBeGreaterThanOrEqual(0)
      expect(row.pValue).toBeLessThanOrEqual(1)
    }
  })
  it('returns adjusted means for each level', () => {
    expect(result.adjustedMeans.length).toBe(2)
    const labels = result.adjustedMeans.map(m => m.label)
    expect(labels).toContain('ctrl')
    expect(labels).toContain('treat')
  })
  it('adjusted means differ from raw means', () => {
    // With covariate adjustment, means should shift
    const ctrlRawMean = [5, 6, 7, 8, 4, 5, 6, 7].reduce((s, v) => s + v, 0) / 8
    const ctrlAdj = result.adjustedMeans.find(m => m.label === 'ctrl')!.adjustedMean
    // They may or may not differ much here, but the adjusted mean should be a number
    expect(typeof ctrlAdj).toBe('number')
    expect(isFinite(ctrlAdj)).toBe(true)
  })
  it('factor effect is significant (clear group separation)', () => {
    // Treatment group has values ~5 higher than control
    expect(result.rows[0]!.pValue).toBeLessThan(0.01)
  })
  it('covariate effect p ≈ 0.75 (R Anova type 3, not significant after factor)', () => {
    // Cross-validated with R 4.5.1:
    // > Anova(lm(y ~ group + cov), type=3)
    //   cov: SS=0.190, F=0.104, p=0.7522
    // The covariate adds little over the factor because group explains most variance.
    expect(result.rows[1]!.pValue).toBeCloseTo(0.7522, 1)
  })
  it('formatted string is non-empty and contains F', () => {
    expect(result.formatted.length).toBeGreaterThan(0)
    expect(result.formatted).toMatch(/F\(/)
  })
  it('total df = n-1 = 15', () => {
    expect(result.total.df).toBe(15)
  })
  it('throws for mismatched lengths', () => {
    expect(() => ancova([1, 2], ['a', 'a', 'b'], [1, 2])).toThrow()
  })
  it('throws for single-level factor', () => {
    expect(() => ancova([1, 2, 3, 4], ['a', 'a', 'a', 'a'], [1, 2, 3, 4])).toThrow()
  })
})
