/**
 * Tests for repeatedMeasuresANOVA, mauchlysTest, epsilonCorrections.
 *
 * Cross-validated with R:
 * > fit <- aov(y ~ cond + Error(subj/cond), data=df)
 * > summary(fit)
 *
 * Mauchly's test and epsilon corrections verified against manual R
 * computation using Helmert orthonormal contrasts and eigen decomposition.
 */
import { describe, it, expect } from 'vitest'
import {
  repeatedMeasuresANOVA,
  mauchlysTest,
  epsilonCorrections,
  tTestPaired,
} from '../../src/stats/comparison.js'

// ── Test Case 1: Simple 3-condition, 5 subjects ────────────────────────
// R reference:
//   y <- c(3,2,1,4,3, 5,4,6,5,7, 7,8,6,9,8)
//   subj <- factor(rep(1:5, 3))
//   cond <- factor(rep(c("C1","C2","C3"), each=5))
//   df <- data.frame(y=y, subj=subj, cond=cond)
//   fit <- aov(y ~ cond + Error(subj/cond), data=df)
//   summary(fit)
//
// Error: subj:cond
//   cond: Df=2, SS=62.8, MS=31.4, F=29.4375, p=0.000205
//   Residuals: Df=8, SS=8.533
//
// Mauchly W = 0.3955, chi-sq = 2.783, df=2, p = 0.249
// GG epsilon = 0.623, HF epsilon = 0.769

const data3cond = [
  [3, 5, 7],
  [2, 4, 8],
  [1, 6, 6],
  [4, 5, 9],
  [3, 7, 8],
]

describe('repeatedMeasuresANOVA — 3 conditions, 5 subjects', () => {
  const result = repeatedMeasuresANOVA(data3cond, { labels: ['C1', 'C2', 'C3'] })

  it('F statistic matches R (tol=0.01)', () => {
    expect(result.statistic).toBeCloseTo(29.4375, 2)
  })

  it('p-value matches R (tol=0.001)', () => {
    expect(result.pValue).toBeCloseTo(0.000205, 3)
  })

  it('SS decomposition matches R', () => {
    expect(result.ssConditions).toBeCloseTo(62.8, 4)
    expect(result.ssSubjects).toBeCloseTo(7.0667, 2)
    expect(result.ssError).toBeCloseTo(8.5333, 2)
    expect(result.ssTotal).toBeCloseTo(78.4, 4)
  })

  it('SS components sum to SS_total', () => {
    expect(result.ssConditions + result.ssSubjects + result.ssError)
      .toBeCloseTo(result.ssTotal, 8)
  })

  it('degrees of freedom correct', () => {
    expect(result.dfConditions).toBe(2)
    expect(result.dfSubjects).toBe(4)
    expect(result.dfError).toBe(8)
  })

  it('mean squares correct', () => {
    expect(result.msConditions).toBeCloseTo(31.4, 4)
    expect(result.msError).toBeCloseTo(1.0667, 2)
  })

  it('partial eta-squared matches R (tol=0.001)', () => {
    expect(result.effectSize.value).toBeCloseTo(0.8804, 3)
    expect(result.effectSize.name).toBe('η²_p')
  })

  it('sphericity test matches R', () => {
    expect(result.sphericity).not.toBeNull()
    expect(result.sphericity!.W).toBeCloseTo(0.3955, 2)
    expect(result.sphericity!.chiSq).toBeCloseTo(2.783, 1)
    expect(result.sphericity!.df).toBe(2)
    expect(result.sphericity!.pValue).toBeCloseTo(0.249, 2)
  })

  it('epsilon corrections match R', () => {
    expect(result.epsilonGG).toBeCloseTo(0.623, 2)
    expect(result.epsilonHF).toBeCloseTo(0.769, 2)
  })

  it('sphericity not violated (p=0.249 > 0.05), no correction', () => {
    expect(result.correction).toBe('none')
    expect(result.correctedDf).toBeUndefined()
  })

  it('condition descriptives populated', () => {
    expect(result.conditions.length).toBe(3)
    expect(result.conditions[0]!.label).toBe('C1')
    expect(result.conditions[0]!.mean).toBeCloseTo(2.6, 4)
    expect(result.conditions[1]!.mean).toBeCloseTo(5.4, 4)
    expect(result.conditions[2]!.mean).toBeCloseTo(7.6, 4)
    for (const c of result.conditions) {
      expect(c.n).toBe(5)
      expect(c.sd).toBeGreaterThan(0)
    }
  })

  it('formatted string is valid APA', () => {
    expect(result.formatted).toMatch(/^F\(/)
    expect(result.formatted).toMatch(/p\s*[<=]/)
    expect(result.formatted).toMatch(/η²_p/)
  })

  it('testName is correct', () => {
    expect(result.testName).toBe('Repeated Measures ANOVA')
  })

  it('n is total observations', () => {
    expect(result.n).toBe(15)
  })
})

// ── Test Case 2: 2 conditions — equivalent to paired t-test ───────────
// R reference:
//   x1 <- c(5,6,7,8,9); x2 <- c(7,5,10,9,13)
//   t.test(x1, x2, paired=TRUE)
//   t = -2.0925, df=4, p = 0.1045
//   t² = 4.3784
//
//   fit <- aov(y ~ cond + Error(subj/cond), ...)
//   F = 4.378, df=(1,4), p = 0.1045

const data2cond = [
  [5, 7],
  [6, 5],
  [7, 10],
  [8, 9],
  [9, 13],
]

describe('repeatedMeasuresANOVA — 2 conditions (= paired t-test)', () => {
  const result = repeatedMeasuresANOVA(data2cond)
  const paired = tTestPaired(
    data2cond.map(r => r[0]!),
    data2cond.map(r => r[1]!)
  )

  it('F = t² (R cross-validated)', () => {
    const tSq = (paired.statistic) ** 2
    expect(result.statistic).toBeCloseTo(tSq, 6)
    expect(result.statistic).toBeCloseTo(4.3784, 2)
  })

  it('p-value matches paired t-test (tol=0.001)', () => {
    expect(result.pValue).toBeCloseTo(paired.pValue, 4)
    expect(result.pValue).toBeCloseTo(0.1045, 3)
  })

  it('sphericity is null for k=2', () => {
    expect(result.sphericity).toBeNull()
  })

  it('no correction applied for k=2', () => {
    expect(result.correction).toBe('none')
    expect(result.epsilonGG).toBe(1)
    expect(result.epsilonHF).toBe(1)
  })

  it('SS decomposition matches R', () => {
    expect(result.ssConditions).toBeCloseTo(8.1, 4)
    expect(result.ssSubjects).toBeCloseTo(39.4, 4)
    expect(result.ssError).toBeCloseTo(7.4, 4)
  })

  it('df correct', () => {
    expect(result.dfConditions).toBe(1)
    expect(result.dfSubjects).toBe(4)
    expect(result.dfError).toBe(4)
  })

  it('default labels used when none provided', () => {
    expect(result.conditions[0]!.label).toBe('C1')
    expect(result.conditions[1]!.label).toBe('C2')
  })
})

// ── Test Case 3: All identical values → F = 0 ─────────────────────────

describe('repeatedMeasuresANOVA — all identical values', () => {
  const dataConst = [
    [5, 5, 5],
    [5, 5, 5],
    [5, 5, 5],
  ]
  const result = repeatedMeasuresANOVA(dataConst)

  it('F = 0 when no condition differences', () => {
    expect(result.statistic).toBe(0)
  })

  it('p = 1 when F = 0', () => {
    expect(result.pValue).toBeCloseTo(1, 4)
  })

  it('all SS = 0', () => {
    expect(result.ssConditions).toBeCloseTo(0, 10)
    expect(result.ssSubjects).toBeCloseTo(0, 10)
    expect(result.ssError).toBeCloseTo(0, 10)
    expect(result.ssTotal).toBeCloseTo(0, 10)
  })

  it('partial eta-squared = 0', () => {
    expect(result.effectSize.value).toBe(0)
  })
})

// ── Test Case 4: 4 conditions with violated sphericity ─────────────────
// R reference:
//   set.seed(42); base <- rnorm(10, sd=2)
//   ... (see R code above)
//   F = 8.6854, p(uncorrected) = 0.000338
//   Mauchly W = 0.1067, chi-sq = 17.279, df=5, p = 0.004
//   GG epsilon = 0.585, HF = 0.715
//   GG-corrected df: (1.755, 15.793), GG-corrected p = 0.00368

const data4cond = [
  [3.394352, 3.588598, 7.108267, 9.771910],
  [0.013926, -1.020051, 3.985116, 3.065317],
  [0.031826, 1.640298, 6.831567, 10.517073],
  [1.126331, 2.873063, 2.438946, 3.632201],
  [0.741876, 2.756133, 5.323402, -0.032869],
  [0.105726, 0.572516, -2.363275, 7.951841],
  [2.880918, 3.894409, 3.669667, 4.966078],
  [-1.517546, -0.070900, 0.257959, 13.031188],
  [2.816614, 5.266896, -0.205776, 7.879616],
  [0.534628, 0.554574, 2.982940, 9.152811],
]

describe('repeatedMeasuresANOVA — 4 conditions, sphericity violated', () => {
  const result = repeatedMeasuresANOVA(data4cond)

  it('F statistic matches R (tol=0.01)', () => {
    expect(result.statistic).toBeCloseTo(8.6854, 2)
  })

  it('SS decomposition matches R', () => {
    expect(result.ssConditions).toBeCloseTo(206.283, 1)
    expect(result.ssError).toBeCloseTo(213.754, 1)
    expect(result.ssSubjects).toBeCloseTo(72.361, 1)
  })

  it('SS components sum to SS_total', () => {
    expect(result.ssConditions + result.ssSubjects + result.ssError)
      .toBeCloseTo(result.ssTotal, 4)
    expect(result.ssTotal).toBeCloseTo(492.399, 1)
  })

  it('Mauchly test detects sphericity violation', () => {
    expect(result.sphericity).not.toBeNull()
    expect(result.sphericity!.W).toBeCloseTo(0.1067, 2)
    expect(result.sphericity!.chiSq).toBeCloseTo(17.279, 0)
    expect(result.sphericity!.df).toBe(5)
    expect(result.sphericity!.pValue).toBeCloseTo(0.004, 2)
    expect(result.sphericity!.pValue).toBeLessThan(0.05)
  })

  it('GG correction applied (sphericity violated)', () => {
    expect(result.correction).toBe('Greenhouse-Geisser')
  })

  it('epsilon corrections match R', () => {
    expect(result.epsilonGG).toBeCloseTo(0.585, 2)
    expect(result.epsilonHF).toBeCloseTo(0.715, 2)
  })

  it('corrected df present and match R', () => {
    expect(result.correctedDf).toBeDefined()
    expect(result.correctedDf![0]).toBeCloseTo(1.755, 1)
    expect(result.correctedDf![1]).toBeCloseTo(15.793, 0)
  })

  it('p-value uses GG-corrected df, matches R', () => {
    // GG-corrected p from R: 0.003684
    expect(result.pValue).toBeCloseTo(0.00368, 3)
  })

  it('partial eta-squared matches R', () => {
    expect(result.effectSize.value).toBeCloseTo(0.4911, 2)
  })

  it('formatted string shows GG prefix', () => {
    expect(result.formatted).toMatch(/^F_GG\(/)
  })
})

// ── mauchlysTest standalone ────────────────────────────────────────────

describe('mauchlysTest', () => {
  it('matches R for 3-condition data', () => {
    const m = mauchlysTest(data3cond)
    expect(m.W).toBeCloseTo(0.3955, 2)
    expect(m.chiSq).toBeCloseTo(2.783, 1)
    expect(m.df).toBe(2)
    expect(m.pValue).toBeCloseTo(0.249, 2)
  })

  it('matches R for 4-condition data', () => {
    const m = mauchlysTest(data4cond)
    expect(m.W).toBeCloseTo(0.1067, 2)
    expect(m.chiSq).toBeCloseTo(17.279, 0)
    expect(m.df).toBe(5)
    expect(m.pValue).toBeCloseTo(0.004, 2)
  })

  it('throws for k < 3', () => {
    expect(() => mauchlysTest(data2cond)).toThrow('at least 3 conditions')
  })

  it('throws for n < 2', () => {
    expect(() => mauchlysTest([[1, 2, 3]])).toThrow('at least 2 subjects')
  })
})

// ── epsilonCorrections standalone ──────────────────────────────────────

describe('epsilonCorrections', () => {
  it('GG and HF match R for 3-condition data', () => {
    const eps = epsilonCorrections(data3cond)
    expect(eps.gg).toBeCloseTo(0.623, 2)
    expect(eps.hf).toBeCloseTo(0.769, 2)
  })

  it('GG and HF match R for 4-condition data', () => {
    const eps = epsilonCorrections(data4cond)
    expect(eps.gg).toBeCloseTo(0.585, 2)
    expect(eps.hf).toBeCloseTo(0.715, 2)
  })

  it('GG ≤ HF (always)', () => {
    const eps3 = epsilonCorrections(data3cond)
    expect(eps3.gg).toBeLessThanOrEqual(eps3.hf + 1e-10)
    const eps4 = epsilonCorrections(data4cond)
    expect(eps4.gg).toBeLessThanOrEqual(eps4.hf + 1e-10)
  })

  it('GG and HF clamped to [1/(k-1), 1]', () => {
    const eps3 = epsilonCorrections(data3cond)
    expect(eps3.gg).toBeGreaterThanOrEqual(0.5)  // 1/(3-1)
    expect(eps3.gg).toBeLessThanOrEqual(1)
    expect(eps3.hf).toBeGreaterThanOrEqual(0.5)
    expect(eps3.hf).toBeLessThanOrEqual(1)
  })

  it('perfect sphericity gives GG ≈ 1', () => {
    // Compound symmetry: same variance, same covariance
    const dataSym = [
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
      [10, 11, 12],
      [13, 14, 15],
    ]
    const eps = epsilonCorrections(dataSym)
    expect(eps.gg).toBeCloseTo(1.0, 1)
  })
})

// ── Edge cases ────────────────────────────────────────────────────────

describe('repeatedMeasuresANOVA — edge cases', () => {
  it('throws for fewer than 2 subjects', () => {
    expect(() => repeatedMeasuresANOVA([[1, 2, 3]])).toThrow('at least 2 subjects')
  })

  it('throws for fewer than 2 conditions', () => {
    expect(() => repeatedMeasuresANOVA([[1], [2]])).toThrow('at least 2 conditions')
  })

  it('throws for mismatched row lengths', () => {
    expect(() => repeatedMeasuresANOVA([[1, 2, 3], [4, 5]])).toThrow('columns, expected 3')
  })

  it('handles large F (zero within-subject error)', () => {
    // All subjects have identical pattern but different levels
    const data = [
      [1, 5, 9],
      [2, 6, 10],
      [3, 7, 11],
    ]
    const result = repeatedMeasuresANOVA(data)
    // SS_error should be 0 (perfect additivity)
    expect(result.ssError).toBeCloseTo(0, 8)
    // F should be Infinity (or very large) since MS_error = 0
    expect(result.statistic).toBe(Infinity)
  })

  it('works with 2 subjects (minimum)', () => {
    const data = [
      [1, 5, 9],
      [2, 3, 7],
    ]
    const result = repeatedMeasuresANOVA(data)
    expect(result.dfSubjects).toBe(1)
    expect(result.dfConditions).toBe(2)
    expect(result.dfError).toBe(2)
    expect(result.statistic).toBeGreaterThan(0)
  })

  it('custom labels are used', () => {
    const result = repeatedMeasuresANOVA(data3cond, { labels: ['Pre', 'Post', 'Follow-up'] })
    expect(result.conditions[0]!.label).toBe('Pre')
    expect(result.conditions[1]!.label).toBe('Post')
    expect(result.conditions[2]!.label).toBe('Follow-up')
  })
})
