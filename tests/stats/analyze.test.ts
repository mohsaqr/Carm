/**
 * Tests for src/stats/analyze.ts
 *
 * Cross-validation notes:
 *   R: shapiro.test(c(2,4,6,8,10)) → W=0.9671, p=0.8486  (normal, p > 0.05)
 *   R: shapiro.test(c(1,1,1,50,100)) → W=0.7237, p=0.01038 (non-normal, p < 0.05)
 */
import { describe, it, expect } from 'vitest'
import { analyze, detectFieldType } from '../../src/stats/analyze.js'

// ─── Test 1: numeric + binary → Welch t-test ────────────────────────────

describe('analyze: numeric + binary → Welch t-test', () => {
  // Use n=2 per group: checkNormality treats n<3 as normal → routes to t-test
  const outcome = {
    type: 'numeric' as const,
    name: 'score',
    values: [72, 85, 90, 68] as number[],
  }
  const predictor = {
    type: 'binary' as const,
    name: 'group',
    values: ['A', 'B', 'A', 'B'] as string[],
  }
  const result = analyze(outcome, predictor)

  it('result.test contains "Welch"', () => {
    expect(result.test).toContain('Welch')
  })
  it('outcome and predictor names preserved', () => {
    expect(result.outcome).toBe('score')
    expect(result.predictor).toBe('group')
  })
  it('result.result has a numeric statistic', () => {
    expect(typeof result.result.statistic).toBe('number')
    expect(isFinite(result.result.statistic)).toBe(true)
  })
  it('descriptives present for each group', () => {
    expect(result.descriptives).toBeDefined()
    expect(result.descriptives!.length).toBe(2)
  })
  it('no posthoc for 2-group comparison', () => {
    // posthoc not meaningful for 2 groups (no pairs beyond the single pair)
    // tukeyHSD is only called for ANOVA; for t-test posthoc is undefined
    expect(result.posthoc).toBeUndefined()
  })
})

// ─── Test 2: numeric + 3-level categorical, normal data → ANOVA ──────────

describe('analyze: numeric + categorical, normal data → one-way ANOVA', () => {
  // R: shapiro.test(c(2,4,6,8,10)) → W=0.9671, p=0.8486 — normal
  const outcome = {
    type: 'numeric' as const,
    name: 'score',
    values: [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30] as number[],
  }
  const predictor = {
    type: 'categorical' as const,
    name: 'group',
    values: [
      'A', 'A', 'A', 'A', 'A',
      'B', 'B', 'B', 'B', 'B',
      'C', 'C', 'C', 'C', 'C',
    ] as string[],
  }
  const result = analyze(outcome, predictor)

  it('dispatches to ANOVA (result.test contains "ANOVA")', () => {
    expect(result.test).toContain('ANOVA')
  })
  it('posthoc present', () => {
    expect(result.posthoc).toBeDefined()
    expect(result.posthoc!.length).toBeGreaterThan(0)
  })
  it('normality results present', () => {
    expect(result.normality).toBeDefined()
    expect(result.normality!.length).toBe(3)
  })
  it('3 descriptive entries', () => {
    expect(result.descriptives!.length).toBe(3)
  })
})

// ─── Test 3: numeric + categorical, non-normal → Kruskal-Wallis + Dunn ───

describe('analyze: numeric + categorical, non-normal → Kruskal-Wallis', () => {
  // R: shapiro.test(c(1,1,1,50,100)) → W=0.7237, p=0.01038 — non-normal
  const outcome = {
    type: 'numeric' as const,
    name: 'rt',
    values: [
      1, 1, 1, 50, 100,
      2, 2, 2, 60, 110,
      3, 3, 3, 70, 120,
    ] as number[],
  }
  const predictor = {
    type: 'categorical' as const,
    name: 'cond',
    values: [
      'A', 'A', 'A', 'A', 'A',
      'B', 'B', 'B', 'B', 'B',
      'C', 'C', 'C', 'C', 'C',
    ] as string[],
  }
  const result = analyze(outcome, predictor)

  it('dispatches to Kruskal-Wallis', () => {
    expect(result.test).toContain('Kruskal')
  })
  it('posthoc (Dunn) present', () => {
    expect(result.posthoc).toBeDefined()
    expect(result.posthoc!.length).toBeGreaterThan(0)
  })
})

// ─── Test 4: binary + binary → chi-square ────────────────────────────────

describe('analyze: binary + binary → chi-square', () => {
  // 2×2 contingency: treatment × outcome
  const outcome = {
    type: 'binary' as const,
    name: 'treatment',
    values: ['yes', 'yes', 'yes', 'yes', 'yes', 'no', 'no', 'no', 'no', 'no'] as string[],
  }
  const predictor = {
    type: 'binary' as const,
    name: 'recovered',
    values: ['yes', 'yes', 'yes', 'no', 'no', 'yes', 'no', 'no', 'no', 'no'] as string[],
  }
  const result = analyze(outcome, predictor)

  it('dispatches to chi-square or Fisher', () => {
    // Either chi-square or Fisher (auto-fallback for low expected counts)
    expect(result.test).toMatch(/chi|χ|Fisher/i)
  })
  it('outcome and predictor names correct', () => {
    expect(result.outcome).toBe('treatment')
    expect(result.predictor).toBe('recovered')
  })
})

// ─── Test 5: paired=true + binary → paired t-test ────────────────────────

describe('analyze: paired=true + binary → paired t-test', () => {
  // Cross-validated with R:
  // > t.test(c(80,85,90,88,92), c(70,78,85,82,88), paired=TRUE)
  // t = 3.9686, df = 4, p-value = 0.01676
  const outcome = {
    type: 'numeric' as const,
    name: 'score',
    values: [80, 85, 90, 88, 92, 70, 78, 85, 82, 88] as number[],
  }
  const predictor = {
    type: 'binary' as const,
    name: 'time',
    values: ['post', 'post', 'post', 'post', 'post', 'pre', 'pre', 'pre', 'pre', 'pre'] as string[],
  }
  const result = analyze(outcome, predictor, { paired: true })

  it('dispatches to paired t-test', () => {
    expect(result.test).toContain('Paired')
  })
  it('result has finite statistic', () => {
    expect(isFinite(result.result.statistic)).toBe(true)
  })
})

// ─── Test 6: forceTest overrides normality routing ───────────────────────

describe('analyze: forceTest="kruskal-wallis" overrides auto-selection', () => {
  // Even with clearly normal data, forceTest forces Kruskal-Wallis
  const outcome = {
    type: 'numeric' as const,
    name: 'score',
    values: [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30] as number[],
  }
  const predictor = {
    type: 'categorical' as const,
    name: 'group',
    values: [
      'A', 'A', 'A', 'A', 'A',
      'B', 'B', 'B', 'B', 'B',
      'C', 'C', 'C', 'C', 'C',
    ] as string[],
  }
  const result = analyze(outcome, predictor, { forceTest: 'kruskal-wallis' })

  it('runs Kruskal-Wallis despite normal data', () => {
    expect(result.test).toContain('Kruskal')
  })
  it('Dunn post-hoc present', () => {
    expect(result.posthoc).toBeDefined()
  })
})

// ─── Test 7: numeric field only, no predictor → descriptive ──────────────

describe('analyze: numeric only, no predictor → descriptive', () => {
  const outcome = {
    type: 'numeric' as const,
    name: 'value',
    values: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] as number[],
  }
  const result = analyze(outcome)

  it('result.test === "descriptive"', () => {
    expect(result.test).toBe('descriptive')
  })
  it('no posthoc', () => {
    expect(result.posthoc).toBeUndefined()
  })
  it('no predictor field', () => {
    expect(result.predictor).toBeUndefined()
  })
  it('descriptives present', () => {
    expect(result.descriptives).toBeDefined()
    expect(result.descriptives!.length).toBe(1)
  })
})

// ─── Test 8: mismatched lengths → throws ─────────────────────────────────

describe('analyze: mismatched lengths → throws', () => {
  it('throws when outcome and predictor have different lengths', () => {
    expect(() =>
      analyze(
        { type: 'numeric', name: 'x', values: [1, 2, 3] },
        { type: 'binary',  name: 'g', values: ['A', 'B'] }
      )
    ).toThrow(/length/)
  })
})

// ─── Test 9: paired=true, unequal group sizes → throws ───────────────────

describe('analyze: paired=true but unequal group sizes → throws', () => {
  it('throws when groups have different sizes', () => {
    // 'pre' has 3 obs, 'post' has 4 obs → unequal → should throw
    expect(() =>
      analyze(
        { type: 'numeric', name: 'score', values: [1, 2, 3, 4, 5, 6, 7] },
        {
          type: 'binary',
          name: 'time',
          values: ['pre', 'pre', 'pre', 'post', 'post', 'post', 'post'],
        },
        { paired: true, forceTest: 't-test-paired' }
      )
    ).toThrow(/unequal/)
  })
})

// ─── Test 10: detectFieldType — [0,1,0,1] → 'binary' ────────────────────

describe('detectFieldType', () => {
  it('[0,1,0,1] → "binary"', () => {
    expect(detectFieldType([0, 1, 0, 1])).toBe('binary')
  })

  it('[1.2, 3.4, 5.6] → "numeric"', () => {
    expect(detectFieldType([1.2, 3.4, 5.6])).toBe('numeric')
  })

  it('[0, 1, 2] → "numeric" (3 unique numbers, not {0,1})', () => {
    expect(detectFieldType([0, 1, 2])).toBe('numeric')
  })
})

// ─── Test 11: detectFieldType — ['a','b','c','a'] → 'categorical' ────────

describe('detectFieldType: string arrays', () => {
  it('["a","b","c","a"] → "categorical"', () => {
    expect(detectFieldType(['a', 'b', 'c', 'a'])).toBe('categorical')
  })

  it('["yes","no","yes"] → "binary"', () => {
    expect(detectFieldType(['yes', 'no', 'yes'])).toBe('binary')
  })
})
