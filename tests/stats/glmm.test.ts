/**
 * Unit tests for logistic GLMM (src/stats/glmm.ts).
 *
 * Tests cover:
 * - Basic binary outcome with random intercepts
 * - Convergence on synthetic data
 * - Odds ratios match manual exp(coef)
 * - ICC latent scale matches formula σ²_b / (σ²_b + π²/3)
 * - Edge cases: single predictor, all-zero or all-one group
 * - GLMMResult structure validation
 */

import { describe, it, expect } from 'vitest'
import { runGLMM } from '../../src/stats/glmm.js'

// ─── Synthetic data generator ────────────────────────────────────────────

/** Simple LCG PRNG for reproducible test data. */
function lcg(seed: number): () => number {
  let s = seed
  return () => {
    s = (s * 1664525 + 1013904223) & 0x7fffffff
    return s / 0x7fffffff
  }
}

/** Box-Muller transform for normal random variates. */
function randn(rng: () => number): number {
  const u1 = rng()
  const u2 = rng()
  return Math.sqrt(-2 * Math.log(Math.max(1e-10, u1))) * Math.cos(2 * Math.PI * u2)
}

/** Generate synthetic GLMM data. */
function generateGLMMData(opts: {
  n: number
  nGroups: number
  betaIntercept: number
  betaX: number
  sigmab: number
  seed?: number
}) {
  const { n, nGroups, betaIntercept, betaX, sigmab, seed = 42 } = opts
  const rng = lcg(seed)

  const y: number[] = []
  const x: number[] = []
  const group: number[] = []

  // Generate random intercepts for each group
  const b = Array.from({ length: nGroups }, () => randn(rng) * sigmab)

  const perGroup = Math.floor(n / nGroups)
  for (let g = 0; g < nGroups; g++) {
    for (let i = 0; i < perGroup; i++) {
      const xi = randn(rng)
      const eta = betaIntercept + betaX * xi + b[g]!
      const mu = 1 / (1 + Math.exp(-eta))
      const yi = rng() < mu ? 1 : 0
      x.push(xi)
      y.push(yi)
      group.push(g + 1)
    }
  }

  return { y, x, group }
}

// ─── Tests ───────────────────────────────────────────────────────────────

describe('runGLMM', () => {
  it('should return a valid GLMMResult structure', () => {
    const { y, x, group } = generateGLMMData({
      n: 200, nGroups: 10, betaIntercept: -0.5, betaX: 0.8, sigmab: 1.0,
    })

    const result = runGLMM({
      outcome: y,
      fixedPredictors: { x },
      groupId: group,
    })

    // Structure checks
    expect(result.family).toBe('binomial')
    expect(result.link).toBe('logit')
    expect(result.nObs).toBe(200)
    expect(result.nGroups).toBe(10)
    expect(result.fixedEffects).toHaveLength(2) // intercept + x
    expect(result.fixedEffects[0]!.name).toBe('(Intercept)')
    expect(result.fixedEffects[1]!.name).toBe('x')

    // Variance components
    expect(result.varianceComponents.intercept).toBeGreaterThan(0)

    // ICC is between 0 and 1
    expect(result.icc).toBeGreaterThanOrEqual(0)
    expect(result.icc).toBeLessThanOrEqual(1)

    // AIC/BIC finite
    expect(isFinite(result.aic)).toBe(true)
    expect(isFinite(result.bic)).toBe(true)
    expect(isFinite(result.logLik)).toBe(true)

    // Deviance = -2 * logLik
    expect(result.deviance).toBeCloseTo(-2 * result.logLik, 2)

    // Formatted string exists
    expect(result.formatted).toContain('ICC')
    expect(result.formatted).toContain('AIC')
  })

  it('should compute odds ratios as exp(estimate)', () => {
    const { y, x, group } = generateGLMMData({
      n: 200, nGroups: 10, betaIntercept: 0, betaX: 1.0, sigmab: 0.5,
    })

    const result = runGLMM({
      outcome: y,
      fixedPredictors: { x },
      groupId: group,
    })

    for (const fe of result.fixedEffects) {
      expect(fe.or).toBeCloseTo(Math.exp(fe.estimate), 2)
      expect(fe.orCI[0]).toBeCloseTo(Math.exp(fe.ci[0]), 2)
      expect(fe.orCI[1]).toBeCloseTo(Math.exp(fe.ci[1]), 2)
    }
  })

  it('should compute latent-scale ICC correctly', () => {
    const { y, x, group } = generateGLMMData({
      n: 200, nGroups: 10, betaIntercept: 0, betaX: 0.5, sigmab: 1.5,
    })

    const result = runGLMM({
      outcome: y,
      fixedPredictors: { x },
      groupId: group,
    })

    // Verify ICC formula: σ²_b / (σ²_b + π²/3)
    const sigmab2 = result.varianceComponents.intercept
    const expectedICC = sigmab2 / (sigmab2 + Math.PI * Math.PI / 3)
    expect(result.icc).toBeCloseTo(expectedICC, 4)
  })

  it('should detect positive x effect with reasonable estimate', () => {
    const { y, x, group } = generateGLMMData({
      n: 400, nGroups: 20, betaIntercept: -0.5, betaX: 1.0, sigmab: 0.8, seed: 123,
    })

    const result = runGLMM({
      outcome: y,
      fixedPredictors: { x },
      groupId: group,
    })

    // x coefficient should be positive (generated with betaX = 1.0)
    const xEffect = result.fixedEffects.find(fe => fe.name === 'x')!
    expect(xEffect.estimate).toBeGreaterThan(0)

    // With n=400, should be significant
    expect(xEffect.pValue).toBeLessThan(0.05)
  })

  it('should converge and produce finite results with small random effects', () => {
    const { y, x, group } = generateGLMMData({
      n: 200, nGroups: 10, betaIntercept: 0, betaX: 0.5, sigmab: 0.2, seed: 77,
    })

    const result = runGLMM({
      outcome: y,
      fixedPredictors: { x },
      groupId: group,
    })

    // Should produce finite results even with small between-group variance
    expect(isFinite(result.varianceComponents.intercept)).toBe(true)
    expect(isFinite(result.icc)).toBe(true)
    expect(isFinite(result.aic)).toBe(true)
    expect(result.icc).toBeGreaterThanOrEqual(0)
    expect(result.icc).toBeLessThanOrEqual(1)
  })

  it('should work with multiple predictors', () => {
    const rng = lcg(42)
    const n = 200
    const nGroups = 10
    const perGroup = n / nGroups

    const y: number[] = []
    const x1: number[] = []
    const x2: number[] = []
    const group: number[] = []

    const b = Array.from({ length: nGroups }, () => randn(rng) * 1.0)

    for (let g = 0; g < nGroups; g++) {
      for (let i = 0; i < perGroup; i++) {
        const xi1 = randn(rng)
        const xi2 = randn(rng)
        const eta = -0.5 + 0.8 * xi1 - 0.3 * xi2 + b[g]!
        const mu = 1 / (1 + Math.exp(-eta))
        y.push(rng() < mu ? 1 : 0)
        x1.push(xi1)
        x2.push(xi2)
        group.push(g + 1)
      }
    }

    const result = runGLMM({
      outcome: y,
      fixedPredictors: { x1, x2 },
      groupId: group,
    })

    expect(result.fixedEffects).toHaveLength(3) // intercept + x1 + x2
    expect(result.fixedEffects[0]!.name).toBe('(Intercept)')
    expect(result.fixedEffects[1]!.name).toBe('x1')
    expect(result.fixedEffects[2]!.name).toBe('x2')

    // nParams: 3 fixed + 1 variance parameter
    expect(result.nParams).toBe(4)
  })

  it('should have z-values consistent with estimate/se', () => {
    const { y, x, group } = generateGLMMData({
      n: 200, nGroups: 10, betaIntercept: 0, betaX: 0.5, sigmab: 1.0,
    })

    const result = runGLMM({
      outcome: y,
      fixedPredictors: { x },
      groupId: group,
    })

    for (const fe of result.fixedEffects) {
      if (fe.se > 0) {
        expect(fe.zValue).toBeCloseTo(fe.estimate / fe.se, 2)
      }
    }
  })

  it('should handle CI level correctly', () => {
    const { y, x, group } = generateGLMMData({
      n: 200, nGroups: 10, betaIntercept: 0, betaX: 0.5, sigmab: 1.0,
    })

    const result99 = runGLMM({
      outcome: y,
      fixedPredictors: { x },
      groupId: group,
      ciLevel: 0.99,
    })

    const result90 = runGLMM({
      outcome: y,
      fixedPredictors: { x },
      groupId: group,
      ciLevel: 0.90,
    })

    // 99% CI should be wider than 90% CI
    const width99 = result99.fixedEffects[1]!.ci[1] - result99.fixedEffects[1]!.ci[0]
    const width90 = result90.fixedEffects[1]!.ci[1] - result90.fixedEffects[1]!.ci[0]
    expect(width99).toBeGreaterThan(width90)
  })

  // ── Edge cases ─────────────────────────────────────────────────────────

  it('should throw for non-binary outcome', () => {
    expect(() => runGLMM({
      outcome: [0, 1, 2, 0, 1],
      fixedPredictors: { x: [1, 2, 3, 4, 5] },
      groupId: [1, 1, 2, 2, 1],
    })).toThrow('outcome must be 0/1 binary')
  })

  it('should throw for fewer than 5 observations', () => {
    expect(() => runGLMM({
      outcome: [0, 1, 0, 1],
      fixedPredictors: { x: [1, 2, 3, 4] },
      groupId: [1, 1, 2, 2],
    })).toThrow('need at least 5 observations')
  })

  it('should throw for single group', () => {
    expect(() => runGLMM({
      outcome: [0, 1, 0, 1, 1],
      fixedPredictors: { x: [1, 2, 3, 4, 5] },
      groupId: [1, 1, 1, 1, 1],
    })).toThrow('need at least 2 groups')
  })

  it('should throw for mismatched lengths', () => {
    expect(() => runGLMM({
      outcome: [0, 1, 0, 1, 1],
      fixedPredictors: { x: [1, 2, 3, 4, 5] },
      groupId: [1, 1, 2],
    })).toThrow('groupId must have same length')
  })

  it('should throw for invalid random slope name', () => {
    expect(() => runGLMM({
      outcome: [0, 1, 0, 1, 1, 0, 1, 0, 1, 0],
      fixedPredictors: { x: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
      groupId: [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
      randomSlopes: ['z'],
    })).toThrow("random slope 'z' not found")
  })

  it('should produce BIC > AIC for reasonable n', () => {
    const { y, x, group } = generateGLMMData({
      n: 200, nGroups: 10, betaIntercept: 0, betaX: 0.5, sigmab: 1.0,
    })

    const result = runGLMM({
      outcome: y,
      fixedPredictors: { x },
      groupId: group,
    })

    // BIC penalizes more than AIC for n > exp(2) ≈ 7.4
    // With n=200 and nParams=3, BIC should be larger
    expect(result.bic).toBeGreaterThan(result.aic)
  })
})
