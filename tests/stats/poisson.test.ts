/**
 * Tests for Poisson regression (src/stats/regression.ts: poissonRegression)
 *
 * Cross-validated with R (verified 2026-02-27):
 *
 * Simple:
 * > mod <- glm(y ~ x, family = poisson, data = data.frame(y=c(0,1,2,3,5,8,13,21,34,55), x=1:10))
 * > coef(mod) => (Intercept)=-0.9198, x=0.4941
 * > deviance(mod) => 1.398
 * > AIC(mod) => 40.946
 * > logLik(mod) => -18.473
 *
 * Multiple predictors:
 * > mod <- glm(y ~ x1 + x2, family = poisson, data = df)
 * > coef(mod) => (Intercept)=1.017, x1=0.554, x2=-0.289
 *
 * Zero-heavy:
 * > mod <- glm(y ~ x, family=poisson, data=data.frame(y=c(0,0,0,0,0,1,0,2,0,1), x=1:10))
 * > coef(mod) => (Intercept)=-3.428, x=0.366
 */
import { describe, it, expect } from 'vitest'
import { poissonRegression } from '../../src/stats/regression.js'
import ref from '../fixtures/lmm-poisson-ref.json'

// ─── Simple Poisson regression ───────────────────────────────────────────

describe('poissonRegression simple', () => {
  const y = [0, 1, 2, 3, 5, 8, 13, 21, 34, 55]
  const x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  const result = poissonRegression(y, [{ name: 'x', values: x }])

  it('intercept matches R (≈ -0.920)', () => {
    const intercept = result.coefficients.find(c => c.name === '(Intercept)')!
    expect(intercept.estimate).toBeCloseTo(ref.poisson_simple.intercept, 3)
  })

  it('slope matches R (≈ 0.494)', () => {
    const slope = result.coefficients.find(c => c.name === 'x')!
    expect(slope.estimate).toBeCloseTo(ref.poisson_simple.slope, 3)
  })

  it('SE of intercept matches R', () => {
    const intercept = result.coefficients.find(c => c.name === '(Intercept)')!
    expect(intercept.se).toBeCloseTo(ref.poisson_simple.se_intercept, 2)
  })

  it('SE of slope matches R', () => {
    const slope = result.coefficients.find(c => c.name === 'x')!
    expect(slope.se).toBeCloseTo(ref.poisson_simple.se_slope, 3)
  })

  it('AIC matches R (≈ 40.95)', () => {
    expect(result.aic).toBeCloseTo(ref.poisson_simple.aic, 0)
  })

  it('BIC matches R (≈ 41.55)', () => {
    expect(result.bic).toBeCloseTo(ref.poisson_simple.bic, 0)
  })

  it('n is correct', () => {
    expect(result.n).toBe(10)
  })

  it('fitted values are positive', () => {
    for (const f of result.fitted) {
      expect(f).toBeGreaterThan(0)
    }
  })

  it('residuals have correct length', () => {
    expect(result.residuals.length).toBe(10)
  })

  it('pseudo-R² is in [0, 1] and high for good fit', () => {
    expect(result.r2).toBeGreaterThan(0.5)
    expect(result.r2).toBeLessThanOrEqual(1)
  })

  it('formatted string contains deviance', () => {
    expect(result.formatted).toMatch(/Deviance/)
  })
})

// ─── Multiple predictor Poisson ──────────────────────────────────────────

describe('poissonRegression multiple predictors', () => {
  const { y, x1, x2 } = ref.poisson_multi
  const result = poissonRegression(
    y as number[],
    [
      { name: 'x1', values: x1 as number[] },
      { name: 'x2', values: x2 as number[] },
    ]
  )

  it('intercept matches R (≈ 1.017)', () => {
    expect(result.coefficients.find(c => c.name === '(Intercept)')!.estimate).toBeCloseTo(ref.poisson_multi.intercept, 2)
  })

  it('x1 coefficient matches R (≈ 0.554)', () => {
    expect(result.coefficients.find(c => c.name === 'x1')!.estimate).toBeCloseTo(ref.poisson_multi.x1_coef, 2)
  })

  it('x2 coefficient matches R (≈ -0.289)', () => {
    expect(result.coefficients.find(c => c.name === 'x2')!.estimate).toBeCloseTo(ref.poisson_multi.x2_coef, 2)
  })

  it('AIC matches R', () => {
    expect(result.aic).toBeCloseTo(ref.poisson_multi.aic, 0)
  })

  it('n is 50', () => {
    expect(result.n).toBe(50)
  })

  it('all fitted values are positive', () => {
    for (const f of result.fitted) {
      expect(f).toBeGreaterThan(0)
    }
  })
})

// ─── Zero-heavy Poisson ─────────────────────────────────────────────────

describe('poissonRegression zero-heavy', () => {
  const y = [0, 0, 0, 0, 0, 1, 0, 2, 0, 1]
  const x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  const result = poissonRegression(y, [{ name: 'x', values: x }])

  it('intercept matches R (≈ -3.43)', () => {
    expect(result.coefficients.find(c => c.name === '(Intercept)')!.estimate).toBeCloseTo(ref.poisson_zeros.intercept, 1)
  })

  it('slope matches R (≈ 0.366)', () => {
    expect(result.coefficients.find(c => c.name === 'x')!.estimate).toBeCloseTo(ref.poisson_zeros.slope, 1)
  })

  it('AIC matches R (≈ 17.40)', () => {
    expect(result.aic).toBeCloseTo(ref.poisson_zeros.aic, 0)
  })
})

// ─── Edge cases ──────────────────────────────────────────────────────────

describe('poissonRegression edge cases', () => {
  it('throws on negative y', () => {
    expect(() =>
      poissonRegression([-1, 0, 1], [{ name: 'x', values: [1, 2, 3] }])
    ).toThrow(/non-negative/)
  })

  it('handles all-zeros y', () => {
    const result = poissonRegression(
      [0, 0, 0, 0, 0],
      [{ name: 'x', values: [1, 2, 3, 4, 5] }]
    )
    // Intercept should be very negative (log of ~0)
    expect(result.coefficients[0]!.estimate).toBeLessThan(-5)
    expect(result.fitted.every(f => f >= 0)).toBe(true)
  })

  it('handles single predictor, n=5', () => {
    const result = poissonRegression(
      [1, 2, 3, 4, 5],
      [{ name: 'x', values: [1, 2, 3, 4, 5] }]
    )
    expect(result.coefficients.length).toBe(2) // intercept + x
    expect(result.n).toBe(5)
    expect(isFinite(result.aic)).toBe(true)
  })

  it('perfect exponential fit recovers parameters', () => {
    // y = exp(0.5 + 0.3*x) approximately
    const x = [1, 2, 3, 4, 5, 6, 7, 8]
    const y = x.map(xi => Math.round(Math.exp(0.5 + 0.3 * xi)))
    const result = poissonRegression(y, [{ name: 'x', values: x }])
    // Should recover intercept ≈ 0.5 and slope ≈ 0.3 approximately
    expect(result.coefficients[0]!.estimate).toBeCloseTo(0.5, 0)
    expect(result.coefficients[1]!.estimate).toBeCloseTo(0.3, 0)
  })

  it('no predictor (intercept-only) works', () => {
    const result = poissonRegression([1, 2, 3, 4, 5], [])
    expect(result.coefficients.length).toBe(1) // just intercept
    // intercept should be log(mean(y)) = log(3)
    expect(result.coefficients[0]!.estimate).toBeCloseTo(Math.log(3), 2)
  })
})
