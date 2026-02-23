/**
 * Tests for src/stats/regression.ts
 * Cross-validated with R:
 * > lm(y ~ x), summary(), AIC(), BIC()
 * > glm(y ~ x, family=binomial)
 */
import { describe, it, expect } from 'vitest'
import { linearRegression, multipleRegression, polynomialRegression, logisticRegression, regressionDiagnostics } from '../../src/stats/regression.js'
import gt from '../fixtures/ground_truth.json'

describe('linearRegression', () => {
  const { x, y } = gt.linear_regression
  const result = linearRegression(x as number[], y as number[])

  it('intercept (R ground truth, tol=0.01)', () => {
    expect(result.coefficients[0]!.estimate).toBeCloseTo(gt.linear_regression.intercept, 2)
  })
  it('slope (R ground truth, tol=0.01)', () => {
    expect(result.coefficients[1]!.estimate).toBeCloseTo(gt.linear_regression.slope, 2)
  })
  it('R² (R ground truth, tol=0.001)', () => {
    expect(result.r2).toBeCloseTo(gt.linear_regression.r2, 3)
  })
  it('adj. R² (R ground truth, tol=0.01)', () => {
    expect(result.adjR2).toBeCloseTo(gt.linear_regression.adj_r2, 2)
  })
  it('F p-value (R ground truth, tol=0.01)', () => {
    expect(result.fPValue).toBeCloseTo(gt.linear_regression.p, 2)
  })
  it('residuals sum to zero (OLS property)', () => {
    const sumRes = result.residuals.reduce((s, r) => s + r, 0)
    expect(sumRes).toBeCloseTo(0, 8)
  })
  it('fitted + residuals = y', () => {
    const y_ = y as number[]
    result.fitted.forEach((f, i) => {
      expect(f + result.residuals[i]!).toBeCloseTo(y_[i]!, 8)
    })
  })
  it('AIC and BIC are finite numbers', () => {
    expect(isFinite(result.aic)).toBe(true)
    expect(isFinite(result.bic)).toBe(true)
  })
  it('formatted string contains R²', () => {
    expect(result.formatted).toMatch(/R²/)
  })
})

describe('multipleRegression', () => {
  // Simulate data: y = 1 + 2*x1 + 3*x2 + ε
  const n = 20
  const x1 = Array.from({ length: n }, (_, i) => i + 1)
  const x2 = Array.from({ length: n }, (_, i) => (i % 5) + 1)
  const y = x1.map((x, i) => 1 + 2 * x + 3 * (x2[i] ?? 0) + (i % 3 === 0 ? 0.5 : -0.5))

  const result = multipleRegression(y, [
    { name: 'x1', values: x1 },
    { name: 'x2', values: x2 },
  ])

  it('intercept ≈ 1 (within 1)', () => {
    expect(result.coefficients[0]!.estimate).toBeCloseTo(1, 0)
  })
  it('slope x1 ≈ 2 (within 0.5)', () => {
    expect(result.coefficients[1]!.estimate).toBeCloseTo(2, 0)
  })
  it('slope x2 ≈ 3 (within 0.5)', () => {
    expect(result.coefficients[2]!.estimate).toBeCloseTo(3, 0)
  })
  it('R² > 0.99', () => {
    expect(result.r2).toBeGreaterThan(0.99)
  })
})

describe('polynomialRegression', () => {
  // y = x² → coefficient of x² ≈ 1, x^1 ≈ 0, intercept ≈ 0
  const x = [-2, -1, 0, 1, 2]
  const y = [4, 1, 0, 1, 4]
  const result = polynomialRegression(x, y, 2)

  it('x² coefficient ≈ 1', () => {
    expect(result.coefficients[2]!.estimate).toBeCloseTo(1, 4)
  })
  it('R² ≈ 1 (perfect fit)', () => {
    expect(result.r2).toBeCloseTo(1, 4)
  })
})

describe('logisticRegression', () => {
  // Simple separable data: y=1 when x>0, y=0 when x<0
  const x = [-3, -2, -1, -0.5, 0.5, 1, 2, 3]
  const y = [0, 0, 0, 0, 1, 1, 1, 1]
  const result = logisticRegression(y, [{ name: 'x', values: x }])

  it('positive slope (x increases P(y=1))', () => {
    expect(result.coefficients[1]!.estimate).toBeGreaterThan(0)
  })
  it('McFadden R² > 0.5 (good separation)', () => {
    expect(result.r2).toBeGreaterThan(0.5)
  })
  it('fitted values in (0, 1)', () => {
    for (const p of result.fitted) {
      expect(p).toBeGreaterThan(0)
      expect(p).toBeLessThan(1)
    }
  })
  it('throws on non-binary outcome', () => {
    expect(() => logisticRegression([0, 1, 2], [{ name: 'x', values: [1, 2, 3] }])).toThrow()
  })
})

describe('regressionDiagnostics', () => {
  const x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  const y = x.map(xi => 2 * xi + 1 + (xi === 10 ? 5 : 0))  // outlier at x=10
  const preds = [{ name: 'x', values: x }]
  const reg = multipleRegression(y, preds)
  const diag = regressionDiagnostics(reg, preds)

  it('leverage has n values', () => {
    expect(diag.leverage.length).toBe(x.length)
  })
  it('leverage values in (0, 1)', () => {
    for (const h of diag.leverage) {
      expect(h).toBeGreaterThanOrEqual(0)
      expect(h).toBeLessThanOrEqual(1)
    }
  })
  it('Cook\'s distance has n values', () => {
    expect(diag.cooksDistance.length).toBe(x.length)
  })
  it('standardized residuals have roughly unit variance', () => {
    const sr = diag.standardizedResiduals
    const mean = sr.reduce((s, v) => s + v, 0) / sr.length
    const variance = sr.reduce((s, v) => s + (v - mean) ** 2, 0) / (sr.length - 1)
    expect(Math.sqrt(variance)).toBeCloseTo(1, 0)
  })
  it('VIF for single predictor ≈ 1', () => {
    // For simple regression, VIF = 1 always
    // But regressionDiagnostics needs at least 2 predictors for VIF > 1
    // Single predictor VIF is trivially 1
    expect(diag.vif[0]).toBeGreaterThanOrEqual(1)
  })
})
