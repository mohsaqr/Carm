/**
 * Regression analysis module.
 * Simple and multiple OLS, logistic regression, polynomial regression,
 * diagnostics (R², AIC, BIC, VIF, residual plots).
 */

import { mean as _mean, variance as _variance, fDistPValue, tDistPValue, tDistQuantile, roundTo } from '../core/math.js'
import { formatRegression } from '../core/apa.js'
import { Matrix } from '../core/matrix.js'
import type { RegressionResult, RegressionCoef } from '../core/types.js'

// ─── OLS Helper ───────────────────────────────────────────────────────────

/**
 * Fit OLS regression: y = Xβ + ε
 * Returns coefficients, SE, t-statistics, p-values, R², AIC, BIC.
 *
 * Cross-validated with R:
 * > lm(y ~ x1 + x2, data = df)
 */
function fitOLS(
  X: Matrix,
  y: readonly number[],
  coefNames: readonly string[],
  ciLevel = 0.95
): RegressionResult {
  const n = y.length
  const p = X.cols  // includes intercept column

  const Xt = X.transpose()
  const XtX = Xt.multiply(X)
  const XtY = Xt.multiply(Matrix.colVec(y))
  const XtXInv = XtX.inverse()
  const betaM = XtXInv.multiply(XtY)
  const beta = Array.from({ length: p }, (_, i) => betaM.get(i, 0))

  // Fitted values and residuals
  const fitted = Array.from({ length: n }, (_, i) => {
    let val = 0
    for (let j = 0; j < p; j++) val += X.get(i, j) * (beta[j] ?? 0)
    return val
  })
  const residuals = y.map((v, i) => v - (fitted[i] ?? 0))

  // RSS, TSS
  const yMean = _mean(y)
  const ss_res = residuals.reduce((s, r) => s + r * r, 0)
  const ss_tot = y.reduce((s, v) => s + (v - yMean) ** 2, 0)
  const r2 = ss_tot > 0 ? Math.max(0, 1 - ss_res / ss_tot) : 0
  const adjR2 = 1 - (1 - r2) * (n - 1) / (n - p)

  // sigma² = RSS / (n - p)
  const dfRes = n - p
  if (dfRes <= 0) throw new Error('fitOLS: not enough degrees of freedom')
  const sigma2 = ss_res / dfRes

  // SE of coefficients: diag(sigma² * (X'X)^-1)
  const covBeta = XtXInv.scale(sigma2)
  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, dfRes)

  const coefficients: RegressionCoef[] = beta.map((b, i) => {
    const se = Math.sqrt(Math.max(0, covBeta.get(i, i)))
    const t = se === 0 ? 0 : b / se
    const pVal = tDistPValue(t, dfRes)
    const ci: readonly [number, number] = [b - tCrit * se, b + tCrit * se]
    return {
      name: coefNames[i] ?? `β${i}`,
      estimate: roundTo(b, 6),
      se: roundTo(se, 6),
      tValue: roundTo(t, 4),
      pValue: roundTo(pVal, 4),
      ci,
    }
  })

  // F-statistic for overall model
  const dfModel = p - 1  // excluding intercept
  const ss_reg = ss_tot - ss_res
  const F = sigma2 === 0 || dfModel === 0 ? 0 : (ss_reg / dfModel) / sigma2
  const fPValue = fDistPValue(F, dfModel, dfRes)

  // AIC and BIC: -2 * logLik + penalty
  // logLik for normal errors: -n/2 * log(2π) - n/2 * log(σ²) - RSS/(2σ²)
  const logLik = -n / 2 * (Math.log(2 * Math.PI) + Math.log(ss_res / n) + 1)
  const aic = -2 * logLik + 2 * (p + 1)
  const bic = -2 * logLik + Math.log(n) * (p + 1)

  const formatted = formatRegression(r2, adjR2, F, dfModel, dfRes, fPValue)

  return {
    coefficients,
    r2: roundTo(r2, 6),
    adjR2: roundTo(adjR2, 6),
    fStatistic: roundTo(F, 4),
    fDf: [dfModel, dfRes],
    fPValue: roundTo(fPValue, 4),
    aic: roundTo(aic, 2),
    bic: roundTo(bic, 2),
    residuals,
    fitted,
    n,
    formatted,
  }
}

// ─── Simple linear regression ─────────────────────────────────────────────

/**
 * Simple linear regression: y = β₀ + β₁·x
 */
export function linearRegression(
  x: readonly number[],
  y: readonly number[],
  ciLevel = 0.95
): RegressionResult {
  if (x.length !== y.length) throw new Error('linearRegression: arrays must have equal length')
  if (x.length < 3) throw new Error('linearRegression: need at least 3 observations')
  const n = x.length
  const X = Matrix.fromArray(Array.from({ length: n }, (_, i) => [1, x[i] ?? 0]))
  return fitOLS(X, y, ['(Intercept)', 'x'], ciLevel)
}

// ─── Multiple linear regression ───────────────────────────────────────────

/**
 * Multiple linear regression: y = β₀ + β₁·x₁ + ... + βₖ·xₖ
 * `predictors`: named columns { name: values[] }
 */
export function multipleRegression(
  y: readonly number[],
  predictors: ReadonlyArray<{ name: string; values: readonly number[] }>,
  ciLevel = 0.95
): RegressionResult {
  if (predictors.length === 0) throw new Error('multipleRegression: need at least 1 predictor')
  const n = y.length
  for (const p of predictors) {
    if (p.values.length !== n) throw new Error(`multipleRegression: predictor '${p.name}' length mismatch`)
  }

  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map(p => p.values[i] ?? 0)])
  )
  const names = ['(Intercept)', ...predictors.map(p => p.name)]
  return fitOLS(X, y, names, ciLevel)
}

// ─── Polynomial regression ────────────────────────────────────────────────

/**
 * Polynomial regression: y = β₀ + β₁·x + β₂·x² + ... + βₖ·xᵏ
 */
export function polynomialRegression(
  x: readonly number[],
  y: readonly number[],
  degree: number,
  ciLevel = 0.95
): RegressionResult {
  if (degree < 1) throw new Error('polynomialRegression: degree must be ≥ 1')
  if (x.length !== y.length) throw new Error('polynomialRegression: arrays must match length')
  const n = x.length
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) =>
      [1, ...Array.from({ length: degree }, (_, d) => (x[i] ?? 0) ** (d + 1))]
    )
  )
  const names = ['(Intercept)', ...Array.from({ length: degree }, (_, d) => `x^${d + 1}`)]
  return fitOLS(X, y, names, ciLevel)
}

// ─── Logistic regression ──────────────────────────────────────────────────

/**
 * Binary logistic regression via IRLS (iteratively reweighted least squares).
 * Outcome y must be 0/1.
 *
 * Cross-validated with R:
 * > glm(y ~ x1 + x2, family = binomial, data = df)
 */
export function logisticRegression(
  y: readonly number[],
  predictors: ReadonlyArray<{ name: string; values: readonly number[] }>,
  ciLevel = 0.95,
  maxIter = 100,
  tol = 1e-8
): RegressionResult {
  for (const v of y) {
    if (v !== 0 && v !== 1) throw new Error('logisticRegression: y must be 0 or 1')
  }
  const n = y.length
  const p = predictors.length + 1  // +intercept

  // Design matrix
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map(pr => pr.values[i] ?? 0)])
  )
  const names = ['(Intercept)', ...predictors.map(pr => pr.name)]

  // IRLS: start with β = 0
  let beta = new Array<number>(p).fill(0)

  for (let iter = 0; iter < maxIter; iter++) {
    // η = Xβ, μ = logistic(η), W = diag(μ(1-μ))
    const eta = Array.from({ length: n }, (_, i) => {
      let v = 0
      for (let j = 0; j < p; j++) v += X.get(i, j) * (beta[j] ?? 0)
      return v
    })
    const mu = eta.map(e => 1 / (1 + Math.exp(-e)))
    const w = mu.map(m => Math.max(1e-10, m * (1 - m)))

    // Weighted least squares: X'WX · Δβ = X'W(y - μ)
    const Xw = Matrix.fromArray(
      Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_, j) => X.get(i, j) * Math.sqrt(w[i]!)))
    )
    const yAdj = Array.from({ length: n }, (_, i) => Math.sqrt(w[i]!) * ((y[i] ?? 0) - (mu[i] ?? 0)))

    try {
      const Xwt = Xw.transpose()
      const XwtXw = Xwt.multiply(Xw)
      const XwtY = Xwt.multiply(Matrix.colVec(yAdj))
      const delta = XwtXw.inverse().multiply(XwtY)
      let maxChange = 0
      for (let j = 0; j < p; j++) {
        const d = delta.get(j, 0)
        beta[j] = (beta[j] ?? 0) + d
        maxChange = Math.max(maxChange, Math.abs(d))
      }
      if (maxChange < tol) break
    } catch {
      break  // singular — stop iterating
    }
  }

  // Compute SEs from Fisher information matrix
  const eta = Array.from({ length: n }, (_, i) => {
    let v = 0
    for (let j = 0; j < p; j++) v += X.get(i, j) * (beta[j] ?? 0)
    return v
  })
  const mu = eta.map(e => 1 / (1 + Math.exp(-e)))
  const w = mu.map(m => Math.max(1e-10, m * (1 - m)))

  const Xw = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_, j) => X.get(i, j) * Math.sqrt(w[i]!)))
  )
  let cov: Matrix
  try {
    cov = Xw.transpose().multiply(Xw).inverse()
  } catch {
    cov = Matrix.identity(p)
  }

  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, n - p)
  const coefficients: RegressionCoef[] = beta.map((b, i) => {
    const se = Math.sqrt(Math.max(0, cov.get(i, i)))
    const z = se === 0 ? 0 : b / se
    const pVal = 2 * (1 - normCDFLocal(Math.abs(z)))
    const ci: readonly [number, number] = [b - tCrit * se, b + tCrit * se]
    return {
      name: names[i] ?? `β${i}`,
      estimate: roundTo(b, 6),
      se: roundTo(se, 6),
      tValue: roundTo(z, 4),
      pValue: roundTo(pVal, 4),
      ci,
    }
  })

  // Log-likelihood
  const logLik = mu.reduce((s, m, i) => {
    const yi = y[i] ?? 0
    return s + yi * Math.log(Math.max(1e-15, m)) + (1 - yi) * Math.log(Math.max(1e-15, 1 - m))
  }, 0)

  // Null log-likelihood (intercept only)
  const pMean = _mean([...y])
  const nullLogLik = n * (pMean * Math.log(Math.max(1e-15, pMean)) + (1 - pMean) * Math.log(Math.max(1e-15, 1 - pMean)))

  // McFadden pseudo-R²
  const r2 = 1 - logLik / nullLogLik

  const aic = -2 * logLik + 2 * p
  const bic = -2 * logLik + Math.log(n) * p
  const residuals = y.map((v, i) => (v ?? 0) - (mu[i] ?? 0))

  return {
    coefficients,
    r2: roundTo(r2, 6),
    adjR2: roundTo(r2, 6),  // McFadden's for logistic
    fStatistic: NaN,
    fDf: [p - 1, n - p],
    fPValue: NaN,
    aic: roundTo(aic, 2),
    bic: roundTo(bic, 2),
    residuals,
    fitted: mu,
    n,
    formatted: `McFadden R² = ${roundTo(r2, 3)}, AIC = ${roundTo(aic, 1)}`,
  }
}

function normCDFLocal(z: number): number {
  const x = Math.abs(z) / Math.SQRT2
  const t = 1 / (1 + 0.3275911 * x)
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
  const erf = 1 - poly * Math.exp(-x * x)
  return 0.5 * (1 + (z >= 0 ? erf : -erf))
}

// ─── Regression diagnostics ───────────────────────────────────────────────

export interface RegressionDiagnostics {
  readonly leverage: readonly number[]       // hat matrix diagonal
  readonly cooksDistance: readonly number[]  // Cook's D
  readonly standardizedResiduals: readonly number[]
  readonly vif: readonly number[]            // variance inflation factors
}

/**
 * Compute regression diagnostics.
 * Returns leverage (hat values), Cook's distance, standardized residuals, VIF.
 */
export function regressionDiagnostics(
  result: RegressionResult,
  predictors: ReadonlyArray<{ name: string; values: readonly number[] }>
): RegressionDiagnostics {
  const n = result.n
  const p = result.coefficients.length

  // Rebuild design matrix
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map(pr => pr.values[i] ?? 0)])
  )

  // Hat matrix diagonal: h_ii = X(X'X)^-1X'
  const Xt = X.transpose()
  let XtXInv: Matrix
  try { XtXInv = Xt.multiply(X).inverse() } catch { XtXInv = Matrix.identity(p) }

  const hat = X.multiply(XtXInv).multiply(Xt)
  const leverage = Array.from({ length: n }, (_, i) => hat.get(i, i))

  // MSE
  const mse = result.residuals.reduce((s, r) => s + r * r, 0) / (n - p)

  // Standardized residuals
  const standardizedResiduals = result.residuals.map((r, i) => {
    const denom = Math.sqrt(mse * (1 - (leverage[i] ?? 0)))
    return denom === 0 ? 0 : r / denom
  })

  // Cook's distance
  const cooksDistance = result.residuals.map((r, i) => {
    const h = leverage[i] ?? 0
    return (r * r * h) / (p * mse * (1 - h) ** 2)
  })

  // VIF: 1 / (1 - R²_j) for each predictor
  const vif = predictors.map((_, j) => {
    const otherPreds = predictors.filter((__, k) => k !== j)
    if (otherPreds.length === 0) return 1

    const xj = predictors[j]!.values
    const others = otherPreds.map(p => ({ name: p.name, values: p.values }))
    try {
      const res = multipleRegression(xj, others)
      return 1 / Math.max(1e-10, 1 - res.r2)
    } catch {
      return NaN
    }
  })

  return { leverage, cooksDistance, standardizedResiduals, vif }
}

