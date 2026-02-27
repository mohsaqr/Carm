/**
 * Regression analysis module.
 * Simple and multiple OLS, logistic regression, polynomial regression,
 * diagnostics (R², AIC, BIC, VIF, residual plots).
 */

import { mean as _mean, variance as _variance, fDistPValue, tDistPValue, tDistQuantile, normalQuantile, roundTo } from '../core/math.js'
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
  // Clamp RSS to avoid log(0) = -Inf on perfect fit
  const rssSafe = Math.max(ss_res, 1e-15)
  const logLik = -n / 2 * (Math.log(2 * Math.PI) + Math.log(rssSafe / n) + 1)
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
 * Implements Fisher scoring matching R's glm.fit exactly:
 *   Working response: z = η + (y - μ) / w   where w = μ(1-μ)
 *   Solve WLS: (X'WX) β_new = X'Wz
 *   Equivalently: (Xw'Xw) β_new = Xw' zw   with Xw = √W·X, zw = √W·z
 *
 * Cross-validated with R:
 * > glm(y ~ x1 + x2, family = binomial, data = df)
 *
 * Cross-validated with Python:
 * > import statsmodels.api as sm
 * > sm.GLM(y, sm.add_constant(X), family=sm.families.Binomial()).fit()
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

  // Helper: compute η = Xβ
  const computeEta = (b: number[]): number[] =>
    Array.from({ length: n }, (_, i) => {
      let v = 0
      for (let j = 0; j < p; j++) v += X.get(i, j) * (b[j] ?? 0)
      return v
    })

  // Helper: logistic with eta clamping to prevent overflow
  // Clamp output to (eps, 1-eps) matching R's family$linkinv boundary behavior
  const EPS_MU = 1e-15
  const logistic = (e: number): number => {
    const mu = 1 / (1 + Math.exp(-Math.min(700, Math.max(-700, e))))
    return Math.min(1 - EPS_MU, Math.max(EPS_MU, mu))
  }

  // Helper: compute deviance = -2 * logLik
  const computeDeviance = (yArr: readonly number[], muArr: number[]): number => {
    let dev = 0
    for (let i = 0; i < n; i++) {
      const yi = yArr[i] ?? 0
      const mi = muArr[i]!
      dev += -2 * (yi * Math.log(Math.max(1e-15, mi)) + (1 - yi) * Math.log(Math.max(1e-15, 1 - mi)))
    }
    return dev
  }

  // IRLS with step-halving (matching R's glm.fit Fisher scoring)
  // Initialize: intercept = log(p̄/(1-p̄)), other coefficients = 0
  let beta = new Array<number>(p).fill(0)
  const pMeanInit = Math.min(1 - 1e-6, Math.max(1e-6, _mean([...y])))
  beta[0] = Math.log(pMeanInit / (1 - pMeanInit))

  let prevDeviance = Infinity

  for (let iter = 0; iter < maxIter; iter++) {
    // η = Xβ, μ = logistic(η), w = μ(1-μ)
    const eta = computeEta(beta)
    const mu = eta.map(logistic)
    const w = mu.map(m => Math.max(1e-10, m * (1 - m)))

    // Deviance = -2 * logLik
    const dev = computeDeviance(y, mu)

    // Convergence: |dev - devold| / (0.1 + |dev|) < tol  (R's formula)
    if (iter > 0 && Math.abs(dev - prevDeviance) / (0.1 + Math.abs(dev)) < tol) break
    prevDeviance = dev

    // Working response: z_i = η_i + (y_i - μ_i) / w_i
    // This is the linearized pseudo-observation from the Fisher scoring derivation.
    // R's glm.fit: z <- eta + (y - mu) / mu.eta.val  where mu.eta.val = dmu/deta = mu*(1-mu)
    const z = Array.from({ length: n }, (_, i) => eta[i]! + ((y[i] ?? 0) - mu[i]!) / w[i]!)

    // Weighted system: Xw = √W·X, zw = √W·z
    // Solve (Xw'Xw) β_new = Xw' zw  ⟹  β_new = (X'WX)⁻¹ X'Wz
    const sqrtW = w.map(Math.sqrt)
    const Xw = Matrix.fromArray(
      Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_, j) => X.get(i, j) * sqrtW[i]!))
    )
    const zw = z.map((zi, i) => zi * sqrtW[i]!)

    try {
      const Xwt = Xw.transpose()
      const XwtXw = Xwt.multiply(Xw)
      const XwtZw = Xwt.multiply(Matrix.colVec(zw))
      const betaNewM = XwtXw.inverse().multiply(XwtZw)
      const betaNew = Array.from({ length: p }, (_, j) => betaNewM.get(j, 0))

      // Step-halving: if deviance increases or is non-finite, halve toward betaNew
      let accepted = false
      let stepScale = 1.0
      for (let half = 0; half < 10; half++) {
        const betaTry = beta.map((b, j) => b + stepScale * (betaNew[j]! - b))

        // Compute trial deviance
        const etaTry = computeEta(betaTry)
        const muTry = etaTry.map(logistic)
        const trialDev = computeDeviance(y, muTry)

        if (isFinite(trialDev) && (trialDev <= dev + 1e-8 || half === 9)) {
          beta = betaTry
          accepted = true
          break
        }
        stepScale *= 0.5
      }
      if (!accepted) break
    } catch {
      break  // singular — stop iterating
    }
  }

  // ── Final quantities at converged β ──────────────────────────────────
  const eta = computeEta(beta)
  const mu = eta.map(logistic)  // clamped via logistic helper
  const w = mu.map(m => Math.max(1e-10, m * (1 - m)))

  // Fisher information: (X'WX)⁻¹ = covariance of β̂
  const sqrtW = w.map(Math.sqrt)
  const Xw = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_, j) => X.get(i, j) * sqrtW[i]!))
  )
  let covBeta: Matrix
  try {
    covBeta = Xw.transpose().multiply(Xw).inverse()
  } catch {
    covBeta = Matrix.identity(p)
  }

  // CIs use normal quantile (MLE asymptotics), not t-distribution
  const zCrit = normalQuantile(1 - (1 - ciLevel) / 2)
  const coefficients: RegressionCoef[] = beta.map((b, i) => {
    const se = Math.sqrt(Math.max(0, covBeta.get(i, i)))
    const zStat = se === 0 ? 0 : b / se
    const pVal = 2 * (1 - normCDFLocal(Math.abs(zStat)))
    const ci: readonly [number, number] = [b - zCrit * se, b + zCrit * se]
    return {
      name: names[i] ?? `β${i}`,
      estimate: roundTo(b, 10),
      se: roundTo(se, 10),
      tValue: roundTo(zStat, 6),
      pValue: roundTo(pVal, 6),
      ci,
    }
  })

  // Log-likelihood at converged estimates
  const logLik = mu.reduce((s, m, i) => {
    const yi = y[i] ?? 0
    return s + yi * Math.log(Math.max(1e-15, m)) + (1 - yi) * Math.log(Math.max(1e-15, 1 - m))
  }, 0)

  // Null log-likelihood (intercept only)
  const pMeanRaw = _mean([...y])
  const pMean = Math.min(1 - 1e-12, Math.max(1e-12, pMeanRaw))
  const nullLogLik = n * (pMean * Math.log(Math.max(1e-15, pMean)) + (1 - pMean) * Math.log(Math.max(1e-15, 1 - pMean)))

  // McFadden pseudo-R²
  const r2 = Math.abs(nullLogLik) < 1e-12 ? NaN : 1 - logLik / nullLogLik

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

