/**
 * Regression analysis module.
 * Simple and multiple OLS, logistic regression, polynomial regression,
 * diagnostics (R², AIC, BIC, VIF, residual plots).
 */

import { mean as _mean, variance as _variance, fDistPValue, tDistPValue, tDistQuantile, normalQuantile, logGamma, digamma, trigamma, roundTo } from '../core/math.js'
import { formatRegression, formatPoisson, formatNegBin } from '../core/apa.js'
import { Matrix } from '../core/matrix.js'
import type { RegressionResult, RegressionCoef, OrdinalRegressionResult } from '../core/types.js'

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

/**
 * Compute residual sum of squares from a design matrix and response.
 * RSS = y'y - y'X(X'X)⁻¹X'y = ||y - X·β̂||²
 * Used by two-way ANOVA and ANCOVA for Type II/III SS via model comparison.
 */
export function computeRSS(X: Matrix, y: readonly number[]): number {
  const Xt = X.transpose()
  const XtX = Xt.multiply(X)
  const XtY = Xt.multiply(Matrix.colVec(y))
  const beta = XtX.inverse().multiply(XtY)
  const n = y.length
  let rss = 0
  for (let i = 0; i < n; i++) {
    let fitted = 0
    for (let j = 0; j < X.cols; j++) fitted += X.get(i, j) * beta.get(j, 0)
    rss += (y[i]! - fitted) ** 2
  }
  return rss
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

// ─── Poisson regression ──────────────────────────────────────────────────

/**
 * Poisson regression via IRLS (iteratively reweighted least squares).
 * Link: log(μ) = Xβ,  Variance: V(μ) = μ
 * Outcome y must be non-negative (ideally integer counts).
 *
 * Implements Fisher scoring matching R's glm(family = poisson):
 *   Working response: z = η + (y - μ) / μ
 *   Weights: W = diag(μ)
 *   Solve WLS: (X'WX) β_new = X'Wz
 *
 * Cross-validated with R:
 * > glm(y ~ x1 + x2, family = poisson, data = df)
 * > coef(mod); summary(mod)$coefficients; AIC(mod); deviance(mod)
 */
export function poissonRegression(
  y: readonly number[],
  predictors: ReadonlyArray<{ name: string; values: readonly number[] }>,
  ciLevel = 0.95,
  maxIter = 100,
  tol = 1e-8
): RegressionResult {
  const n = y.length
  for (let i = 0; i < n; i++) {
    if ((y[i] ?? 0) < 0) throw new Error('poissonRegression: y must be non-negative')
  }
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

  // Helper: exp with overflow protection
  const safeExp = (x: number): number => Math.exp(Math.min(700, Math.max(-700, x)))
  const EPS_MU = 1e-10

  // Helper: compute Poisson deviance = 2·Σ[y·log(y/μ) − (y − μ)]
  // Convention: 0·log(0) = 0
  const computeDeviance = (yArr: readonly number[], muArr: number[]): number => {
    let dev = 0
    for (let i = 0; i < n; i++) {
      const yi = yArr[i] ?? 0
      const mi = muArr[i]!
      if (yi > 0) dev += yi * Math.log(yi / Math.max(EPS_MU, mi))
      dev -= (yi - mi)
    }
    return 2 * dev
  }

  // IRLS with step-halving
  // Initialize: intercept = log(mean(y)), other coefficients = 0
  let beta = new Array<number>(p).fill(0)
  const yMeanRaw = _mean([...y])
  const yMeanSafe = Math.max(EPS_MU, yMeanRaw)
  beta[0] = Math.log(yMeanSafe)

  let prevDeviance = Infinity

  for (let iter = 0; iter < maxIter; iter++) {
    const eta = computeEta(beta)
    const mu = eta.map(e => Math.max(EPS_MU, safeExp(e)))
    const dev = computeDeviance(y, mu)

    // Convergence check (R's formula)
    if (iter > 0 && Math.abs(dev - prevDeviance) / (0.1 + Math.abs(dev)) < tol) break
    prevDeviance = dev

    // Working response: z_i = η_i + (y_i − μ_i) / μ_i
    // Weights: w_i = μ_i
    const z = Array.from({ length: n }, (_, i) => eta[i]! + ((y[i] ?? 0) - mu[i]!) / mu[i]!)
    const w = mu.map(m => Math.max(EPS_MU, m))

    // Weighted LS: (X'WX)β = X'Wz via √W transformation
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

      // Step-halving: if deviance increases, halve step
      let accepted = false
      let stepScale = 1.0
      for (let half = 0; half < 10; half++) {
        const betaTry = beta.map((b, j) => b + stepScale * (betaNew[j]! - b))
        const etaTry = computeEta(betaTry)
        const muTry = etaTry.map(e => Math.max(EPS_MU, safeExp(e)))
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
  const mu = eta.map(e => Math.max(EPS_MU, safeExp(e)))

  // Fisher information: (X'WX)⁻¹  where W = diag(μ)
  const sqrtW = mu.map(m => Math.sqrt(Math.max(EPS_MU, m)))
  const Xw = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_, j) => X.get(i, j) * sqrtW[i]!))
  )
  let covBeta: Matrix
  try {
    covBeta = Xw.transpose().multiply(Xw).inverse()
  } catch {
    covBeta = Matrix.identity(p)
  }

  // Wald z-tests (normal approximation)
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
      tValue: roundTo(zStat, 6),  // actually z, but same field
      pValue: roundTo(pVal, 6),
      ci,
    }
  })

  // Log-likelihood: Σ[y·log(μ) − μ − log(y!)]
  let logLik = 0
  for (let i = 0; i < n; i++) {
    const yi = y[i] ?? 0
    const mi = mu[i]!
    logLik += yi * Math.log(Math.max(EPS_MU, mi)) - mi - logGamma(yi + 1)
  }

  // Null log-likelihood (intercept-only model: μ = mean(y))
  let nullLogLik = 0
  for (let i = 0; i < n; i++) {
    const yi = y[i] ?? 0
    nullLogLik += yi * Math.log(Math.max(EPS_MU, yMeanSafe)) - yMeanSafe - logGamma(yi + 1)
  }

  // Deviance and null deviance
  const deviance = computeDeviance(y, mu)
  const nullDeviance = computeDeviance(y, new Array(n).fill(yMeanSafe))

  // McFadden pseudo-R² = 1 − deviance/null_deviance
  const r2 = nullDeviance > 1e-12 ? 1 - deviance / nullDeviance : 0

  const aic = -2 * logLik + 2 * p
  const bic = -2 * logLik + Math.log(n) * p
  const residuals = y.map((v, i) => (v ?? 0) - (mu[i] ?? 0))

  return {
    coefficients,
    r2: roundTo(r2, 6),
    adjR2: roundTo(r2, 6),  // McFadden pseudo-R² for Poisson
    fStatistic: NaN,
    fDf: [p - 1, n - p],
    fPValue: NaN,
    aic: roundTo(aic, 2),
    bic: roundTo(bic, 2),
    residuals,
    fitted: mu,
    n,
    formatted: formatPoisson(deviance, nullDeviance, aic),
  }
}

// ─── Quasi-Poisson regression ─────────────────────────────────────────────

/**
 * Quasi-Poisson regression.
 * Thin wrapper around poissonRegression that estimates the dispersion parameter
 * phi = Pearson chi^2/(n-p) and scales standard errors by sqrt(phi).
 * AIC/BIC are NaN (not defined for quasi-likelihood).
 *
 * Cross-validated with R:
 * > glm(y ~ x, family = quasipoisson, data = df)
 * > summary(mod)$dispersion
 */
export function quasiPoissonRegression(
  y: readonly number[],
  predictors: ReadonlyArray<{ name: string; values: readonly number[] }>,
  ciLevel = 0.95,
  maxIter = 100,
  tol = 1e-8
): RegressionResult & { readonly dispersion: number } {
  const base = poissonRegression(y, predictors, ciLevel, maxIter, tol)
  const n = y.length
  const p = predictors.length + 1

  // Pearson chi^2 = sum (y_i - mu_i)^2 / mu_i
  let pearsonChi2 = 0
  for (let i = 0; i < n; i++) {
    const yi = y[i]!
    const mi = base.fitted[i]!
    pearsonChi2 += (yi - mi) ** 2 / Math.max(1e-10, mi)
  }
  const dispersion = pearsonChi2 / (n - p)

  // Scale SEs by sqrt(phi), recompute z and p
  const sqrtPhi = Math.sqrt(dispersion)
  const zCrit = normalQuantile(1 - (1 - ciLevel) / 2)
  const coefficients: RegressionCoef[] = base.coefficients.map(c => {
    const newSE = c.se * sqrtPhi
    const newZ = newSE === 0 ? 0 : c.estimate / newSE
    const newP = 2 * (1 - normCDFLocal(Math.abs(newZ)))
    const newCI: readonly [number, number] = [c.estimate - zCrit * newSE, c.estimate + zCrit * newSE]
    return { ...c, se: roundTo(newSE, 10), tValue: roundTo(newZ, 6), pValue: roundTo(newP, 6), ci: newCI }
  })

  return {
    ...base,
    coefficients,
    aic: NaN,  // not defined for quasi-likelihood
    bic: NaN,
    dispersion,
    formatted: `Quasi-Poisson (phi = ${roundTo(dispersion, 2)}), ${base.formatted}`,
  }
}

// ─── Negative binomial regression ─────────────────────────────────────────

/**
 * Negative binomial regression via IRLS.
 * Link: log(mu) = X*beta, Variance: V(mu) = mu + mu^2/theta
 * Uses outer loop for theta: method-of-moments init, Newton steps on profile logLik.
 *
 * Cross-validated with R:
 * > library(MASS)
 * > glm.nb(y ~ x1 + x2, data = df)
 */
export function negativeBinomialRegression(
  y: readonly number[],
  predictors: ReadonlyArray<{ name: string; values: readonly number[] }>,
  ciLevel = 0.95,
  maxIter = 100,
  tol = 1e-8
): RegressionResult & { readonly theta: number } {
  const n = y.length
  for (let i = 0; i < n; i++) {
    if ((y[i] ?? 0) < 0) throw new Error('negativeBinomialRegression: y must be non-negative')
  }
  const p = predictors.length + 1
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map(pr => pr.values[i] ?? 0)])
  )
  const names = ['(Intercept)', ...predictors.map(pr => pr.name)]

  const computeEta = (b: number[]): number[] =>
    Array.from({ length: n }, (_, i) => {
      let v = 0
      for (let j = 0; j < p; j++) v += X.get(i, j) * (b[j] ?? 0)
      return v
    })

  const safeExp = (x: number): number => Math.exp(Math.min(700, Math.max(-700, x)))
  const EPS_MU = 1e-10

  // Initialize beta from Poisson fit
  const poisFit = poissonRegression(y, predictors, ciLevel, maxIter, tol)
  let beta = poisFit.coefficients.map(c => c.estimate)
  let eta = computeEta(beta)
  let mu = eta.map(e => Math.max(EPS_MU, safeExp(e)))

  // Initialize theta via method of moments: Var(Y) = mu + mu^2/theta -> theta = mu^2/(Var-mu)
  const yMean = y.reduce((s, v) => s + v, 0) / n
  const yVar = y.reduce((s, v) => s + (v - yMean) ** 2, 0) / (n - 1)
  let theta = yVar > yMean ? yMean * yMean / (yVar - yMean) : 10

  // Outer loop for theta
  for (let outer = 0; outer < 25; outer++) {
    // Inner IRLS with current theta: V(mu) = mu + mu^2/theta, w = mu^2/(mu+mu^2/theta) = mu*theta/(theta+mu)
    for (let iter = 0; iter < maxIter; iter++) {
      eta = computeEta(beta)
      mu = eta.map(e => Math.max(EPS_MU, safeExp(e)))
      const w = mu.map(m => m * theta / (theta + m))
      const z = Array.from({ length: n }, (_, i) => eta[i]! + (y[i]! - mu[i]!) / mu[i]!)

      const sqrtW = w.map(wi => Math.sqrt(Math.max(EPS_MU, wi)))
      const Xw = Matrix.fromArray(
        Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_, j) => X.get(i, j) * sqrtW[i]!))
      )
      const zw = z.map((zi, i) => zi * sqrtW[i]!)

      try {
        const Xwt = Xw.transpose()
        const betaNewM = Xwt.multiply(Xw).inverse().multiply(Xwt.multiply(Matrix.colVec(zw)))
        const betaNew = Array.from({ length: p }, (_, j) => betaNewM.get(j, 0))

        const maxDelta = Math.max(...betaNew.map((b, i) => Math.abs(b - beta[i]!)))
        beta = betaNew
        if (maxDelta < tol) break
      } catch { break }
    }

    // Update theta via Newton step on profile log-likelihood
    // score = d(logLik)/d(theta), info = d^2(logLik)/d(theta)^2
    eta = computeEta(beta)
    mu = eta.map(e => Math.max(EPS_MU, safeExp(e)))

    let score = 0, info = 0
    for (let i = 0; i < n; i++) {
      const yi = y[i]!
      const mi = mu[i]!
      score += digamma(yi + theta) - digamma(theta) + Math.log(theta) + 1 - Math.log(theta + mi) - (yi + theta) / (theta + mi)
      info += trigamma(yi + theta) - trigamma(theta) + 1 / theta - 2 / (theta + mi) + (yi + theta) / ((theta + mi) * (theta + mi))
    }

    if (Math.abs(info) > EPS_MU) {
      const thetaNew = theta - score / info
      if (thetaNew > 0.001) {
        const thetaDelta = Math.abs(thetaNew - theta)
        theta = thetaNew
        if (thetaDelta < tol) break
      }
    }
  }

  // Final quantities
  eta = computeEta(beta)
  mu = eta.map(e => Math.max(EPS_MU, safeExp(e)))
  const w = mu.map(m => m * theta / (theta + m))
  const sqrtW = w.map(wi => Math.sqrt(Math.max(EPS_MU, wi)))
  const Xw = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_, j) => X.get(i, j) * sqrtW[i]!))
  )
  let covBeta: Matrix
  try { covBeta = Xw.transpose().multiply(Xw).inverse() } catch { covBeta = Matrix.identity(p) }

  const zCrit = normalQuantile(1 - (1 - ciLevel) / 2)
  const coefficients: RegressionCoef[] = beta.map((b, i) => {
    const se = Math.sqrt(Math.max(0, covBeta.get(i, i)))
    const zStat = se === 0 ? 0 : b / se
    const pVal = 2 * (1 - normCDFLocal(Math.abs(zStat)))
    const ci: readonly [number, number] = [b - zCrit * se, b + zCrit * se]
    return { name: names[i] ?? `b${i}`, estimate: roundTo(b, 10), se: roundTo(se, 10), tValue: roundTo(zStat, 6), pValue: roundTo(pVal, 6), ci }
  })

  // Log-likelihood
  let logLik = 0
  for (let i = 0; i < n; i++) {
    const yi = y[i]!
    const mi = mu[i]!
    logLik += logGamma(yi + theta) - logGamma(theta) - logGamma(yi + 1) + theta * Math.log(theta / (theta + mi)) + yi * Math.log(mi / (theta + mi))
  }

  // Deviance: 2 * sum[y*log(y/mu) - (y+theta)*log((y+theta)/(mu+theta))]
  let deviance = 0
  for (let i = 0; i < n; i++) {
    const yi = y[i]!
    const mi = mu[i]!
    deviance += yi > 0 ? yi * Math.log(yi / mi) : 0
    deviance -= (yi + theta) * Math.log((yi + theta) / (mi + theta))
  }
  deviance *= 2

  const aic = -2 * logLik + 2 * (p + 1)  // +1 for theta
  const bic = -2 * logLik + Math.log(n) * (p + 1)
  const residuals = y.map((v, i) => v - mu[i]!)

  return {
    coefficients,
    r2: NaN, adjR2: NaN,
    fStatistic: NaN, fDf: [p - 1, n - p], fPValue: NaN,
    aic: roundTo(aic, 2), bic: roundTo(bic, 2),
    residuals, fitted: mu, n,
    formatted: formatNegBin(deviance, theta, aic),
    theta,
  }
}

// ─── Ordinal logistic regression ──────────────────────────────────────────

/**
 * Ordinal logistic regression (proportional odds / cumulative logit model).
 * P(Y <= j | x) = logistic(alpha_j - x'beta), j = 1,...,J-1
 * Uses Fisher scoring (Newton-Raphson on the logLik).
 *
 * Cross-validated with R:
 * > library(MASS)
 * > polr(factor(y) ~ x1 + x2, method = "logistic", data = df)
 */
export function ordinalLogisticRegression(
  y: readonly number[],  // integer categories 1..J
  predictors: ReadonlyArray<{ name: string; values: readonly number[] }>,
  ciLevel = 0.95,
  maxIter = 100,
  tol = 1e-8
): OrdinalRegressionResult {
  const n = y.length
  if (n === 0) throw new Error('ordinalLogisticRegression: empty data')
  for (const pr of predictors) {
    if (pr.values.length !== n) throw new Error(`ordinalLogisticRegression: predictor '${pr.name}' length mismatch`)
  }

  // Get unique sorted categories
  const categories = [...new Set(y)].sort((a, b) => a - b)
  const J = categories.length
  if (J < 2) throw new Error('ordinalLogisticRegression: need at least 2 categories')
  const K = J - 1  // number of thresholds

  // Map y to category indices 0..J-1
  const catMap = new Map<number, number>()
  for (let j = 0; j < J; j++) catMap.set(categories[j]!, j)
  const yIdx = y.map(v => catMap.get(v)!)

  const nPred = predictors.length
  const nParams = K + nPred  // thresholds + coefficients

  // Build predictor matrix (no intercept column - thresholds act as intercepts)
  const xMat: number[][] = Array.from({ length: n }, (_, i) =>
    predictors.map(pr => pr.values[i] ?? 0)
  )

  // Logistic function with clamping
  const logistic = (e: number): number => {
    const v = 1 / (1 + Math.exp(-Math.min(700, Math.max(-700, e))))
    return Math.min(1 - 1e-15, Math.max(1e-15, v))
  }

  // Initialize thresholds from empirical cumulative logits
  const cumCounts = new Array<number>(J).fill(0)
  for (const idx of yIdx) cumCounts[idx]!++
  const thresholds = new Array<number>(K)
  let cumSum = 0
  for (let j = 0; j < K; j++) {
    cumSum += cumCounts[j]!
    const p = Math.min(1 - 1e-6, Math.max(1e-6, cumSum / n))
    thresholds[j] = Math.log(p / (1 - p))
  }

  // Initialize coefficients to 0
  let betaCoefs = new Array<number>(nPred).fill(0)

  // Pack/unpack parameter vector: [alpha_1, ..., alpha_{K}, beta_1, ..., beta_p]
  const pack = (): number[] => [...thresholds, ...betaCoefs]
  const unpack = (params: number[]): void => {
    for (let j = 0; j < K; j++) thresholds[j] = params[j]!
    for (let j = 0; j < nPred; j++) betaCoefs[j] = params[K + j]!
  }

  // Compute cumulative probabilities P(Y <= j | x_i) for each observation and category
  const computeGammas = (): number[][] => {
    // gammas[i][j] = P(Y <= j | x_i), j = 0..K-1
    // P(Y <= j) = logistic(alpha_j - x'beta)
    return Array.from({ length: n }, (_, i) => {
      let xBeta = 0
      for (let k = 0; k < nPred; k++) xBeta += xMat[i]![k]! * betaCoefs[k]!
      return Array.from({ length: K }, (_, j) => logistic(thresholds[j]! - xBeta))
    })
  }

  // Compute category probabilities P(Y = j | x_i) from cumulative
  const computeProbs = (gammas: number[][]): number[][] => {
    return Array.from({ length: n }, (_, i) => {
      return Array.from({ length: J }, (_, j) => {
        let prob: number
        if (j === 0) {
          prob = gammas[i]![0]!
        } else if (j === J - 1) {
          prob = 1 - gammas[i]![K - 1]!
        } else {
          prob = gammas[i]![j]! - gammas[i]![j - 1]!
        }
        return Math.max(1e-15, prob)
      })
    })
  }

  // Fisher scoring iterations
  for (let iter = 0; iter < maxIter; iter++) {
    const gammas = computeGammas()
    const probs = computeProbs(gammas)

    // Compute gradient and Hessian
    const grad = new Array<number>(nParams).fill(0)
    const hess: number[][] = Array.from({ length: nParams }, () => new Array<number>(nParams).fill(0))

    for (let i = 0; i < n; i++) {
      const yi = yIdx[i]!
      const pi = probs[i]!
      const gi = gammas[i]!

      // For each threshold alpha_j:
      // d(logLik_i)/d(alpha_j) depends on whether yi == j, yi == j+1, or neither
      // gamma_j * (1 - gamma_j) is the derivative of logistic
      for (let j = 0; j < K; j++) {
        const gj = gi[j]!
        const dgj = gj * (1 - gj)  // d(gamma_j)/d(alpha_j)

        // d(log p_yi)/d(alpha_j)
        let dLogP = 0
        if (yi === j) {
          // p_j = gamma_j - gamma_{j-1}, d(p_j)/d(alpha_j) = dgj
          dLogP = dgj / pi[yi]!
        } else if (yi === j + 1) {
          // p_{j+1} = gamma_{j+1} - gamma_j, d(p_{j+1})/d(alpha_j) = -dgj
          dLogP = -dgj / pi[yi]!
        }
        // All other cases: d(p_yi)/d(alpha_j) = 0

        grad[j]! += dLogP

        // For the Hessian (expected Fisher information):
        // E[-d^2 logLik / d(alpha_j) d(alpha_k)] contributions
        // Use expected information: I_{jk} = sum_i sum_m (1/p_im) * dp_im/d(theta_j) * dp_im/d(theta_k)
      }

      // For coefficients beta_k:
      // d(gamma_j)/d(beta_k) = -gamma_j*(1-gamma_j)*x_{ik}
      for (let k = 0; k < nPred; k++) {
        const xik = xMat[i]![k]!
        let dLogP = 0
        if (yi === 0) {
          // p_0 = gamma_0, dp_0/dbeta_k = -dg0 * xik
          const dg0 = gi[0]! * (1 - gi[0]!)
          dLogP = -dg0 * xik / pi[0]!
        } else if (yi === J - 1) {
          // p_{J-1} = 1 - gamma_{K-1}, dp/dbeta_k = dg_{K-1} * xik
          const dgK = gi[K - 1]! * (1 - gi[K - 1]!)
          dLogP = dgK * xik / pi[J - 1]!
        } else {
          // p_yi = gamma_yi - gamma_{yi-1}
          const dgY = gi[yi]! * (1 - gi[yi]!)
          const dgY1 = gi[yi - 1]! * (1 - gi[yi - 1]!)
          dLogP = (-dgY + dgY1) * xik / pi[yi]!
        }
        grad[K + k]! += dLogP
      }

      // Expected Fisher information (using observed per-category contributions)
      // I_{ab} = sum_i sum_m (1/p_im) * dp_im/da * dp_im/db
      // dp_im/d(alpha_j): nonzero only for m=j (positive) and m=j+1 (negative)
      // dp_im/d(beta_k): nonzero for m=0 (neg), m in [1..K-1] (diff), m=J-1 (pos)

      // Precompute dp/d(param) for each category m
      const dpda: number[][] = Array.from({ length: J }, (_, m) => {
        const result = new Array<number>(nParams).fill(0)
        // Threshold derivatives
        if (m === 0) {
          // dp_0/d(alpha_0) = dg_0
          result[0] = gi[0]! * (1 - gi[0]!)
        } else if (m === J - 1) {
          // dp_{J-1}/d(alpha_{K-1}) = -dg_{K-1}
          result[K - 1] = -(gi[K - 1]! * (1 - gi[K - 1]!))
        } else {
          // dp_m/d(alpha_m) = dg_m
          result[m] = gi[m]! * (1 - gi[m]!)
          // dp_m/d(alpha_{m-1}) = -dg_{m-1}
          result[m - 1] = -(gi[m - 1]! * (1 - gi[m - 1]!))
        }
        // Beta derivatives: dp_m/d(beta_k) = -x_ik * [dg_m (if m<J-1) - dg_{m-1} (if m>0)]
        for (let k = 0; k < nPred; k++) {
          const xik = xMat[i]![k]!
          let val = 0
          if (m < J - 1) val -= gi[m]! * (1 - gi[m]!) * xik
          if (m > 0) val += gi[m - 1]! * (1 - gi[m - 1]!) * xik
          result[K + k] = val
        }
        return result
      })

      // Accumulate Fisher information
      for (let m = 0; m < J; m++) {
        const invP = 1 / pi[m]!
        for (let a = 0; a < nParams; a++) {
          for (let b = a; b < nParams; b++) {
            const val = invP * dpda[m]![a]! * dpda[m]![b]!
            hess[a]![b]! += val
            if (b !== a) hess[b]![a]! += val
          }
        }
      }
    }

    // Newton step: delta = H^{-1} * g
    const H = Matrix.fromArray(hess)
    const g = Matrix.colVec(grad)
    let delta: number[]
    try {
      const deltaM = H.inverse().multiply(g)
      delta = Array.from({ length: nParams }, (_, i) => deltaM.get(i, 0))
    } catch {
      break  // singular Hessian
    }

    // Update parameters
    const params = pack()
    const newParams = params.map((p, i) => p + delta[i]!)
    unpack(newParams)

    // Enforce threshold ordering: alpha_1 < alpha_2 < ... < alpha_{K}
    for (let j = 1; j < K; j++) {
      if (thresholds[j]! <= thresholds[j - 1]!) {
        thresholds[j] = thresholds[j - 1]! + 0.01
      }
    }

    // Check convergence
    const maxDelta = Math.max(...delta.map(Math.abs))
    if (maxDelta < tol) break
  }

  // Final quantities
  const gammas = computeGammas()
  const probs = computeProbs(gammas)

  // Log-likelihood
  let logLik = 0
  for (let i = 0; i < n; i++) {
    logLik += Math.log(probs[i]![yIdx[i]!]!)
  }

  // Compute final Hessian for SEs
  const finalHess: number[][] = Array.from({ length: nParams }, () => new Array<number>(nParams).fill(0))
  for (let i = 0; i < n; i++) {
    const pi = probs[i]!
    const gi = gammas[i]!

    const dpda: number[][] = Array.from({ length: J }, (_, m) => {
      const result = new Array<number>(nParams).fill(0)
      if (m === 0) {
        result[0] = gi[0]! * (1 - gi[0]!)
      } else if (m === J - 1) {
        result[K - 1] = -(gi[K - 1]! * (1 - gi[K - 1]!))
      } else {
        result[m] = gi[m]! * (1 - gi[m]!)
        result[m - 1] = -(gi[m - 1]! * (1 - gi[m - 1]!))
      }
      for (let k = 0; k < nPred; k++) {
        const xik = xMat[i]![k]!
        let val = 0
        if (m < J - 1) val -= gi[m]! * (1 - gi[m]!) * xik
        if (m > 0) val += gi[m - 1]! * (1 - gi[m - 1]!) * xik
        result[K + k] = val
      }
      return result
    })

    for (let m = 0; m < J; m++) {
      const invP = 1 / pi[m]!
      for (let a = 0; a < nParams; a++) {
        for (let b = a; b < nParams; b++) {
          const val = invP * dpda[m]![a]! * dpda[m]![b]!
          finalHess[a]![b]! += val
          if (b !== a) finalHess[b]![a]! += val
        }
      }
    }
  }

  let covMatrix: Matrix
  try { covMatrix = Matrix.fromArray(finalHess).inverse() } catch { covMatrix = Matrix.identity(nParams) }

  const zCritVal = normalQuantile(1 - (1 - ciLevel) / 2)

  // Threshold results
  const thresholdResults = thresholds.map((alpha, j) => {
    const se = Math.sqrt(Math.max(0, covMatrix.get(j, j)))
    const z = se === 0 ? 0 : alpha / se
    const pVal = 2 * (1 - normCDFLocal(Math.abs(z)))
    return {
      name: `${categories[j]!}|${categories[j + 1]!}`,
      estimate: roundTo(alpha, 6),
      se: roundTo(se, 6),
      z: roundTo(z, 4),
      pValue: roundTo(pVal, 6),
    }
  })

  // Coefficient results
  const coefficientResults = betaCoefs.map((b, k) => {
    const se = Math.sqrt(Math.max(0, covMatrix.get(K + k, K + k)))
    const z = se === 0 ? 0 : b / se
    const pVal = 2 * (1 - normCDFLocal(Math.abs(z)))
    const ciR: readonly [number, number] = [b - zCritVal * se, b + zCritVal * se]
    return {
      name: predictors[k]!.name,
      estimate: roundTo(b, 6),
      se: roundTo(se, 6),
      z: roundTo(z, 4),
      pValue: roundTo(pVal, 6),
      ci: ciR,
    }
  })

  const aic = -2 * logLik + 2 * nParams
  const bic = -2 * logLik + Math.log(n) * nParams

  const formatted = `Ordinal logistic: logLik = ${roundTo(logLik, 1)}, AIC = ${roundTo(aic, 1)}, BIC = ${roundTo(bic, 1)}, J = ${J}`

  return {
    thresholds: thresholdResults,
    coefficients: coefficientResults,
    logLik: roundTo(logLik, 4),
    aic: roundTo(aic, 2),
    bic: roundTo(bic, 2),
    n,
    nCategories: J,
    formatted,
  }
}

