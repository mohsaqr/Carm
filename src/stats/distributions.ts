/**
 * Distribution Fitting Module — PDFs, CDFs, quantiles, MLE, and goodness-of-fit.
 *
 * Exports ~38 functions covering:
 *   - Continuous PDFs: normal, t, chi-square, F, exponential, uniform, gamma, beta
 *   - CDFs & quantiles for distributions NOT in core/math.ts: exponential, uniform,
 *     gamma, beta, F quantile, binomial, Poisson
 *   - MLE fitting via fitDistribution()
 *   - Anderson-Darling and Kolmogorov-Smirnov goodness-of-fit tests
 *
 * Naming convention: core/math.ts already exports normalCDF, chiSqCDF, tDistCDF, etc.
 * This module does NOT re-export those. It only exports genuinely new functions.
 */

import type { StatResult } from '../core/types.js'

// ─── Types ───────────────────────────────────────────────────────────────────

/** Supported distribution names for MLE fitting. */
export type FitDistName = 'normal' | 'exponential' | 'uniform' | 'poisson' | 'gamma' | 'beta'

/** Result from maximum-likelihood distribution fitting. */
export interface FitDistributionResult {
  readonly distribution: FitDistName
  readonly params: Readonly<Record<string, number>>
  readonly logLik: number
  readonly aic: number
  readonly bic: number
  readonly n: number
  readonly formatted: string
}
import {
  logGamma,
  logBeta,
  incompleteBeta,
  incompleteGamma,
  normalCDF as _normalCDF,
  normalQuantile as _normalQuantile,
  fDistCDF as _fDistCDF,
  digamma,
  trigamma,
  mean as _mean,
  sortAsc,
} from '../core/math.js'

// ─── Section A: Continuous PDFs ──────────────────────────────────────────────

/**
 * Normal (Gaussian) PDF.
 * f(x) = (1 / (σ√(2π))) exp(-(x-μ)² / (2σ²))
 */
export function normalPDF(x: number, mu = 0, sigma = 1): number {
  if (sigma <= 0) throw new Error(`normalPDF: sigma must be positive, got ${sigma}`)
  const z = (x - mu) / sigma
  return Math.exp(-0.5 * z * z) / (sigma * Math.sqrt(2 * Math.PI))
}

/**
 * Student's t PDF.
 * f(x | ν) = Γ((ν+1)/2) / (√(νπ) Γ(ν/2)) (1 + x²/ν)^(-(ν+1)/2)
 * Reference: NIST Handbook §6.8
 */
export function tPDF(x: number, df: number): number {
  if (df <= 0) throw new Error(`tPDF: df must be positive, got ${df}`)
  const logCoeff = logGamma((df + 1) / 2) - logGamma(df / 2) - 0.5 * Math.log(df * Math.PI)
  return Math.exp(logCoeff + (-(df + 1) / 2) * Math.log(1 + x * x / df))
}

/**
 * Chi-square PDF.
 * f(x | k) = x^(k/2-1) exp(-x/2) / (2^(k/2) Γ(k/2))
 * Reference: NIST Handbook §6.9
 */
export function chiSqPDF(x: number, df: number): number {
  if (df <= 0) throw new Error(`chiSqPDF: df must be positive, got ${df}`)
  if (x < 0) return 0
  if (x === 0) return df === 2 ? 0.5 : (df < 2 ? Infinity : 0)
  const halfDf = df / 2
  const logP = (halfDf - 1) * Math.log(x) - x / 2 - halfDf * Math.log(2) - logGamma(halfDf)
  return Math.exp(logP)
}

/**
 * F-distribution PDF.
 * f(x | d1, d2) = √((d1·x)^d1 · d2^d2 / (d1·x+d2)^(d1+d2)) / (x · B(d1/2, d2/2))
 * Reference: NIST Handbook §6.10
 */
export function fPDF(x: number, df1: number, df2: number): number {
  if (df1 <= 0 || df2 <= 0) throw new Error(`fPDF: df must be positive, got df1=${df1}, df2=${df2}`)
  if (x < 0) return 0
  if (x === 0) return df1 === 2 ? 1 / Math.exp(logBeta(df1 / 2, df2 / 2)) : (df1 < 2 ? Infinity : 0)
  const logP = (df1 / 2) * Math.log(df1 / df2) + (df1 / 2 - 1) * Math.log(x)
    - ((df1 + df2) / 2) * Math.log(1 + df1 * x / df2) - logBeta(df1 / 2, df2 / 2)
  return Math.exp(logP)
}

/**
 * Exponential PDF.
 * f(x | λ) = λ exp(-λx), x ≥ 0
 */
export function exponentialPDF(x: number, rate = 1): number {
  if (rate <= 0) throw new Error(`exponentialPDF: rate must be positive, got ${rate}`)
  if (x < 0) return 0
  return rate * Math.exp(-rate * x)
}

/**
 * Uniform PDF.
 * f(x | a, b) = 1/(b-a), a ≤ x ≤ b
 */
export function uniformPDF(x: number, min = 0, max = 1): number {
  if (max <= min) throw new Error(`uniformPDF: max must exceed min, got min=${min}, max=${max}`)
  if (x < min || x > max) return 0
  return 1 / (max - min)
}

/**
 * Gamma PDF.
 * f(x | α, β) = β^α / Γ(α) · x^(α-1) · exp(-βx), x > 0
 * Where α = shape, β = rate (= 1/scale).
 * Reference: NIST Handbook §6.14
 */
export function gammaPDF(x: number, shape: number, rate = 1): number {
  if (shape <= 0) throw new Error(`gammaPDF: shape must be positive, got ${shape}`)
  if (rate <= 0) throw new Error(`gammaPDF: rate must be positive, got ${rate}`)
  if (x < 0) return 0
  if (x === 0) return shape === 1 ? rate : (shape < 1 ? Infinity : 0)
  const logP = shape * Math.log(rate) - logGamma(shape) + (shape - 1) * Math.log(x) - rate * x
  return Math.exp(logP)
}

/**
 * Beta PDF.
 * f(x | α, β) = x^(α-1) (1-x)^(β-1) / B(α, β), 0 < x < 1
 * Reference: NIST Handbook §6.16
 */
export function betaPDF(x: number, shape1: number, shape2: number): number {
  if (shape1 <= 0 || shape2 <= 0)
    throw new Error(`betaPDF: shapes must be positive, got shape1=${shape1}, shape2=${shape2}`)
  if (x < 0 || x > 1) return 0
  if (x === 0) return shape1 === 1 ? 1 / Math.exp(logBeta(shape1, shape2)) : (shape1 < 1 ? Infinity : 0)
  if (x === 1) return shape2 === 1 ? 1 / Math.exp(logBeta(shape1, shape2)) : (shape2 < 1 ? Infinity : 0)
  const logP = (shape1 - 1) * Math.log(x) + (shape2 - 1) * Math.log(1 - x) - logBeta(shape1, shape2)
  return Math.exp(logP)
}

// ─── Section B: CDFs and Quantiles (new distributions) ───────────────────────

/**
 * Exponential CDF: P(X ≤ x | λ) = 1 - exp(-λx)
 */
export function exponentialCDF(x: number, rate = 1): number {
  if (rate <= 0) throw new Error(`exponentialCDF: rate must be positive, got ${rate}`)
  if (x < 0) return 0
  return 1 - Math.exp(-rate * x)
}

/**
 * Exponential quantile: Q(p | λ) = -ln(1-p) / λ
 */
export function exponentialQuantile(p: number, rate = 1): number {
  if (rate <= 0) throw new Error(`exponentialQuantile: rate must be positive, got ${rate}`)
  if (p < 0 || p > 1) throw new Error(`exponentialQuantile: p must be in [0,1], got ${p}`)
  if (p === 1) return Infinity
  return -Math.log(1 - p) / rate
}

/**
 * Uniform CDF: P(X ≤ x | a, b) = (x - a) / (b - a), clamped to [0, 1]
 */
export function uniformCDF(x: number, min = 0, max = 1): number {
  if (max <= min) throw new Error(`uniformCDF: max must exceed min, got min=${min}, max=${max}`)
  if (x <= min) return 0
  if (x >= max) return 1
  return (x - min) / (max - min)
}

/**
 * Uniform quantile: Q(p | a, b) = a + p(b - a)
 */
export function uniformQuantile(p: number, min = 0, max = 1): number {
  if (max <= min) throw new Error(`uniformQuantile: max must exceed min, got min=${min}, max=${max}`)
  if (p < 0 || p > 1) throw new Error(`uniformQuantile: p must be in [0,1], got ${p}`)
  return min + p * (max - min)
}

/**
 * Gamma CDF: P(X ≤ x | shape, rate) = P(shape, rate·x) where P is regularized incomplete gamma.
 * Uses: incompleteGamma(shape, rate * x) from core/math.ts
 */
export function gammaCDF(x: number, shape: number, rate = 1): number {
  if (shape <= 0) throw new Error(`gammaCDF: shape must be positive, got ${shape}`)
  if (rate <= 0) throw new Error(`gammaCDF: rate must be positive, got ${rate}`)
  if (x <= 0) return 0
  return incompleteGamma(shape, rate * x)
}

/**
 * Gamma quantile via bisection.
 */
export function gammaQuantile(p: number, shape: number, rate = 1): number {
  if (shape <= 0) throw new Error(`gammaQuantile: shape must be positive, got ${shape}`)
  if (rate <= 0) throw new Error(`gammaQuantile: rate must be positive, got ${rate}`)
  if (p < 0 || p > 1) throw new Error(`gammaQuantile: p must be in [0,1], got ${p}`)
  if (p === 0) return 0
  if (p === 1) return Infinity

  // Initial bracket: mean ± 10 SD of Gamma(shape, rate)
  const mean = shape / rate
  const sd = Math.sqrt(shape) / rate
  let lo = Math.max(1e-15, mean - 10 * sd)
  let hi = mean + 20 * sd

  // Expand hi until CDF(hi) > p
  while (gammaCDF(hi, shape, rate) < p) hi *= 2

  // Bisection
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2
    if (gammaCDF(mid, shape, rate) < p) lo = mid
    else hi = mid
    if ((hi - lo) / (lo + 1e-15) < 1e-12) break
  }
  return (lo + hi) / 2
}

/**
 * Beta CDF: P(X ≤ x | α, β) = I_x(α, β) where I is the regularized incomplete beta.
 */
export function betaCDF(x: number, shape1: number, shape2: number): number {
  if (shape1 <= 0 || shape2 <= 0)
    throw new Error(`betaCDF: shapes must be positive, got shape1=${shape1}, shape2=${shape2}`)
  if (x <= 0) return 0
  if (x >= 1) return 1
  return incompleteBeta(x, shape1, shape2)
}

/**
 * Beta quantile via bisection.
 */
export function betaQuantile(p: number, shape1: number, shape2: number): number {
  if (shape1 <= 0 || shape2 <= 0)
    throw new Error(`betaQuantile: shapes must be positive, got shape1=${shape1}, shape2=${shape2}`)
  if (p < 0 || p > 1) throw new Error(`betaQuantile: p must be in [0,1], got ${p}`)
  if (p === 0) return 0
  if (p === 1) return 1

  // Bisection on [0, 1]
  let lo = 0, hi = 1
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2
    if (betaCDF(mid, shape1, shape2) < p) lo = mid
    else hi = mid
    if (hi - lo < 1e-12) break
  }
  return (lo + hi) / 2
}

/**
 * F-distribution quantile via bisection.
 * MISSING from core/math.ts — needed for this module and other tests.
 */
export function fQuantile(p: number, df1: number, df2: number): number {
  if (df1 <= 0 || df2 <= 0) throw new Error(`fQuantile: df must be positive, got df1=${df1}, df2=${df2}`)
  if (p < 0 || p > 1) throw new Error(`fQuantile: p must be in [0,1], got ${p}`)
  if (p === 0) return 0
  if (p === 1) return Infinity

  // Initial bracket
  let lo = 0, hi = df2 > 4 ? df2 / (df2 - 2) * 5 : 100
  while (_fDistCDF(hi, df1, df2) < p) hi *= 2

  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2
    if (_fDistCDF(mid, df1, df2) < p) lo = mid
    else hi = mid
    if ((hi - lo) / (lo + 1e-15) < 1e-12) break
  }
  return (lo + hi) / 2
}

// ─── Section C: Discrete PMFs, CDFs, and Quantiles ───────────────────────────

/**
 * Binomial PMF: P(X = k | n, p) = C(n,k) p^k (1-p)^(n-k)
 * Computed in log-space for numerical stability.
 */
export function binomialPMF(k: number, n: number, prob: number): number {
  if (n < 0 || !Number.isInteger(n)) throw new Error(`binomialPMF: n must be non-negative integer, got ${n}`)
  if (prob < 0 || prob > 1) throw new Error(`binomialPMF: prob must be in [0,1], got ${prob}`)
  k = Math.floor(k)
  if (k < 0 || k > n) return 0
  if (prob === 0) return k === 0 ? 1 : 0
  if (prob === 1) return k === n ? 1 : 0
  // log C(n,k) = logGamma(n+1) - logGamma(k+1) - logGamma(n-k+1)
  const logC = logGamma(n + 1) - logGamma(k + 1) - logGamma(n - k + 1)
  return Math.exp(logC + k * Math.log(prob) + (n - k) * Math.log(1 - prob))
}

/**
 * Binomial CDF: P(X ≤ k | n, p) = I_{1-p}(n-k, k+1)
 * Uses regularized incomplete beta.
 * Reference: NIST Handbook §6.3
 */
export function binomialCDF(k: number, n: number, prob: number): number {
  if (n < 0 || !Number.isInteger(n)) throw new Error(`binomialCDF: n must be non-negative integer, got ${n}`)
  if (prob < 0 || prob > 1) throw new Error(`binomialCDF: prob must be in [0,1], got ${prob}`)
  k = Math.floor(k)
  if (k < 0) return 0
  if (k >= n) return 1
  // P(X ≤ k) = I_{1-p}(n-k, k+1)
  return incompleteBeta(1 - prob, n - k, k + 1)
}

/**
 * Binomial quantile via integer search.
 */
export function binomialQuantile(p: number, n: number, prob: number): number {
  if (p < 0 || p > 1) throw new Error(`binomialQuantile: p must be in [0,1], got ${p}`)
  if (n < 0 || !Number.isInteger(n)) throw new Error(`binomialQuantile: n must be non-negative integer, got ${n}`)
  if (prob < 0 || prob > 1) throw new Error(`binomialQuantile: prob must be in [0,1], got ${prob}`)
  if (p === 0) return 0
  if (p === 1) return n

  // Start near the mean and search
  const mu = n * prob
  let lo = Math.max(0, Math.floor(mu - 4 * Math.sqrt(mu * (1 - prob))))
  // Scan upward from lo
  for (let k = lo; k <= n; k++) {
    if (binomialCDF(k, n, prob) >= p) return k
  }
  return n
}

/**
 * Poisson PMF: P(X = k | λ) = λ^k exp(-λ) / k!
 * Computed in log-space for stability.
 */
export function poissonPMF(k: number, lambda: number): number {
  if (lambda < 0) throw new Error(`poissonPMF: lambda must be non-negative, got ${lambda}`)
  k = Math.floor(k)
  if (k < 0) return 0
  if (lambda === 0) return k === 0 ? 1 : 0
  const logP = k * Math.log(lambda) - lambda - logGamma(k + 1)
  return Math.exp(logP)
}

/**
 * Poisson CDF: P(X ≤ k | λ) = 1 - P(k+1, λ) where P is regularized incomplete gamma.
 * Uses: 1 - incompleteGamma(k+1, lambda)
 * Reference: NIST Handbook §6.7
 */
export function poissonCDF(k: number, lambda: number): number {
  if (lambda < 0) throw new Error(`poissonCDF: lambda must be non-negative, got ${lambda}`)
  k = Math.floor(k)
  if (k < 0) return 0
  if (lambda === 0) return 1
  return 1 - incompleteGamma(k + 1, lambda)
}

/**
 * Poisson quantile via integer search.
 */
export function poissonQuantile(p: number, lambda: number): number {
  if (p < 0 || p > 1) throw new Error(`poissonQuantile: p must be in [0,1], got ${p}`)
  if (lambda < 0) throw new Error(`poissonQuantile: lambda must be non-negative, got ${lambda}`)
  if (p === 0) return 0

  // Start near the mean and search upward
  const start = Math.max(0, Math.floor(lambda - 3 * Math.sqrt(lambda)))
  for (let k = start; ; k++) {
    if (poissonCDF(k, lambda) >= p) return k
  }
}

// ─── Section D: MLE Fitting ──────────────────────────────────────────────────

/**
 * Fit a distribution to data via maximum likelihood estimation.
 *
 * | Distribution | Algorithm |
 * |---|---|
 * | Normal       | Closed-form: μ=mean, σ=√(Σ(x-μ)²/n) [MLE denominator n] |
 * | Exponential  | Closed-form: rate = 1/mean |
 * | Uniform      | Closed-form: min=min(x), max=max(x) |
 * | Poisson      | Closed-form: λ = mean |
 * | Gamma        | Choi-Wette (1969) init + Newton-Raphson |
 * | Beta         | Method of moments init + Newton-Raphson |
 *
 * @param data - Sample observations
 * @param distribution - Which distribution to fit
 * @param options - maxIter and tol for iterative methods
 * @returns FitDistributionResult with params, logLik, AIC, BIC
 */
export function fitDistribution(
  data: readonly number[],
  distribution: FitDistName,
  options?: { maxIter?: number; tol?: number }
): FitDistributionResult {
  const n = data.length
  if (n === 0) throw new Error('fitDistribution: data must be non-empty')
  const maxIter = options?.maxIter ?? 1000
  const tol = options?.tol ?? 1e-8

  switch (distribution) {
    case 'normal':
      return fitNormal(data, n)
    case 'exponential':
      return fitExponential(data, n)
    case 'uniform':
      return fitUniform(data, n)
    case 'poisson':
      return fitPoisson(data, n)
    case 'gamma':
      return fitGamma(data, n, maxIter, tol)
    case 'beta':
      return fitBeta(data, n, maxIter, tol)
    default:
      throw new Error(`fitDistribution: unknown distribution '${distribution as string}'`)
  }
}

function fitNormal(data: readonly number[], n: number): FitDistributionResult {
  const mu = _mean(data)
  // MLE uses n denominator, not n-1
  const sigma = Math.sqrt(data.reduce((s, x) => s + (x - mu) ** 2, 0) / n)
  const logLik = data.reduce((s, x) => s + Math.log(normalPDF(x, mu, sigma)), 0)
  const k = 2 // 2 params: mu, sigma
  return {
    distribution: 'normal',
    params: { mu, sigma },
    logLik,
    aic: -2 * logLik + 2 * k,
    bic: -2 * logLik + k * Math.log(n),
    n,
    formatted: `Normal(μ = ${mu.toFixed(4)}, σ = ${sigma.toFixed(4)}), log-lik = ${logLik.toFixed(2)}, AIC = ${(-2 * logLik + 2 * k).toFixed(2)}`,
  }
}

function fitExponential(data: readonly number[], n: number): FitDistributionResult {
  if (data.some(x => x < 0)) throw new Error('fitDistribution: exponential requires non-negative data')
  const mu = _mean(data)
  if (mu <= 0) throw new Error('fitDistribution: exponential requires positive mean')
  const rate = 1 / mu
  const logLik = data.reduce((s, x) => s + Math.log(exponentialPDF(x, rate)), 0)
  const k = 1
  return {
    distribution: 'exponential',
    params: { rate },
    logLik,
    aic: -2 * logLik + 2 * k,
    bic: -2 * logLik + k * Math.log(n),
    n,
    formatted: `Exponential(rate = ${rate.toFixed(4)}), log-lik = ${logLik.toFixed(2)}, AIC = ${(-2 * logLik + 2 * k).toFixed(2)}`,
  }
}

function fitUniform(data: readonly number[], n: number): FitDistributionResult {
  const minVal = Math.min(...data)
  const maxVal = Math.max(...data)
  if (minVal === maxVal) throw new Error('fitDistribution: uniform requires non-degenerate data')
  const logLik = -n * Math.log(maxVal - minVal)
  const k = 2
  return {
    distribution: 'uniform',
    params: { min: minVal, max: maxVal },
    logLik,
    aic: -2 * logLik + 2 * k,
    bic: -2 * logLik + k * Math.log(n),
    n,
    formatted: `Uniform(min = ${minVal.toFixed(4)}, max = ${maxVal.toFixed(4)}), log-lik = ${logLik.toFixed(2)}, AIC = ${(-2 * logLik + 2 * k).toFixed(2)}`,
  }
}

function fitPoisson(data: readonly number[], n: number): FitDistributionResult {
  if (data.some(x => x < 0 || !Number.isInteger(x)))
    throw new Error('fitDistribution: poisson requires non-negative integer data')
  const lambda = _mean(data)
  const logLik = data.reduce((s, x) => s + Math.log(poissonPMF(x, lambda)), 0)
  const k = 1
  return {
    distribution: 'poisson',
    params: { lambda },
    logLik,
    aic: -2 * logLik + 2 * k,
    bic: -2 * logLik + k * Math.log(n),
    n,
    formatted: `Poisson(λ = ${lambda.toFixed(4)}), log-lik = ${logLik.toFixed(2)}, AIC = ${(-2 * logLik + 2 * k).toFixed(2)}`,
  }
}

/**
 * Gamma MLE using Choi-Wette (1969) initialization + Newton-Raphson.
 * The profile log-likelihood for shape α given rate β=α/x̄:
 *   ∂ℓ/∂α = n(log(α) - digamma(α) - log(x̄) + log_mean_x) = 0
 * where s = log(x̄) - mean(log(x)).
 * Newton update: α_{k+1} = α_k - (log(α_k) - digamma(α_k) - s) / (1/α_k - trigamma(α_k))
 * Reference: Choi & Wette (1969), Technometrics 11(4):683-690
 */
function fitGamma(data: readonly number[], n: number, maxIter: number, tol: number): FitDistributionResult {
  if (data.some(x => x <= 0)) throw new Error('fitDistribution: gamma requires positive data')

  const xBar = _mean(data)
  const logXBar = Math.log(xBar)
  const meanLogX = data.reduce((s, x) => s + Math.log(x), 0) / n
  const s = logXBar - meanLogX // s > 0 for typical data

  // Choi-Wette initial estimate: α₀ = (0.5000876 + 0.1648852 s - 0.0544274 s²) / s
  let alpha: number
  if (s < 0.5672) {
    alpha = (0.5000876 + 0.1648852 * s - 0.0544274 * s * s) / s
  } else {
    alpha = (8.898919 + 9.05995 * s + 0.9775373 * s * s) / (s * (17.79728 + 11.968477 * s + s * s))
  }

  // Newton-Raphson
  for (let iter = 0; iter < maxIter; iter++) {
    const fVal = Math.log(alpha) - digamma(alpha) - s
    const fPrime = 1 / alpha - trigamma(alpha)
    const delta = fVal / fPrime
    alpha -= delta
    if (alpha <= 0) alpha = 0.001 // guard against negative
    if (Math.abs(delta) < tol) break
  }

  const rate = alpha / xBar
  const logLik = data.reduce((s2, x) => s2 + Math.log(gammaPDF(x, alpha, rate)), 0)
  const k = 2
  return {
    distribution: 'gamma',
    params: { shape: alpha, rate },
    logLik,
    aic: -2 * logLik + 2 * k,
    bic: -2 * logLik + k * Math.log(n),
    n,
    formatted: `Gamma(shape = ${alpha.toFixed(4)}, rate = ${rate.toFixed(4)}), log-lik = ${logLik.toFixed(2)}, AIC = ${(-2 * logLik + 2 * k).toFixed(2)}`,
  }
}

/**
 * Beta MLE using method of moments initialization + Newton-Raphson.
 * Profile log-likelihood gradient:
 *   ∂ℓ/∂α = Σ log(x_i) - n(digamma(α) - digamma(α+β))
 *   ∂ℓ/∂β = Σ log(1-x_i) - n(digamma(β) - digamma(α+β))
 * Method of moments init:
 *   α₀ = x̄((x̄(1-x̄)/s² - 1)
 *   β₀ = (1-x̄)(x̄(1-x̄)/s² - 1)
 * Reference: Minka (2000), "Estimating a Dirichlet distribution"
 */
function fitBeta(data: readonly number[], n: number, maxIter: number, tol: number): FitDistributionResult {
  if (data.some(x => x <= 0 || x >= 1))
    throw new Error('fitDistribution: beta requires data strictly in (0, 1)')

  const xBar = _mean(data)
  const s2 = data.reduce((s, x) => s + (x - xBar) ** 2, 0) / n

  // Method of moments initial estimates
  const common = xBar * (1 - xBar) / s2 - 1
  let alpha = xBar * common
  let beta = (1 - xBar) * common
  if (alpha <= 0) alpha = 0.5
  if (beta <= 0) beta = 0.5

  const sumLogX = data.reduce((s, x) => s + Math.log(x), 0)
  const sumLog1mX = data.reduce((s, x) => s + Math.log(1 - x), 0)

  // Newton-Raphson with 2x2 Hessian
  for (let iter = 0; iter < maxIter; iter++) {
    const psiAB = digamma(alpha + beta)
    const triAB = trigamma(alpha + beta)

    // Gradient
    const gAlpha = sumLogX - n * (digamma(alpha) - psiAB)
    const gBeta = sumLog1mX - n * (digamma(beta) - psiAB)

    // Hessian
    const hAA = -n * (trigamma(alpha) - triAB)
    const hBB = -n * (trigamma(beta) - triAB)
    const hAB = n * triAB

    // 2x2 inverse: [hBB, -hAB; -hAB, hAA] / det
    const det = hAA * hBB - hAB * hAB
    if (Math.abs(det) < 1e-30) break

    const dAlpha = (hBB * gAlpha - hAB * gBeta) / det
    const dBeta = (hAA * gBeta - hAB * gAlpha) / det

    alpha -= dAlpha
    beta -= dBeta

    if (alpha <= 0) alpha = 0.001
    if (beta <= 0) beta = 0.001

    if (Math.abs(dAlpha) + Math.abs(dBeta) < tol) break
  }

  const logLik = data.reduce((s, x) => s + Math.log(betaPDF(x, alpha, beta)), 0)
  const k = 2
  return {
    distribution: 'beta',
    params: { shape1: alpha, shape2: beta },
    logLik,
    aic: -2 * logLik + 2 * k,
    bic: -2 * logLik + k * Math.log(n),
    n,
    formatted: `Beta(α = ${alpha.toFixed(4)}, β = ${beta.toFixed(4)}), log-lik = ${logLik.toFixed(2)}, AIC = ${(-2 * logLik + 2 * k).toFixed(2)}`,
  }
}

// ─── Section E: Goodness-of-Fit Tests ────────────────────────────────────────

/**
 * Anderson-Darling goodness-of-fit test.
 *
 * Tests whether data follow a specified distribution. If distribution/params are
 * omitted, defaults to normal with MLE parameters.
 *
 * Test statistic: A² = -n - (1/n) Σ_{i=1}^{n} (2i - 1)[ln(z_i) + ln(1 - z_{n+1-i})]
 * where z_i = CDF(x_{(i)}) using the hypothesized distribution.
 *
 * For normal case, p-value uses Stephens (1974) modification and polynomial approximation.
 * Reference: Stephens (1974), "EDF Statistics for Goodness of Fit and Some Comparisons", JASA
 *            D'Agostino & Stephens (1986), "Goodness-of-Fit Techniques"
 *
 * @param data - Sample observations
 * @param distribution - Distribution to test against (default: 'normal')
 * @param params - Distribution parameters. If omitted, estimated via MLE.
 */
export function andersonDarling(
  data: readonly number[],
  distribution?: FitDistName,
  params?: Readonly<Record<string, number>>
): StatResult {
  const n = data.length
  if (n < 3) throw new Error('andersonDarling: need at least 3 observations')

  const dist = distribution ?? 'normal'
  const sorted = sortAsc([...data])

  // Get CDF function and fit params if needed.
  // For normal case without explicit params, use sample mean/sd (n-1 denominator)
  // to match R's nortest::ad.test behavior, NOT MLE sigma (n denominator).
  let fitParams: Readonly<Record<string, number>>
  if (params) {
    fitParams = params
  } else if (dist === 'normal') {
    const mu = _mean(data)
    const s = Math.sqrt(data.reduce((acc, x) => acc + (x - mu) ** 2, 0) / (n - 1))
    fitParams = { mu, sigma: s }
  } else {
    fitParams = fitDistribution(data, dist).params
  }
  const cdfFn = getCDF(dist, fitParams)

  // Compute z_i = CDF(x_(i))
  const z = sorted.map(x => cdfFn(x))

  // Clamp z to (epsilon, 1-epsilon) to avoid log(0)
  const eps = 1e-15
  const zClamp = z.map(v => Math.max(eps, Math.min(1 - eps, v)))

  // A² = -n - (1/n) Σ (2i-1)[ln(z_i) + ln(1 - z_{n+1-i})]
  let sum = 0
  for (let i = 0; i < n; i++) {
    sum += (2 * (i + 1) - 1) * (Math.log(zClamp[i]!) + Math.log(1 - zClamp[n - 1 - i]!))
  }
  const A2 = -n - sum / n

  // p-value depends on whether params were estimated
  const paramsEstimated = !params
  const pValue = adPValue(A2, n, dist, paramsEstimated)

  return {
    testName: 'Anderson-Darling',
    statistic: A2,
    df: NaN,
    pValue,
    effectSize: { value: A2, name: 'A²', interpretation: 'negligible' },
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted: `Anderson-Darling A² = ${A2.toFixed(4)}, p ${pValue < 0.001 ? '< .001' : `= ${pValue.toFixed(3).replace(/^0/, '')}`}`,
  }
}

/**
 * Kolmogorov-Smirnov goodness-of-fit test.
 *
 * Tests whether data follow a specified distribution. If distribution/params are
 * omitted, defaults to normal with MLE parameters.
 *
 * Test statistic: D = max_i |F_n(x_i) - F_0(x_i)| (two-sided)
 * P-value: asymptotic Kolmogorov distribution.
 *
 * Reference: Conover (1999), "Practical Nonparametric Statistics"
 *
 * @param data - Sample observations
 * @param distribution - Distribution to test against (default: 'normal')
 * @param params - Distribution parameters. If omitted, estimated via MLE.
 */
export function kolmogorovSmirnov(
  data: readonly number[],
  distribution?: FitDistName,
  params?: Readonly<Record<string, number>>
): StatResult {
  const n = data.length
  if (n < 1) throw new Error('kolmogorovSmirnov: need at least 1 observation')

  const dist = distribution ?? 'normal'
  const sorted = sortAsc([...data])

  // Get CDF function and fit params if needed
  const fitParams = params ?? fitDistribution(data, dist).params
  const cdfFn = getCDF(dist, fitParams)

  // D = max(D+, D-)
  let dPlus = 0, dMinus = 0
  for (let i = 0; i < n; i++) {
    const fi = cdfFn(sorted[i]!)
    const dp = (i + 1) / n - fi      // F_n(x_i) - F_0(x_i) using upper step
    const dm = fi - i / n            // F_0(x_i) - F_n(x_{i-1}) using lower step
    dPlus = Math.max(dPlus, dp)
    dMinus = Math.max(dMinus, dm)
  }
  const D = Math.max(dPlus, dMinus)

  // Asymptotic p-value: P(D_n > d) ≈ 2 Σ_{k=1}^{∞} (-1)^{k+1} exp(-2k²n·d²)
  // For params estimated from data, this is conservative (true p is smaller).
  const pValue = ksPValue(D, n)

  return {
    testName: 'Kolmogorov-Smirnov',
    statistic: D,
    df: NaN,
    pValue,
    effectSize: { value: D, name: 'D', interpretation: 'negligible' },
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted: `Kolmogorov-Smirnov D = ${D.toFixed(4)}, p ${pValue < 0.001 ? '< .001' : `= ${pValue.toFixed(3).replace(/^0/, '')}`}`,
  }
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/** Get a CDF function for a named distribution with given parameters. */
function getCDF(dist: FitDistName, params: Readonly<Record<string, number>>): (x: number) => number {
  switch (dist) {
    case 'normal': {
      const mu = params['mu'] ?? 0
      const sigma = params['sigma'] ?? 1
      return (x: number) => _normalCDF((x - mu) / sigma)
    }
    case 'exponential': {
      const rate = params['rate'] ?? 1
      return (x: number) => exponentialCDF(x, rate)
    }
    case 'uniform': {
      const min = params['min'] ?? 0
      const max = params['max'] ?? 1
      return (x: number) => uniformCDF(x, min, max)
    }
    case 'gamma': {
      const shape = params['shape'] ?? 1
      const rate = params['rate'] ?? 1
      return (x: number) => gammaCDF(x, shape, rate)
    }
    case 'beta': {
      const shape1 = params['shape1'] ?? 1
      const shape2 = params['shape2'] ?? 1
      return (x: number) => betaCDF(x, shape1, shape2)
    }
    case 'poisson':
      throw new Error('andersonDarling/kolmogorovSmirnov: Poisson distribution not supported (discrete)')
    default:
      throw new Error(`getCDF: unknown distribution '${dist as string}'`)
  }
}

/**
 * Anderson-Darling p-value computation.
 * For normal distribution with estimated parameters, uses the Stephens (1986) modification:
 *   A*² = A²(1 + 0.75/n + 2.25/n²)
 * Then polynomial approximation from D'Agostino & Stephens (1986), Table 4.9:
 *   if A*² < 0.2: p = 1 - exp(-13.436 + 101.14·A*² - 223.73·A*²²)
 *   if 0.2 ≤ A*² < 0.34: p = 1 - exp(-8.318 + 42.796·A*² - 59.938·A*²²)
 *   if 0.34 ≤ A*² < 0.6: p = exp(0.9177 - 4.279·A*² - 1.38·A*²²)
 *   if 0.6 ≤ A*² ≤ 13: p = exp(1.2937 - 5.709·A*² + 0.0186·A*²²)
 *   if A*² > 13: p ≈ 0
 *
 * Reference: D'Agostino & Stephens (1986), pp. 122-123
 */
function adPValue(A2: number, n: number, dist: FitDistName, paramsEstimated: boolean): number {
  if (dist === 'normal' && paramsEstimated) {
    // Stephens modification
    const Astar = A2 * (1 + 0.75 / n + 2.25 / (n * n))
    if (Astar < 0.2) return 1 - Math.exp(-13.436 + 101.14 * Astar - 223.73 * Astar * Astar)
    if (Astar < 0.34) return 1 - Math.exp(-8.318 + 42.796 * Astar - 59.938 * Astar * Astar)
    if (Astar < 0.6) return Math.exp(0.9177 - 4.279 * Astar - 1.38 * Astar * Astar)
    if (Astar <= 13) return Math.exp(1.2937 - 5.709 * Astar + 0.0186 * Astar * Astar)
    return 0
  }

  // For other distributions or known parameters: asymptotic approximation
  // Marsaglia & Marsaglia (2004) give critical values; use Lewis (1961) series
  // as a reasonable approximation for the general case.
  // Approximate via asymptotic: P(A² > z) ≈ as per Marsaglia (2004)
  return adPValueAsymptotic(A2)
}

/**
 * Asymptotic p-value for Anderson-Darling when CDF is fully specified (no estimated params).
 * Uses Marsaglia's (2004) rational approximation.
 * Reference: Marsaglia & Marsaglia (2004), "Evaluating the Anderson-Darling Distribution", JCAM
 */
function adPValueAsymptotic(A2: number): number {
  if (A2 <= 0) return 1
  // Approximation from Marsaglia (2004)
  if (A2 < 2) {
    const z = A2
    const p = Math.exp(-1.2337141 / z) / Math.sqrt(z) * (
      2.00012 + (0.247105 - (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * z) * z) * z) * z) * z
    )
    return Math.max(0, Math.min(1, 1 - p))
  }
  // For A² ≥ 2
  const z = A2
  const p = Math.exp(-Math.exp(1.0776 - (2.30695 - (0.43424 - (0.082433 - (0.008056 - 0.0003146 * z) * z) * z) * z) * z))
  return Math.max(0, Math.min(1, p))
}

/**
 * Kolmogorov-Smirnov p-value using asymptotic formula.
 * P(D_n > d) ≈ 2 Σ_{k=1}^{∞} (-1)^{k+1} exp(-2k² (√n · d + 0.12 + 0.11/√n)²)
 * The correction terms (0.12, 0.11) are from Stephens (1970) for improved small-sample accuracy.
 * Reference: R's ks.test implementation, Marsaglia et al. (2003)
 */
function ksPValue(D: number, n: number): number {
  if (D <= 0) return 1
  if (D >= 1) return 0

  // Use the classic asymptotic: 2 Σ (-1)^{k+1} exp(-2k²t²) where t = √n · D
  const t = Math.sqrt(n) * D
  let sum = 0
  for (let k = 1; k <= 100; k++) {
    const sign = k % 2 === 0 ? -1 : 1
    const term = sign * Math.exp(-2 * k * k * t * t)
    sum += term
    if (Math.abs(term) < 1e-15) break
  }
  return Math.max(0, Math.min(1, 2 * sum))
}
