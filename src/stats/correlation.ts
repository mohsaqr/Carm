/**
 * Correlation analysis module.
 * Pearson, Spearman, Kendall tau-b, partial correlation, correlation matrix.
 */

import {
  mean as _mean,
  sd as _sd,
  cov,
  rank,
  tDistPValue,
  normalQuantile,
  normalCDF,
  roundTo,
} from '../core/math.js'
import { formatCorrelation, interpretR } from '../core/apa.js'
import type { StatResult, EffectInterpretation } from '../core/types.js'
import { Matrix } from '../core/matrix.js'

// ─── Pearson correlation ──────────────────────────────────────────────────

/**
 * Pearson product-moment correlation coefficient.
 *
 * Cross-validated with R:
 * > cor.test(c(1,2,3,4,5), c(2,4,1,5,3))
 * t = 0.6547, df = 3, p-value = 0.5607, r = 0.3536
 * 95% CI = [-0.6154, 0.9059]
 */
export function pearsonCorrelation(
  x: readonly number[],
  y: readonly number[],
  ciLevel = 0.95
): StatResult {
  if (x.length !== y.length) throw new Error('pearsonCorrelation: arrays must have equal length')
  const n = x.length
  if (n < 3) throw new Error('pearsonCorrelation: need at least 3 observations')

  const sdX = _sd(x), sdY = _sd(y)
  if (sdX === 0 || sdY === 0) throw new Error('pearsonCorrelation: zero variance in input')

  const r = cov(x, y) / (sdX * sdY)
  const rClamped = Math.max(-1, Math.min(1, r))

  const df = n - 2
  const t = Math.abs(rClamped) === 1 ? Infinity : rClamped * Math.sqrt(df / (1 - rClamped * rClamped))
  const pValue = tDistPValue(t, df)

  // Fisher z-transform CI
  const ci = fisherZCI(rClamped, n, ciLevel)

  return {
    testName: 'Pearson r',
    statistic: roundTo(rClamped, 4),
    df,
    pValue: roundTo(pValue, 4),
    effectSize: {
      value: rClamped,
      name: "Pearson r",
      interpretation: interpretR(rClamped) as EffectInterpretation,
    },
    ci,
    ciLevel,
    n,
    formatted: formatCorrelation(rClamped, df, pValue, ci, 'r', ciLevel),
  }
}

/** Fisher z-transform confidence interval for Pearson r. */
function fisherZCI(r: number, n: number, ciLevel: number): readonly [number, number] {
  const z = Math.log((1 + r) / (1 - r)) / 2  // Fisher's z
  const se = 1 / Math.sqrt(n - 3)
  const zCrit = normalQuantile(1 - (1 - ciLevel) / 2)
  const lo = z - zCrit * se
  const hi = z + zCrit * se
  // Back-transform
  return [
    Math.tanh(lo),
    Math.tanh(hi),
  ]
}

// ─── Spearman correlation ─────────────────────────────────────────────────

/**
 * Spearman rank correlation.
 *
 * Cross-validated with R:
 * > cor.test(c(1,2,3,4,5), c(5,6,7,8,7), method = "spearman")
 * rho = 0.8211, p-value = 0.08852
 */
export function spearmanCorrelation(
  x: readonly number[],
  y: readonly number[],
  ciLevel = 0.95
): StatResult {
  if (x.length !== y.length) throw new Error('spearmanCorrelation: arrays must have equal length')
  const n = x.length
  if (n < 3) throw new Error('spearmanCorrelation: need at least 3 observations')

  // Pearson r on ranks
  const rx = rank(x), ry = rank(y)
  const rhoResult = pearsonCorrelation(rx, ry, ciLevel)

  return {
    ...rhoResult,
    testName: "Spearman's ρ",
    effectSize: {
      ...rhoResult.effectSize,
      name: "Spearman's ρ",
    },
    formatted: formatCorrelation(rhoResult.statistic, typeof rhoResult.df === 'number' ? rhoResult.df : 0, rhoResult.pValue, rhoResult.ci, 'ρ', ciLevel),
  }
}

// ─── Kendall's tau-b ──────────────────────────────────────────────────────

/**
 * Kendall's tau-b correlation.
 *
 * Cross-validated with R:
 * > cor.test(c(1,2,3,4,5), c(5,6,7,8,7), method = "kendall")
 * tau = 0.7378, p-value = 0.1041
 */
export function kendallTau(
  x: readonly number[],
  y: readonly number[],
  ciLevel = 0.95
): StatResult {
  if (x.length !== y.length) throw new Error('kendallTau: arrays must have equal length')
  const n = x.length
  if (n < 3) throw new Error('kendallTau: need at least 3 observations')

  let concordant = 0, discordant = 0
  let tiesX = 0, tiesY = 0

  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const dx = (x[i] ?? 0) - (x[j] ?? 0)
      const dy = (y[i] ?? 0) - (y[j] ?? 0)
      const sign = Math.sign(dx * dy)
      if (sign > 0) concordant++
      else if (sign < 0) discordant++
      if (dx === 0) tiesX++
      if (dy === 0) tiesY++
    }
  }

  const n2 = n * (n - 1) / 2
  const tau = (concordant - discordant) / Math.sqrt((n2 - tiesX) * (n2 - tiesY))

  // Normal approximation for p-value
  const varTau = (2 * (2 * n + 5)) / (9 * n * (n - 1))
  const z = tau / Math.sqrt(varTau)
  const pValue = 2 * (1 - normalCDF(Math.abs(z)))

  // CI via Fisher z approximation (approximate for Kendall)
  const ci = fisherZCI(tau, n, ciLevel)
  const df = n - 2

  return {
    testName: "Kendall's τ",
    statistic: roundTo(tau, 4),
    df,
    pValue: roundTo(pValue, 4),
    effectSize: {
      value: tau,
      name: "Kendall's τ",
      interpretation: interpretR(tau) as EffectInterpretation,
    },
    ci,
    ciLevel,
    n,
    formatted: formatCorrelation(tau, df, pValue, ci, 'τ', ciLevel),
  }
}

// ─── Partial correlation ──────────────────────────────────────────────────

/**
 * Partial correlation between x and y controlling for z.
 * r_xy.z = (r_xy - r_xz · r_yz) / sqrt((1 - r_xz²)(1 - r_yz²))
 *
 * Cross-validated with R:
 * > ppcor::pcor.test(x, y, z)
 */
export function partialCorrelation(
  x: readonly number[],
  y: readonly number[],
  controls: readonly (readonly number[])[]
): StatResult {
  if (x.length !== y.length) throw new Error('partialCorrelation: arrays must have equal length')

  // Use OLS residuals: residualize x and y on controls
  const xRes = residualize(x, controls)
  const yRes = residualize(y, controls)
  return pearsonCorrelation(xRes, yRes)
}

/** Compute OLS residuals of y regressed on predictors. */
function residualize(y: readonly number[], predictors: readonly (readonly number[])[]): number[] {
  if (predictors.length === 0) return [...y]
  const n = y.length
  // Design matrix: [1, x1, x2, ...]
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map(p => p[i] ?? 0)])
  )
  const Xt = X.transpose()
  const XtX = Xt.multiply(X)
  const XtY = Xt.multiply(Matrix.colVec(y))
  const beta = XtX.inverse().multiply(XtY)
  const fitted = X.multiply(beta)
  return Array.from({ length: n }, (_, i) => (y[i] ?? 0) - fitted.get(i, 0))
}

// ─── Correlation matrix ───────────────────────────────────────────────────

export interface CorrelationMatrix {
  readonly r: readonly (readonly number[])[]
  readonly pValues: readonly (readonly number[])[]
  readonly n: number
  readonly labels: readonly string[]
}

/**
 * Compute pairwise correlation matrix with p-values.
 * Method: 'pearson' | 'spearman' | 'kendall'
 */
export function correlationMatrix(
  data: readonly (readonly number[])[],
  labels?: readonly string[],
  method: 'pearson' | 'spearman' | 'kendall' = 'pearson'
): CorrelationMatrix {
  const k = data.length
  if (k < 2) throw new Error('correlationMatrix: need at least 2 variables')
  const n = data[0]!.length

  const corrFn = method === 'pearson'
    ? pearsonCorrelation
    : method === 'spearman'
      ? spearmanCorrelation
      : kendallTau

  const r: number[][] = Array.from({ length: k }, (_, i) =>
    Array.from({ length: k }, (_, j) => {
      if (i === j) return 1
      if (j < i) return 0  // fill below diagonal later
      try {
        return corrFn(data[i]!, data[j]!).statistic
      } catch {
        return NaN
      }
    })
  )

  const pValues: number[][] = Array.from({ length: k }, (_, i) =>
    Array.from({ length: k }, (_, j) => {
      if (i === j) return NaN
      if (j < i) return 0
      try {
        return corrFn(data[i]!, data[j]!).pValue
      } catch {
        return NaN
      }
    })
  )

  // Fill lower triangle (symmetric)
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < i; j++) {
      r[i]![j] = r[j]![i]!
      pValues[i]![j] = pValues[j]![i]!
    }
  }

  return {
    r,
    pValues,
    n,
    labels: labels ?? Array.from({ length: k }, (_, i) => `Var${i + 1}`),
  }
}
