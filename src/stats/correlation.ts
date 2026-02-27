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
  normalSurvival,
  pKendallExact,
  pSpearmanExact,
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
  const rho = rhoResult.statistic

  // Check for ties in ranks
  const hasTies = new Set(rx).size < n || new Set(ry).size < n

  // Exact/Edgeworth p-value for n ≤ 1290 without ties (matching R's cor.test)
  // R uses AS89 (prho.c): exact enumeration for n ≤ 9, Edgeworth series for n ≥ 10
  let pValue = rhoResult.pValue
  if (!hasTies && n <= 1290) {
    // Compute D = Σ(d_i²) directly from ranks (avoid rounding through rho)
    let D = 0
    for (let i = 0; i < n; i++) {
      const d = rx[i]! - ry[i]!
      D += d * d
    }
    const q = (n * n * n - n) / 6 // expected D under H0

    // Two-sided test matching R's cor.test logic:
    // R calls pspearman(q, n, lower.tail) which passes round(q) + 2*lower.tail to C_pRho
    // prho computes P(S >= is) [upper] or P(S < is) [lower]
    let p: number
    if (D > q) {
      // Negative correlation: upper tail P(S >= D)
      p = pSpearmanExact(Math.round(D), n, false)
    } else {
      // Positive correlation: lower tail P(S < D+2) = P(S <= D)
      p = pSpearmanExact(Math.round(D) + 2, n, true)
    }
    pValue = Math.min(2 * p, 1)
  }

  return {
    ...rhoResult,
    pValue,
    testName: "Spearman's ρ",
    effectSize: {
      ...rhoResult.effectSize,
      name: "Spearman's ρ",
    },
    formatted: formatCorrelation(rho, typeof rhoResult.df === 'number' ? rhoResult.df : 0, pValue, rhoResult.ci, 'ρ', ciLevel),
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

  // Tie-corrected variance formula — Kendall (1970)
  // Compute tie group sizes for X and Y
  const xCounts = new Map<number, number>()
  const yCounts = new Map<number, number>()
  for (let i = 0; i < n; i++) {
    const xv = x[i]!
    const yv = y[i]!
    xCounts.set(xv, (xCounts.get(xv) ?? 0) + 1)
    yCounts.set(yv, (yCounts.get(yv) ?? 0) + 1)
  }

  const v0 = n * (n - 1) * (2 * n + 5)
  let vt = 0  // tie correction for X
  let vu = 0  // tie correction for Y
  let v1t = 0 // t_i * (t_i - 1)
  let v1u = 0 // u_j * (u_j - 1)
  let v2t = 0 // t_i * (t_i - 1) * (t_i - 2)
  let v2u = 0 // u_j * (u_j - 1) * (u_j - 2)

  for (const t of xCounts.values()) {
    vt += t * (t - 1) * (2 * t + 5)
    v1t += t * (t - 1)
    v2t += t * (t - 1) * (t - 2)
  }
  for (const u of yCounts.values()) {
    vu += u * (u - 1) * (2 * u + 5)
    v1u += u * (u - 1)
    v2u += u * (u - 1) * (u - 2)
  }

  const sigma2 = (v0 - vt - vu) / 18
    + (v1t * v1u) / (2 * n * (n - 1))
    + (v2t * v2u) / (9 * n * (n - 1) * (n - 2))

  const z = (concordant - discordant) / Math.sqrt(sigma2)

  // Check for ties
  const hasTies = tiesX > 0 || tiesY > 0

  // Exact p-value for n < 50 without ties (matching R's cor.test)
  let pValue: number
  if (!hasTies && n < 50) {
    // T = number of concordant pairs = (S + n2) / 2 where S = concordant - discordant
    const S = concordant - discordant
    const T = (S + n2) / 2
    // Two-sided: P(tau >= |observed|) = P(T >= T_obs) + P(T <= n2 - T_obs)
    // By symmetry: p = 2 * min(P(T <= T_obs), P(T >= T_obs))
    const pLower = pKendallExact(Math.round(T), n)
    const pUpper = 1 - pKendallExact(Math.round(T) - 1, n)
    pValue = Math.min(1, 2 * Math.min(pLower, pUpper))
  } else {
    pValue = 2 * normalSurvival(Math.abs(z))
  }

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
  const result = pearsonCorrelation(xRes, yRes)

  // Correct degrees of freedom: df = n - 2 - q (q = number of control variables)
  const n = x.length
  const q = controls.length
  const correctedDf = n - 2 - q
  if (correctedDf < 1) return result  // can't correct if df < 1

  const r = result.effectSize.value
  const t = Math.abs(r) === 1 ? Infinity : r * Math.sqrt(correctedDf / (1 - r * r))
  const pValue = tDistPValue(t, correctedDf)

  return {
    ...result,
    df: correctedDf,
    pValue: roundTo(pValue, 4),
    formatted: formatCorrelation(result.statistic, correctedDf, pValue, result.ci, 'r', result.ciLevel),
  }
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
        return corrFn(data[i]!, data[j]!).effectSize.value
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

// ─── Point-Biserial correlation ──────────────────────────────────────────

/**
 * Point-biserial correlation between a binary variable and a continuous variable.
 * Validates binary is 0/1, delegates to Pearson, renames to r_pb.
 *
 * Cross-validated with R:
 * > cor.test(binary, continuous)
 * # Point-biserial is just Pearson r when one var is binary
 */
export function pointBiserialCorrelation(
  binary: readonly number[],
  continuous: readonly number[],
  ciLevel = 0.95
): StatResult {
  // validate binary is all 0 or 1
  for (const v of binary) {
    if (v !== 0 && v !== 1) throw new Error('pointBiserialCorrelation: binary must contain only 0 and 1')
  }
  const result = pearsonCorrelation(binary, continuous, ciLevel)
  return {
    ...result,
    testName: 'Point-biserial r',
    effectSize: {
      ...result.effectSize,
      name: 'r_pb',
    },
    formatted: formatCorrelation(result.statistic, typeof result.df === 'number' ? result.df : 0, result.pValue, result.ci, 'r_pb', ciLevel),
  }
}
