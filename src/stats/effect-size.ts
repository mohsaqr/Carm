/**
 * Effect size calculations for Carm.
 * Cohen's d, Hedges' g, eta-squared, omega-squared, rank-biserial correlation.
 * All functions return structured EffectSize objects.
 */

import { mean as _mean, variance as _variance, sd as _sd } from '../core/math.js'
import {
  interpretCohensD,
  interpretEtaSq,
  interpretR,
} from '../core/apa.js'
import type { EffectSize, EffectInterpretation } from '../core/types.js'

// ─── Cohen's d ────────────────────────────────────────────────────────────

/**
 * Cohen's d for independent samples.
 * Uses pooled SD (equal variances assumed).
 * Formula: d = (M₁ - M₂) / SD_pooled
 * Reference: Cohen (1988), "Statistical Power Analysis for the Behavioral Sciences"
 *
 * Cross-validated with R:
 * > library(effsize)
 * > cohen.d(c(1,2,3,4,5), c(3,4,5,6,7))
 * d = -1.2649...  (negative because group2 > group1)
 */
export function cohensD(x1: readonly number[], x2: readonly number[]): EffectSize {
  if (x1.length < 2 || x2.length < 2) throw new Error('cohensD: need at least 2 observations per group')
  const n1 = x1.length, n2 = x2.length
  const m1 = _mean(x1), m2 = _mean(x2)
  const v1 = _variance(x1), v2 = _variance(x2)
  // Pooled SD
  const sdPooled = Math.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  if (sdPooled === 0) return { value: 0, name: "Cohen's d", interpretation: 'negligible' }
  const d = (m1 - m2) / sdPooled
  return {
    value: d,
    name: "Cohen's d",
    interpretation: interpretCohensD(d) as EffectInterpretation,
  }
}

/**
 * Cohen's d for paired samples (dependent t-test).
 * Formula: d = M_diff / SD_diff
 */
export function cohensDPaired(diffs: readonly number[]): EffectSize {
  if (diffs.length < 2) throw new Error('cohensDPaired: need at least 2 differences')
  const m = _mean(diffs)
  const s = _sd(diffs)
  const d = s === 0 ? 0 : m / s
  return {
    value: d,
    name: "Cohen's d",
    interpretation: interpretCohensD(d) as EffectInterpretation,
  }
}

// ─── Hedges' g ────────────────────────────────────────────────────────────

/**
 * Hedges' g — bias-corrected version of Cohen's d.
 * Correction factor J = 1 - 3/(4·df - 1), df = n1 + n2 - 2.
 * Reference: Hedges (1981), Journal of Educational Statistics
 */
export function hedgesG(x1: readonly number[], x2: readonly number[]): EffectSize {
  const d = cohensD(x1, x2)
  const df = x1.length + x2.length - 2
  const J = 1 - 3 / (4 * df - 1)
  const g = d.value * J
  return {
    value: g,
    name: "Hedges' g",
    interpretation: interpretCohensD(g) as EffectInterpretation,
  }
}

// ─── Eta-squared ──────────────────────────────────────────────────────────

/**
 * Eta-squared: η² = SS_between / SS_total
 * For one-way ANOVA. Between 0 and 1.
 * Reference: Cohen (1973)
 */
export function etaSquared(ssBetween: number, ssTotal: number): EffectSize {
  if (ssTotal <= 0) return { value: 0, name: 'η²', interpretation: 'negligible' }
  const eta2 = Math.max(0, Math.min(1, ssBetween / ssTotal))
  return {
    value: eta2,
    name: 'η²',
    interpretation: interpretEtaSq(eta2) as EffectInterpretation,
  }
}

/**
 * Omega-squared: ω² = (SS_between - df_between · MS_within) / (SS_total + MS_within)
 * Less biased than eta-squared.
 * Reference: Hays (1963)
 */
export function omegaSquared(
  ssBetween: number,
  ssTotal: number,
  dfBetween: number,
  msWithin: number
): EffectSize {
  const denom = ssTotal + msWithin
  if (denom <= 0) return { value: 0, name: 'ω²', interpretation: 'negligible' }
  const omega2 = Math.max(0, (ssBetween - dfBetween * msWithin) / denom)
  return {
    value: omega2,
    name: 'ω²',
    interpretation: interpretEtaSq(omega2) as EffectInterpretation,
  }
}

// ─── Rank-biserial correlation ────────────────────────────────────────────

/**
 * Rank-biserial correlation for Mann-Whitney U.
 * r = 1 - (2U) / (n1 * n2)
 * Reference: Wendt (1972)
 */
export function rankBiserial(U: number, n1: number, n2: number): EffectSize {
  const r = 1 - 2 * U / (n1 * n2)
  return {
    value: r,
    name: 'r (rank-biserial)',
    interpretation: interpretR(r) as EffectInterpretation,
  }
}

/**
 * Rank-biserial correlation for Wilcoxon signed-rank (paired).
 * r = T / (n(n+1)/2) where T = sum of positive ranks (or negative).
 */
export function rankBiserialWilcoxon(T: number, n: number): EffectSize {
  const maxT = n * (n + 1) / 2
  const r = maxT > 0 ? T / maxT * 2 - 1 : 0
  return {
    value: r,
    name: 'r (rank-biserial)',
    interpretation: interpretR(r) as EffectInterpretation,
  }
}

// ─── Eta² for Kruskal-Wallis ──────────────────────────────────────────────

/**
 * Eta-squared for Kruskal-Wallis: η²_H = (H - k + 1) / (n - k)
 * Reference: Tomczak & Tomczak (2014)
 */
export function etaSquaredKW(H: number, k: number, n: number): EffectSize {
  const eta2 = Math.max(0, (H - k + 1) / (n - k))
  return {
    value: eta2,
    name: 'η²_H',
    interpretation: interpretEtaSq(eta2) as EffectInterpretation,
  }
}

// ─── CI for Cohen's d ─────────────────────────────────────────────────────

/**
 * Normal approximation CI for Cohen's d.
 * Reference: Hedges & Olkin (1985), "Statistical Methods for Meta-Analysis"
 */
export function cohensDCI(
  d: number,
  n1: number,
  n2: number,
  ciLevel = 0.95
): readonly [number, number] {
  const sampleSE = Math.sqrt((n1 + n2) / (n1 * n2) + d * d / (2 * (n1 + n2)))
  const z = normalQuantileInline(1 - (1 - ciLevel) / 2)
  return [d - z * sampleSE, d + z * sampleSE]
}

/** Inline normal quantile to avoid circular import. */
function normalQuantileInline(p: number): number {
  const a = [2.515517, 0.802853, 0.010328]
  const b = [1.432788, 0.189269, 0.001308]
  const t = Math.sqrt(-2 * Math.log(p <= 0.5 ? p : 1 - p))
  const num = a[0]! + a[1]! * t + a[2]! * t * t
  const den = 1 + b[0]! * t + b[1]! * t * t + b[2]! * t * t * t
  const x = t - num / den
  return p <= 0.5 ? -x : x
}
