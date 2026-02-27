/**
 * Descriptive statistics module.
 * Computes mean, median, mode, variance, SD, SE, skewness, kurtosis,
 * percentiles, confidence intervals, and the Shapiro-Wilk normality test.
 */

import {
  mean as _mean,
  median as _median,
  variance as _variance,
  sd as _sd,
  se as _se,
  quantile,
  sortAsc,
  tDistQuantile,
  normalCDF,
  normalQuantile,
  roundTo,
} from '../core/math.js'
import { formatP } from '../core/apa.js'
import type { DescriptiveResult } from '../core/types.js'

// ─── Mode ────────────────────────────────────────────────────────────────

function mode(x: readonly number[]): number[] {
  const counts = new Map<number, number>()
  for (const v of x) counts.set(v, (counts.get(v) ?? 0) + 1)
  const maxCount = Math.max(...counts.values())
  const modes: number[] = []
  for (const [v, c] of counts) {
    if (c === maxCount) modes.push(v)
  }
  return modes.sort((a, b) => a - b)
}

// ─── Trimmed mean ────────────────────────────────────────────────────────

/** α-trimmed mean: removes α proportion from each tail. */
export function trimmedMean(x: readonly number[], alpha = 0.05): number {
  const sorted = sortAsc(x)
  const n = sorted.length
  const trim = Math.floor(n * alpha)
  const trimmed = sorted.slice(trim, n - trim)
  if (trimmed.length === 0) throw new Error('trimmedMean: too much trimming, no data remains')
  return _mean(trimmed)
}

// ─── Skewness & Kurtosis ─────────────────────────────────────────────────

/** Sample skewness (adjusted Fisher-Pearson). */
export function skewness(x: readonly number[]): number {
  const n = x.length
  if (n < 3) throw new Error('skewness: need at least 3 observations')
  const m = _mean(x)
  const s = _sd(x)
  if (s === 0) return 0
  const sum = x.reduce((acc, v) => acc + ((v - m) / s) ** 3, 0)
  return (n / ((n - 1) * (n - 2))) * sum
}

/** Sample excess kurtosis (adjusted Fisher-Pearson). */
export function kurtosis(x: readonly number[]): number {
  const n = x.length
  if (n < 4) throw new Error('kurtosis: need at least 4 observations')
  const m = _mean(x)
  const s = _sd(x)
  if (s === 0) return 0
  const sum = x.reduce((acc, v) => acc + ((v - m) / s) ** 4, 0)
  const k1 = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3)) * sum
  const k2 = 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))
  return k1 - k2
}

// ─── Confidence interval for mean ────────────────────────────────────────

/** t-based CI for the mean. Returns [lower, upper]. */
export function ciMean(
  x: readonly number[],
  ciLevel = 0.95
): readonly [number, number] {
  const n = x.length
  if (n < 2) throw new Error('ciMean: need at least 2 observations')
  const m = _mean(x)
  const s = _se(x)
  const t = tDistQuantile(1 - (1 - ciLevel) / 2, n - 1)
  return [m - t * s, m + t * s]
}

// ─── Shapiro-Wilk normality test ──────────────────────────────────────────

// AS R94 polynomial coefficient arrays (ascending degree: c[0] + c[1]*x + c[2]*x² + ...)
// Source: Royston (1995) Applied Statistics 44(4):547-551, Table 1
const SW_C1 = [0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056] as const
const SW_C2 = [0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633] as const

/** Polynomial evaluation, ascending degree: c[0] + c[1]*x + ... (Horner's method) */
function swPoly(c: readonly number[], x: number): number {
  let r = 0
  for (let i = c.length - 1; i >= 0; i--) r = r * x + (c[i] ?? 0)
  return r
}

/**
 * Shapiro-Wilk W test for normality.
 * Full AS R94 algorithm (Royston 1995) — valid for n=3..5000.
 * Reference: Royston (1995), Applied Statistics, 44(4):547-551.
 *
 * Cross-validated with R:
 * > shapiro.test(c(1,2,3,4,5,6,7,8,9,10))
 * W = 0.9728, p-value = 0.9177
 */
export function shapiroWilk(x: readonly number[]): { statistic: number; pValue: number } {
  const n = x.length
  if (n < 3) throw new Error('shapiroWilk: need at least 3 observations')
  if (n > 5000) throw new Error('shapiroWilk: n > 5000 not supported')

  const sorted = sortAsc(x)
  const nn2 = Math.floor(n / 2)

  // a[0..nn2-1]: upper-half Shapiro-Wilk coefficients, all positive after computation.
  // a[0] = outermost (largest), a[nn2-1] = innermost (smallest).
  const a: number[] = new Array(nn2).fill(0)

  if (n === 3) {
    a[0] = Math.SQRT1_2  // 1/sqrt(2) ≈ 0.7071
  } else if (n === 4) {
    a[0] = 0.6872; a[1] = 0.1677
  } else if (n === 5) {
    a[0] = 0.6646; a[1] = 0.2413
  } else {
    // n >= 6: full AS R94 algorithm — port of R's swilk.c
    //
    // Step 1: expected normal quantiles (1-indexed in R: a[i] = Φ⁻¹((i−0.375)/(n+0.25))).
    // These are negative for i=1..nn2 (small-probability quantiles).
    for (let i = 0; i < nn2; i++) {
      a[i] = normalQuantile((i + 1 - 0.375) / (n + 0.25))
    }

    // Step 2: sum of squares of the full antisymmetric array (2 × upper half).
    let summ2 = 0
    for (let i = 0; i < nn2; i++) summ2 += a[i]! * a[i]!
    summ2 *= 2
    const ssumm2 = Math.sqrt(summ2)
    const rsn = 1 / Math.sqrt(n)

    // Save originals before overwriting — needed to compute fac after a[0]/a[1] are replaced.
    const a0orig = a[0]!  // most-negative quantile → becomes outermost positive coefficient
    const a1orig = a[1]!  // second-most-negative → becomes second positive coefficient

    // Step 3: polynomial correction for outermost coefficient (R's a[1]).
    // Subtracting a0orig (negative) / ssumm2 effectively adds a positive term.
    const a1corr = swPoly(SW_C1, rsn) - a0orig / ssumm2

    if (n > 5) {
      // Step 4: polynomial correction for second coefficient (R's a[2]).
      const a2corr = -a1orig / ssumm2 + swPoly(SW_C2, rsn)

      // Step 5: scale factor — ensures the full antisymmetric array has unit sum of squares.
      // num/den: variance accounted for by the two corrected outer coefficients vs. their originals.
      const num = summ2 - 2 * a0orig * a0orig - 2 * a1orig * a1orig
      const den = 1 - 2 * a1corr * a1corr - 2 * a2corr * a2corr
      const fac = num > 0 && den > 0 ? Math.sqrt(num / den) : 1

      a[0] = a1corr; a[1] = a2corr
      // Remaining inner quantiles (still negative) divided by −fac become positive.
      for (let i = 2; i < nn2; i++) a[i] = (a[i] ?? 0) / -fac

    } else {
      // n === 6: correct only the outermost coefficient.
      const num = summ2 - 2 * a0orig * a0orig
      const den = 1 - 2 * a1corr * a1corr
      const fac = num > 0 && den > 0 ? Math.sqrt(num / den) : 1

      a[0] = a1corr
      for (let i = 1; i < nn2; i++) a[i] = (a[i] ?? 0) / -fac
    }
  }

  // W = (Σᵢ a[i] · (x[n−1−i] − x[i]))² / SST
  // Equivalent to the standard Σ aᵢ·x₍ᵢ₎ with a full antisymmetric array, but using only the
  // upper half: each pair contributes a[i]·x[n−1−i] − (−a[i])·x[i] = a[i]·(x[n−1−i] − x[i]).
  let w1 = 0
  for (let i = 0; i < nn2; i++) {
    w1 += (a[i] ?? 0) * ((sorted[n - 1 - i] ?? 0) - (sorted[i] ?? 0))
  }
  const sst = _variance(sorted) * (n - 1)
  const W = sst > 0 ? Math.min(1, (w1 * w1) / sst) : 1

  const pValue = shapiroWilkPValue(W, n)
  return { statistic: roundTo(W, 4), pValue: roundTo(pValue, 4) }
}

/**
 * Ascending-degree polynomial: c[0] + c[1]*x + c[2]*x² + ...  (Horner's method)
 * Matches R's poly() in swilk.c: zero-order coefficient first.
 */
function polyAsc(c: readonly number[], x: number): number {
  let r = c[c.length - 1]!
  for (let i = c.length - 2; i >= 0; i--) r = r * x + c[i]!
  return r
}

// R swilk.c coefficient arrays (Royston 1995, Applied Statistics 44(4):547-551)
// Ascending degree: c[0] + c[1]*x + c[2]*x² + c[3]*x³
const SW_G  = [-2.273, 0.459] as const                             // gamma (n ≤ 11), var = n
const SW_C3 = [0.544, -0.39978, 0.025054, -6.714e-4] as const     // mu (n ≤ 11), var = n
const SW_C4 = [1.3822, -0.77857, 0.062767, -0.0020322] as const   // log(sigma) (n ≤ 11), var = n
const SW_C5 = [-1.5861, -0.31082, -0.083751, 0.0038915] as const  // mu (n ≥ 12), var = log(n)
const SW_C6 = [-0.4803, -0.082676, 0.0030302] as const            // log(sigma) (n ≥ 12), var = log(n)

function shapiroWilkPValue(W: number, n: number): number {
  // n = 3: exact p-value (R's swilk.c special case)
  if (n === 3) {
    const pi6 = 6 / Math.PI
    const stqr = normalQuantile(0.75) // ≈ 0.6744898
    const pw = pi6 * (Math.asin(Math.sqrt(W)) - stqr)
    return Math.max(0, Math.min(1, pw))
  }

  const y0 = Math.log(1 - W)
  let z: number

  if (n <= 11) {
    // Royston (1995) small-sample approximation, variable = n
    const gamma = polyAsc(SW_G, n)
    if (y0 >= gamma) return 0 // effectively zero (R returns 1e-99)
    const y = -Math.log(gamma - y0)
    const mu = polyAsc(SW_C3, n)
    const sigma = Math.exp(polyAsc(SW_C4, n))
    z = (y - mu) / sigma
  } else {
    // Royston (1995) large-sample approximation, variable = log(n)
    const xx = Math.log(n)
    const mu = polyAsc(SW_C5, xx)
    const sigma = Math.exp(polyAsc(SW_C6, xx))
    z = (y0 - mu) / sigma
  }

  // Upper tail: large W → normal → large p
  return Math.max(0, Math.min(1, 1 - normalCDF(z)))
}

// ─── Main descriptive function ────────────────────────────────────────────

/**
 * Compute full descriptive statistics for a numeric vector.
 *
 * Cross-validated with R:
 * > x <- c(2,4,4,4,5,5,7,9)
 * > mean(x)  # 5
 * > sd(x)    # 2
 * > e1071::skewness(x, type=2)  # 0.4895
 */
export function describe(x: readonly number[], ciLevel = 0.95): DescriptiveResult {
  if (x.length === 0) throw new Error('describe: empty array')
  const n = x.length
  const m = _mean(x)
  const med = _median(x)
  const modes = mode(x)
  const tm = trimmedMean(x, 0.05)
  const s = _sd(x)
  const sem = _se(x)
  const v = _variance(x)
  const sorted = sortAsc(x)
  const mn = sorted[0]!
  const mx = sorted[n - 1]!
  const q1 = quantile(x, 0.25)
  const q3 = quantile(x, 0.75)
  const iqr = q3 - q1
  const skew = n >= 3 ? skewness(x) : 0
  const kurt = n >= 4 ? kurtosis(x) : 0
  const ci = ciMean(x, ciLevel)
  const sw = n >= 3 ? shapiroWilk(x) : { statistic: NaN, pValue: NaN }

  const ciPct = Math.round(ciLevel * 100)
  const formatted = [
    `n = ${n}`,
    `M = ${roundTo(m, 2)}, SD = ${roundTo(s, 2)}, SE = ${roundTo(sem, 2)}`,
    `Mdn = ${roundTo(med, 2)}, IQR = ${roundTo(iqr, 2)}`,
    `Skew = ${roundTo(skew, 2)}, Kurt = ${roundTo(kurt, 2)}`,
    `${ciPct}% CI [${roundTo(ci[0], 2)}, ${roundTo(ci[1], 2)}]`,
    n >= 3 ? `Shapiro-Wilk W = ${roundTo(sw.statistic, 3)}, ${formatP(sw.pValue)}` : '',
  ].filter(Boolean).join('; ')

  return {
    n,
    mean: roundTo(m, 6),
    median: roundTo(med, 6),
    mode: modes,
    trimmedMean: roundTo(tm, 6),
    sd: roundTo(s, 6),
    se: roundTo(sem, 6),
    variance: roundTo(v, 6),
    min: mn,
    max: mx,
    range: mx - mn,
    iqr: roundTo(iqr, 6),
    q1: roundTo(q1, 6),
    q3: roundTo(q3, 6),
    skewness: roundTo(skew, 6),
    kurtosis: roundTo(kurt, 6),
    ci,
    ciLevel,
    shapiroWilk: sw,
    formatted,
  }
}

// Named re-export of math utilities for external use
export { _mean as mean, _median as median, _sd as sd, _se as se, _variance as variance, quantile }
