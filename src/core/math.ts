/**
 * Core mathematical utilities for Carm.
 * Beta/Gamma functions, Nelder-Mead optimizer, p-value adjustment.
 * All pure functions — no DOM, no side effects.
 */

import type { PAdjMethod } from './types.js'

// ─── Gamma / Beta functions ────────────────────────────────────────────────

/**
 * Natural log of the Gamma function via Lanczos approximation.
 * Reference: Press et al., "Numerical Recipes in C", 3rd ed., §6.1
 */
export function logGamma(z: number): number {
  if (z <= 0) throw new Error(`logGamma: z must be positive, got ${z}`)
  if (z < 0.5) {
    // Reflection formula: Γ(z)Γ(1-z) = π/sin(πz)
    return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z)
  }
  const g = 7
  const c = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
  ]
  const x = z - 1
  let sum = c[0]!
  for (let i = 1; i < g + 2; i++) sum += (c[i] ?? 0) / (x + i)
  const t = x + g + 0.5
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(sum)
}

/** Gamma function. */
export function gamma(z: number): number {
  return Math.exp(logGamma(z))
}

/** Log of Beta function: log B(a, b) = logΓ(a) + logΓ(b) - logΓ(a+b). */
export function logBeta(a: number, b: number): number {
  return logGamma(a) + logGamma(b) - logGamma(a + b)
}

/** Beta function. */
export function betaFn(a: number, b: number): number {
  return Math.exp(logBeta(a, b))
}

/**
 * Regularized incomplete beta function I_x(a, b) via continued fraction.
 * Used for t, F, and beta distribution CDFs.
 * Reference: Press et al., §6.4 (Lentz continued fraction algorithm)
 */
export function incompleteBeta(x: number, a: number, b: number): number {
  if (x < 0 || x > 1) throw new Error(`incompleteBeta: x must be in [0,1], got ${x}`)
  if (x === 0) return 0
  if (x === 1) return 1

  // Use symmetry when x > (a+1)/(a+b+2) for faster convergence
  if (x > (a + 1) / (a + b + 2)) {
    return 1 - incompleteBeta(1 - x, b, a)
  }

  const lbeta_ab = logBeta(a, b)
  const front = Math.exp(Math.log(x) * a + Math.log(1 - x) * b - lbeta_ab) / a

  // Lentz continued fraction
  const MAX_ITER = 200
  const EPS = 3e-7
  const FPMIN = 1e-30

  let f = 1
  let C = 1
  let D = 1 - (a + b) * x / (a + 1)
  if (Math.abs(D) < FPMIN) D = FPMIN
  D = 1 / D
  f = D

  for (let m = 1; m <= MAX_ITER; m++) {
    // Even step
    let num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
    D = 1 + num * D
    if (Math.abs(D) < FPMIN) D = FPMIN
    C = 1 + num / C
    if (Math.abs(C) < FPMIN) C = FPMIN
    D = 1 / D
    f *= C * D

    // Odd step
    num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
    D = 1 + num * D
    if (Math.abs(D) < FPMIN) D = FPMIN
    C = 1 + num / C
    if (Math.abs(C) < FPMIN) C = FPMIN
    D = 1 / D
    const delta = C * D
    f *= delta

    if (Math.abs(delta - 1) < EPS) break
  }
  return front * f
}

/**
 * Regularized incomplete gamma function P(a, x) via series expansion.
 * Used for chi-square and Poisson CDFs.
 * Reference: Press et al., §6.2
 */
export function incompleteGamma(a: number, x: number): number {
  if (x < 0) throw new Error(`incompleteGamma: x must be ≥ 0, got ${x}`)
  if (x === 0) return 0

  if (x < a + 1) {
    // Series representation
    let term = 1 / a
    let sum = term
    for (let n = 1; n < 200; n++) {
      term *= x / (a + n)
      sum += term
      if (Math.abs(term) < Math.abs(sum) * 3e-7) break
    }
    return sum * Math.exp(-x + a * Math.log(x) - logGamma(a))
  } else {
    // Continued fraction (complement)
    return 1 - incompleteGammaComplement(a, x)
  }
}

function incompleteGammaComplement(a: number, x: number): number {
  // Lentz continued fraction for Q(a,x)
  const FPMIN = 1e-30
  let f = x + 1 - a
  if (Math.abs(f) < FPMIN) f = FPMIN
  let C = f, D = 0
  for (let i = 1; i <= 200; i++) {
    const an = -i * (i - a)
    const bn = x + 2 * i + 1 - a
    D = bn + an * D
    if (Math.abs(D) < FPMIN) D = FPMIN
    C = bn + an / C
    if (Math.abs(C) < FPMIN) C = FPMIN
    D = 1 / D
    const delta = C * D
    f *= delta
    if (Math.abs(delta - 1) < 3e-7) break
  }
  return Math.exp(-x + a * Math.log(x) - logGamma(a)) / f
}

// ─── t-distribution CDF ───────────────────────────────────────────────────

/** Two-tailed p-value from t statistic with df degrees of freedom. */
export function tDistPValue(t: number, df: number): number {
  const x = df / (df + t * t)
  const p = incompleteBeta(x, df / 2, 0.5)
  return p  // already two-tailed
}

/** CDF of t-distribution P(T ≤ t | df). */
export function tDistCDF(t: number, df: number): number {
  const x = df / (df + t * t)
  const p = 0.5 * incompleteBeta(x, df / 2, 0.5)
  return t >= 0 ? 1 - p : p
}

/** Quantile (inverse CDF) of t-distribution. Uses bisection. */
export function tDistQuantile(p: number, df: number): number {
  if (p <= 0 || p >= 1) throw new Error(`tDistQuantile: p must be in (0,1), got ${p}`)
  // Bisection on CDF
  let lo = -50, hi = 50
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2
    if (tDistCDF(mid, df) < p) lo = mid
    else hi = mid
    if (hi - lo < 1e-10) break
  }
  return (lo + hi) / 2
}

// ─── F-distribution CDF ──────────────────────────────────────────────────

/** P(F ≤ f | df1, df2) via regularized incomplete beta. */
export function fDistCDF(f: number, df1: number, df2: number): number {
  if (f <= 0) return 0
  const x = df1 * f / (df1 * f + df2)
  return incompleteBeta(x, df1 / 2, df2 / 2)
}

/** Upper tail p-value for F-distribution. */
export function fDistPValue(f: number, df1: number, df2: number): number {
  return 1 - fDistCDF(f, df1, df2)
}

// ─── Chi-square distribution ─────────────────────────────────────────────

/** P(χ² ≤ x | df). */
export function chiSqCDF(x: number, df: number): number {
  if (x <= 0) return 0
  return incompleteGamma(df / 2, x / 2)
}

/** Upper tail p-value for chi-square. */
export function chiSqPValue(x: number, df: number): number {
  return 1 - chiSqCDF(x, df)
}

/** Quantile of chi-square distribution via bisection. */
export function chiSqQuantile(p: number, df: number): number {
  if (p <= 0) return 0
  if (p >= 1) return Infinity
  let lo = 0, hi = df + 10 * Math.sqrt(2 * df)
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2
    if (chiSqCDF(mid, df) < p) lo = mid
    else hi = mid
    if (hi - lo < 1e-10) break
  }
  return (lo + hi) / 2
}

// ─── Normal distribution ──────────────────────────────────────────────────

/** Standard normal CDF using complementary error function. */
export function normalCDF(z: number): number {
  return 0.5 * (1 + erf(z / Math.SQRT2))
}

/** Inverse standard normal CDF via rational approximation (Beasley-Springer-Moro). */
export function normalQuantile(p: number): number {
  if (p <= 0 || p >= 1) throw new Error(`normalQuantile: p must be in (0,1), got ${p}`)
  // Rational approximation from Abramowitz & Stegun §26.2.16
  const a = [2.515517, 0.802853, 0.010328]
  const b = [1.432788, 0.189269, 0.001308]
  const t = Math.sqrt(-2 * Math.log(p <= 0.5 ? p : 1 - p))
  const num = a[0]! + a[1]! * t + a[2]! * t * t
  const den = 1 + b[0]! * t + b[1]! * t * t + b[2]! * t * t * t
  const x = t - num / den
  return p <= 0.5 ? -x : x
}

/** Error function approximation (Horner's method). */
function erf(x: number): number {
  // Abramowitz & Stegun formula 7.1.26
  const t = 1 / (1 + 0.3275911 * Math.abs(x))
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
  const result = 1 - poly * Math.exp(-x * x)
  return x >= 0 ? result : -result
}

// ─── Nelder-Mead Optimizer ────────────────────────────────────────────────

export interface NelderMeadOptions {
  readonly maxIter?: number
  readonly tol?: number
  readonly alpha?: number   // reflection
  readonly beta?: number    // contraction
  readonly gamma?: number   // expansion
  readonly delta?: number   // shrink
}

export interface NelderMeadResult {
  readonly x: readonly number[]
  readonly fval: number
  readonly iterations: number
  readonly converged: boolean
}

/**
 * Nelder-Mead simplex optimizer for unconstrained minimization.
 * Reference: Nelder & Mead (1965), The Computer Journal 7(4):308-313
 */
export function nelderMead(
  fn: (x: readonly number[]) => number,
  x0: readonly number[],
  opts: NelderMeadOptions = {}
): NelderMeadResult {
  const n = x0.length
  const maxIter = opts.maxIter ?? 1000 * n
  const tol = opts.tol ?? 1e-8
  const alpha = opts.alpha ?? 1.0
  const beta = opts.beta ?? 0.5
  const gam = opts.gamma ?? 2.0
  const delta = opts.delta ?? 0.5

  // Initialize simplex: n+1 vertices
  const simplex: number[][] = Array.from({ length: n + 1 }, (_, i) => {
    const v = [...x0]
    if (i > 0) v[i - 1] = (v[i - 1] ?? 0) + ((Math.abs(v[i - 1] ?? 0) > 1e-8) ? 0.05 * (v[i - 1] ?? 0) : 0.00025)
    return v
  })
  let fvals = simplex.map(v => fn(v))

  let iter = 0
  let converged = false

  for (; iter < maxIter; iter++) {
    // Sort by function value and keep simplex in sync
    const order = fvals.map((_, i) => i).sort((a, b) => fvals[a]! - fvals[b]!)
    const sorted = order.map(i => simplex[i]!)
    fvals = order.map(i => fvals[i]!)
    // Critical: sync simplex so simplex[i] corresponds to fvals[i]
    for (let i = 0; i <= n; i++) simplex[i] = sorted[i]!

    // Check convergence
    const fRange = (fvals[n] ?? 0) - (fvals[0] ?? 0)
    if (fRange < tol) { converged = true; break }

    // Centroid of all but worst (simplex[n] is now the worst)
    const centroid = Array.from({ length: n }, (_, j) =>
      simplex.slice(0, n).reduce((s, v) => s + (v[j] ?? 0), 0) / n
    )

    // Reflect worst
    const worst = simplex[n]!
    const reflected = centroid.map((c, j) => c + alpha * (c - (worst[j] ?? 0)))
    const fr = fn(reflected)

    if (fr < (fvals[0] ?? 0)) {
      // Expansion
      const expanded = centroid.map((c, j) => c + gam * (reflected[j]! - c))
      const fe = fn(expanded)
      simplex[n] = fe < fr ? expanded : reflected
      fvals[n] = fe < fr ? fe : fr
    } else if (fr < (fvals[n - 1] ?? 0)) {
      simplex[n] = reflected
      fvals[n] = fr
    } else {
      // Contraction
      const contracted = centroid.map((c, j) => c + beta * ((worst[j] ?? 0) - c))
      const fc = fn(contracted)
      if (fc < (fvals[n] ?? 0)) {
        simplex[n] = contracted
        fvals[n] = fc
      } else {
        // Shrink toward best (simplex[0])
        const best = simplex[0]!
        for (let i = 1; i <= n; i++) {
          simplex[i] = simplex[i]!.map((v, j) => (best[j] ?? 0) + delta * (v - (best[j] ?? 0)))
          fvals[i] = fn(simplex[i]!)
        }
      }
    }
  }

  const order0 = fvals.map((_, i) => i).sort((a, b) => fvals[a]! - fvals[b]!)
  return {
    x: simplex[order0[0]!]!,
    fval: fvals[order0[0]!]!,
    iterations: iter,
    converged,
  }
}

// ─── P-value adjustment ──────────────────────────────────────────────────

/**
 * Adjust p-values for multiple comparisons.
 * Methods: 'bonferroni', 'holm', 'BH' (Benjamini-Hochberg), 'BY', 'none'.
 * Reference: R stats::p.adjust
 */
export function adjustPValues(pValues: readonly number[], method: PAdjMethod): number[] {
  const n = pValues.length
  if (n === 0) return []
  if (method === 'none') return [...pValues]

  if (method === 'bonferroni') {
    return pValues.map(p => Math.min(1, p * n))
  }

  if (method === 'holm') {
    const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => a.p - b.p)
    const adj = new Array<number>(n)
    let running = 0
    for (let k = 0; k < n; k++) {
      const { p, i } = order[k]!
      running = Math.max(running, p * (n - k))
      adj[i] = Math.min(1, running)
    }
    return adj
  }

  if (method === 'BH') {
    // Benjamini-Hochberg
    const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => b.p - a.p) // descending
    const adj = new Array<number>(n)
    let minSoFar = Infinity
    for (let k = 0; k < n; k++) {
      const { p, i } = order[k]!
      const rank = n - k  // rank from smallest (1-indexed)
      const adjusted = p * n / rank
      minSoFar = Math.min(minSoFar, adjusted)
      adj[i] = Math.min(1, minSoFar)
    }
    return adj
  }

  if (method === 'BY') {
    // Benjamini-Yekutieli
    const cm = pValues.reduce((s, _, k) => s + 1 / (k + 1), 0)
    const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => b.p - a.p)
    const adj = new Array<number>(n)
    let minSoFar = Infinity
    for (let k = 0; k < n; k++) {
      const { p, i } = order[k]!
      const rank = n - k
      const adjusted = p * n * cm / rank
      minSoFar = Math.min(minSoFar, adjusted)
      adj[i] = Math.min(1, minSoFar)
    }
    return adj
  }

  throw new Error(`Unknown p-adjust method: ${method}`)
}

// ─── Descriptive utilities ────────────────────────────────────────────────

/** Sample mean. */
export function mean(x: readonly number[]): number {
  if (x.length === 0) throw new Error('mean: empty array')
  return x.reduce((s, v) => s + v, 0) / x.length
}

/** Sample variance (n-1 denominator). */
export function variance(x: readonly number[]): number {
  if (x.length < 2) throw new Error('variance: need at least 2 observations')
  const m = mean(x)
  return x.reduce((s, v) => s + (v - m) ** 2, 0) / (x.length - 1)
}

/** Sample standard deviation. */
export function sd(x: readonly number[]): number {
  return Math.sqrt(variance(x))
}

/** Standard error of the mean. */
export function se(x: readonly number[]): number {
  return sd(x) / Math.sqrt(x.length)
}

/** Median. */
export function median(x: readonly number[]): number {
  if (x.length === 0) throw new Error('median: empty array')
  const sorted = [...x].sort((a, b) => a - b)
  const mid = Math.floor(sorted.length / 2)
  return sorted.length % 2 === 0
    ? ((sorted[mid - 1]! + sorted[mid]!) / 2)
    : sorted[mid]!
}

/** Quantile (using linear interpolation, same as R type=7). */
export function quantile(x: readonly number[], p: number): number {
  if (p < 0 || p > 1) throw new Error(`quantile: p must be in [0,1], got ${p}`)
  const sorted = [...x].sort((a, b) => a - b)
  const n = sorted.length
  const h = (n - 1) * p
  const lo = Math.floor(h)
  const hi = Math.ceil(h)
  if (lo === hi) return sorted[lo]!
  return sorted[lo]! + (h - lo) * (sorted[hi]! - sorted[lo]!)
}

/** Sort array ascending (returns new array). */
export function sortAsc(x: readonly number[]): number[] {
  return [...x].sort((a, b) => a - b)
}

/** Sum of squared deviations from mean. */
export function ss(x: readonly number[]): number {
  const m = mean(x)
  return x.reduce((s, v) => s + (v - m) ** 2, 0)
}

/** Ranks (average ties), 1-indexed. */
export function rank(x: readonly number[]): number[] {
  const n = x.length
  const order = x.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v)
  const ranks = new Array<number>(n)

  let i = 0
  while (i < n) {
    let j = i
    while (j < n && order[j]!.v === order[i]!.v) j++
    const avgRank = (i + j - 1) / 2 + 1
    for (let k = i; k < j; k++) ranks[order[k]!.i] = avgRank
    i = j
  }
  return ranks
}

/** Covariance between two arrays (n-1 denominator). */
export function cov(x: readonly number[], y: readonly number[]): number {
  if (x.length !== y.length) throw new Error('cov: arrays must have equal length')
  if (x.length < 2) throw new Error('cov: need at least 2 observations')
  const mx = mean(x), my = mean(y)
  return x.reduce((s, v, i) => s + (v - mx) * ((y[i] ?? 0) - my), 0) / (x.length - 1)
}

/** Clamp a number to [lo, hi]. */
export function clamp(v: number, lo: number, hi: number): number {
  return Math.max(lo, Math.min(hi, v))
}

/** Round to n decimal places. */
export function roundTo(v: number, n: number): number {
  const factor = 10 ** n
  return Math.round(v * factor) / factor
}
