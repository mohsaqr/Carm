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

/**
 * Upper-tail probability P(Z > z) — numerically stable for large z.
 * Uses erfc directly to avoid catastrophic cancellation in `1 - normalCDF(z)`.
 * Matches R's `pnorm(z, lower.tail = FALSE)`.
 */
export function normalSurvival(z: number): number {
  if (z >= 0) {
    // P(Z > z) = 0.5 * erfc(z / sqrt(2))
    return 0.5 * erfc(z / Math.SQRT2)
  }
  // P(Z > z) for z < 0: use symmetry = 1 - 0.5*erfc(-z/sqrt(2))
  return 1 - 0.5 * erfc(-z / Math.SQRT2)
}

// ─── Exact distribution functions for rank correlations ───────────────────

/**
 * Exact CDF of Kendall's tau concordance count: P(T ≤ q).
 * T = number of concordant pairs in a permutation of n elements.
 * Uses recursive DP with memoization (R's kendall.c algorithm, Best & Gipps 1974).
 * Only valid for n < 50 (no ties). Falls back to normal approximation otherwise.
 */
export function pKendallExact(q: number, n: number): number {
  q = Math.floor(q + 1e-7)
  const maxT = n * (n - 1) / 2
  if (q < 0) return 0
  if (q > maxT) return 1

  // ckendall(k, nn): number of permutations of nn elements with exactly k concordant pairs
  const cache = new Map<number, Float64Array>()
  function ckendall(k: number, nn: number): number {
    const u = nn * (nn - 1) / 2
    if (k < 0 || k > u) return 0
    if (nn === 1) return k === 0 ? 1 : 0

    if (!cache.has(nn)) {
      const arr = new Float64Array(u + 1)
      for (let i = 0; i <= u; i++) arr[i] = -1
      cache.set(nn, arr)
    }
    const w = cache.get(nn)!
    if (w[k]! < 0) {
      let s = 0
      for (let i = 0; i < nn; i++) s += ckendall(k - i, nn - 1)
      w[k] = s
    }
    return w[k]!
  }

  let p = 0
  for (let j = 0; j <= q; j++) p += ckendall(j, n)

  // Divide by n!
  let nfact = 1
  for (let i = 2; i <= n; i++) nfact *= i
  return p / nfact
}

/**
 * Exact/Edgeworth CDF for Spearman's rho statistic D = Σ(d_i²).
 * For n ≤ 9: exact enumeration of all n! permutations (AS89).
 * For n > 9: Edgeworth series expansion (AS89 large-n approximation).
 * Reference: Best & Roberts (1975), R's prho.c
 *
 * @param D - Sum of squared rank differences Σ(rank_x_i - rank_y_i)²
 * @param n - Sample size
 * @param lowerTail - If true, returns P(D' ≤ D); if false, P(D' ≥ D)
 */
export function pSpearmanExact(D: number, n: number, lowerTail: boolean): number {
  const n3 = n * (n * n - 1) / 3 // max possible D

  if (n <= 1) return lowerTail ? 0 : 1
  if (D <= 0) return lowerTail ? 0 : 1
  if (D > n3) return lowerTail ? 1 : 0

  if (n <= 9) {
    // Exact enumeration — generate all n! permutations
    const perm = new Array<number>(n)
    for (let i = 0; i < n; i++) perm[i] = i + 1

    let nfac = 1
    for (let i = 1; i <= n; i++) nfac *= i

    let ifr = 0 // count of perms where D' ≥ D
    // Use the same permutation algorithm as R's prho.c
    for (let m = 0; m < nfac; m++) {
      // Compute sum of squared differences for current permutation
      let ise = 0
      for (let i = 0; i < n; i++) {
        const d = i + 1 - perm[i]!
        ise += d * d
      }
      if (D <= ise) ifr++

      // Advance to next permutation (R's cyclic rotation method)
      let n1 = n
      do {
        const mt = perm[0]!
        for (let i = 1; i < n1; i++) perm[i - 1] = perm[i]!
        n1--
        perm[n1] = mt
        if (mt !== n1 + 1 || n1 <= 1) break
      } while (true)
    }

    return (lowerTail ? nfac - ifr : ifr) / nfac
  }

  // Edgeworth series expansion for n > 9 (from R's prho.c, AS89)
  const c1 = 0.2274, c2 = 0.2531, c3 = 0.1745, c4 = 0.0758
  const c5 = 0.1033, c6 = 0.3932, c7 = 0.0879, c8 = 0.0151
  const c9 = 0.0072, c10 = 0.0831, c11 = 0.0131, c12 = 4.6e-4

  const yn = n
  const b = 1 / yn
  const x = (6 * (D - 1) * b / (yn * yn - 1) - 1) * Math.sqrt(yn - 1)
  let y = x * x
  const u = x * b * (c1 + b * (c2 + c3 * b) +
    y * (-c4 + b * (c5 + c6 * b) -
      y * b * (c7 + c8 * b -
        y * (c9 - c10 * b + y * b * (c11 - c12 * y)))))
  y = x * x // restore y
  const correction = u / Math.exp(y / 2)

  // pnorm(x, lower_tail) + correction
  let pv: number
  if (lowerTail) {
    pv = -correction + normalCDF(x)
  } else {
    pv = correction + normalSurvival(x)
  }
  if (pv < 0) pv = 0
  if (pv > 1) pv = 1
  return pv
}

/**
 * Inverse standard normal CDF (quantile function).
 * Peter Acklam's rational approximation — max absolute error ~1.15e-9.
 * Reference: https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/
 */
export function normalQuantile(p: number): number {
  if (p <= 0 || p >= 1) throw new Error(`normalQuantile: p must be in (0,1), got ${p}`)

  // Coefficients for the rational approximations
  const a1 = -3.969683028665376e+01, a2 =  2.209460984245205e+02
  const a3 = -2.759285104469687e+02, a4 =  1.383577518672690e+02
  const a5 = -3.066479806614716e+01, a6 =  2.506628277459239e+00
  const b1 = -5.447609879822406e+01, b2 =  1.615858368580409e+02
  const b3 = -1.556989798598866e+02, b4 =  6.680131188771972e+01
  const b5 = -1.328068155288572e+01
  const c1 = -7.784894002430293e-03, c2 = -3.223964580411365e-01
  const c3 = -2.400758277161838e+00, c4 = -2.549732539343734e+00
  const c5 =  4.374664141464968e+00, c6 =  2.938163982698783e+00
  const d1 =  7.784695709041462e-03, d2 =  3.224671290700398e-01
  const d3 =  2.445134137142996e+00, d4 =  3.754408661907416e+00

  const pLow = 0.02425, pHigh = 1 - pLow

  if (p < pLow) {
    // Lower tail region
    const q = Math.sqrt(-2 * Math.log(p))
    return (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
           ((((d1*q+d2)*q+d3)*q+d4)*q+1)
  } else if (p <= pHigh) {
    // Central region
    const q = p - 0.5, r = q * q
    return (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q /
           (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
  } else {
    // Upper tail region (mirror of lower)
    const q = Math.sqrt(-2 * Math.log(1 - p))
    return -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
            ((((d1*q+d2)*q+d3)*q+d4)*q+1)
  }
}

/** Error function approximation (Horner's method). */
function erf(x: number): number {
  // Abramowitz & Stegun formula 7.1.26
  const t = 1 / (1 + 0.3275911 * Math.abs(x))
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
  const result = 1 - poly * Math.exp(-x * x)
  return x >= 0 ? result : -result
}

/** Complementary error function erfc(x) = 1 - erf(x), stable for large x. */
function erfc(x: number): number {
  // For x >= 0: erfc(x) = poly * exp(-x²), directly from A&S 7.1.26
  const ax = Math.abs(x)
  const t = 1 / (1 + 0.3275911 * ax)
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
  const val = poly * Math.exp(-ax * ax)
  return x >= 0 ? val : 2 - val
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
  // Explicit Number() coercion guards against JS string-concatenation if TypeScript
  // safety is bypassed (e.g. CSV-sourced data passed without type-checking).
  return x.reduce((s, v) => s + Number(v), 0) / x.length
}

/** Sample variance (n-1 denominator). */
export function variance(x: readonly number[]): number {
  if (x.length < 2) throw new Error('variance: need at least 2 observations')
  const m = mean(x)
  // If the mean is non-finite (e.g. any value is Infinity), return NaN explicitly
  // rather than letting Infinity - Infinity produce an opaque NaN mid-loop.
  if (!isFinite(m)) return NaN
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

// ─── Studentized Range Distribution (ptukey) ─────────────────────────────

/**
 * Standard normal PDF.
 */
function normalPDF(z: number): number {
  return Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI)
}

/**
 * CDF of the studentized range distribution P(Q ≤ q | k, df).
 *
 * Uses multi-interval Gauss-Legendre quadrature over the variable
 * v = 2χ²/df which has density Gamma(df/2, rate=df/4) with mean=2.
 * The integral: P(Q ≤ q | k, df) = ∫_0^∞ f_v(v) · P_∞(q·√(v/2)) dv
 *
 * Follows R's Copenhaver & Holland (1988) approach: adaptive subinterval
 * widths based on df, 16-point Gauss-Legendre per subinterval.
 *
 * Reference: Copenhaver & Holland (1988), JRSS-B 50:36-45;
 * R source code src/nmath/ptukey.c
 */
export function ptukeyApprox(q: number, k: number, df: number): number {
  if (q <= 0) return 0
  if (!isFinite(q) || isNaN(q)) return NaN
  if (k < 2) return NaN

  // For very large df, use the normal (df=∞) case directly
  if (df > 5000) return ptukeyInfDf(q, k)

  // Determine subinterval width based on df (narrower for larger df
  // so the chi-square peak near v=2 is adequately resolved)
  let wStep: number
  if (df <= 100) wStep = 1.0
  else if (df <= 800) wStep = 0.5
  else if (df <= 5000) wStep = 0.25
  else wStep = 0.125

  // Density of v = 2χ²/df: v ~ Gamma(shape=df/2, rate=df/4)
  // f_v(v) = (df/4)^(df/2) / Gamma(df/2) · v^(df/2-1) · exp(-df·v/4)
  // Mean=2, Variance=8/df
  const halfDf = df / 2
  const rate = df / 4
  const logNorm = halfDf * Math.log(rate) - logGamma(halfDf)

  // 16-point Gauss-Legendre per subinterval (matches R)
  const glPts = 16
  const gl = gaussLegendreRef(glPts)

  let result = 0
  let vLo = 0

  // Integrate subintervals until contribution is negligible
  for (let interval = 0; interval < 200; interval++) {
    const vHi = vLo + wStep
    const mid = (vLo + vHi) / 2
    const half = (vHi - vLo) / 2

    let subResult = 0
    for (let i = 0; i < glPts; i++) {
      const v = mid + half * gl.nodes[i]!
      const w = half * gl.weights[i]!
      if (v <= 0) continue

      // Gamma(df/2, rate=df/4) density
      const logDens = logNorm + (halfDf - 1) * Math.log(v) - rate * v
      const dens = Math.exp(logDens)
      if (dens < 1e-300) continue

      // P_∞(q · √(v/2))  where √(v/2) = √(χ²/df)
      const scaledQ = q * Math.sqrt(v / 2)
      const pInf = ptukeyInfDf(scaledQ, k)

      subResult += w * dens * pInf
    }

    result += subResult
    vLo = vHi

    // Stop when this subinterval contributed negligibly AND we're past the peak
    if (vLo > 2 && Math.abs(subResult) < 1e-10) break
  }

  return Math.max(0, Math.min(1, result))
}

/**
 * CDF of studentized range for df = ∞ (normal case).
 * P(Q ≤ q | k, ∞) = k * ∫_{-∞}^{∞} φ(z) * [Φ(z + q) - Φ(z)]^{k-1} dz
 *
 * Uses Gauss-Hermite quadrature for the infinite-range integral.
 */
function ptukeyInfDf(q: number, k: number): number {
  if (q <= 0) return 0

  // Use Gauss-Legendre on a truncated range. φ(z) is negligible for |z| > 8
  const lo = -8
  const hi = 8
  const npts = 64

  const { nodes, weights } = gaussLegendreAB(npts, lo, hi)

  let result = 0
  for (let i = 0; i < npts; i++) {
    const z = nodes[i]!
    const w = weights[i]!
    const phi = normalPDF(z)
    const diff = normalCDF(z + q) - normalCDF(z)
    if (diff <= 0) continue
    result += w * phi * Math.pow(diff, k - 1)
  }

  return Math.max(0, Math.min(1, k * result))
}

/**
 * Upper-tail p-value for studentized range distribution.
 * P(Q > q | k, df)
 */
export function pValueStudentizedRangeApprox(q: number, k: number, df: number): number {
  return Math.max(0, Math.min(1, 1 - ptukeyApprox(q, k, df)))
}

// ─── Gauss-Legendre Quadrature Nodes & Weights ──────────────────────────

/**
 * Gauss-Legendre nodes and weights on [a, b].
 * Computes from [-1, 1] reference nodes via affine transformation.
 */
function gaussLegendreAB(n: number, a: number, b: number): { nodes: number[]; weights: number[] } {
  const ref = gaussLegendreRef(n)
  const mid = (a + b) / 2
  const half = (b - a) / 2
  return {
    nodes: ref.nodes.map(t => mid + half * t),
    weights: ref.weights.map(w => half * w),
  }
}

/**
 * Gauss-Legendre nodes and weights on [-1, 1].
 * Uses Newton iteration on Legendre polynomials.
 * Reference: Golub & Welsch (1969)
 */
function gaussLegendreRef(n: number): { nodes: number[]; weights: number[] } {
  const nodes: number[] = []
  const weights: number[] = []
  const m = Math.ceil(n / 2)

  for (let i = 0; i < m; i++) {
    // Initial guess from Tricomi approximation
    let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5))

    let pp = 0
    for (let iter = 0; iter < 100; iter++) {
      let p0 = 1
      let p1 = z
      for (let j = 2; j <= n; j++) {
        const p2 = ((2 * j - 1) * z * p1 - (j - 1) * p0) / j
        p0 = p1
        p1 = p2
      }
      pp = n * (z * p1 - p0) / (z * z - 1)
      const zOld = z
      z = z - p1 / pp
      if (Math.abs(z - zOld) < 1e-15) break
    }

    nodes.push(z)
    nodes.push(-z)
    const w = 2 / ((1 - z * z) * pp * pp)
    weights.push(w)
    weights.push(w)
  }

  // If n is odd, the center node (z=0) was added twice — fix
  if (n % 2 === 1) {
    nodes.pop()
    weights.pop()
  }

  return { nodes, weights }
}
