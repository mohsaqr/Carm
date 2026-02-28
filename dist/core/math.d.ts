/**
 * Core mathematical utilities for Carm.
 * Beta/Gamma functions, Nelder-Mead optimizer, p-value adjustment.
 * All pure functions — no DOM, no side effects.
 */
import type { PAdjMethod } from './types.js';
/**
 * Natural log of the Gamma function via Lanczos approximation.
 * Reference: Press et al., "Numerical Recipes in C", 3rd ed., §6.1
 */
export declare function logGamma(z: number): number;
/** Gamma function. */
export declare function gamma(z: number): number;
/**
 * Digamma function ψ(x) = d/dx ln Γ(x).
 * Asymptotic series for x ≥ 8, recurrence ψ(x) = ψ(x+1) - 1/x for small x.
 * Reference: Abramowitz & Stegun §6.3.18
 */
export declare function digamma(x: number): number;
/**
 * Trigamma function ψ'(x) = d²/dx² ln Γ(x).
 * Asymptotic series for x ≥ 8, recurrence ψ'(x) = ψ'(x+1) + 1/x² for small x.
 * Reference: Abramowitz & Stegun §6.4.12
 *
 * Asymptotic: ψ'(z) = 1/z + 1/(2z²) + 1/(6z³) - 1/(30z⁵) + 1/(42z⁷) - 1/(30z⁹) + ...
 *           = 1/z + 1/(2z²) + (1/z³)·[1/6 - (1/z²)·(1/30 - (1/z²)·(1/42 - ...))]
 */
export declare function trigamma(x: number): number;
/** Log of Beta function: log B(a, b) = logΓ(a) + logΓ(b) - logΓ(a+b). */
export declare function logBeta(a: number, b: number): number;
/** Beta function. */
export declare function betaFn(a: number, b: number): number;
/**
 * Regularized incomplete beta function I_x(a, b) via continued fraction.
 * Used for t, F, and beta distribution CDFs.
 * Reference: Press et al., §6.4 (Lentz continued fraction algorithm)
 */
export declare function incompleteBeta(x: number, a: number, b: number): number;
/**
 * Regularized incomplete gamma function P(a, x) via series expansion.
 * Used for chi-square and Poisson CDFs.
 * Reference: Press et al., §6.2
 */
export declare function incompleteGamma(a: number, x: number): number;
/** Two-tailed p-value from t statistic with df degrees of freedom. */
export declare function tDistPValue(t: number, df: number): number;
/** CDF of t-distribution P(T ≤ t | df). */
export declare function tDistCDF(t: number, df: number): number;
/** Quantile (inverse CDF) of t-distribution. Uses bisection. */
export declare function tDistQuantile(p: number, df: number): number;
/** P(F ≤ f | df1, df2) via regularized incomplete beta. */
export declare function fDistCDF(f: number, df1: number, df2: number): number;
/** Upper tail p-value for F-distribution. */
export declare function fDistPValue(f: number, df1: number, df2: number): number;
/** P(χ² ≤ x | df). */
export declare function chiSqCDF(x: number, df: number): number;
/** Upper tail p-value for chi-square. */
export declare function chiSqPValue(x: number, df: number): number;
/** Quantile of chi-square distribution via bisection. */
export declare function chiSqQuantile(p: number, df: number): number;
/** Standard normal CDF using complementary error function. */
export declare function normalCDF(z: number): number;
/**
 * Upper-tail probability P(Z > z) — numerically stable for large z.
 * Uses erfc directly to avoid catastrophic cancellation in `1 - normalCDF(z)`.
 * Matches R's `pnorm(z, lower.tail = FALSE)`.
 */
export declare function normalSurvival(z: number): number;
/**
 * Exact CDF of Kendall's tau concordance count: P(T ≤ q).
 * T = number of concordant pairs in a permutation of n elements.
 * Uses recursive DP with memoization (R's kendall.c algorithm, Best & Gipps 1974).
 * Only valid for n < 50 (no ties). Falls back to normal approximation otherwise.
 */
export declare function pKendallExact(q: number, n: number): number;
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
export declare function pSpearmanExact(D: number, n: number, lowerTail: boolean): number;
/**
 * Inverse standard normal CDF (quantile function).
 * Peter Acklam's rational approximation — max absolute error ~1.15e-9.
 * Reference: https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/
 */
export declare function normalQuantile(p: number): number;
export interface NelderMeadOptions {
    readonly maxIter?: number;
    readonly tol?: number;
    readonly alpha?: number;
    readonly beta?: number;
    readonly gamma?: number;
    readonly delta?: number;
}
export interface NelderMeadResult {
    readonly x: readonly number[];
    readonly fval: number;
    readonly iterations: number;
    readonly converged: boolean;
}
/**
 * Nelder-Mead simplex optimizer for unconstrained minimization.
 * Reference: Nelder & Mead (1965), The Computer Journal 7(4):308-313
 */
export declare function nelderMead(fn: (x: readonly number[]) => number, x0: readonly number[], opts?: NelderMeadOptions): NelderMeadResult;
/**
 * Adjust p-values for multiple comparisons.
 * Methods: 'bonferroni', 'holm', 'BH' (Benjamini-Hochberg), 'BY', 'none'.
 * Reference: R stats::p.adjust
 */
export declare function adjustPValues(pValues: readonly number[], method: PAdjMethod): number[];
/** Sample mean. */
export declare function mean(x: readonly number[]): number;
/** Sample variance (n-1 denominator). */
export declare function variance(x: readonly number[]): number;
/** Sample standard deviation. */
export declare function sd(x: readonly number[]): number;
/** Standard error of the mean. */
export declare function se(x: readonly number[]): number;
/** Median. */
export declare function median(x: readonly number[]): number;
/** Quantile (using linear interpolation, same as R type=7). */
export declare function quantile(x: readonly number[], p: number): number;
/** Sort array ascending (returns new array). */
export declare function sortAsc(x: readonly number[]): number[];
/** Sum of squared deviations from mean. */
export declare function ss(x: readonly number[]): number;
/** Ranks (average ties), 1-indexed. */
export declare function rank(x: readonly number[]): number[];
/** Covariance between two arrays (n-1 denominator). */
export declare function cov(x: readonly number[], y: readonly number[]): number;
/** Clamp a number to [lo, hi]. */
export declare function clamp(v: number, lo: number, hi: number): number;
/** Round to n decimal places. */
export declare function roundTo(v: number, n: number): number;
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
export declare function ptukeyApprox(q: number, k: number, df: number): number;
/**
 * Upper-tail p-value for studentized range distribution.
 * P(Q > q | k, df)
 */
export declare function pValueStudentizedRangeApprox(q: number, k: number, df: number): number;
//# sourceMappingURL=math.d.ts.map