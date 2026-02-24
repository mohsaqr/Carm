import { P as PAdjMethod } from './types-DC8rlZlK.cjs';

/**
 * Matrix class for Carm.
 * Implements multiply, transpose, inverse (via Cholesky or LU), log-determinant, and SVD.
 * Pure computation — no DOM, no D3, no side effects.
 */
declare class Matrix {
    readonly rows: number;
    readonly cols: number;
    private readonly _data;
    constructor(rows: number, cols: number, data?: readonly number[]);
    /** Build from 2-D array (row-major). */
    static fromArray(arr: readonly (readonly number[])[]): Matrix;
    /** Identity matrix of size n. */
    static identity(n: number): Matrix;
    /** Zero matrix. */
    static zeros(rows: number, cols: number): Matrix;
    /** Get element at (i, j) — 0-indexed. */
    get(i: number, j: number): number;
    /** Return 2-D array representation. */
    toArray(): number[][];
    /** Return flat row-major copy. */
    toFlat(): number[];
    /** Matrix transpose. */
    transpose(): Matrix;
    /** Matrix multiplication: this × other. */
    multiply(other: Matrix): Matrix;
    /** Scalar multiplication. */
    scale(s: number): Matrix;
    /** Element-wise add. */
    add(other: Matrix): Matrix;
    /** Element-wise subtract. */
    subtract(other: Matrix): Matrix;
    /**
     * Cholesky decomposition for symmetric positive-definite matrices.
     * Returns lower-triangular L such that this = L * L^T.
     * Throws if matrix is not SPD.
     */
    cholesky(): Matrix;
    /**
     * Log-determinant via Cholesky: log|A| = 2 * Σ log(L_ii).
     * Only valid for symmetric positive-definite matrices.
     */
    logDet(): number;
    /**
     * Inverse via LU decomposition with partial pivoting.
     * Works for any non-singular square matrix.
     * Formula: Doolittle LU, then forward/back substitution for each column of I.
     */
    inverse(): Matrix;
    /**
     * Singular Value Decomposition: A = U · S · V^T
     * Returns { U, S (diagonal values), V }.
     * Algorithm: Golub-Reinsch (one-sided Jacobi for small matrices).
     * Reference: Golub & Van Loan, "Matrix Computations", 4th ed., Algorithm 8.6.2
     */
    svd(): {
        U: Matrix;
        S: number[];
        V: Matrix;
    };
    /**
     * Pseudo-inverse via SVD: A+ = V · S^{-1} · U^T
     */
    pseudoInverse(tol?: number): Matrix;
    /** Trace (sum of diagonal elements). */
    trace(): number;
    /** Extract diagonal as array. */
    diagonal(): number[];
    /** Column vector as Matrix from array. */
    static colVec(data: readonly number[]): Matrix;
    /** Row vector as Matrix from array. */
    static rowVec(data: readonly number[]): Matrix;
    /**
     * Eigendecomposition for symmetric matrices via Jacobi iterations.
     * Returns { values, vectors } where vectors are columns of the eigenvector matrix.
     * Reference: Golub & Van Loan, Algorithm 8.4.1
     */
    eigen(): {
        values: number[];
        vectors: Matrix;
    };
}
/** Solve linear system A·x = b using the (already computed) inverse.  */
declare function solveLinear(A: Matrix, b: readonly number[]): number[];

/**
 * Core mathematical utilities for Carm.
 * Beta/Gamma functions, Nelder-Mead optimizer, p-value adjustment.
 * All pure functions — no DOM, no side effects.
 */

/**
 * Natural log of the Gamma function via Lanczos approximation.
 * Reference: Press et al., "Numerical Recipes in C", 3rd ed., §6.1
 */
declare function logGamma(z: number): number;
/** Gamma function. */
declare function gamma(z: number): number;
/** Log of Beta function: log B(a, b) = logΓ(a) + logΓ(b) - logΓ(a+b). */
declare function logBeta(a: number, b: number): number;
/** Beta function. */
declare function betaFn(a: number, b: number): number;
/**
 * Regularized incomplete beta function I_x(a, b) via continued fraction.
 * Used for t, F, and beta distribution CDFs.
 * Reference: Press et al., §6.4 (Lentz continued fraction algorithm)
 */
declare function incompleteBeta(x: number, a: number, b: number): number;
/**
 * Regularized incomplete gamma function P(a, x) via series expansion.
 * Used for chi-square and Poisson CDFs.
 * Reference: Press et al., §6.2
 */
declare function incompleteGamma(a: number, x: number): number;
/** Two-tailed p-value from t statistic with df degrees of freedom. */
declare function tDistPValue(t: number, df: number): number;
/** CDF of t-distribution P(T ≤ t | df). */
declare function tDistCDF(t: number, df: number): number;
/** Quantile (inverse CDF) of t-distribution. Uses bisection. */
declare function tDistQuantile(p: number, df: number): number;
/** P(F ≤ f | df1, df2) via regularized incomplete beta. */
declare function fDistCDF(f: number, df1: number, df2: number): number;
/** Upper tail p-value for F-distribution. */
declare function fDistPValue(f: number, df1: number, df2: number): number;
/** P(χ² ≤ x | df). */
declare function chiSqCDF(x: number, df: number): number;
/** Upper tail p-value for chi-square. */
declare function chiSqPValue(x: number, df: number): number;
/** Quantile of chi-square distribution via bisection. */
declare function chiSqQuantile(p: number, df: number): number;
/** Standard normal CDF using complementary error function. */
declare function normalCDF(z: number): number;
/**
 * Inverse standard normal CDF (quantile function).
 * Peter Acklam's rational approximation — max absolute error ~1.15e-9.
 * Reference: https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/
 */
declare function normalQuantile(p: number): number;
interface NelderMeadOptions {
    readonly maxIter?: number;
    readonly tol?: number;
    readonly alpha?: number;
    readonly beta?: number;
    readonly gamma?: number;
    readonly delta?: number;
}
interface NelderMeadResult {
    readonly x: readonly number[];
    readonly fval: number;
    readonly iterations: number;
    readonly converged: boolean;
}
/**
 * Nelder-Mead simplex optimizer for unconstrained minimization.
 * Reference: Nelder & Mead (1965), The Computer Journal 7(4):308-313
 */
declare function nelderMead(fn: (x: readonly number[]) => number, x0: readonly number[], opts?: NelderMeadOptions): NelderMeadResult;
/**
 * Adjust p-values for multiple comparisons.
 * Methods: 'bonferroni', 'holm', 'BH' (Benjamini-Hochberg), 'BY', 'none'.
 * Reference: R stats::p.adjust
 */
declare function adjustPValues(pValues: readonly number[], method: PAdjMethod): number[];
/** Sample mean. */
declare function mean(x: readonly number[]): number;
/** Sample variance (n-1 denominator). */
declare function variance(x: readonly number[]): number;
/** Sample standard deviation. */
declare function sd(x: readonly number[]): number;
/** Standard error of the mean. */
declare function se(x: readonly number[]): number;
/** Median. */
declare function median(x: readonly number[]): number;
/** Quantile (using linear interpolation, same as R type=7). */
declare function quantile(x: readonly number[], p: number): number;
/** Sort array ascending (returns new array). */
declare function sortAsc(x: readonly number[]): number[];
/** Sum of squared deviations from mean. */
declare function ss(x: readonly number[]): number;
/** Ranks (average ties), 1-indexed. */
declare function rank(x: readonly number[]): number[];
/** Covariance between two arrays (n-1 denominator). */
declare function cov(x: readonly number[], y: readonly number[]): number;
/** Clamp a number to [lo, hi]. */
declare function clamp(v: number, lo: number, hi: number): number;
/** Round to n decimal places. */
declare function roundTo(v: number, n: number): number;

export { sortAsc as A, ss as B, tDistCDF as C, tDistPValue as D, tDistQuantile as E, variance as F, Matrix as M, type NelderMeadOptions as N, type NelderMeadResult as a, adjustPValues as b, betaFn as c, chiSqCDF as d, chiSqPValue as e, chiSqQuantile as f, clamp as g, cov as h, fDistCDF as i, fDistPValue as j, gamma as k, incompleteBeta as l, incompleteGamma as m, logBeta as n, logGamma as o, mean as p, median as q, nelderMead as r, normalCDF as s, normalQuantile as t, quantile as u, rank as v, roundTo as w, sd as x, se as y, solveLinear as z };
