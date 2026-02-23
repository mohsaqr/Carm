import { E as EffectInterpretation } from '../types-DC8rlZlK.cjs';
export { A as AnalysisResult, a as AnalyzeOptions, D as DataMatrix, b as DescriptiveResult, c as EffectSize, F as Field, d as FieldType, e as FixedEffect, f as FrequencyRow, g as FrequencyTestResult, G as GroupData, h as GroupField, L as LMMResult, N as NumericField, P as PAdjMethod, i as PCAResult, j as PairwiseResult, R as RegressionCoef, k as RegressionResult, S as StatResult } from '../types-DC8rlZlK.cjs';
export { N as NelderMeadOptions, a as NelderMeadResult, b as adjustPValues, c as betaFn, d as chiSqCDF, e as chiSqPValue, f as chiSqQuantile, g as clamp, h as cov, i as fDistCDF, j as fDistPValue, k as gamma, l as incompleteBeta, m as incompleteGamma, n as logBeta, o as logGamma, p as mean, q as median, r as nelderMead, s as normalCDF, t as normalQuantile, u as quantile, v as rank, w as roundTo, x as sd, y as se, z as sortAsc, A as ss, B as tDistCDF, C as tDistPValue, D as tDistQuantile, E as variance } from '../math-g4nrtyHp.cjs';

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
 * APA 7th edition string formatting for statistical results.
 * Every stat function uses these formatters to produce the `formatted` string
 * that gets embedded directly into plots as subtitles/captions.
 */
/**
 * Format a p-value in APA style.
 * < .001 → "p < .001"
 * otherwise → "p = .025" (leading zero dropped, 3 decimal places)
 */
declare function formatP(p: number): string;
/**
 * Format a statistic value: 2 decimal places, no leading zero for absolute values.
 * Used for t, F, z, χ², etc.
 */
declare function formatStat(v: number, decimals?: number): string;
/**
 * Format a confidence interval: [lo, hi] → "[0.08, 1.22]"
 */
declare function formatCI(ci: readonly [number, number], decimals?: number): string;
/**
 * Format degrees of freedom — integer if whole number, else 2 decimal places.
 */
declare function formatDF(df: number | readonly [number, number]): string;
/**
 * Full APA string for a t-test result.
 * e.g. "t(48) = 2.31, p = .025, d = 0.65, 95% CI [0.08, 1.22]"
 */
declare function formatTTest(t: number, df: number, pValue: number, effectSize: number, effectName: string, ci: readonly [number, number], ciLevel?: number): string;
/**
 * Full APA string for ANOVA.
 * e.g. "F(2, 87) = 4.21, p = .018, η² = 0.09"
 */
declare function formatANOVA(F: number, df1: number, df2: number, pValue: number, effectSize: number, effectName?: string): string;
/**
 * Full APA string for chi-square test.
 * e.g. "χ²(3) = 8.46, p = .037, V = 0.29"
 */
declare function formatChiSq(chiSq: number, df: number, pValue: number, effectSize: number, effectName?: string): string;
/**
 * APA string for Pearson/Spearman/Kendall correlation.
 * e.g. "r(48) = 0.54, p = .003, 95% CI [0.28, 0.73]"
 */
declare function formatCorrelation(r: number, df: number, pValue: number, ci: readonly [number, number], method?: string, ciLevel?: number): string;
/**
 * APA string for linear regression.
 * e.g. "R² = 0.42, F(3, 96) = 23.2, p < .001"
 */
declare function formatRegression(r2: number, adjR2: number, F: number, df1: number, df2: number, pValue: number): string;
/**
 * APA string for Mann-Whitney / Wilcoxon.
 * e.g. "W = 423, p = .031, r = 0.27"
 */
declare function formatMannWhitney(W: number, pValue: number, r: number): string;
/**
 * APA string for Kruskal-Wallis.
 * e.g. "H(2) = 12.3, p = .002, η²_H = 0.18"
 */
declare function formatKruskalWallis(H: number, df: number, pValue: number, etaSq: number): string;
/**
 * APA string for LMM.
 * e.g. "ICC = 0.42, AIC = 1234.5"
 */
declare function formatLMM(icc: number, aic: number, bic: number, logLik: number): string;
/**
 * Interpret Cohen's d effect size.
 */
declare function interpretCohensD(d: number): string;
/**
 * Interpret eta-squared / omega-squared.
 */
declare function interpretEtaSq(eta: number): string;
/**
 * Interpret Pearson r.
 */
declare function interpretR(r: number): string;
/**
 * Interpret Cramér's V.
 */
declare function interpretCramerV(v: number, df: number): string;

/** Generic effect size interpretation helper. */
declare function interpretEffect(value: number, thresholds: readonly [number, number, number]): EffectInterpretation;

export { EffectInterpretation, Matrix, formatANOVA, formatCI, formatChiSq, formatCorrelation, formatDF, formatKruskalWallis, formatLMM, formatMannWhitney, formatP, formatRegression, formatStat, formatTTest, interpretCohensD, interpretCramerV, interpretEffect, interpretEtaSq, interpretR, solveLinear };
