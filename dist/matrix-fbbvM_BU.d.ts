/**
 * Core type definitions for Carm.
 * All shared interfaces live here — nothing imports from other src/ files.
 */
/** Interpretation threshold labels for effect sizes. */
type EffectInterpretation = 'negligible' | 'small' | 'medium' | 'large' | 'very large';
/** A named effect size with its value and qualitative label. */
interface EffectSize {
    readonly value: number;
    readonly name: string;
    readonly interpretation: EffectInterpretation;
}
/**
 * Universal statistical result — every test returns this shape.
 * `formatted` is the APA-style string for embedding in plots.
 */
interface StatResult {
    readonly testName: string;
    readonly statistic: number;
    readonly df: number | readonly [number, number];
    readonly pValue: number;
    readonly effectSize: EffectSize;
    readonly ci: readonly [number, number];
    readonly ciLevel: number;
    readonly n: number;
    readonly formatted: string;
}
/** Descriptive statistics summary for a single numeric vector. */
interface DescriptiveResult {
    readonly n: number;
    readonly mean: number;
    readonly median: number;
    readonly mode: readonly number[];
    readonly trimmedMean: number;
    readonly sd: number;
    readonly se: number;
    readonly variance: number;
    readonly min: number;
    readonly max: number;
    readonly range: number;
    readonly iqr: number;
    readonly q1: number;
    readonly q3: number;
    readonly skewness: number;
    readonly kurtosis: number;
    readonly ci: readonly [number, number];
    readonly ciLevel: number;
    readonly shapiroWilk: {
        statistic: number;
        pValue: number;
    };
    readonly formatted: string;
}
/** A single fixed effect from an LMM. */
interface FixedEffect {
    readonly name: string;
    readonly estimate: number;
    readonly se: number;
    readonly tValue: number;
    readonly pValue: number;
    readonly ci: readonly [number, number];
}
/** Full result from a linear mixed model. */
interface LMMResult {
    readonly fixedEffects: readonly FixedEffect[];
    readonly varianceComponents: {
        readonly intercept: number;
        readonly residual: number;
        readonly slopes?: Readonly<Record<string, number>>;
    };
    readonly icc: number;
    readonly logLik: number;
    readonly aic: number;
    readonly bic: number;
    readonly nObs: number;
    readonly nGroups: number;
    readonly formatted: string;
}
/** PCA result with loadings, scores, and variance explained. */
interface PCAResult {
    readonly loadings: readonly (readonly number[])[];
    readonly scores: readonly (readonly number[])[];
    readonly eigenvalues: readonly number[];
    readonly varianceExplained: readonly number[];
    readonly cumulativeVariance: readonly number[];
    readonly nComponents: number;
}
/** Regression coefficient with CI. */
interface RegressionCoef {
    readonly name: string;
    readonly estimate: number;
    readonly se: number;
    readonly tValue: number;
    readonly pValue: number;
    readonly ci: readonly [number, number];
}
/** Result from a linear regression (OLS). */
interface RegressionResult {
    readonly coefficients: readonly RegressionCoef[];
    readonly r2: number;
    readonly adjR2: number;
    readonly fStatistic: number;
    readonly fDf: readonly [number, number];
    readonly fPValue: number;
    readonly aic: number;
    readonly bic: number;
    readonly residuals: readonly number[];
    readonly fitted: readonly number[];
    readonly n: number;
    readonly formatted: string;
}
/** Frequency table row. */
interface FrequencyRow {
    readonly value: string | number;
    readonly count: number;
    readonly relative: number;
    readonly cumulative: number;
}
/** Result from chi-square or Fisher test. */
interface FrequencyTestResult extends StatResult {
    readonly table: readonly FrequencyRow[];
    readonly expectedCounts?: readonly (readonly number[])[];
}
/** Pairwise comparison result (post-hoc). */
interface PairwiseResult {
    readonly group1: string;
    readonly group2: string;
    readonly meanDiff: number;
    readonly se: number;
    readonly statistic: number;
    readonly pValue: number;
    readonly pValueAdj: number;
    readonly ci: readonly [number, number];
    readonly significant: boolean;
}
/** Configuration for p-value adjustment. */
type PAdjMethod = 'bonferroni' | 'holm' | 'BH' | 'BY' | 'none';
/** Generic data table: rows × columns. */
type DataMatrix = readonly (readonly number[])[];
/** Named group data for comparison functions. */
interface GroupData {
    readonly values: readonly number[];
    readonly label: string;
}
/** Fit indices for factor models (EFA and CFA). */
interface FactorFit {
    readonly chiSq: number;
    readonly df: number;
    readonly pValue: number;
    readonly rmsea: number;
    readonly rmseaCI: readonly [number, number];
    readonly cfi: number;
    readonly tli: number;
    readonly srmr: number;
    readonly aic: number;
    readonly bic: number;
}
/** A single parameter estimate with SE, z, p, and standardized value. */
interface ParameterEstimate {
    readonly estimate: number;
    readonly se: number;
    readonly z: number;
    readonly pValue: number;
    readonly stdAll: number;
}
/** Diagnostics for factor analysis adequacy: KMO, Bartlett, MAP, parallel analysis. */
interface FADiagnostics {
    readonly kmo: number;
    readonly kmoPerItem: readonly number[];
    readonly bartlett: {
        readonly chiSq: number;
        readonly df: number;
        readonly pValue: number;
    };
    readonly mapSuggested: number;
    readonly parallelEigenvalues: readonly number[];
    readonly parallelSimulated: readonly number[];
    readonly parallelSuggested: number;
}
/** Result from exploratory factor analysis (EFA). */
interface FAResult {
    readonly loadings: readonly (readonly number[])[];
    readonly standardizedLoadings: readonly (readonly number[])[];
    readonly uniqueness: readonly number[];
    readonly communalities: readonly number[];
    readonly factorCorrelations: readonly (readonly number[])[];
    readonly fit: FactorFit;
    readonly eigenvalues: readonly number[];
    readonly nFactors: number;
    readonly rotation: string;
    readonly extraction: string;
    readonly variableNames: readonly string[];
    readonly factorNames: readonly string[];
    readonly formatted: string;
}
/** Result from confirmatory factor analysis (CFA). */
interface CFAResult extends FAResult {
    readonly parameterEstimates: {
        readonly loadings: readonly (readonly ParameterEstimate[])[];
        readonly uniquenesses: readonly ParameterEstimate[];
        readonly factorCovariances: readonly (readonly ParameterEstimate[])[];
    };
    readonly model: Readonly<Record<string, readonly number[]>>;
}
/** Field type descriptor for the analyze() dispatch layer. */
type FieldType = 'numeric' | 'binary' | 'categorical' | 'ordinal';
/**
 * A numeric outcome field.
 * @field type    - Always 'numeric'
 * @field name    - Column / variable name, used in result labels
 * @field values  - The raw numeric observations
 */
interface NumericField {
    readonly type: 'numeric';
    readonly name: string;
    readonly values: readonly number[];
}
/**
 * A grouping field whose values are labels (string or 0/1 int).
 * Declared as 'binary' (exactly 2 unique values) or 'categorical' (3+).
 * @field type    - 'binary' | 'categorical'
 * @field name    - Column / variable name
 * @field values  - One label per observation, parallel to the outcome field
 */
interface GroupField {
    readonly type: 'binary' | 'categorical';
    readonly name: string;
    readonly values: readonly (string | number)[];
}
/** Union of all field kinds accepted by analyze(). */
type Field = NumericField | GroupField;
/**
 * Options bag for analyze().
 * @field ciLevel          - Confidence level for all CIs (default 0.95)
 * @field paired           - If true and predictor is binary, uses paired t-test
 *                           instead of independent. Requires both groups to be
 *                           the same length.
 * @field pAdjMethod       - Multiple-comparison correction for post-hoc tests
 *                           (default 'holm'). Passed to tukeyHSD / dunnTest.
 * @field forceTest        - Skip auto-selection and run a specific test by name:
 *                           't-test-independent' | 't-test-paired' |
 *                           'one-way-anova' | 'kruskal-wallis' |
 *                           'chi-square' | 'fisher' | 'mann-whitney' |
 *                           'wilcoxon'
 * @field equalVariances   - Passed to tTestIndependent when resolving to
 *                           't-test-independent'. Default false (Welch's).
 * @field normalityAlpha   - Shapiro-Wilk p-value threshold for auto-routing
 *                           parametric vs non-parametric. Default 0.05.
 */
interface AnalyzeOptions {
    readonly ciLevel?: number;
    readonly paired?: boolean;
    readonly pAdjMethod?: PAdjMethod;
    readonly forceTest?: string;
    readonly equalVariances?: boolean;
    readonly normalityAlpha?: number;
}
/**
 * Result envelope returned by analyze().
 * @field test          - Human-readable name of the test that was run
 * @field outcome       - Name of the outcome field
 * @field predictor     - Name of the predictor/grouping field (if any)
 * @field result        - The full statistical result (StatResult or
 *                        FrequencyTestResult, depending on test)
 * @field descriptives  - Per-group or overall DescriptiveResult (always
 *                        included for numeric outcomes)
 * @field posthoc       - Pairwise comparisons (only when 3+ groups)
 * @field normality     - Shapiro-Wilk results per group used for routing
 *                        (informational — helps explain why parametric/
 *                        non-parametric was chosen)
 */
interface AnalysisResult {
    readonly test: string;
    readonly outcome: string;
    readonly predictor?: string;
    readonly result: StatResult | FrequencyTestResult;
    readonly descriptives?: readonly DescriptiveResult[];
    readonly posthoc?: readonly PairwiseResult[];
    readonly normality?: ReadonlyArray<{
        group: string;
        W: number;
        p: number;
    }>;
}

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

export { type AnalysisResult as A, type CFAResult as C, type DataMatrix as D, type EffectInterpretation as E, type FADiagnostics as F, type GroupData as G, type LMMResult as L, Matrix as M, type NumericField as N, type PAdjMethod as P, type RegressionCoef as R, type StatResult as S, type AnalyzeOptions as a, type DescriptiveResult as b, type EffectSize as c, type FAResult as d, type FactorFit as e, type Field as f, type FieldType as g, type FixedEffect as h, type FrequencyRow as i, type FrequencyTestResult as j, type GroupField as k, type PCAResult as l, type PairwiseResult as m, type ParameterEstimate as n, type RegressionResult as o, solveLinear as s };
