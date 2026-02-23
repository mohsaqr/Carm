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

export type { AnalysisResult as A, DataMatrix as D, EffectInterpretation as E, Field as F, GroupData as G, LMMResult as L, NumericField as N, PAdjMethod as P, RegressionCoef as R, StatResult as S, AnalyzeOptions as a, DescriptiveResult as b, EffectSize as c, FieldType as d, FixedEffect as e, FrequencyRow as f, FrequencyTestResult as g, GroupField as h, PCAResult as i, PairwiseResult as j, RegressionResult as k };
