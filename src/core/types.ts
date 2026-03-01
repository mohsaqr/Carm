/**
 * Core type definitions for Carm.
 * All shared interfaces live here — nothing imports from other src/ files.
 */

/** Interpretation threshold labels for effect sizes. */
export type EffectInterpretation = 'negligible' | 'small' | 'medium' | 'large' | 'very large'

/** A named effect size with its value and qualitative label. */
export interface EffectSize {
  readonly value: number
  readonly name: string          // "Cohen's d", "eta²", "Cramér's V", etc.
  readonly interpretation: EffectInterpretation
}

/**
 * Universal statistical result — every test returns this shape.
 * `formatted` is the APA-style string for embedding in plots.
 */
export interface StatResult {
  readonly testName: string
  readonly statistic: number
  readonly df: number | readonly [number, number]
  readonly pValue: number
  readonly effectSize: EffectSize
  readonly ci: readonly [number, number]
  readonly ciLevel: number       // e.g. 0.95
  readonly n: number
  readonly formatted: string     // "t(48) = 2.31, p = .025, d = 0.65, 95% CI [0.08, 1.22]"
}

/** Descriptive statistics summary for a single numeric vector. */
export interface DescriptiveResult {
  readonly n: number
  readonly mean: number
  readonly median: number
  readonly mode: readonly number[]
  readonly trimmedMean: number   // 5% trimmed
  readonly sd: number
  readonly se: number
  readonly variance: number
  readonly min: number
  readonly max: number
  readonly range: number
  readonly iqr: number
  readonly q1: number
  readonly q3: number
  readonly skewness: number
  readonly kurtosis: number      // excess kurtosis
  readonly ci: readonly [number, number]
  readonly ciLevel: number
  readonly shapiroWilk: { statistic: number; pValue: number }
  readonly andersonDarling?: { statistic: number; pValue: number }
  readonly formatted: string
}

/** A single fixed effect from an LMM. */
export interface FixedEffect {
  readonly name: string
  readonly estimate: number
  readonly se: number
  readonly tValue: number
  readonly pValue: number
  readonly ci: readonly [number, number]
}

/** Full result from a linear mixed model. */
export interface LMMResult {
  readonly fixedEffects: readonly FixedEffect[]
  readonly varianceComponents: {
    readonly intercept: number   // σ²_b (random intercept variance)
    readonly residual: number    // σ²_e
    readonly slopes?: Readonly<Record<string, number>>
  }
  readonly randomCorrelations?: Readonly<Record<string, number>>  // correlations between random effects
  readonly icc: number           // intraclass correlation
  readonly logLik: number
  readonly aic: number
  readonly bic: number
  readonly nObs: number
  readonly nGroups: number
  readonly nParams: number       // total model parameters (for LRT)
  readonly method: 'REML' | 'ML'
  readonly r2Marginal: number    // Nakagawa R² — fixed effects only
  readonly r2Conditional: number // Nakagawa R² — fixed + random
  readonly formatted: string
}

/** PCA result with loadings, scores, and variance explained. */
export interface PCAResult {
  readonly loadings: readonly (readonly number[])[]   // [nComponents][nVars]
  readonly scores: readonly (readonly number[])[]     // [nObs][nComponents]
  readonly eigenvalues: readonly number[]
  readonly varianceExplained: readonly number[]       // proportions
  readonly cumulativeVariance: readonly number[]
  readonly nComponents: number
}

/** Regression coefficient with CI. */
export interface RegressionCoef {
  readonly name: string
  readonly estimate: number
  readonly se: number
  readonly tValue: number
  readonly pValue: number
  readonly ci: readonly [number, number]
}

/** Result from a linear regression (OLS). */
export interface RegressionResult {
  readonly coefficients: readonly RegressionCoef[]
  readonly r2: number
  readonly adjR2: number
  readonly fStatistic: number
  readonly fDf: readonly [number, number]
  readonly fPValue: number
  readonly aic: number
  readonly bic: number
  readonly residuals: readonly number[]
  readonly fitted: readonly number[]
  readonly residualNormality?: { test: string; statistic: number; pValue: number }
  readonly n: number
  readonly formatted: string
}

/** Frequency table row. */
export interface FrequencyRow {
  readonly value: string | number
  readonly count: number
  readonly relative: number     // proportion
  readonly cumulative: number   // cumulative proportion
}

/** Result from chi-square or Fisher test. */
export interface FrequencyTestResult extends StatResult {
  readonly table: readonly FrequencyRow[]
  readonly expectedCounts?: readonly (readonly number[])[]
}

/** Pairwise comparison result (post-hoc). */
export interface PairwiseResult {
  readonly group1: string
  readonly group2: string
  readonly meanDiff: number
  readonly se: number
  readonly statistic: number
  readonly pValue: number
  readonly pValueAdj: number
  readonly ci: readonly [number, number]
  readonly significant: boolean
}

/** Configuration for p-value adjustment. */
export type PAdjMethod = 'bonferroni' | 'holm' | 'BH' | 'BY' | 'none'

/** Generic data table: rows × columns. */
export type DataMatrix = readonly (readonly number[])[]

/** Named group data for comparison functions. */
export interface GroupData {
  readonly values: readonly number[]
  readonly label: string
}

// ─── Factor Analysis types ────────────────────────────────────────────────

/** Fit indices for factor models (EFA and CFA). */
export interface FactorFit {
  readonly chiSq: number
  readonly df: number
  readonly pValue: number
  readonly rmsea: number
  readonly rmseaCI: readonly [number, number]   // 90% CI
  readonly cfi: number
  readonly tli: number
  readonly srmr: number
  readonly aic: number
  readonly bic: number
}

/** A single parameter estimate with SE, z, p, and standardized value. */
export interface ParameterEstimate {
  readonly estimate: number
  readonly se: number
  readonly z: number
  readonly pValue: number
  readonly stdAll: number  // STDYX standardized
}

/** Diagnostics for factor analysis adequacy: KMO, Bartlett, MAP, parallel analysis. */
export interface FADiagnostics {
  readonly kmo: number
  readonly kmoPerItem: readonly number[]
  readonly bartlett: { readonly chiSq: number; readonly df: number; readonly pValue: number }
  readonly mapSuggested: number
  readonly parallelEigenvalues: readonly number[]       // observed eigenvalues
  readonly parallelSimulated: readonly number[]         // 95th percentile simulated thresholds
  readonly parallelSuggested: number
}

/** Result from exploratory factor analysis (EFA). */
export interface FAResult {
  readonly loadings: readonly (readonly number[])[]            // [items × factors] rotated
  readonly standardizedLoadings: readonly (readonly number[])[]
  readonly uniqueness: readonly number[]
  readonly communalities: readonly number[]
  readonly factorCorrelations: readonly (readonly number[])[]  // Phi matrix
  readonly fit: FactorFit
  readonly eigenvalues: readonly number[]
  readonly nFactors: number
  readonly rotation: string
  readonly extraction: string
  readonly variableNames: readonly string[]
  readonly factorNames: readonly string[]
  readonly formatted: string  // APA string
}

/** Result from confirmatory factor analysis (CFA). */
export interface CFAResult extends FAResult {
  readonly parameterEstimates: {
    readonly loadings: readonly (readonly ParameterEstimate[])[]   // [factor][item]
    readonly uniquenesses: readonly ParameterEstimate[]
    readonly factorCovariances: readonly (readonly ParameterEstimate[])[]
  }
  readonly model: Readonly<Record<string, readonly number[]>>
}

// ─── Repeated Measures ANOVA ─────────────────────────────────────────────

/** Result from repeated measures one-way ANOVA. */
export interface RMANOVAResult extends StatResult {
  readonly conditions: readonly {
    readonly label: string
    readonly mean: number
    readonly sd: number
    readonly n: number
  }[]
  readonly ssConditions: number
  readonly ssSubjects: number
  readonly ssError: number
  readonly ssTotal: number
  readonly dfConditions: number
  readonly dfSubjects: number
  readonly dfError: number
  readonly msConditions: number
  readonly msError: number
  readonly sphericity: {
    readonly W: number
    readonly chiSq: number
    readonly df: number
    readonly pValue: number
  } | null
  readonly epsilonGG: number
  readonly epsilonHF: number
  readonly correctedDf?: readonly [number, number]
  readonly correction: 'none' | 'Greenhouse-Geisser' | 'Huynh-Feldt'
}

// ─── Two-Way ANOVA ───────────────────────────────────────────────────────

/** A single row in the Two-Way ANOVA table (factor A, B, interaction, residual, total). */
export interface ANOVATableRow {
  readonly source: string
  readonly ss: number
  readonly df: number
  readonly ms: number
  readonly F: number
  readonly pValue: number
  readonly etaSq: number
}

/** Result from a two-way factorial ANOVA. */
export interface TwoWayANOVAResult {
  readonly rows: readonly ANOVATableRow[]
  readonly residual: { readonly ss: number; readonly df: number; readonly ms: number }
  readonly total: { readonly ss: number; readonly df: number }
  readonly n: number
  readonly formatted: string
}

// ─── ANCOVA ─────────────────────────────────────────────────────────────

/** Result from analysis of covariance. */
export interface ANCOVAResult {
  readonly rows: readonly ANOVATableRow[]
  readonly residual: { readonly ss: number; readonly df: number; readonly ms: number }
  readonly total: { readonly ss: number; readonly df: number }
  readonly adjustedMeans: readonly { readonly label: string; readonly adjustedMean: number }[]
  readonly n: number
  readonly formatted: string
}

// ─── Ordinal Logistic Regression ────────────────────────────────────────

/** Result from ordinal logistic (proportional odds) regression. */
export interface OrdinalRegressionResult {
  readonly thresholds: readonly { readonly name: string; readonly estimate: number; readonly se: number; readonly z: number; readonly pValue: number }[]
  readonly coefficients: readonly { readonly name: string; readonly estimate: number; readonly se: number; readonly z: number; readonly pValue: number; readonly ci: readonly [number, number] }[]
  readonly logLik: number
  readonly aic: number
  readonly bic: number
  readonly n: number
  readonly nCategories: number
  readonly formatted: string
}

// ─── Equivalence Testing (TOST) ─────────────────────────────────────────

/** Result from TOST equivalence testing. */
export interface EquivalenceResult {
  readonly testName: string
  readonly estimate: number              // observed value (mean diff, r, d)
  readonly bounds: readonly [number, number]  // [-delta, +delta]
  readonly t1: number                    // t-stat for H0: diff ≤ -delta
  readonly t2: number                    // t-stat for H0: diff ≥ +delta
  readonly p1: number                    // p-value for lower bound test
  readonly p2: number                    // p-value for upper bound test
  readonly pValue: number                // max(p1, p2)
  readonly df: number
  readonly ci: readonly [number, number] // (1-2α) CI
  readonly ciLevel: number               // actual CI level (e.g. 0.90 for α=0.05)
  readonly equivalent: boolean           // pValue < α
  readonly n: number
  readonly effectSize?: EffectSize
  readonly formatted: string
}

// ─── Logistic GLMM ─────────────────────────────────────────────────────

/** A single fixed effect from a logistic GLMM (z-tests, not t-tests). */
export interface GLMMFixedEffect {
  readonly name: string
  readonly estimate: number       // on log-odds scale
  readonly se: number
  readonly zValue: number         // Wald z
  readonly pValue: number
  readonly ci: readonly [number, number]
  readonly or: number             // exp(estimate) = odds ratio
  readonly orCI: readonly [number, number]  // exp(CI bounds)
}

/** Full result from a logistic GLMM. */
export interface GLMMResult {
  readonly fixedEffects: readonly GLMMFixedEffect[]
  readonly varianceComponents: {
    readonly intercept: number    // σ²_b (random intercept variance)
    readonly slopes?: Readonly<Record<string, number>>
  }
  readonly randomCorrelations?: Readonly<Record<string, number>>
  readonly icc: number            // latent-scale ICC = σ²_b / (σ²_b + π²/3)
  readonly logLik: number         // Laplace log-likelihood
  readonly deviance: number       // -2 * logLik
  readonly aic: number
  readonly bic: number
  readonly nObs: number
  readonly nGroups: number
  readonly nParams: number
  readonly family: 'binomial'
  readonly link: 'logit'
  readonly formatted: string
}

// ─── Bootstrap CI ───────────────────────────────────────────────────────

/** Result from bootstrap confidence interval estimation. */
export interface BootstrapCIResult {
  readonly estimate: number
  readonly ci: readonly [number, number]
  readonly se: number
  readonly ciLevel: number
  readonly nBoot: number
  readonly method: 'percentile' | 'bca'
  readonly formatted: string
}

// ─── analyze() dispatch layer ─────────────────────────────────────────────

/** Field type descriptor for the analyze() dispatch layer. */
export type FieldType = 'numeric' | 'binary' | 'categorical' | 'ordinal'

/**
 * A numeric outcome field.
 * @field type    - Always 'numeric'
 * @field name    - Column / variable name, used in result labels
 * @field values  - The raw numeric observations
 */
export interface NumericField {
  readonly type: 'numeric'
  readonly name: string
  readonly values: readonly number[]
}

/**
 * A grouping field whose values are labels (string or 0/1 int).
 * Declared as 'binary' (exactly 2 unique values) or 'categorical' (3+).
 * @field type    - 'binary' | 'categorical'
 * @field name    - Column / variable name
 * @field values  - One label per observation, parallel to the outcome field
 */
export interface GroupField {
  readonly type: 'binary' | 'categorical'
  readonly name: string
  readonly values: readonly (string | number)[]
}

/** Union of all field kinds accepted by analyze(). */
export type Field = NumericField | GroupField

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
 * @field equivalenceDelta - Equivalence bound for TOST tests. Required when
 *                           forceTest is a TOST variant.
 */
export interface AnalyzeOptions {
  readonly ciLevel?: number
  readonly paired?: boolean
  readonly pAdjMethod?: PAdjMethod
  readonly forceTest?: string
  readonly equalVariances?: boolean
  readonly normalityAlpha?: number
  readonly equivalenceDelta?: number
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
export interface AnalysisResult {
  readonly test: string
  readonly outcome: string
  readonly predictor?: string
  readonly result: StatResult | FrequencyTestResult | EquivalenceResult
  readonly descriptives?: readonly DescriptiveResult[]
  readonly posthoc?: readonly PairwiseResult[]
  readonly normality?: ReadonlyArray<{ group: string; W: number; p: number }>
}
