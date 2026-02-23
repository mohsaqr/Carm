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
  readonly icc: number           // intraclass correlation
  readonly logLik: number
  readonly aic: number
  readonly bic: number
  readonly nObs: number
  readonly nGroups: number
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
