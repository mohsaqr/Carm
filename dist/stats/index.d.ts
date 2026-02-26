import { b as DescriptiveResult, c as EffectSize, j as FrequencyTestResult, S as StatResult, i as FrequencyRow, G as GroupData, P as PAdjMethod, m as PairwiseResult, o as RegressionResult, l as PCAResult, L as LMMResult, f as Field, a as AnalyzeOptions, A as AnalysisResult, g as FieldType, C as CFAResult, d as FAResult, F as FADiagnostics } from '../matrix-fbbvM_BU.js';
export { C as ClusterDiagnostics, a as CorrelationMatrix, b as CovarianceModel, D as DBSCANOptions, c as DBSCANResult, G as GMMOptions, d as GMMRangeEntry, e as GMMResult, H as HACMerge, f as HACOptions, g as HACResult, K as KMeansOptions, h as KMeansRangeEntry, i as KMeansResult, L as LCAOptions, j as LCAResult, k as LTAOptions, l as LTAResult, m as LinkageMethod, P as PointType, n as correlationMatrix, o as cutTree, p as cutTreeHeight, q as euclideanDistMatrix, r as findBestGMM, s as fitGMM, t as fitGMMRange, u as fitKMeansRange, v as fitLCA, w as fitLTA, x as kDistancePlot, y as kendallTau, z as partialCorrelation, A as pearsonCorrelation, B as predictGMM, E as predictKMeans, F as runDBSCAN, I as runHierarchical, J as runKMeans, M as silhouetteScores, N as spearmanCorrelation } from '../clustering-28_Mu8Yg.js';
export { p as mean, q as median, u as quantile, x as sd, y as se, E as variance } from '../math-BE2LWmnL.js';

/**
 * Descriptive statistics module.
 * Computes mean, median, mode, variance, SD, SE, skewness, kurtosis,
 * percentiles, confidence intervals, and the Shapiro-Wilk normality test.
 */

/** α-trimmed mean: removes α proportion from each tail. */
declare function trimmedMean(x: readonly number[], alpha?: number): number;
/** Sample skewness (adjusted Fisher-Pearson). */
declare function skewness(x: readonly number[]): number;
/** Sample excess kurtosis (adjusted Fisher-Pearson). */
declare function kurtosis(x: readonly number[]): number;
/** t-based CI for the mean. Returns [lower, upper]. */
declare function ciMean(x: readonly number[], ciLevel?: number): readonly [number, number];
/**
 * Shapiro-Wilk W test for normality.
 * Full AS R94 algorithm (Royston 1995) — valid for n=3..5000.
 * Reference: Royston (1995), Applied Statistics, 44(4):547-551.
 *
 * Cross-validated with R:
 * > shapiro.test(c(1,2,3,4,5,6,7,8,9,10))
 * W = 0.9728, p-value = 0.9177
 */
declare function shapiroWilk(x: readonly number[]): {
    statistic: number;
    pValue: number;
};
/**
 * Compute full descriptive statistics for a numeric vector.
 *
 * Cross-validated with R:
 * > x <- c(2,4,4,4,5,5,7,9)
 * > mean(x)  # 5
 * > sd(x)    # 2
 * > e1071::skewness(x, type=2)  # 0.4895
 */
declare function describe(x: readonly number[], ciLevel?: number): DescriptiveResult;

/**
 * Effect size calculations for Carm.
 * Cohen's d, Hedges' g, eta-squared, omega-squared, rank-biserial correlation.
 * All functions return structured EffectSize objects.
 */

/**
 * Cohen's d for independent samples.
 * Uses pooled SD (equal variances assumed).
 * Formula: d = (M₁ - M₂) / SD_pooled
 * Reference: Cohen (1988), "Statistical Power Analysis for the Behavioral Sciences"
 *
 * Cross-validated with R:
 * > library(effsize)
 * > cohen.d(c(1,2,3,4,5), c(3,4,5,6,7))
 * d = -1.2649...  (negative because group2 > group1)
 */
declare function cohensD(x1: readonly number[], x2: readonly number[]): EffectSize;
/**
 * Cohen's d for paired samples (dependent t-test).
 * Formula: d = M_diff / SD_diff
 */
declare function cohensDPaired(diffs: readonly number[]): EffectSize;
/**
 * Hedges' g — bias-corrected version of Cohen's d.
 * Correction factor J = 1 - 3/(4·df - 1), df = n1 + n2 - 2.
 * Reference: Hedges (1981), Journal of Educational Statistics
 */
declare function hedgesG(x1: readonly number[], x2: readonly number[]): EffectSize;
/**
 * Eta-squared: η² = SS_between / SS_total
 * For one-way ANOVA. Between 0 and 1.
 * Reference: Cohen (1973)
 */
declare function etaSquared(ssBetween: number, ssTotal: number): EffectSize;
/**
 * Omega-squared: ω² = (SS_between - df_between · MS_within) / (SS_total + MS_within)
 * Less biased than eta-squared.
 * Reference: Hays (1963)
 */
declare function omegaSquared(ssBetween: number, ssTotal: number, dfBetween: number, msWithin: number): EffectSize;
/**
 * Rank-biserial correlation for Mann-Whitney U.
 * r = 1 - (2U) / (n1 * n2)
 * Reference: Wendt (1972)
 */
declare function rankBiserial(U: number, n1: number, n2: number): EffectSize;
/**
 * Rank-biserial correlation for Wilcoxon signed-rank (paired).
 * r = T / (n(n+1)/2) where T = sum of positive ranks (or negative).
 */
declare function rankBiserialWilcoxon(T: number, n: number): EffectSize;
/**
 * Eta-squared for Kruskal-Wallis: η²_H = (H - k + 1) / (n - k)
 * Reference: Tomczak & Tomczak (2014)
 */
declare function etaSquaredKW(H: number, k: number, n: number): EffectSize;
/**
 * Normal approximation CI for Cohen's d.
 * Reference: Hedges & Olkin (1985), "Statistical Methods for Meta-Analysis"
 */
declare function cohensDCI(d: number, n1: number, n2: number, ciLevel?: number): readonly [number, number];

/**
 * Frequency analysis module.
 * Frequency tables, chi-square test of independence, Fisher's exact test,
 * Cramér's V, phi coefficient, odds ratio, goodness-of-fit test.
 */

/**
 * Build a frequency table from an array of values.
 * Returns rows sorted by value with absolute, relative, and cumulative frequencies.
 */
declare function frequencyTable(data: readonly (string | number)[]): FrequencyRow[];
/** Convert grouped data to a contingency table (rows = group1, cols = group2). */
declare function contingencyTable(group1: readonly (string | number)[], group2: readonly (string | number)[]): {
    table: number[][];
    rowLabels: (string | number)[];
    colLabels: (string | number)[];
};
/**
 * Pearson chi-square test of independence for a contingency table.
 *
 * Cross-validated with R:
 * > chisq.test(matrix(c(10,20,30,40), nrow=2))
 * X-squared = 0.6061, df = 1, p-value = 0.4363
 * Cramér's V = 0.0875
 */
declare function chiSquareTest(observed: readonly (readonly number[])[], yatesCorrection?: boolean): FrequencyTestResult;
/**
 * Fisher's exact test for 2×2 contingency tables.
 * Uses the hypergeometric distribution.
 *
 * Cross-validated with R:
 * > fisher.test(matrix(c(11, 5, 3, 6), nrow=2))
 * p-value = 0.2684, odds ratio = 4.0
 */
declare function fisherExactTest(a: number, b: number, c: number, d: number): StatResult;
/**
 * Phi coefficient for 2×2 tables: φ = (ad - bc) / sqrt((a+b)(c+d)(a+c)(b+d))
 */
declare function phiCoefficient(a: number, b: number, c: number, d: number): number;
/**
 * Chi-square goodness-of-fit test.
 * Tests whether observed frequencies match expected proportions.
 *
 * Cross-validated with R:
 * > chisq.test(c(30, 20, 50), p = c(1/3, 1/3, 1/3))
 * X-squared = 10, df = 2, p-value = 0.006738
 */
declare function goodnessOfFit(observed: readonly number[], expected?: readonly number[]): StatResult;

/**
 * Group comparison module.
 * t-tests, one-way ANOVA, Mann-Whitney U, Wilcoxon signed-rank,
 * Kruskal-Wallis, Friedman test.
 */

/**
 * Welch's t-test for independent samples (unequal variances, default).
 * Set `equalVariances = true` for Student's t-test.
 *
 * Cross-validated with R:
 * > t.test(c(2,3,5,6,8), c(1,4,7,9,10), var.equal = FALSE)
 * t = -0.4851, df = 7.6, p-value = 0.6408
 * > t.test(c(2,3,5,6,8), c(1,4,7,9,10), var.equal = TRUE)
 * t = -0.4851, df = 8, p-value = 0.6403
 */
declare function tTestIndependent(x1: readonly number[], x2: readonly number[], equalVariances?: boolean, ciLevel?: number, alternative?: 'two.sided' | 'less' | 'greater'): StatResult;
/**
 * Paired samples t-test.
 *
 * Cross-validated with R:
 * > t.test(c(1,2,3,4,5), c(2,4,6,8,10), paired = TRUE)
 * t = -3.873, df = 4, p-value = 0.01789
 */
declare function tTestPaired(x1: readonly number[], x2: readonly number[], ciLevel?: number): StatResult;
interface ANOVAResult extends StatResult {
    readonly groups: readonly {
        readonly label: string;
        readonly n: number;
        readonly mean: number;
        readonly sd: number;
    }[];
    readonly ssBetween: number;
    readonly ssWithin: number;
    readonly ssTotal: number;
    readonly msBetween: number;
    readonly msWithin: number;
    readonly dfBetween: number;
    readonly dfWithin: number;
}
/**
 * One-way ANOVA.
 *
 * Cross-validated with R:
 * > oneway.test(value ~ group, data = df, var.equal = TRUE)
 * > aov(value ~ group, data = df)
 */
declare function oneWayANOVA(groups: readonly GroupData[]): ANOVAResult;
/**
 * Mann-Whitney U test (Wilcoxon rank-sum test).
 * Uses normal approximation for large n.
 *
 * Cross-validated with R:
 * > wilcox.test(c(1,2,3,4,5), c(3,4,5,6,7))
 * W = 5, p-value = 0.09502
 */
declare function mannWhitneyU(x1: readonly number[], x2: readonly number[], alternative?: 'two.sided' | 'less' | 'greater'): StatResult;
/**
 * Wilcoxon signed-rank test for paired data.
 *
 * Cross-validated with R:
 * > wilcox.test(c(1,2,3,4,5), c(2,4,6,8,10), paired = TRUE)
 * V = 0, p-value = 0.0625
 */
declare function wilcoxonSignedRank(x1: readonly number[], x2: readonly number[]): StatResult;
/**
 * Kruskal-Wallis H test (non-parametric one-way ANOVA).
 *
 * Cross-validated with R:
 * > kruskal.test(list(c(1,2,3), c(4,5,6), c(7,8,9)))
 * Kruskal-Wallis chi-squared = 7.2, df = 2, p-value = 0.02732
 */
declare function kruskalWallis(groups: readonly GroupData[]): StatResult;
/**
 * Friedman test for repeated measures (non-parametric ANOVA).
 * Data is a matrix: rows = subjects, columns = conditions.
 *
 * Cross-validated with R:
 * > friedman.test(matrix(c(1,2,3, 4,5,6, 7,8,9), nrow=3))
 */
declare function friedmanTest(data: readonly (readonly number[])[]): StatResult;

/**
 * Post-hoc tests for group comparisons.
 * Tukey HSD, Games-Howell, Dunn's test.
 */

/**
 * Tukey's Honest Significant Difference test.
 * Assumes equal variances (uses pooled MSE from ANOVA).
 * Uses Studentized range distribution approximation.
 *
 * Cross-validated with R:
 * > TukeyHSD(aov(value ~ group, data = df))
 */
declare function tukeyHSD(groups: readonly GroupData[], msWithin: number, dfWithin: number, ciLevel?: number): PairwiseResult[];
/**
 * Games-Howell test — does not assume equal variances.
 * Use when Levene's test is significant or group sizes differ substantially.
 *
 * Cross-validated with R:
 * > rstatix::games_howell_test(df, value ~ group)
 */
declare function gamesHowell(groups: readonly GroupData[], ciLevel?: number): PairwiseResult[];
/**
 * Dunn's post-hoc test following Kruskal-Wallis.
 * Uses rank sums and z-scores with adjustable p-value correction.
 *
 * Cross-validated with R:
 * > dunn.test::dunn.test(values, groups, method = "bonferroni")
 */
declare function dunnTest(groups: readonly GroupData[], method?: PAdjMethod): PairwiseResult[];

/**
 * Regression analysis module.
 * Simple and multiple OLS, logistic regression, polynomial regression,
 * diagnostics (R², AIC, BIC, VIF, residual plots).
 */

/**
 * Simple linear regression: y = β₀ + β₁·x
 */
declare function linearRegression(x: readonly number[], y: readonly number[], ciLevel?: number): RegressionResult;
/**
 * Multiple linear regression: y = β₀ + β₁·x₁ + ... + βₖ·xₖ
 * `predictors`: named columns { name: values[] }
 */
declare function multipleRegression(y: readonly number[], predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>, ciLevel?: number): RegressionResult;
/**
 * Polynomial regression: y = β₀ + β₁·x + β₂·x² + ... + βₖ·xᵏ
 */
declare function polynomialRegression(x: readonly number[], y: readonly number[], degree: number, ciLevel?: number): RegressionResult;
/**
 * Binary logistic regression via IRLS (iteratively reweighted least squares).
 * Outcome y must be 0/1.
 *
 * Cross-validated with R:
 * > glm(y ~ x1 + x2, family = binomial, data = df)
 */
declare function logisticRegression(y: readonly number[], predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>, ciLevel?: number, maxIter?: number, tol?: number): RegressionResult;
interface RegressionDiagnostics {
    readonly leverage: readonly number[];
    readonly cooksDistance: readonly number[];
    readonly standardizedResiduals: readonly number[];
    readonly vif: readonly number[];
}
/**
 * Compute regression diagnostics.
 * Returns leverage (hat values), Cook's distance, standardized residuals, VIF.
 */
declare function regressionDiagnostics(result: RegressionResult, predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>): RegressionDiagnostics;

/**
 * PCA (Principal Component Analysis) via SVD.
 * Also provides varimax rotation and factor loading computation.
 */

/**
 * PCA via SVD on the standardized data matrix.
 * Equivalent to eigen-decomposition of the correlation matrix.
 *
 * Cross-validated with R:
 * > prcomp(data, scale. = TRUE)
 * > summary(pca)  # check proportion of variance explained
 */
declare function runPCA(data: readonly (readonly number[])[], nComponents?: number, scale?: boolean): PCAResult;
/**
 * Varimax rotation of PCA loadings.
 * Maximizes the variance of squared loadings within each factor.
 * Reference: Kaiser (1958), Psychometrika 23:187-200
 *
 * Cross-validated with R:
 * > varimax(pca$rotation[, 1:3])
 */
declare function varimaxRotation(loadings: readonly (readonly number[])[], maxIter?: number, tol?: number): {
    rotatedLoadings: number[][];
    rotationMatrix: number[][];
};
interface ScreeData {
    readonly components: readonly number[];
    readonly eigenvalues: readonly number[];
    readonly varianceExplained: readonly number[];
    readonly cumulativeVariance: readonly number[];
}
/** Extract scree plot data from a PCA result. */
declare function screeData(pca: PCAResult): ScreeData;

/**
 * Linear Mixed Models (LMM) via REML.
 * Model: y = Xβ + Zb + ε
 *   b ~ N(0, σ²_b · I), ε ~ N(0, σ²_e · I)
 * Random intercepts + optional random slopes.
 * Optimization: Nelder-Mead on REML profile log-likelihood.
 *
 * Cross-validated with R:
 * > lme4::lmer(y ~ x + (1|group), data = df, REML = TRUE)
 */

interface LMMInput {
    readonly outcome: readonly number[];
    readonly fixedPredictors: Readonly<Record<string, readonly number[]>>;
    readonly groupId: readonly (string | number)[];
    readonly randomSlopes?: readonly string[];
    readonly ciLevel?: number;
}
/**
 * Fit a linear mixed model with random intercepts (and optionally random slopes).
 *
 * Cross-validated with R lme4:
 * > mod <- lmer(y ~ x + (1|group), data = df, REML = TRUE)
 * > fixef(mod)
 * > VarCorr(mod)
 * > icc(mod)
 */
declare function runLMM(input: LMMInput): LMMResult;
/**
 * Compute BLUPs (Best Linear Unbiased Predictors) — the random intercepts.
 * b_hat = σ²_b Z'V^{-1}(y - Xβ)
 */
declare function computeBLUPs(input: LMMInput, result: LMMResult): ReadonlyArray<{
    group: string | number;
    blup: number;
}>;

/**
 * Field-based analysis dispatch.
 * Pass named fields with declared types; analyze() selects and runs the
 * right statistical test automatically.
 */

/**
 * Infer the FieldType of an array of raw values.
 *
 * Rules (in order):
 *   1. If every value is a finite number and exactly 2 unique values both
 *      in {0, 1} → 'binary'
 *   2. If every value is a finite number → 'numeric'
 *   3. If there are exactly 2 unique string/mixed values → 'binary'
 *   4. Otherwise → 'categorical'
 *
 * @param values - Raw column values (string | number), length ≥ 1
 * @returns Inferred FieldType
 *
 * @example
 * detectFieldType([0, 1, 0, 1])      // → 'binary'
 * detectFieldType([1.2, 3.4, 5.6])   // → 'numeric'
 * detectFieldType(['A','B','A','C'])  // → 'categorical'
 * detectFieldType(['yes','no'])       // → 'binary'
 */
declare function detectFieldType(values: readonly (string | number)[]): FieldType;
/**
 * High-level statistical dispatch: pass fields, get the right test result.
 *
 * Automatically:
 *   - Detects field types from the declared .type property
 *   - Splits numeric outcome by group labels
 *   - Selects parametric vs non-parametric via Shapiro-Wilk
 *   - Runs the selected test with remaining options forwarded
 *   - Computes descriptive statistics for numeric outcomes
 *   - Runs post-hoc tests when 3+ groups are present
 *
 * @param outcome   - The outcome/dependent variable field
 * @param predictor - The grouping/independent variable field (optional)
 * @param opts      - Tuning options (all optional)
 *
 * @returns AnalysisResult with test name, StatResult, optional descriptives,
 *          optional posthoc, and the normality check used for routing.
 *
 * @throws Error if outcome and predictor have different lengths
 * @throws Error if paired=true but group sizes are unequal
 * @throws Error if forceTest names an unknown test
 *
 * @example — independent t-test (auto-detected)
 * analyze(
 *   { type: 'numeric', name: 'score', values: [72, 85, 90, 68, 77] },
 *   { type: 'binary',  name: 'group', values: ['A','B','A','B','A'] }
 * )
 *
 * @example — force Kruskal-Wallis
 * analyze(
 *   { type: 'numeric',     name: 'rt',   values: [...] },
 *   { type: 'categorical', name: 'cond', values: [...] },
 *   { forceTest: 'kruskal-wallis' }
 * )
 *
 * @example — paired t-test
 * analyze(
 *   { type: 'numeric', name: 'post', values: [80, 85, 90] },
 *   { type: 'binary',  name: 'time', values: ['pre','post','pre'] },
 *   { paired: true }
 * )
 */
declare function analyze(outcome: Field, predictor?: Field, opts?: AnalyzeOptions): AnalysisResult;

/**
 * Data preprocessing for clustering and PCA.
 *
 * Provides center, standardize, log, and sqrt transforms with inverse.
 * Uses mean/sd from core/math.ts.
 *
 * Cross-validate with R:
 * > scale(data)                      # standardize
 * > scale(data, scale = FALSE)       # center only
 * > log(data)                        # log transform
 * > sqrt(data)                       # sqrt transform
 */
type PreprocessMethod = 'none' | 'center' | 'standardize' | 'log' | 'sqrt';
interface PreprocessOptions {
    readonly method?: PreprocessMethod;
}
interface PreprocessResult {
    readonly data: readonly (readonly number[])[];
    readonly colMeans: readonly number[];
    readonly colSDs: readonly number[];
    readonly method: PreprocessMethod;
    readonly centered: boolean;
    readonly scaled: boolean;
}
/**
 * Preprocess a numeric data matrix.
 *
 * - 'none':        pass-through (colMeans/colSDs still computed for reference)
 * - 'center':      subtract column mean (R: scale(x, scale=FALSE))
 * - 'standardize': subtract mean, divide by SD (R: scale(x))
 * - 'log':         natural log (requires all values > 0)
 * - 'sqrt':        square root (requires all values >= 0)
 *
 * Zero-variance columns get SD = 1 to avoid division by zero.
 *
 * @param data - N × D numeric matrix
 * @param options - preprocessing configuration
 * @returns PreprocessResult with transformed data and parameters
 */
declare function preprocessData(data: readonly (readonly number[])[], options?: PreprocessOptions): PreprocessResult;
/**
 * Inverse transform preprocessed data back to the original scale.
 *
 * @param data - N × D preprocessed matrix
 * @param params - the PreprocessResult containing transform parameters
 * @returns data in original scale
 */
declare function inverseTransform(data: readonly (readonly number[])[], params: PreprocessResult): readonly (readonly number[])[];

/**
 * Factor Analysis module for Carm.
 * Exploratory Factor Analysis (EFA), Confirmatory Factor Analysis (CFA),
 * and psychometric diagnostics (KMO, Bartlett, MAP, Parallel Analysis).
 *
 * All functions are pure — no DOM, no D3, no side effects.
 * Uses splitmix32 PRNG for deterministic reproducibility.
 *
 * Algorithms adapted from:
 * - Iterated PAF: Gorsuch (1983), Factor Analysis (2nd ed.)
 * - ML extraction: Jöreskog (1967), Psychometrika 32(4)
 * - Varimax: Kaiser (1958), Psychometrika 23:187-200
 * - Oblimin/Promax: Jennrich (2002), Psychometrika 67(1)
 * - CFA: Bollen (1989), Structural Equations with Latent Variables
 * - Fit indices: Browne & Cudeck (1993), Hu & Bentler (1999)
 */

interface EFAOptions {
    readonly nFactors?: number;
    readonly extraction?: 'paf' | 'ml';
    readonly rotation?: 'varimax' | 'oblimin' | 'promax' | 'quartimin' | 'none';
    readonly seed?: number;
    readonly maxIter?: number;
    readonly tol?: number;
    readonly variableNames?: readonly string[];
}
interface CFAOptions {
    readonly maxIter?: number;
    readonly tol?: number;
    readonly variableNames?: readonly string[];
    readonly factorNames?: readonly string[];
}
interface FADiagnosticsOptions {
    readonly seed?: number;
    readonly parallelIterations?: number;
}
/**
 * Compute factor analysis diagnostics: KMO, Bartlett's test,
 * Velicer's MAP, and parallel analysis.
 *
 * Cross-validated with R:
 * > library(psych)
 * > KMO(data)
 * > cortest.bartlett(cor(data), n = nrow(data))
 * > fa.parallel(data, fm = "ml", fa = "fa")
 */
declare function runFADiagnostics(data: readonly (readonly number[])[], options?: FADiagnosticsOptions): FADiagnostics;
/**
 * Exploratory Factor Analysis with extraction (PAF or ML) and rotation.
 *
 * Cross-validated with R:
 * > library(psych)
 * > fa(data, nfactors = 3, fm = "ml", rotate = "promax")
 * > fa(data, nfactors = 3, fm = "minres", rotate = "varimax")
 */
declare function runEFA(data: readonly (readonly number[])[], options?: EFAOptions): FAResult;
/**
 * Confirmatory Factor Analysis via ML estimation.
 *
 * Model is specified as a record: { F1: [0, 1, 2], F2: [3, 4, 5] }
 * where keys are factor names and values are 0-indexed item indices.
 *
 * Cross-validated with R:
 * > library(lavaan)
 * > model <- 'F1 =~ x1 + x2 + x3\n F2 =~ x4 + x5 + x6'
 * > fit <- cfa(model, data)
 * > standardizedSolution(fit)
 * > fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli", "srmr"))
 */
declare function runCFA(data: readonly (readonly number[])[], model: Readonly<Record<string, readonly number[]>>, options?: CFAOptions): CFAResult;

export { type ANOVAResult, type CFAOptions, type EFAOptions, type FADiagnosticsOptions, type LMMInput, type PreprocessMethod, type PreprocessOptions, type PreprocessResult, type RegressionDiagnostics, type ScreeData, analyze, chiSquareTest, ciMean, cohensD, cohensDCI, cohensDPaired, computeBLUPs, contingencyTable, describe, detectFieldType, dunnTest, etaSquared, etaSquaredKW, fisherExactTest, frequencyTable, friedmanTest, gamesHowell, goodnessOfFit, hedgesG, inverseTransform, kruskalWallis, kurtosis, linearRegression, logisticRegression, mannWhitneyU, multipleRegression, omegaSquared, oneWayANOVA, phiCoefficient, polynomialRegression, preprocessData, rankBiserial, rankBiserialWilcoxon, regressionDiagnostics, runCFA, runEFA, runFADiagnostics, runLMM, runPCA, screeData, shapiroWilk, skewness, tTestIndependent, tTestPaired, trimmedMean, tukeyHSD, varimaxRotation, wilcoxonSignedRank };
