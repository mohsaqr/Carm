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
export declare function formatP(p: number): string;
/**
 * Format a statistic value: 2 decimal places, no leading zero for absolute values.
 * Used for t, F, z, χ², etc.
 */
export declare function formatStat(v: number, decimals?: number): string;
/**
 * Format a confidence interval: [lo, hi] → "[0.08, 1.22]"
 */
export declare function formatCI(ci: readonly [number, number], decimals?: number): string;
/**
 * Format degrees of freedom — integer if whole number, else 2 decimal places.
 */
export declare function formatDF(df: number | readonly [number, number]): string;
/**
 * Full APA string for a t-test result.
 * e.g. "t(48) = 2.31, p = .025, d = 0.65, 95% CI [0.08, 1.22]"
 */
export declare function formatTTest(t: number, df: number, pValue: number, effectSize: number, effectName: string, ci: readonly [number, number], ciLevel?: number): string;
/**
 * Full APA string for ANOVA.
 * e.g. "F(2, 87) = 4.21, p = .018, η² = 0.09"
 */
export declare function formatANOVA(F: number, df1: number, df2: number, pValue: number, effectSize: number, effectName?: string): string;
/**
 * Full APA string for chi-square test.
 * e.g. "χ²(3) = 8.46, p = .037, V = 0.29"
 */
export declare function formatChiSq(chiSq: number, df: number, pValue: number, effectSize: number, effectName?: string): string;
/**
 * APA string for Pearson/Spearman/Kendall correlation.
 * e.g. "r(48) = 0.54, p = .003, 95% CI [0.28, 0.73]"
 */
export declare function formatCorrelation(r: number, df: number, pValue: number, ci: readonly [number, number], method?: string, ciLevel?: number): string;
/**
 * APA string for linear regression.
 * e.g. "R² = 0.42, F(3, 96) = 23.2, p < .001"
 */
export declare function formatRegression(r2: number, adjR2: number, F: number, df1: number, df2: number, pValue: number): string;
/**
 * APA string for Mann-Whitney / Wilcoxon.
 * e.g. "W = 423, p = .031, r = 0.27"
 */
export declare function formatMannWhitney(W: number, pValue: number, r: number): string;
/**
 * APA string for Kruskal-Wallis.
 * e.g. "H(2) = 12.3, p = .002, η²_H = 0.18"
 */
export declare function formatKruskalWallis(H: number, df: number, pValue: number, etaSq: number): string;
/**
 * APA string for LMM.
 * e.g. "ICC = 0.42, AIC = 1234.5"
 */
export declare function formatLMM(icc: number, aic: number, bic: number, logLik: number): string;
/**
 * Interpret Cohen's d effect size.
 */
export declare function interpretCohensD(d: number): string;
/**
 * Interpret eta-squared / omega-squared.
 */
export declare function interpretEtaSq(eta: number): string;
/**
 * Interpret Pearson r.
 */
export declare function interpretR(r: number): string;
/**
 * Interpret Cramér's V.
 */
export declare function interpretCramerV(v: number, df: number): string;
import type { EffectInterpretation, FactorFit } from './types.js';
/**
 * APA string for CFA/EFA fit indices.
 * e.g. "χ²(24) = 28.42, p = .241, RMSEA = 0.042 [0.000, 0.085], CFI = 0.987, TLI = 0.982, SRMR = 0.038"
 */
export declare function formatCFAFit(fit: FactorFit): string;
/** Generic effect size interpretation helper. */
export declare function interpretEffect(value: number, thresholds: readonly [number, number, number]): EffectInterpretation;
//# sourceMappingURL=apa.d.ts.map