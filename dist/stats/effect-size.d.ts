/**
 * Effect size calculations for Carm.
 * Cohen's d, Hedges' g, eta-squared, omega-squared, rank-biserial correlation.
 * All functions return structured EffectSize objects.
 */
import type { EffectSize } from '../core/types.js';
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
export declare function cohensD(x1: readonly number[], x2: readonly number[]): EffectSize;
/**
 * Cohen's d for paired samples (dependent t-test).
 * Formula: d = M_diff / SD_diff
 */
export declare function cohensDPaired(diffs: readonly number[]): EffectSize;
/**
 * Hedges' g — bias-corrected version of Cohen's d.
 * Correction factor J = 1 - 3/(4·df - 1), df = n1 + n2 - 2.
 * Reference: Hedges (1981), Journal of Educational Statistics
 */
export declare function hedgesG(x1: readonly number[], x2: readonly number[]): EffectSize;
/**
 * Eta-squared: η² = SS_between / SS_total
 * For one-way ANOVA. Between 0 and 1.
 * Reference: Cohen (1973)
 */
export declare function etaSquared(ssBetween: number, ssTotal: number): EffectSize;
/**
 * Omega-squared: ω² = (SS_between - df_between · MS_within) / (SS_total + MS_within)
 * Less biased than eta-squared.
 * Reference: Hays (1963)
 */
export declare function omegaSquared(ssBetween: number, ssTotal: number, dfBetween: number, msWithin: number): EffectSize;
/**
 * Rank-biserial correlation for Mann-Whitney U.
 * r = 1 - (2U) / (n1 * n2)
 * Reference: Wendt (1972)
 */
export declare function rankBiserial(U: number, n1: number, n2: number): EffectSize;
/**
 * Rank-biserial correlation for Wilcoxon signed-rank (paired).
 * r = T / (n(n+1)/2) where T = sum of positive ranks (or negative).
 */
export declare function rankBiserialWilcoxon(T: number, n: number): EffectSize;
/**
 * Eta-squared for Kruskal-Wallis: η²_H = (H - k + 1) / (n - k)
 * Reference: Tomczak & Tomczak (2014)
 */
export declare function etaSquaredKW(H: number, k: number, n: number): EffectSize;
/**
 * Normal approximation CI for Cohen's d.
 * Reference: Hedges & Olkin (1985), "Statistical Methods for Meta-Analysis"
 */
export declare function cohensDCI(d: number, n1: number, n2: number, ciLevel?: number): readonly [number, number];
//# sourceMappingURL=effect-size.d.ts.map