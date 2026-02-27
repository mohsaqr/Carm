/**
 * Group comparison module.
 * t-tests, one-way ANOVA, Mann-Whitney U, Wilcoxon signed-rank,
 * Kruskal-Wallis, Friedman test.
 */
import type { StatResult, GroupData, RMANOVAResult } from '../core/types.js';
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
export declare function tTestIndependent(x1: readonly number[], x2: readonly number[], equalVariances?: boolean, ciLevel?: number, alternative?: 'two.sided' | 'less' | 'greater'): StatResult;
/**
 * Paired samples t-test.
 *
 * Cross-validated with R:
 * > t.test(c(1,2,3,4,5), c(2,4,6,8,10), paired = TRUE)
 * t = -3.873, df = 4, p-value = 0.01789
 */
export declare function tTestPaired(x1: readonly number[], x2: readonly number[], ciLevel?: number): StatResult;
/**
 * Levene's test for homogeneity of variances.
 *
 * Computes absolute deviations from each group's center, then runs
 * a one-way ANOVA on those deviations. A significant result (p < .05)
 * indicates unequal variances across groups.
 *
 * @param groups  Two or more groups with labels and values.
 * @param center  'median' (default, Brown-Forsythe variant — more robust)
 *                or 'mean' (original Levene).
 *
 * Cross-validated with R:
 * > car::leveneTest(y ~ group, center = median)
 * > # or manually:
 * > medians <- tapply(y, group, median)
 * > abs_dev <- abs(y - medians[group])
 * > anova(lm(abs_dev ~ group))
 */
export declare function leveneTest(groups: readonly GroupData[], center?: 'median' | 'mean'): StatResult;
export interface ANOVAResult extends StatResult {
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
    readonly levene: {
        readonly statistic: number;
        readonly pValue: number;
        readonly df: readonly [number, number];
        readonly homogeneous: boolean;
    };
}
/**
 * One-way ANOVA.
 *
 * Cross-validated with R:
 * > oneway.test(value ~ group, data = df, var.equal = TRUE)
 * > aov(value ~ group, data = df)
 */
export declare function oneWayANOVA(groups: readonly GroupData[]): ANOVAResult;
/**
 * Mann-Whitney U test (Wilcoxon rank-sum test).
 * Uses normal approximation for large n.
 *
 * Cross-validated with R:
 * > wilcox.test(c(1,2,3,4,5), c(3,4,5,6,7))
 * W = 5, p-value = 0.09502
 */
export declare function mannWhitneyU(x1: readonly number[], x2: readonly number[], alternative?: 'two.sided' | 'less' | 'greater'): StatResult;
/**
 * Wilcoxon signed-rank test for paired data.
 *
 * Cross-validated with R:
 * > wilcox.test(c(1,2,3,4,5), c(2,4,6,8,10), paired = TRUE)
 * V = 0, p-value = 0.0625
 */
export declare function wilcoxonSignedRank(x1: readonly number[], x2: readonly number[]): StatResult;
/**
 * Kruskal-Wallis H test (non-parametric one-way ANOVA).
 *
 * Cross-validated with R:
 * > kruskal.test(list(c(1,2,3), c(4,5,6), c(7,8,9)))
 * Kruskal-Wallis chi-squared = 7.2, df = 2, p-value = 0.02732
 */
export declare function kruskalWallis(groups: readonly GroupData[]): StatResult;
/**
 * Friedman test for repeated measures (non-parametric ANOVA).
 * Data is a matrix: rows = subjects, columns = conditions.
 *
 * Cross-validated with R:
 * > friedman.test(matrix(c(1,2,3, 4,5,6, 7,8,9), nrow=3))
 */
export declare function friedmanTest(data: readonly (readonly number[])[]): StatResult;
/**
 * Mauchly's test of sphericity for repeated measures designs.
 * Tests whether the covariance matrix of orthonormalized contrasts
 * is proportional to an identity matrix.
 *
 * @param data  Subjects × conditions matrix (n × k, k ≥ 3).
 * @returns { W, chiSq, df, pValue } where W is Mauchly's statistic.
 *
 * Reference: Mauchly (1940), "Significance test for sphericity of a
 * normal n-variate distribution"
 *
 * Cross-validated with R:
 * > library(ez)
 * > ezANOVA(data=df, dv=y, wid=subj, within=cond)$`Mauchly's Test for Sphericity`
 */
export declare function mauchlysTest(data: readonly (readonly number[])[]): {
    readonly W: number;
    readonly chiSq: number;
    readonly df: number;
    readonly pValue: number;
};
/**
 * Greenhouse-Geisser and Huynh-Feldt epsilon corrections for sphericity.
 *
 * @param data  Subjects × conditions matrix (n × k).
 * @returns { gg, hf } clamped to [1/(k-1), 1].
 *
 * Reference: Greenhouse & Geisser (1959); Huynh & Feldt (1976)
 */
export declare function epsilonCorrections(data: readonly (readonly number[])[]): {
    readonly gg: number;
    readonly hf: number;
};
/**
 * Repeated measures one-way ANOVA.
 *
 * Decomposes total variance into SS_subjects + SS_conditions + SS_error.
 * Tests whether condition means differ, accounting for within-subject correlation.
 * Automatically runs Mauchly's sphericity test (k ≥ 3) and applies
 * Greenhouse-Geisser correction when sphericity is violated (p < .05).
 *
 * @param data     n × k matrix: data[i][j] = score for subject i in condition j.
 * @param options  { labels?: string[], ciLevel?: number }
 *
 * Cross-validated with R:
 * > fit <- aov(y ~ cond + Error(subj/cond), data=df)
 * > summary(fit)
 * > library(ez); ezANOVA(data=df, dv=y, wid=subj, within=cond)
 */
export declare function repeatedMeasuresANOVA(data: readonly (readonly number[])[], options?: {
    readonly labels?: readonly string[];
    readonly ciLevel?: number;
}): RMANOVAResult;
//# sourceMappingURL=comparison.d.ts.map