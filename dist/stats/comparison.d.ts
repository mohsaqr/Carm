/**
 * Group comparison module.
 * t-tests, one-way ANOVA, Mann-Whitney U, Wilcoxon signed-rank,
 * Kruskal-Wallis, Friedman test.
 */
import type { StatResult, GroupData } from '../core/types.js';
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
//# sourceMappingURL=comparison.d.ts.map