/**
 * Frequency analysis module.
 * Frequency tables, chi-square test of independence, Fisher's exact test,
 * Cramér's V, phi coefficient, odds ratio, goodness-of-fit test.
 */
import type { FrequencyRow, FrequencyTestResult, StatResult } from '../core/types.js';
/**
 * Build a frequency table from an array of values.
 * Returns rows sorted by value with absolute, relative, and cumulative frequencies.
 */
export declare function frequencyTable(data: readonly (string | number)[]): FrequencyRow[];
/** Convert grouped data to a contingency table (rows = group1, cols = group2). */
export declare function contingencyTable(group1: readonly (string | number)[], group2: readonly (string | number)[]): {
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
export declare function chiSquareTest(observed: readonly (readonly number[])[], yatesCorrection?: boolean): FrequencyTestResult;
/**
 * Fisher's exact test for 2×2 contingency tables.
 * Uses the hypergeometric distribution.
 *
 * Cross-validated with R:
 * > fisher.test(matrix(c(11, 5, 3, 6), nrow=2))
 * p-value = 0.2684, odds ratio = 4.0
 */
export declare function fisherExactTest(a: number, b: number, c: number, d: number): StatResult;
/**
 * Phi coefficient for 2×2 tables: φ = (ad - bc) / sqrt((a+b)(c+d)(a+c)(b+d))
 */
export declare function phiCoefficient(a: number, b: number, c: number, d: number): number;
/**
 * Chi-square goodness-of-fit test.
 * Tests whether observed frequencies match expected proportions.
 *
 * Cross-validated with R:
 * > chisq.test(c(30, 20, 50), p = c(1/3, 1/3, 1/3))
 * X-squared = 10, df = 2, p-value = 0.006738
 */
export declare function goodnessOfFit(observed: readonly number[], expected?: readonly number[]): StatResult;
/**
 * Cramér's V with bootstrap confidence interval via multinomial resampling.
 *
 * Cross-validated with R:
 * > library(rcompanion)
 * > m <- matrix(c(10,20,30,40,50,60), nrow=2)
 * > cramerV(m, ci=TRUE, R=2000)
 */
export declare function cramersVWithCI(observed: readonly (readonly number[])[], ciLevel?: number, nBoot?: number, seed?: number): StatResult;
/**
 * McNemar's test for paired nominal data (2×2 table).
 * Inputs b and c are off-diagonal counts: b = 0→1 switches, c = 1→0 switches.
 *
 * Cross-validated with R:
 * > mcnemar.test(matrix(c(20, 5, 10, 15), nrow=2), correct=FALSE)
 * McNemar's chi-squared = 1.6667, df = 1, p-value = 0.1967
 */
export declare function mcnemarsTest(b: number, c: number, correction?: boolean): StatResult;
/**
 * Exact binomial test with Clopper-Pearson CI.
 * Uses regularized incomplete beta function for exact p-values and CI.
 *
 * Cross-validated with R:
 * > binom.test(6, 10, p = 0.5)
 * p-value = 0.7539, 95% CI [0.2624, 0.8784]
 */
export declare function binomialTest(successes: number, trials: number, p0?: number, alternative?: 'two.sided' | 'less' | 'greater', ciLevel?: number): StatResult;
/**
 * Z-test for proportions (1-sample or 2-sample).
 *
 * Cross-validated with R:
 * > prop.test(c(30, 40), c(100, 100), correct=FALSE)
 * X-squared = 2.2222, p-value = 0.1360  (z = -1.4907)
 *
 * > prop.test(60, 100, p=0.5, correct=FALSE)
 * X-squared = 4, p-value = 0.0455  (z = 2.0)
 */
export declare function proportionsZTest(x1: number, n1: number, x2?: number, n2?: number, p0?: number, alternative?: 'two.sided' | 'less' | 'greater', ciLevel?: number, yates?: boolean): StatResult;
//# sourceMappingURL=frequency.d.ts.map