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
//# sourceMappingURL=frequency.d.ts.map