/**
 * Descriptive statistics module.
 * Computes mean, median, mode, variance, SD, SE, skewness, kurtosis,
 * percentiles, confidence intervals, and the Shapiro-Wilk normality test.
 */
import { mean as _mean, median as _median, variance as _variance, sd as _sd, se as _se, quantile } from '../core/math.js';
import type { DescriptiveResult } from '../core/types.js';
/** α-trimmed mean: removes α proportion from each tail. */
export declare function trimmedMean(x: readonly number[], alpha?: number): number;
/** Sample skewness (adjusted Fisher-Pearson). */
export declare function skewness(x: readonly number[]): number;
/** Sample excess kurtosis (adjusted Fisher-Pearson). */
export declare function kurtosis(x: readonly number[]): number;
/** t-based CI for the mean. Returns [lower, upper]. */
export declare function ciMean(x: readonly number[], ciLevel?: number): readonly [number, number];
/**
 * Shapiro-Wilk W test for normality.
 * Full AS R94 algorithm (Royston 1995) — valid for n=3..5000.
 * Reference: Royston (1995), Applied Statistics, 44(4):547-551.
 *
 * Cross-validated with R:
 * > shapiro.test(c(1,2,3,4,5,6,7,8,9,10))
 * W = 0.9728, p-value = 0.9177
 */
export declare function shapiroWilk(x: readonly number[]): {
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
export declare function describe(x: readonly number[], ciLevel?: number): DescriptiveResult;
export { _mean as mean, _median as median, _sd as sd, _se as se, _variance as variance, quantile };
//# sourceMappingURL=descriptive.d.ts.map