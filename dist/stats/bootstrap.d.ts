/**
 * Bootstrap confidence interval estimation.
 * Percentile and BCa methods with seeded PRNG for reproducibility.
 */
import type { BootstrapCIResult } from '../core/types.js';
/**
 * Bootstrap CI for a single-sample statistic.
 * @param data       Input data array
 * @param statFn     Function that computes the statistic from a data array
 * @param options    { nBoot, ciLevel, method, seed }
 *
 * Cross-validated with R:
 * > library(boot)
 * > b <- boot(data, function(d,i) mean(d[i]), R=2000)
 * > boot.ci(b, type="perc")
 * > boot.ci(b, type="bca")
 */
export declare function bootstrapCI(data: readonly number[], statFn: (d: readonly number[]) => number, options?: {
    readonly nBoot?: number;
    readonly ciLevel?: number;
    readonly method?: 'percentile' | 'bca';
    readonly seed?: number;
}): BootstrapCIResult;
/**
 * Bootstrap CI for a two-sample statistic.
 * @param x1       First sample
 * @param x2       Second sample
 * @param statFn   Function that computes the statistic from two data arrays
 * @param options  { nBoot, ciLevel, seed }
 *
 * Cross-validated with R:
 * > library(boot)
 * > b <- boot(cbind(x1,x2), function(d,i) mean(d[i,1])-mean(d[i,2]), R=2000)
 * > boot.ci(b, type="perc")
 */
export declare function bootstrapCITwoSample(x1: readonly number[], x2: readonly number[], statFn: (a: readonly number[], b: readonly number[]) => number, options?: {
    readonly nBoot?: number;
    readonly ciLevel?: number;
    readonly seed?: number;
}): BootstrapCIResult;
//# sourceMappingURL=bootstrap.d.ts.map