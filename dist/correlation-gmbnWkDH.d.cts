import { S as StatResult } from './types-DC8rlZlK.cjs';

/**
 * Correlation analysis module.
 * Pearson, Spearman, Kendall tau-b, partial correlation, correlation matrix.
 */

/**
 * Pearson product-moment correlation coefficient.
 *
 * Cross-validated with R:
 * > cor.test(c(1,2,3,4,5), c(2,4,1,5,3))
 * t = 0.6547, df = 3, p-value = 0.5607, r = 0.3536
 * 95% CI = [-0.6154, 0.9059]
 */
declare function pearsonCorrelation(x: readonly number[], y: readonly number[], ciLevel?: number): StatResult;
/**
 * Spearman rank correlation.
 *
 * Cross-validated with R:
 * > cor.test(c(1,2,3,4,5), c(5,6,7,8,7), method = "spearman")
 * rho = 0.8211, p-value = 0.08852
 */
declare function spearmanCorrelation(x: readonly number[], y: readonly number[], ciLevel?: number): StatResult;
/**
 * Kendall's tau-b correlation.
 *
 * Cross-validated with R:
 * > cor.test(c(1,2,3,4,5), c(5,6,7,8,7), method = "kendall")
 * tau = 0.7378, p-value = 0.1041
 */
declare function kendallTau(x: readonly number[], y: readonly number[], ciLevel?: number): StatResult;
/**
 * Partial correlation between x and y controlling for z.
 * r_xy.z = (r_xy - r_xz · r_yz) / sqrt((1 - r_xz²)(1 - r_yz²))
 *
 * Cross-validated with R:
 * > ppcor::pcor.test(x, y, z)
 */
declare function partialCorrelation(x: readonly number[], y: readonly number[], controls: readonly (readonly number[])[]): StatResult;
interface CorrelationMatrix {
    readonly r: readonly (readonly number[])[];
    readonly pValues: readonly (readonly number[])[];
    readonly n: number;
    readonly labels: readonly string[];
}
/**
 * Compute pairwise correlation matrix with p-values.
 * Method: 'pearson' | 'spearman' | 'kendall'
 */
declare function correlationMatrix(data: readonly (readonly number[])[], labels?: readonly string[], method?: 'pearson' | 'spearman' | 'kendall'): CorrelationMatrix;

export { type CorrelationMatrix as C, pearsonCorrelation as a, correlationMatrix as c, kendallTau as k, partialCorrelation as p, spearmanCorrelation as s };
