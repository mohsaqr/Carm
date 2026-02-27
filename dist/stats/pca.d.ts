/**
 * PCA (Principal Component Analysis) via SVD.
 * Also provides varimax rotation and factor loading computation.
 */
import type { PCAResult } from '../core/types.js';
/**
 * PCA via SVD on the standardized data matrix.
 * Equivalent to eigen-decomposition of the correlation matrix.
 *
 * Cross-validated with R:
 * > prcomp(data, scale. = TRUE)
 * > summary(pca)  # check proportion of variance explained
 */
export declare function runPCA(data: readonly (readonly number[])[], nComponents?: number, scale?: boolean): PCAResult;
/**
 * Varimax rotation of a loading/eigenvector matrix.
 * Maximizes the variance of squared loadings within each factor.
 * Uses SVD-based update (same algorithm as R's stats::varimax).
 *
 * Kaiser normalization (normalize=true, the default) rescales each row
 * to unit length before rotation and restores afterwards, matching R's
 * varimax(x, normalize=TRUE) default.
 *
 * Reference: Kaiser (1958), Psychometrika 23:187-200
 *
 * Cross-validated with R:
 * > varimax(pca$rotation[, 1:3])
 */
export declare function varimaxRotation(loadings: readonly (readonly number[])[], maxIter?: number, tol?: number, normalize?: boolean): {
    rotatedLoadings: number[][];
    rotationMatrix: number[][];
};
export interface ScreeData {
    readonly components: readonly number[];
    readonly eigenvalues: readonly number[];
    readonly varianceExplained: readonly number[];
    readonly cumulativeVariance: readonly number[];
}
/** Extract scree plot data from a PCA result. */
export declare function screeData(pca: PCAResult): ScreeData;
//# sourceMappingURL=pca.d.ts.map