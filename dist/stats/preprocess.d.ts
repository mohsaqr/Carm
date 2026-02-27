/**
 * Data preprocessing for clustering and PCA.
 *
 * Provides center, standardize, log, and sqrt transforms with inverse.
 * Uses mean/sd from core/math.ts.
 *
 * Cross-validate with R:
 * > scale(data)                      # standardize
 * > scale(data, scale = FALSE)       # center only
 * > log(data)                        # log transform
 * > sqrt(data)                       # sqrt transform
 */
export type PreprocessMethod = 'none' | 'center' | 'standardize' | 'log' | 'sqrt';
export interface PreprocessOptions {
    readonly method?: PreprocessMethod;
}
export interface PreprocessResult {
    readonly data: readonly (readonly number[])[];
    readonly colMeans: readonly number[];
    readonly colSDs: readonly number[];
    readonly method: PreprocessMethod;
    readonly centered: boolean;
    readonly scaled: boolean;
}
/**
 * Preprocess a numeric data matrix.
 *
 * - 'none':        pass-through (colMeans/colSDs still computed for reference)
 * - 'center':      subtract column mean (R: scale(x, scale=FALSE))
 * - 'standardize': subtract mean, divide by SD (R: scale(x))
 * - 'log':         natural log (requires all values > 0)
 * - 'sqrt':        square root (requires all values >= 0)
 *
 * Zero-variance columns get SD = 1 to avoid division by zero.
 *
 * @param data - N × D numeric matrix
 * @param options - preprocessing configuration
 * @returns PreprocessResult with transformed data and parameters
 */
export declare function preprocessData(data: readonly (readonly number[])[], options?: PreprocessOptions): PreprocessResult;
/**
 * Inverse transform preprocessed data back to the original scale.
 *
 * @param data - N × D preprocessed matrix
 * @param params - the PreprocessResult containing transform parameters
 * @returns data in original scale
 */
export declare function inverseTransform(data: readonly (readonly number[])[], params: PreprocessResult): readonly (readonly number[])[];
//# sourceMappingURL=preprocess.d.ts.map