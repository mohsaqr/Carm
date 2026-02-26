import { S as StatResult, M as Matrix } from './matrix-fbbvM_BU.js';

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

/**
 * Clustering & Mixture Models: GMM, LCA, LTA, K-Means.
 *
 * - GMM: Gaussian Mixture with EM, K-Means++ init, mclust-style covariance constraints
 * - LCA: Latent Class Analysis for binary data (MLE, matches poLCA)
 * - LTA: Latent Transition Analysis (Hidden Markov LCA) with Baum-Welch in log-space
 * - K-Means: Lloyd's algorithm with K-Means++ init and empty-cluster re-seeding
 *
 * All functions are deterministic via a seeded PRNG (default seed: 42).
 * Cross-validate against: mclust (GMM), poLCA (LCA), seqHMM (LTA), stats::kmeans (K-Means).
 */

type CovarianceModel = 'VVV' | 'EEE' | 'VVI' | 'EEI' | 'VII' | 'EII';
interface ClusterDiagnostics {
    readonly converged: boolean;
    readonly iterations: number;
    readonly logLikelihood: number;
    readonly df: number;
    readonly aic: number;
    readonly bic: number;
    readonly icl: number;
    readonly entropy: number;
    readonly avepp: readonly number[];
    readonly formatted: string;
}
interface GMMOptions {
    readonly k: number;
    readonly model?: CovarianceModel;
    readonly seed?: number;
    readonly tol?: number;
    readonly maxIter?: number;
    readonly regCovar?: number;
}
interface GMMResult {
    readonly weights: readonly number[];
    readonly means: readonly number[][];
    readonly covariances: readonly Matrix[];
    readonly posteriors: readonly (readonly number[])[];
    readonly labels: readonly number[];
    readonly diagnostics: ClusterDiagnostics;
}
interface LCAOptions {
    readonly k: number;
    readonly seed?: number;
    readonly tol?: number;
    readonly maxIter?: number;
}
interface LCAResult {
    readonly rho: readonly (readonly number[])[];
    readonly priorWeights: readonly number[];
    readonly posteriors: readonly (readonly number[])[];
    readonly labels: readonly number[];
    readonly diagnostics: ClusterDiagnostics;
}
interface LTAOptions {
    readonly k: number;
    readonly seed?: number;
    readonly tol?: number;
    readonly maxIter?: number;
}
interface LTAResult {
    readonly pi: readonly number[];
    readonly tau: readonly (readonly number[])[];
    readonly rho: readonly (readonly number[])[];
    readonly trajectories: readonly (readonly number[])[];
    readonly posteriors: readonly (readonly (readonly number[])[])[];
    readonly diagnostics: ClusterDiagnostics;
}
interface KMeansOptions {
    readonly k: number;
    readonly seed?: number;
    readonly maxIter?: number;
    readonly tol?: number;
}
interface KMeansResult {
    readonly centroids: readonly (readonly number[])[];
    readonly labels: readonly number[];
    readonly inertia: number;
    readonly converged: boolean;
    readonly iterations: number;
}
/**
 * Fit a Gaussian Mixture Model via Expectation-Maximization.
 *
 * Supports mclust-style covariance constraints:
 * - VVV: Variable volume, variable shape, variable orientation (full covariance per component)
 * - EEE: Equal volume, equal shape, equal orientation (single pooled covariance)
 * - VVI: Variable volume, variable shape, identity orientation (diagonal, per component)
 * - EEI: Equal volume, equal shape, identity orientation (single shared diagonal)
 * - VII: Variable volume, identity shape, identity orientation (spherical, per component)
 * - EII: Equal volume, identity shape, identity orientation (single shared scalar × I)
 *
 * @param data - N × D numeric data matrix (array of observation arrays)
 * @param options - GMM configuration
 * @returns GMMResult with weights, means, covariances, posteriors, labels, diagnostics
 *
 * Cross-validate with R:
 * > library(mclust)
 * > fit <- Mclust(data, G=3, modelNames="VVV")
 * > fit$parameters$mean
 * > fit$parameters$variance$sigma
 * > fit$BIC
 */
declare function fitGMM(data: readonly (readonly number[])[], options: GMMOptions): GMMResult;
/**
 * Predict cluster assignments for new data given a fitted GMM.
 */
declare function predictGMM(data: readonly (readonly number[])[], result: GMMResult): {
    readonly labels: readonly number[];
    readonly posteriors: readonly (readonly number[])[];
};
/**
 * Automatic model selection: fit GMM across a grid of K and covariance models,
 * return the model with the lowest BIC.
 *
 * @param data - N × D data matrix
 * @param kRange - Array of K values to try (default [1,2,3,4,5])
 * @param models - Array of covariance models to try (default all 6)
 * @returns The GMMResult with the lowest BIC
 */
declare function findBestGMM(data: readonly (readonly number[])[], kRange?: readonly number[], models?: readonly CovarianceModel[]): GMMResult;
/**
 * Fit a Latent Class Analysis model for binary data.
 *
 * Uses EM with Bernoulli emission model and Beta(1,1) (uniform) prior smoothing.
 *
 * @param data - N × M binary matrix (0/1 values)
 * @param options - LCA configuration
 * @returns LCAResult with rho (item-response probabilities), priorWeights, posteriors, labels
 *
 * Cross-validate with R:
 * > library(poLCA)
 * > f <- cbind(V1, V2, V3, ...) ~ 1
 * > fit <- poLCA(f, data, nclass=3, nrep=1, probs.start=...)
 * > fit$probs
 */
declare function fitLCA(data: readonly (readonly number[])[], options: LCAOptions): LCAResult;
/**
 * Fit a Latent Transition Analysis model (categorical Hidden Markov Model).
 *
 * Uses Baum-Welch (EM) in log-space for numerical stability.
 * Measurement model is time-invariant (measurement invariance assumption).
 * Viterbi decoding provides most-likely state trajectories.
 *
 * @param data - N × T × M binary tensor (subjects × timepoints × items)
 * @param options - LTA configuration
 * @returns LTAResult with pi, tau, rho, trajectories, posteriors, diagnostics
 *
 * Cross-validate with R:
 * > library(seqHMM)
 * > # or manual forward-backward on small synthetic example
 */
declare function fitLTA(data: readonly (readonly (readonly number[])[])[], options: LTAOptions): LTAResult;
/**
 * K-Means clustering with K-Means++ initialization and empty-cluster re-seeding.
 *
 * @param data - N × D numeric data matrix
 * @param options - K-Means configuration
 * @returns KMeansResult with centroids, labels, inertia
 *
 * Cross-validate with R:
 * > km <- kmeans(data, centers=3, nstart=1, algorithm="Lloyd")
 * > km$centers; km$cluster; km$tot.withinss
 */
declare function runKMeans(data: readonly (readonly number[])[], options: KMeansOptions): KMeansResult;
/**
 * Predict cluster assignments for new data given fitted K-Means centroids.
 */
/** A single entry from fitting GMM at one K value. */
interface GMMRangeEntry {
    readonly k: number;
    readonly model: CovarianceModel;
    readonly result: GMMResult;
}
/** A single entry from fitting KMeans at one K value. */
interface KMeansRangeEntry {
    readonly k: number;
    readonly result: KMeansResult;
}
/**
 * Fit GMM for each K in kRange and return results sorted by K.
 * Skips failed fits (singular covariance etc). Throws only if ALL fail.
 *
 * @param data - N × D numeric data matrix
 * @param kRange - Array of K values to try, e.g. [2,3,4,5,6,7,8,9,10]
 * @param model - Covariance model (default 'VVV')
 * @returns GMMRangeEntry[] sorted by K
 */
declare function fitGMMRange(data: readonly (readonly number[])[], kRange: readonly number[], model?: CovarianceModel): readonly GMMRangeEntry[];
/**
 * Fit KMeans for each K in kRange and return results sorted by K.
 * Skips failed fits. Throws only if ALL fail.
 *
 * @param data - N × D numeric data matrix
 * @param kRange - Array of K values to try
 * @returns KMeansRangeEntry[] sorted by K
 */
declare function fitKMeansRange(data: readonly (readonly number[])[], kRange: readonly number[]): readonly KMeansRangeEntry[];
declare function predictKMeans(data: readonly (readonly number[])[], centroids: readonly (readonly number[])[]): readonly number[];
/**
 * Compute Euclidean distance matrix (N × N, row-major Float64Array).
 *
 * @param data - N × D numeric data matrix
 * @returns flat Float64Array of size N*N with dist[i*N+j] = euclidean(i,j)
 *
 * Cross-validate with R:
 * > as.matrix(dist(data))
 */
declare function euclideanDistMatrix(data: readonly (readonly number[])[]): Float64Array;
/**
 * Compute silhouette scores for each point (excluding noise labels = -1).
 *
 * For each non-noise point i:
 *   a(i) = mean dist to points in same cluster
 *   b(i) = min over other clusters of mean dist to that cluster
 *   s(i) = (b(i) - a(i)) / max(a(i), b(i))
 *
 * @param data - N × D data matrix
 * @param labels - cluster assignments (0-indexed, -1 = noise)
 * @returns scores per point (NaN for noise) and mean silhouette (excluding noise)
 *
 * Cross-validate with R:
 * > library(cluster)
 * > silhouette(labels, dist(data))
 */
declare function silhouetteScores(data: readonly (readonly number[])[], labels: readonly number[]): {
    readonly scores: readonly number[];
    readonly mean: number;
};
type PointType = 'core' | 'border' | 'noise';
interface DBSCANOptions {
    readonly eps: number;
    readonly minPts: number;
    readonly seed?: number;
    readonly preprocess?: 'none' | 'center' | 'standardize' | 'log' | 'sqrt';
}
interface DBSCANResult {
    readonly labels: readonly number[];
    readonly pointTypes: readonly PointType[];
    readonly nClusters: number;
    readonly nNoise: number;
    readonly clusterSizes: readonly number[];
    readonly silhouette: {
        readonly scores: readonly number[];
        readonly mean: number;
    };
    readonly formatted: string;
}
/**
 * DBSCAN clustering (Ester et al. 1996).
 *
 * Algorithm:
 * 1. For each point, compute eps-neighborhood via distance scan.
 * 2. Core points: |neighbors| >= minPts. BFS expansion from cores.
 * 3. Border points: not core, but within eps of a core point.
 * 4. Noise: neither core nor border.
 *
 * Labels: -1 = noise, 0, 1, 2, ... = cluster IDs (0-indexed).
 *
 * @param data - N × D numeric data matrix
 * @param options - DBSCAN configuration (eps, minPts)
 * @returns DBSCANResult with labels, point types, silhouette
 *
 * Cross-validate with R:
 * > library(dbscan)
 * > db <- dbscan(data, eps = 1.5, minPts = 5)
 * > db$cluster  # R: 0=noise, 1-indexed → Carm: -1=noise, 0-indexed
 */
declare function runDBSCAN(data: readonly (readonly number[])[], options: DBSCANOptions): DBSCANResult;
/**
 * Compute k-distance plot data for epsilon estimation.
 *
 * For each point, compute the distance to its k-th nearest neighbor,
 * then return these distances sorted in ascending order.
 * The "elbow" in the sorted plot suggests a good epsilon.
 *
 * @param data - N × D numeric data matrix
 * @param k - neighbor index (typically minPts)
 * @returns sorted k-NN distances (ascending)
 *
 * Cross-validate with R:
 * > library(dbscan)
 * > kNNdist(data, k = 5)  # sorted externally
 */
declare function kDistancePlot(data: readonly (readonly number[])[], k: number): readonly number[];
type LinkageMethod = 'single' | 'complete' | 'average' | 'ward';
interface HACOptions {
    readonly linkage?: LinkageMethod;
    readonly preprocess?: 'none' | 'center' | 'standardize' | 'log' | 'sqrt';
}
interface HACMerge {
    readonly a: number;
    readonly b: number;
    readonly height: number;
}
interface HACResult {
    readonly merges: readonly HACMerge[];
    readonly heights: readonly number[];
    readonly order: readonly number[];
    readonly copheneticCorrelation: number;
    readonly formatted: string;
}
/**
 * Hierarchical agglomerative clustering using Lance-Williams recurrence.
 *
 * Linkage methods and their Lance-Williams coefficients:
 * | Method   | α_i           | α_j           | β        | γ    |
 * |----------|---------------|---------------|----------|------|
 * | single   | 0.5           | 0.5           | 0        | -0.5 |
 * | complete | 0.5           | 0.5           | 0        | 0.5  |
 * | average  | n_i/(n_i+n_j) | n_j/(n_i+n_j) | 0        | 0    |
 * | ward     | (n_i+n_k)/N_t | (n_j+n_k)/N_t | -n_k/N_t | 0    |
 *
 * Ward's method uses squared Euclidean distances internally.
 * Final merge heights are sqrt(distance) to match R's hclust(method="ward.D2").
 *
 * @param data - N × D numeric data matrix
 * @param options - linkage method (default: ward)
 * @returns HACResult with merges, heights, leaf order, cophenetic correlation
 *
 * Cross-validate with R:
 * > hc <- hclust(dist(data), method = "ward.D2")
 * > hc$merge; hc$height; hc$order
 * > cor(cophenetic(hc), dist(data))
 */
declare function runHierarchical(data: readonly (readonly number[])[], options?: HACOptions): HACResult;
/**
 * Cut a dendrogram at K clusters.
 *
 * @param result - HACResult from runHierarchical
 * @param k - number of clusters desired
 * @returns 0-indexed cluster labels
 *
 * Cross-validate with R:
 * > cutree(hc, k = 3)  # R is 1-indexed → Carm 0-indexed
 */
declare function cutTree(result: HACResult, k: number): readonly number[];
/**
 * Cut a dendrogram at a specific height.
 *
 * @param result - HACResult from runHierarchical
 * @param h - height threshold
 * @returns 0-indexed cluster labels
 *
 * Cross-validate with R:
 * > cutree(hc, h = 5.0)
 */
declare function cutTreeHeight(result: HACResult, h: number): readonly number[];

export { pearsonCorrelation as A, predictGMM as B, type ClusterDiagnostics as C, type DBSCANOptions as D, predictKMeans as E, runDBSCAN as F, type GMMOptions as G, type HACMerge as H, runHierarchical as I, runKMeans as J, type KMeansOptions as K, type LCAOptions as L, silhouetteScores as M, spearmanCorrelation as N, type PointType as P, type CorrelationMatrix as a, type CovarianceModel as b, type DBSCANResult as c, type GMMRangeEntry as d, type GMMResult as e, type HACOptions as f, type HACResult as g, type KMeansRangeEntry as h, type KMeansResult as i, type LCAResult as j, type LTAOptions as k, type LTAResult as l, type LinkageMethod as m, correlationMatrix as n, cutTree as o, cutTreeHeight as p, euclideanDistMatrix as q, findBestGMM as r, fitGMM as s, fitGMMRange as t, fitKMeansRange as u, fitLCA as v, fitLTA as w, kDistancePlot as x, kendallTau as y, partialCorrelation as z };
