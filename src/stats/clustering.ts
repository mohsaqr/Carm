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

import { roundTo } from '../core/math.js'
import { Matrix } from '../core/matrix.js'
import { PRNG } from '../core/prng.js'

// ─── Constants ───────────────────────────────────────────────────────────

const LOG_2PI = Math.log(2 * Math.PI)
const MIN_PROB = 1e-300

// ─── Shared Utilities ────────────────────────────────────────────────────

function logSumExp(arr: ArrayLike<number>): number {
  let max = -Infinity
  for (let i = 0; i < arr.length; i++) {
    if (arr[i]! > max) max = arr[i]!
  }
  if (max === -Infinity) return -Infinity
  let sum = 0
  for (let i = 0; i < arr.length; i++) sum += Math.exp(arr[i]! - max)
  return max + Math.log(sum)
}

/** Raw entropy E = -sum_i sum_k (z_ik * log(z_ik)). Used for ICL = BIC + 2E. */
function computeRawEntropy(resp: ArrayLike<number>[], k: number): number {
  let ent = 0
  for (let i = 0; i < resp.length; i++) {
    for (let j = 0; j < k; j++) {
      const z = resp[i]![j]!
      if (z > MIN_PROB) ent -= z * Math.log(z)
    }
  }
  return ent
}

/** Normalized entropy = 1 - E / (N * log(K)), in [0,1]. 1 = perfect separation. */
function computeNormalizedEntropy(resp: ArrayLike<number>[], k: number): number {
  if (k <= 1) return 1  // trivial case: one cluster
  const rawE = computeRawEntropy(resp, k)
  const n = resp.length
  const denom = n * Math.log(k)
  return denom > 0 ? 1 - rawE / denom : 1
}

function computeAvePP(resp: ArrayLike<number>[], k: number): number[] {
  const sums = new Float64Array(k)
  const counts = new Float64Array(k)
  for (let i = 0; i < resp.length; i++) {
    let maxP = -1, best = 0
    for (let j = 0; j < k; j++) {
      if (resp[i]![j]! > maxP) { maxP = resp[i]![j]!; best = j }
    }
    sums[best]! += maxP
    counts[best]!++
  }
  return Array.from(sums).map((s, i) => counts[i]! > 0 ? s / counts[i]! : 0)
}

// ─── Types ───────────────────────────────────────────────────────────────

export type CovarianceModel = 'VVV' | 'EEE' | 'VVI' | 'EEI' | 'VII' | 'EII'

export interface ClusterDiagnostics {
  readonly converged: boolean
  readonly iterations: number
  readonly logLikelihood: number
  readonly df: number
  readonly aic: number
  readonly bic: number
  readonly icl: number
  readonly entropy: number
  readonly avepp: readonly number[]
  readonly formatted: string
}

export interface GMMOptions {
  readonly k: number
  readonly model?: CovarianceModel
  readonly seed?: number
  readonly tol?: number
  readonly maxIter?: number
  readonly regCovar?: number
}

export interface GMMResult {
  readonly weights: readonly number[]
  readonly means: readonly number[][]
  readonly covariances: readonly Matrix[]
  readonly posteriors: readonly (readonly number[])[]
  readonly labels: readonly number[]
  readonly diagnostics: ClusterDiagnostics
}

export interface LCAOptions {
  readonly k: number
  readonly seed?: number
  readonly tol?: number
  readonly maxIter?: number
}

export interface LCAResult {
  readonly rho: readonly (readonly number[])[]
  readonly priorWeights: readonly number[]
  readonly posteriors: readonly (readonly number[])[]
  readonly labels: readonly number[]
  readonly diagnostics: ClusterDiagnostics
}

export interface LTAOptions {
  readonly k: number
  readonly seed?: number
  readonly tol?: number
  readonly maxIter?: number
}

export interface LTAResult {
  readonly pi: readonly number[]
  readonly tau: readonly (readonly number[])[]
  readonly rho: readonly (readonly number[])[]
  readonly trajectories: readonly (readonly number[])[]
  readonly posteriors: readonly (readonly (readonly number[])[])[]
  readonly diagnostics: ClusterDiagnostics
}

export interface KMeansOptions {
  readonly k: number
  readonly seed?: number
  readonly maxIter?: number
  readonly tol?: number
}

export interface KMeansResult {
  readonly centroids: readonly (readonly number[])[]
  readonly labels: readonly number[]
  readonly inertia: number
  readonly converged: boolean
  readonly iterations: number
}

// ─── K-Means++ Initialization ────────────────────────────────────────────

function kMeansPlusPlus(
  data: readonly (readonly number[])[],
  k: number,
  rng: PRNG
): number[][] {
  const n = data.length
  const d = data[0]!.length
  const means: number[][] = [[...data[Math.floor(rng.next() * n)]!]]
  const dists = new Float64Array(n).fill(Infinity)

  for (let j = 1; j < k; j++) {
    const lastMean = means[j - 1]!
    let sumSqDist = 0

    for (let i = 0; i < n; i++) {
      const pt = data[i]!
      let dSq = 0
      for (let dim = 0; dim < d; dim++) {
        const diff = pt[dim]! - lastMean[dim]!
        dSq += diff * diff
      }
      if (dSq < dists[i]!) dists[i] = dSq
      sumSqDist += dists[i]!
    }

    let target = rng.next() * sumSqDist
    let cumulative = 0
    for (let i = 0; i < n; i++) {
      cumulative += dists[i]!
      if (cumulative >= target) {
        means.push([...data[i]!])
        break
      }
    }
    // Safety: if rounding prevents selection, pick last point
    if (means.length <= j) means.push([...data[n - 1]!])
  }
  return means
}

// ─── Multivariate Normal Log-PDF ─────────────────────────────────────────

/**
 * Log-PDF of multivariate normal using eigendecomposition.
 * Σ = U · diag(S) · U^T, so Σ^{-1} = U · diag(1/S) · U^T.
 * Mahalanobis distance = Σ (U^T · (x-μ))_i² / S_i.
 */
function mvnLogPdf(
  x: readonly number[],
  mu: readonly number[],
  uFlat: readonly number[],  // U eigenvector matrix, row-major
  eigenvals: readonly number[],
  d: number
): number {
  let mahal = 0
  let logDet = 0
  for (let i = 0; i < d; i++) {
    // Project (x - mu) onto i-th eigenvector (column i of U)
    let yi = 0
    for (let j = 0; j < d; j++) {
      yi += uFlat[j * d + i]! * (x[j]! - mu[j]!)
    }
    const eig = eigenvals[i]!
    mahal += (yi * yi) / eig
    logDet += Math.log(eig)
  }
  return -0.5 * (d * LOG_2PI + logDet + mahal)
}

// ═════════════════════════════════════════════════════════════════════════
// GMM
// ═════════════════════════════════════════════════════════════════════════

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
export function fitGMM(
  data: readonly (readonly number[])[],
  options: GMMOptions
): GMMResult {
  const n = data.length
  if (n === 0) throw new Error('fitGMM: data cannot be empty')
  const d = data[0]!.length
  const k = options.k
  if (k < 1) throw new Error('fitGMM: k must be >= 1')
  if (k > n) throw new Error('fitGMM: k cannot exceed n')

  const modelType = options.model ?? 'VVV'
  const rng = new PRNG(options.seed ?? 42)
  const tol = options.tol ?? 1e-6
  const maxIter = options.maxIter ?? 200
  const regCovar = options.regCovar ?? 1e-6

  // Initialize with K-Means++
  const means = kMeansPlusPlus(data, k, rng)
  const weights = new Float64Array(k).fill(1 / k)
  const resp = Array.from({ length: n }, () => new Float64Array(k))

  // Eigenvalues (S) and eigenvector flats (U) per component
  const I = Matrix.identity(d)
  const iFlat = I.toFlat()
  const uFlats: number[][] = Array.from({ length: k }, () => [...iFlat])
  const eigenvals: number[][] = Array.from({ length: k }, () => new Array<number>(d).fill(1))

  // Initialize eigenvalues from global variance
  let globalVar = regCovar
  const colMeans: number[] = Array.from({ length: d }, (_, j) => {
    let s = 0
    for (let i = 0; i < n; i++) s += data[i]![j]!
    return s / n
  })
  for (let i = 0; i < n; i++) {
    for (let dim = 0; dim < d; dim++) {
      globalVar += ((data[i]![dim]! - colMeans[dim]!) ** 2) / (n * d)
    }
  }
  for (let j = 0; j < k; j++) eigenvals[j]!.fill(globalVar)

  // EM loop
  let prevLogL = -Infinity
  let converged = false
  let iter = 0

  for (; iter < maxIter; iter++) {
    // E-step: compute responsibilities
    let currentLogL = 0
    const logLiks = new Float64Array(k)

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < k; j++) {
        logLiks[j] = Math.log(Math.max(weights[j]!, MIN_PROB)) +
          mvnLogPdf(data[i]!, means[j]!, uFlats[j]!, eigenvals[j]!, d)
      }
      const marg = logSumExp(logLiks)
      currentLogL += marg
      for (let j = 0; j < k; j++) {
        resp[i]![j] = Math.exp(logLiks[j]! - marg)
      }
    }

    if (Math.abs(currentLogL - prevLogL) < tol) {
      converged = true
      break
    }
    prevLogL = currentLogL

    // M-step: update weights, means, covariances
    const Nk = new Float64Array(k)
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < k; j++) Nk[j]! += resp[i]![j]!
    }

    // Pooled covariance accumulator for EEE/EEI/EII models
    const needPool = modelType === 'EEE' || modelType === 'EEI' || modelType === 'EII'
    const pool: number[][] | null = needPool
      ? Array.from({ length: d }, () => new Array<number>(d).fill(0))
      : null

    for (let j = 0; j < k; j++) {
      const nk = Math.max(Nk[j]!, MIN_PROB)
      weights[j] = Nk[j]! / n

      // Update means
      const mu = means[j]!
      mu.fill(0)
      for (let i = 0; i < n; i++) {
        const r = resp[i]![j]!
        for (let dim = 0; dim < d; dim++) mu[dim]! += r * data[i]![dim]!
      }
      for (let dim = 0; dim < d; dim++) mu[dim]! /= nk

      // Compute empirical covariance
      const cov: number[][] = Array.from({ length: d }, () => new Array<number>(d).fill(0))
      for (let i = 0; i < n; i++) {
        const r = resp[i]![j]!
        for (let r_idx = 0; r_idx < d; r_idx++) {
          const dr = data[i]![r_idx]! - mu[r_idx]!
          for (let c_idx = r_idx; c_idx < d; c_idx++) {
            const val = r * dr * (data[i]![c_idx]! - mu[c_idx]!)
            cov[r_idx]![c_idx]! += val
          }
        }
      }
      // Symmetrize and regularize
      for (let r_idx = 0; r_idx < d; r_idx++) {
        for (let c_idx = r_idx; c_idx < d; c_idx++) {
          cov[r_idx]![c_idx]! /= nk
          cov[c_idx]![r_idx] = cov[r_idx]![c_idx]!
        }
        cov[r_idx]![r_idx]! += regCovar
      }

      // Accumulate pool for EEE/EEI/EII
      if (pool) {
        const w = weights[j]!
        for (let r_idx = 0; r_idx < d; r_idx++) {
          for (let c_idx = 0; c_idx < d; c_idx++) {
            pool[r_idx]![c_idx]! += w * cov[r_idx]![c_idx]!
          }
        }
      }

      // Apply per-component constraints (VVV, VVI, VII)
      if (modelType === 'VVV') {
        const { values, vectors } = Matrix.fromArray(cov).eigen()
        eigenvals[j] = values.map(v => Math.max(v, regCovar))
        uFlats[j] = vectors.toFlat()
      } else if (modelType === 'VVI') {
        // Diagonal covariance — eigenvalues are diagonal entries, U = I
        for (let dim = 0; dim < d; dim++) eigenvals[j]![dim] = cov[dim]![dim]!
        uFlats[j] = [...iFlat]
      } else if (modelType === 'VII') {
        // Spherical — single variance = trace(Σ) / d
        let trace = 0
        for (let dim = 0; dim < d; dim++) trace += cov[dim]![dim]!
        eigenvals[j]!.fill(trace / d)
        uFlats[j] = [...iFlat]
      }
      // EEE/EEI/EII applied after pool is complete (below)
    }

    // Apply pooled constraints
    if (pool) {
      if (modelType === 'EEE') {
        const { values, vectors } = Matrix.fromArray(pool).eigen()
        const sharedEig = values.map(v => Math.max(v, regCovar))
        const sharedU = vectors.toFlat()
        for (let j = 0; j < k; j++) {
          eigenvals[j] = [...sharedEig]
          uFlats[j] = [...sharedU]
        }
      } else if (modelType === 'EEI') {
        // Shared diagonal
        const sharedDiag = new Array<number>(d)
        for (let dim = 0; dim < d; dim++) sharedDiag[dim] = Math.max(pool[dim]![dim]!, regCovar)
        for (let j = 0; j < k; j++) {
          eigenvals[j] = [...sharedDiag]
          uFlats[j] = [...iFlat]
        }
      } else if (modelType === 'EII') {
        // Shared scalar × I
        let trace = 0
        for (let dim = 0; dim < d; dim++) trace += pool[dim]![dim]!
        const scalar = Math.max(trace / d, regCovar)
        for (let j = 0; j < k; j++) {
          eigenvals[j]!.fill(scalar)
          uFlats[j] = [...iFlat]
        }
      }
    }
  }

  // Hard labels
  const labels: number[] = new Array(n)
  for (let i = 0; i < n; i++) {
    let maxP = -1, best = 0
    for (let j = 0; j < k; j++) {
      if (resp[i]![j]! > maxP) { maxP = resp[i]![j]!; best = j }
    }
    labels[i] = best
  }

  // Reconstruct covariance matrices: Σ = U · diag(S) · U^T
  const covariances = eigenvals.map((eig, j) => {
    const U = new Matrix(d, d, uFlats[j]!)
    const D = Matrix.fromArray(Array.from({ length: d }, (_, r) =>
      Array.from({ length: d }, (_, c) => r === c ? eig[r]! : 0)
    ))
    return U.multiply(D).multiply(U.transpose())
  })

  // Diagnostics
  const logL = prevLogL === -Infinity ? 0 : prevLogL
  const dfMap: Record<CovarianceModel, number> = {
    'VVV': k * d * (d + 1) / 2,
    'EEE': d * (d + 1) / 2,
    'VVI': k * d,
    'EEI': d,
    'VII': k,
    'EII': 1,
  }
  const df = (k - 1) + (k * d) + dfMap[modelType]
  const rawEntropy = computeRawEntropy(resp, k)
  const entropy = computeNormalizedEntropy(resp, k)
  const bic = df * Math.log(n) - 2 * logL
  const aic = 2 * df - 2 * logL
  const icl = bic + 2 * rawEntropy

  return {
    weights: Array.from(weights),
    means: means.map(m => [...m]),
    covariances,
    posteriors: resp.map(r => Array.from(r)),
    labels,
    diagnostics: {
      converged,
      iterations: iter,
      logLikelihood: logL,
      df,
      aic,
      bic,
      icl,
      entropy,
      avepp: computeAvePP(resp, k),
      formatted: `GMM (K = ${k}, ${modelType}): BIC = ${roundTo(bic, 1)}, AIC = ${roundTo(aic, 1)}, LL = ${roundTo(logL, 1)}, AvePP = [${computeAvePP(resp, k).map(v => roundTo(v, 2)).join(', ')}]`,
    },
  }
}

/**
 * Predict cluster assignments for new data given a fitted GMM.
 */
export function predictGMM(
  data: readonly (readonly number[])[],
  result: GMMResult
): { readonly labels: readonly number[]; readonly posteriors: readonly (readonly number[])[] } {
  const k = result.weights.length
  const d = result.means[0]!.length
  const uFlats = result.covariances.map(cov => {
    const { vectors } = cov.eigen()
    return vectors.toFlat()
  })
  const eigVals = result.covariances.map(cov => {
    const { values } = cov.eigen()
    return values.map(v => Math.max(v, 1e-12))
  })

  const n = data.length
  const labels: number[] = new Array(n)
  const posteriors: number[][] = new Array(n)
  const logLiks = new Float64Array(k)

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < k; j++) {
      logLiks[j] = Math.log(Math.max(result.weights[j]!, MIN_PROB)) +
        mvnLogPdf(data[i]!, result.means[j]!, uFlats[j]!, eigVals[j]!, d)
    }
    const marg = logSumExp(logLiks)
    const post = new Array<number>(k)
    let maxP = -1, best = 0
    for (let j = 0; j < k; j++) {
      post[j] = Math.exp(logLiks[j]! - marg)
      if (post[j]! > maxP) { maxP = post[j]!; best = j }
    }
    labels[i] = best
    posteriors[i] = post
  }
  return { labels, posteriors }
}

/**
 * Automatic model selection: fit GMM across a grid of K and covariance models,
 * return the model with the lowest BIC.
 *
 * @param data - N × D data matrix
 * @param kRange - Array of K values to try (default [1,2,3,4,5])
 * @param models - Array of covariance models to try (default all 6)
 * @returns The GMMResult with the lowest BIC
 */
export function findBestGMM(
  data: readonly (readonly number[])[],
  kRange: readonly number[] = [1, 2, 3, 4, 5],
  models: readonly CovarianceModel[] = ['VVV', 'EEE', 'VVI', 'EEI', 'VII', 'EII']
): GMMResult {
  let best: GMMResult | null = null
  for (const k of kRange) {
    for (const model of models) {
      try {
        const res = fitGMM(data, { k, model })
        if (!best || res.diagnostics.bic < best.diagnostics.bic) best = res
      } catch {
        // Some model/k combos may fail (e.g., singular covariance) — skip
      }
    }
  }
  if (!best) throw new Error('findBestGMM: all model fits failed')
  return best
}

// ═════════════════════════════════════════════════════════════════════════
// LCA
// ═════════════════════════════════════════════════════════════════════════

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
export function fitLCA(
  data: readonly (readonly number[])[],
  options: LCAOptions
): LCAResult {
  const n = data.length
  if (n === 0) throw new Error('fitLCA: data cannot be empty')
  const m = data[0]!.length
  const k = options.k
  if (k < 1) throw new Error('fitLCA: k must be >= 1')

  const rng = new PRNG(options.seed ?? 42)
  const tol = options.tol ?? 1e-6
  const maxIter = options.maxIter ?? 200

  // Initialize rho randomly in (0.1, 0.9)
  const rho: number[][] = Array.from({ length: k }, () =>
    Array.from({ length: m }, () => 0.1 + rng.next() * 0.8)
  )
  const weights = new Float64Array(k).fill(1 / k)
  const resp = Array.from({ length: n }, () => new Float64Array(k))

  let prevLogL = -Infinity
  let converged = false
  let iter = 0

  for (; iter < maxIter; iter++) {
    // E-step
    let currentLogL = 0
    const logLiks = new Float64Array(k)

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < k; j++) {
        let ll = Math.log(Math.max(weights[j]!, MIN_PROB))
        for (let d = 0; d < m; d++) {
          const r = Math.max(Math.min(rho[j]![d]!, 1 - 1e-12), 1e-12)
          ll += data[i]![d] === 1 ? Math.log(r) : Math.log(1 - r)
        }
        logLiks[j] = ll
      }
      const marg = logSumExp(logLiks)
      currentLogL += marg
      for (let j = 0; j < k; j++) resp[i]![j] = Math.exp(logLiks[j]! - marg)
    }

    if (Math.abs(currentLogL - prevLogL) < tol) {
      converged = true
      break
    }
    prevLogL = currentLogL

    // M-step
    const Nk = new Float64Array(k)
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < k; j++) Nk[j]! += resp[i]![j]!
    }

    for (let j = 0; j < k; j++) {
      weights[j] = Nk[j]! / n
      const sumW = Nk[j]!
      for (let d = 0; d < m; d++) {
        let sumX = 0
        for (let i = 0; i < n; i++) {
          if (data[i]![d] === 1) sumX += resp[i]![j]!
        }
        // MLE with minimal floor to prevent log(0)
        rho[j]![d] = Math.max(Math.min(sumX / sumW, 1 - 1e-10), 1e-10)
      }
    }
  }

  // Hard labels
  const labels: number[] = new Array(n)
  for (let i = 0; i < n; i++) {
    let maxP = -1, best = 0
    for (let j = 0; j < k; j++) {
      if (resp[i]![j]! > maxP) { maxP = resp[i]![j]!; best = j }
    }
    labels[i] = best
  }

  const logL = prevLogL === -Infinity ? 0 : prevLogL
  const df = (k - 1) + (k * m)
  const rawEntropy = computeRawEntropy(resp, k)
  const entropy = computeNormalizedEntropy(resp, k)
  const bic = df * Math.log(n) - 2 * logL
  const aic = 2 * df - 2 * logL
  const icl = bic + 2 * rawEntropy

  return {
    rho: rho.map(r => [...r]),
    priorWeights: Array.from(weights),
    posteriors: resp.map(r => Array.from(r)),
    labels,
    diagnostics: {
      converged,
      iterations: iter,
      logLikelihood: logL,
      df,
      aic,
      bic,
      icl,
      entropy,
      avepp: computeAvePP(resp, k),
      formatted: `LCA (K = ${k}): BIC = ${roundTo(bic, 1)}, AIC = ${roundTo(aic, 1)}, LL = ${roundTo(logL, 1)}`,
    },
  }
}

// ═════════════════════════════════════════════════════════════════════════
// LTA (Latent Transition Analysis — Hidden Markov LCA)
// ═════════════════════════════════════════════════════════════════════════

/**
 * Compute log emission probability for observation x given state s.
 * Uses Bernoulli model in log-space to prevent underflow.
 */
function logEmission(
  x: readonly number[],
  rho_s: readonly number[],
  m: number
): number {
  let ll = 0
  for (let d = 0; d < m; d++) {
    const r = Math.max(Math.min(rho_s[d]!, 1 - 1e-12), 1e-12)
    ll += x[d] === 1 ? Math.log(r) : Math.log(1 - r)
  }
  return ll
}

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
export function fitLTA(
  data: readonly (readonly (readonly number[])[])[],
  options: LTAOptions
): LTAResult {
  const n = data.length
  if (n === 0) throw new Error('fitLTA: data cannot be empty')
  const T = data[0]!.length
  if (T < 2) throw new Error('fitLTA: requires at least 2 timepoints')
  const m = data[0]![0]!.length
  const k = options.k
  if (k < 2) throw new Error('fitLTA: k must be >= 2')

  const rng = new PRNG(options.seed ?? 42)
  const tol = options.tol ?? 1e-6
  const maxIter = options.maxIter ?? 200

  // Initialize parameters
  const pi = new Float64Array(k).fill(1 / k)
  const tau: number[][] = Array.from({ length: k }, (_, j) =>
    Array.from({ length: k }, (_, l) => j === l ? 0.8 : 0.2 / (k - 1))
  )
  const rho: number[][] = Array.from({ length: k }, () =>
    Array.from({ length: m }, () => 0.1 + rng.next() * 0.8)
  )

  let prevLogL = -Infinity
  let converged = false
  let iter = 0

  // Store gamma from last iteration for diagnostics
  let finalGamma: number[][][] = []

  for (; iter < maxIter; iter++) {
    let currentLogL = 0

    // Accumulators for M-step
    const piAcc = new Float64Array(k)
    const tauNum: number[][] = Array.from({ length: k }, () => new Array<number>(k).fill(0))
    const tauDen = new Float64Array(k)
    const rhoNum: number[][] = Array.from({ length: k }, () => new Array<number>(m).fill(0))
    const rhoDen = new Float64Array(k)

    const iterGamma: number[][][] = new Array(n)

    for (let i = 0; i < n; i++) {
      // Pre-compute log emissions: T × K
      const logB: number[][] = Array.from({ length: T }, (_, t) =>
        Array.from({ length: k }, (_, s) => logEmission(data[i]![t]!, rho[s]!, m))
      )

      // Forward pass (log-alpha)
      const logAlpha: number[][] = Array.from({ length: T }, () => new Array<number>(k))
      for (let s = 0; s < k; s++) {
        logAlpha[0]![s] = Math.log(Math.max(pi[s]!, MIN_PROB)) + logB[0]![s]!
      }
      for (let t = 1; t < T; t++) {
        for (let s = 0; s < k; s++) {
          const trans = new Float64Array(k)
          for (let p = 0; p < k; p++) {
            trans[p] = logAlpha[t - 1]![p]! + Math.log(Math.max(tau[p]![s]!, MIN_PROB))
          }
          logAlpha[t]![s] = logB[t]![s]! + logSumExp(trans)
        }
      }

      // Subject log-likelihood
      const subjLL = logSumExp(logAlpha[T - 1]!)
      currentLogL += subjLL

      // Backward pass (log-beta)
      const logBeta: number[][] = Array.from({ length: T }, () => new Array<number>(k))
      for (let s = 0; s < k; s++) logBeta[T - 1]![s] = 0
      for (let t = T - 2; t >= 0; t--) {
        for (let s = 0; s < k; s++) {
          const combined = new Float64Array(k)
          for (let next = 0; next < k; next++) {
            combined[next] = Math.log(Math.max(tau[s]![next]!, MIN_PROB)) +
              logB[t + 1]![next]! + logBeta[t + 1]![next]!
          }
          logBeta[t]![s] = logSumExp(combined)
        }
      }

      // Posterior marginals (gamma)
      const gamma: number[][] = Array.from({ length: T }, () => new Array<number>(k))
      for (let t = 0; t < T; t++) {
        const logRow = new Float64Array(k)
        for (let s = 0; s < k; s++) logRow[s] = logAlpha[t]![s]! + logBeta[t]![s]!
        const den = logSumExp(logRow)
        for (let s = 0; s < k; s++) gamma[t]![s] = Math.exp(logRow[s]! - den)
      }
      iterGamma[i] = gamma

      // Xi: transition posteriors
      for (let t = 0; t < T - 1; t++) {
        const logXi = new Float64Array(k * k)
        for (let j = 0; j < k; j++) {
          for (let l = 0; l < k; l++) {
            logXi[j * k + l] = logAlpha[t]![j]! +
              Math.log(Math.max(tau[j]![l]!, MIN_PROB)) +
              logB[t + 1]![l]! + logBeta[t + 1]![l]!
          }
        }
        const xiDen = logSumExp(logXi)
        for (let j = 0; j < k; j++) {
          for (let l = 0; l < k; l++) {
            tauNum[j]![l]! += Math.exp(logXi[j * k + l]! - xiDen)
          }
        }
      }

      // Accumulate for M-step
      for (let s = 0; s < k; s++) piAcc[s]! += gamma[0]![s]!
      for (let t = 0; t < T - 1; t++) {
        for (let s = 0; s < k; s++) tauDen[s]! += gamma[t]![s]!
      }
      for (let t = 0; t < T; t++) {
        for (let s = 0; s < k; s++) {
          rhoDen[s]! += gamma[t]![s]!
          for (let d = 0; d < m; d++) {
            if (data[i]![t]![d] === 1) rhoNum[s]![d]! += gamma[t]![s]!
          }
        }
      }
    }

    finalGamma = iterGamma

    if (Math.abs(currentLogL - prevLogL) < tol) {
      converged = true
      break
    }
    prevLogL = currentLogL

    // M-step: update pi, tau, rho
    for (let s = 0; s < k; s++) {
      pi[s] = (piAcc[s]! + 1) / (n + k)  // Dirichlet(1) smoothing
    }
    for (let j = 0; j < k; j++) {
      const den = tauDen[j]!
      for (let l = 0; l < k; l++) {
        tau[j]![l] = (tauNum[j]![l]! + 0.1) / (den + 0.1 * k)  // Light smoothing
      }
    }
    for (let s = 0; s < k; s++) {
      const den = rhoDen[s]!
      for (let d = 0; d < m; d++) {
        rho[s]![d] = Math.max(Math.min(rhoNum[s]![d]! / den, 1 - 1e-10), 1e-10)
      }
    }
  }

  // Viterbi decoding
  const trajectories: number[][] = data.map(subj => {
    const vt: number[][] = Array.from({ length: T }, () => new Array<number>(k))
    const ptr: number[][] = Array.from({ length: T }, () => new Array<number>(k))

    for (let s = 0; s < k; s++) {
      vt[0]![s] = Math.log(Math.max(pi[s]!, MIN_PROB)) + logEmission(subj[0]!, rho[s]!, m)
    }
    for (let t = 1; t < T; t++) {
      for (let s = 0; s < k; s++) {
        const emit = logEmission(subj[t]!, rho[s]!, m)
        let maxVal = -Infinity, bestP = 0
        for (let p = 0; p < k; p++) {
          const sc = vt[t - 1]![p]! + Math.log(Math.max(tau[p]![s]!, MIN_PROB))
          if (sc > maxVal) { maxVal = sc; bestP = p }
        }
        vt[t]![s] = emit + maxVal
        ptr[t]![s] = bestP
      }
    }

    const path = new Array<number>(T)
    let maxFinal = -Infinity, bestFinal = 0
    for (let s = 0; s < k; s++) {
      if (vt[T - 1]![s]! > maxFinal) { maxFinal = vt[T - 1]![s]!; bestFinal = s }
    }
    path[T - 1] = bestFinal
    for (let t = T - 2; t >= 0; t--) {
      path[t] = ptr[t + 1]![path[t + 1]!]!
    }
    return path
  })

  // Diagnostics
  const logL = prevLogL === -Infinity ? 0 : prevLogL
  const df = (k - 1) + k * (k - 1) + k * m

  // Flatten gamma for entropy/avepp computation
  const flatGamma: number[][] = []
  for (let i = 0; i < n; i++) {
    for (let t = 0; t < T; t++) {
      flatGamma.push(finalGamma[i]![t]!)
    }
  }
  const rawEntropy = computeRawEntropy(flatGamma, k)
  const entropy = computeNormalizedEntropy(flatGamma, k)
  const bic = df * Math.log(n) - 2 * logL
  const aic = 2 * df - 2 * logL
  const icl = bic + 2 * rawEntropy

  return {
    pi: Array.from(pi),
    tau: tau.map(row => [...row]),
    rho: rho.map(r => [...r]),
    trajectories,
    posteriors: finalGamma.map(subj => subj.map(t => [...t])),
    diagnostics: {
      converged,
      iterations: iter,
      logLikelihood: logL,
      df,
      aic,
      bic,
      icl,
      entropy,
      avepp: computeAvePP(flatGamma, k),
      formatted: `LTA (K = ${k}, T = ${T}): BIC = ${roundTo(bic, 1)}, AIC = ${roundTo(aic, 1)}, LL = ${roundTo(logL, 1)}`,
    },
  }
}

// ═════════════════════════════════════════════════════════════════════════
// K-Means
// ═════════════════════════════════════════════════════════════════════════

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
export function runKMeans(
  data: readonly (readonly number[])[],
  options: KMeansOptions
): KMeansResult {
  const n = data.length
  if (n === 0) throw new Error('runKMeans: data cannot be empty')
  const d = data[0]!.length
  const k = options.k
  if (k < 1) throw new Error('runKMeans: k must be >= 1')
  if (k > n) throw new Error('runKMeans: k cannot exceed n')

  const rng = new PRNG(options.seed ?? 42)
  const maxIter = options.maxIter ?? 300
  const tol = options.tol ?? 1e-6

  const centroids = kMeansPlusPlus(data, k, rng)
  const labels = new Array<number>(n).fill(0)

  let converged = false
  let iter = 0
  let inertia = 0

  for (; iter < maxIter; iter++) {
    const nextCentroids: number[][] = Array.from({ length: k }, () => new Array<number>(d).fill(0))
    const counts = new Array<number>(k).fill(0)
    inertia = 0

    // Assignment step
    for (let i = 0; i < n; i++) {
      let minD = Infinity, bestK = 0
      for (let j = 0; j < k; j++) {
        let dist = 0
        for (let dim = 0; dim < d; dim++) {
          const diff = data[i]![dim]! - centroids[j]![dim]!
          dist += diff * diff
        }
        if (dist < minD) { minD = dist; bestK = j }
      }
      labels[i] = bestK
      counts[bestK]!++
      inertia += minD
      for (let dim = 0; dim < d; dim++) nextCentroids[bestK]![dim]! += data[i]![dim]!
    }

    // Update step with empty-cluster re-seeding
    let shiftSq = 0
    for (let j = 0; j < k; j++) {
      if (counts[j] === 0) {
        // Re-seed: assign the farthest point from its current centroid
        let maxDist = -1, farIdx = 0
        for (let i = 0; i < n; i++) {
          let di = 0
          for (let dim = 0; dim < d; dim++) {
            const diff = data[i]![dim]! - centroids[labels[i]!]![dim]!
            di += diff * diff
          }
          if (di > maxDist) { maxDist = di; farIdx = i }
        }
        for (let dim = 0; dim < d; dim++) {
          const updated = data[farIdx]![dim]!
          const diff = updated - centroids[j]![dim]!
          shiftSq += diff * diff
          centroids[j]![dim] = updated
        }
      } else {
        for (let dim = 0; dim < d; dim++) {
          const updated = nextCentroids[j]![dim]! / counts[j]!
          const diff = updated - centroids[j]![dim]!
          shiftSq += diff * diff
          centroids[j]![dim] = updated
        }
      }
    }

    if (shiftSq < tol) {
      converged = true
      break
    }
  }

  return {
    centroids: centroids.map(c => [...c]),
    labels,
    inertia,
    converged,
    iterations: iter,
  }
}

/**
 * Predict cluster assignments for new data given fitted K-Means centroids.
 */
// ═════════════════════════════════════════════════════════════════════════
// Range-fitting (auto-K selection)
// ═════════════════════════════════════════════════════════════════════════

/** A single entry from fitting GMM at one K value. */
export interface GMMRangeEntry {
  readonly k: number
  readonly model: CovarianceModel
  readonly result: GMMResult
}

/** A single entry from fitting KMeans at one K value. */
export interface KMeansRangeEntry {
  readonly k: number
  readonly result: KMeansResult
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
export function fitGMMRange(
  data: readonly (readonly number[])[],
  kRange: readonly number[],
  model: CovarianceModel = 'VVV',
): readonly GMMRangeEntry[] {
  const entries: GMMRangeEntry[] = []
  for (const k of kRange) {
    try {
      const result = fitGMM(data, { k, model, seed: 42 })
      entries.push({ k, model, result })
    } catch {
      // Skip failed fits (e.g. k > n, singular covariance)
    }
  }
  if (entries.length === 0) throw new Error('fitGMMRange: all fits failed')
  return entries.sort((a, b) => a.k - b.k)
}

/**
 * Fit KMeans for each K in kRange and return results sorted by K.
 * Skips failed fits. Throws only if ALL fail.
 *
 * @param data - N × D numeric data matrix
 * @param kRange - Array of K values to try
 * @returns KMeansRangeEntry[] sorted by K
 */
export function fitKMeansRange(
  data: readonly (readonly number[])[],
  kRange: readonly number[],
): readonly KMeansRangeEntry[] {
  const entries: KMeansRangeEntry[] = []
  for (const k of kRange) {
    try {
      const result = runKMeans(data, { k, seed: 42 })
      entries.push({ k, result })
    } catch {
      // Skip failed fits
    }
  }
  if (entries.length === 0) throw new Error('fitKMeansRange: all fits failed')
  return entries.sort((a, b) => a.k - b.k)
}

export function predictKMeans(
  data: readonly (readonly number[])[],
  centroids: readonly (readonly number[])[]
): readonly number[] {
  const k = centroids.length
  const d = centroids[0]!.length
  return data.map(pt => {
    let minD = Infinity, best = 0
    for (let j = 0; j < k; j++) {
      let dist = 0
      for (let dim = 0; dim < d; dim++) {
        const diff = pt[dim]! - centroids[j]![dim]!
        dist += diff * diff
      }
      if (dist < minD) { minD = dist; best = j }
    }
    return best
  })
}

// ═════════════════════════════════════════════════════════════════════════
// Shared Utilities: Distance Matrix & Silhouette
// ═════════════════════════════════════════════════════════════════════════

/**
 * Compute Euclidean distance matrix (N × N, row-major Float64Array).
 *
 * @param data - N × D numeric data matrix
 * @returns flat Float64Array of size N*N with dist[i*N+j] = euclidean(i,j)
 *
 * Cross-validate with R:
 * > as.matrix(dist(data))
 */
export function euclideanDistMatrix(
  data: readonly (readonly number[])[]
): Float64Array {
  const n = data.length
  const d = data[0]?.length ?? 0
  const dist = new Float64Array(n * n)

  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      let sq = 0
      for (let dim = 0; dim < d; dim++) {
        const diff = data[i]![dim]! - data[j]![dim]!
        sq += diff * diff
      }
      const val = Math.sqrt(sq)
      dist[i * n + j] = val
      dist[j * n + i] = val
    }
  }
  return dist
}

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
export function silhouetteScores(
  data: readonly (readonly number[])[],
  labels: readonly number[]
): { readonly scores: readonly number[]; readonly mean: number } {
  const n = data.length
  const dist = euclideanDistMatrix(data)

  // Find unique non-noise labels
  const labelSet = new Set<number>()
  for (let i = 0; i < n; i++) {
    if (labels[i]! >= 0) labelSet.add(labels[i]!)
  }
  const uniqueLabels = [...labelSet].sort((a, b) => a - b)
  const nClusters = uniqueLabels.length

  if (nClusters < 2) {
    // Silhouette undefined for < 2 clusters
    return { scores: new Array(n).fill(0) as number[], mean: 0 }
  }

  const scores: number[] = new Array(n).fill(NaN) as number[]
  let sumSil = 0
  let countNonNoise = 0

  for (let i = 0; i < n; i++) {
    const ci = labels[i]!
    if (ci < 0) continue // noise

    // a(i) = mean distance to same cluster
    let aSum = 0, aCount = 0
    for (let j = 0; j < n; j++) {
      if (j === i || labels[j] !== ci) continue
      aSum += dist[i * n + j]!
      aCount++
    }
    const ai = aCount > 0 ? aSum / aCount : 0

    // b(i) = min mean-distance to other clusters
    let bi = Infinity
    for (const cl of uniqueLabels) {
      if (cl === ci) continue
      let bSum = 0, bCount = 0
      for (let j = 0; j < n; j++) {
        if (labels[j] !== cl) continue
        bSum += dist[i * n + j]!
        bCount++
      }
      if (bCount > 0) {
        const meanDist = bSum / bCount
        if (meanDist < bi) bi = meanDist
      }
    }

    const maxAB = Math.max(ai, bi)
    scores[i] = maxAB > 0 ? (bi - ai) / maxAB : 0
    sumSil += scores[i]!
    countNonNoise++
  }

  return {
    scores,
    mean: countNonNoise > 0 ? sumSil / countNonNoise : 0,
  }
}

// ═════════════════════════════════════════════════════════════════════════
// DBSCAN
// ═════════════════════════════════════════════════════════════════════════

export type PointType = 'core' | 'border' | 'noise'

export interface DBSCANOptions {
  readonly eps: number
  readonly minPts: number
  readonly seed?: number
  readonly preprocess?: 'none' | 'center' | 'standardize' | 'log' | 'sqrt'
}

export interface DBSCANResult {
  readonly labels: readonly number[]           // -1 = noise, 0-indexed clusters
  readonly pointTypes: readonly PointType[]
  readonly nClusters: number
  readonly nNoise: number
  readonly clusterSizes: readonly number[]
  readonly silhouette: { readonly scores: readonly number[]; readonly mean: number }
  readonly formatted: string
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
export function runDBSCAN(
  data: readonly (readonly number[])[],
  options: DBSCANOptions
): DBSCANResult {
  const n = data.length
  if (n === 0) throw new Error('runDBSCAN: data cannot be empty')
  const { eps, minPts } = options
  if (eps <= 0) throw new Error('runDBSCAN: eps must be > 0')
  if (minPts < 1) throw new Error('runDBSCAN: minPts must be >= 1')

  // Compute pairwise distance matrix
  const dist = euclideanDistMatrix(data)

  // Find eps-neighborhoods
  const neighbors: number[][] = new Array(n)
  for (let i = 0; i < n; i++) {
    const nb: number[] = []
    for (let j = 0; j < n; j++) {
      if (dist[i * n + j]! <= eps) nb.push(j)
    }
    neighbors[i] = nb
  }

  // Identify core points
  const isCore = new Uint8Array(n)
  for (let i = 0; i < n; i++) {
    if (neighbors[i]!.length >= minPts) isCore[i] = 1
  }

  // BFS cluster expansion
  const labels = new Int32Array(n).fill(-1)  // -1 = unvisited/noise
  let clusterId = 0

  for (let i = 0; i < n; i++) {
    if (!isCore[i] || labels[i] !== -1) continue

    // Start new cluster from this core point
    const queue: number[] = [i]
    labels[i] = clusterId

    let head = 0
    while (head < queue.length) {
      const current = queue[head]!
      head++

      for (const nb of neighbors[current]!) {
        if (labels[nb] === -1) {
          labels[nb] = clusterId
          // If neighbor is also core, expand from it
          if (isCore[nb]) queue.push(nb)
        }
      }
    }
    clusterId++
  }

  // Classify point types
  const pointTypes: PointType[] = new Array(n)
  for (let i = 0; i < n; i++) {
    if (isCore[i]) {
      pointTypes[i] = 'core'
    } else if (labels[i]! >= 0) {
      pointTypes[i] = 'border'
    } else {
      pointTypes[i] = 'noise'
    }
  }

  // Cluster sizes
  const nClusters = clusterId
  const clusterSizes = new Array<number>(nClusters).fill(0)
  let nNoise = 0
  for (let i = 0; i < n; i++) {
    if (labels[i]! >= 0) clusterSizes[labels[i]!]!++
    else nNoise++
  }

  // Silhouette (only if >= 2 clusters and some non-noise points)
  const sil = nClusters >= 2
    ? silhouetteScores(data, Array.from(labels))
    : { scores: new Array(n).fill(NaN) as number[], mean: NaN }

  const formatted = `DBSCAN (eps = ${roundTo(eps, 3)}, minPts = ${minPts}): ${nClusters} clusters, ${nNoise} noise points, silhouette = ${isNaN(sil.mean) ? 'N/A' : roundTo(sil.mean, 3)}`

  return {
    labels: Array.from(labels),
    pointTypes,
    nClusters,
    nNoise,
    clusterSizes,
    silhouette: sil,
    formatted,
  }
}

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
export function kDistancePlot(
  data: readonly (readonly number[])[],
  k: number
): readonly number[] {
  const n = data.length
  if (k < 1 || k >= n) throw new Error(`kDistancePlot: k must be in [1, n-1], got ${k}`)

  const dist = euclideanDistMatrix(data)
  const kDists: number[] = new Array(n)

  for (let i = 0; i < n; i++) {
    // Collect distances from point i to all others
    const dists: number[] = new Array(n - 1)
    let idx = 0
    for (let j = 0; j < n; j++) {
      if (j !== i) dists[idx++] = dist[i * n + j]!
    }
    dists.sort((a, b) => a - b)
    kDists[i] = dists[k - 1]!  // k-th nearest (0-indexed, so k-1)
  }

  return kDists.sort((a, b) => a - b)
}

// ═════════════════════════════════════════════════════════════════════════
// Hierarchical Agglomerative Clustering (HAC)
// ═════════════════════════════════════════════════════════════════════════

export type LinkageMethod = 'single' | 'complete' | 'average' | 'ward'

export interface HACOptions {
  readonly linkage?: LinkageMethod
  readonly preprocess?: 'none' | 'center' | 'standardize' | 'log' | 'sqrt'
}

export interface HACMerge {
  readonly a: number   // index of first cluster merged (0..n-1 = obs, n+ = intermediate)
  readonly b: number   // index of second cluster merged
  readonly height: number  // merge height
}

export interface HACResult {
  readonly merges: readonly HACMerge[]
  readonly heights: readonly number[]
  readonly order: readonly number[]   // leaf order for dendrogram (DFS)
  readonly copheneticCorrelation: number
  readonly formatted: string
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
export function runHierarchical(
  data: readonly (readonly number[])[],
  options?: HACOptions
): HACResult {
  const n = data.length
  if (n < 2) throw new Error('runHierarchical: need at least 2 observations')

  const linkage = options?.linkage ?? 'ward'

  // Compute initial pairwise distance matrix (Euclidean)
  const origDist = euclideanDistMatrix(data)

  // For ward, we work with squared distances internally
  const useSquared = linkage === 'ward'
  // Active distance matrix — condensed upper triangle
  // We use a full N×N mutable copy for simplicity during merges
  const D = new Float64Array(n * n)
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const val = useSquared
        ? origDist[i * n + j]! * origDist[i * n + j]!
        : origDist[i * n + j]!
      D[i * n + j] = val
      D[j * n + i] = val
    }
  }

  // Cluster sizes
  const sizes = new Float64Array(2 * n - 1)
  for (let i = 0; i < n; i++) sizes[i] = 1

  // Active cluster set
  const active = new Set<number>()
  for (let i = 0; i < n; i++) active.add(i)

  // Merge tracking
  const merges: HACMerge[] = []
  const mergeHeights: number[] = []

  // Current dimension for distance lookups (grows as we add merged clusters)
  // We'll extend D as needed using a map for merged cluster distances
  const distMap = new Map<string, number>()

  function getDist(a: number, b: number): number {
    if (a < n && b < n) return D[a * n + b]!
    const key = a < b ? `${a},${b}` : `${b},${a}`
    return distMap.get(key) ?? Infinity
  }

  function setDist(a: number, b: number, val: number): void {
    if (a < n && b < n) {
      D[a * n + b] = val
      D[b * n + a] = val
    } else {
      const key = a < b ? `${a},${b}` : `${b},${a}`
      distMap.set(key, val)
    }
  }

  let nextCluster = n  // merged clusters start at index n

  for (let step = 0; step < n - 1; step++) {
    // Find closest pair among active clusters
    let minDist = Infinity
    let mergeA = -1, mergeB = -1
    const activeArr = [...active]
    for (let ii = 0; ii < activeArr.length; ii++) {
      const ci = activeArr[ii]!
      for (let jj = ii + 1; jj < activeArr.length; jj++) {
        const cj = activeArr[jj]!
        const d = getDist(ci, cj)
        if (d < minDist) {
          minDist = d
          mergeA = ci
          mergeB = cj
        }
      }
    }

    // Record merge height
    const height = useSquared ? Math.sqrt(minDist) : minDist
    merges.push({ a: mergeA, b: mergeB, height })
    mergeHeights.push(height)

    // New merged cluster
    const newCluster = nextCluster++
    sizes[newCluster] = sizes[mergeA]! + sizes[mergeB]!
    active.delete(mergeA)
    active.delete(mergeB)

    // Update distances to new cluster via Lance-Williams
    const ni = sizes[mergeA]!
    const nj = sizes[mergeB]!

    for (const ck of active) {
      const nk = sizes[ck]!
      const dik = getDist(mergeA, ck)
      const djk = getDist(mergeB, ck)
      const dij = getDist(mergeA, mergeB)

      let newDist: number
      if (linkage === 'single') {
        newDist = 0.5 * dik + 0.5 * djk - 0.5 * Math.abs(dik - djk)
      } else if (linkage === 'complete') {
        newDist = 0.5 * dik + 0.5 * djk + 0.5 * Math.abs(dik - djk)
      } else if (linkage === 'average') {
        const ai = ni / (ni + nj)
        const aj = nj / (ni + nj)
        newDist = ai * dik + aj * djk
      } else {
        // ward
        const nt = ni + nj + nk
        newDist = ((ni + nk) / nt) * dik + ((nj + nk) / nt) * djk - (nk / nt) * dij
      }
      setDist(newCluster, ck, newDist)
    }

    active.add(newCluster)
  }

  // Build leaf order via DFS of the dendrogram tree
  const order = buildLeafOrder(merges, n)

  // Cophenetic correlation
  const cophCorr = computeCopheneticCorrelation(merges, n, origDist)

  const formatted = `HAC (${linkage}): ${n} observations, cophenetic r = ${roundTo(cophCorr, 3)}`

  return {
    merges,
    heights: mergeHeights,
    order,
    copheneticCorrelation: cophCorr,
    formatted,
  }
}

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
export function cutTree(
  result: HACResult,
  k: number
): readonly number[] {
  const nMerges = result.merges.length
  const n = nMerges + 1
  if (k < 1 || k > n) throw new Error(`cutTree: k must be in [1, ${n}], got ${k}`)

  // Apply the first (n - k) merges using union-find
  const parent = new Int32Array(2 * n - 1)
  for (let i = 0; i < 2 * n - 1; i++) parent[i] = i

  function find(x: number): number {
    while (parent[x] !== x) {
      parent[x] = parent[parent[x]!]!
      x = parent[x]!
    }
    return x
  }

  const nMergesToApply = n - k
  for (let s = 0; s < nMergesToApply; s++) {
    const merge = result.merges[s]!
    const newCluster = n + s
    parent[find(merge.a)] = newCluster
    parent[find(merge.b)] = newCluster
  }

  // Assign labels: find root of each original observation
  const rootToLabel = new Map<number, number>()
  let nextLabel = 0
  const labels = new Array<number>(n)

  for (let i = 0; i < n; i++) {
    const root = find(i)
    let label = rootToLabel.get(root)
    if (label === undefined) {
      label = nextLabel++
      rootToLabel.set(root, label)
    }
    labels[i] = label
  }

  return labels
}

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
export function cutTreeHeight(
  result: HACResult,
  h: number
): readonly number[] {
  const nMerges = result.merges.length
  const n = nMerges + 1

  // Apply all merges with height <= h
  const parent = new Int32Array(2 * n - 1)
  for (let i = 0; i < 2 * n - 1; i++) parent[i] = i

  function find(x: number): number {
    while (parent[x] !== x) {
      parent[x] = parent[parent[x]!]!
      x = parent[x]!
    }
    return x
  }

  for (let s = 0; s < nMerges; s++) {
    const merge = result.merges[s]!
    if (merge.height > h) break
    const newCluster = n + s
    parent[find(merge.a)] = newCluster
    parent[find(merge.b)] = newCluster
  }

  const rootToLabel = new Map<number, number>()
  let nextLabel = 0
  const labels = new Array<number>(n)

  for (let i = 0; i < n; i++) {
    const root = find(i)
    let label = rootToLabel.get(root)
    if (label === undefined) {
      label = nextLabel++
      rootToLabel.set(root, label)
    }
    labels[i] = label
  }

  return labels
}

// ─── HAC helpers ─────────────────────────────────────────────────────────

/** Build leaf order via DFS traversal of the dendrogram. */
function buildLeafOrder(merges: readonly HACMerge[], n: number): number[] {
  // Build tree structure: children[nodeId] = [left, right]
  const children = new Map<number, [number, number]>()
  for (let s = 0; s < merges.length; s++) {
    const merge = merges[s]!
    const nodeId = n + s
    children.set(nodeId, [merge.a, merge.b])
  }

  const root = n + merges.length - 1
  const order: number[] = []

  function dfs(node: number): void {
    const ch = children.get(node)
    if (ch) {
      dfs(ch[0])
      dfs(ch[1])
    } else {
      order.push(node)
    }
  }

  dfs(root)
  return order
}

/**
 * Compute cophenetic correlation: Pearson(original distances, cophenetic distances).
 * The cophenetic distance between two points is the merge height at which
 * they first join the same cluster.
 */
function computeCopheneticCorrelation(
  merges: readonly HACMerge[],
  n: number,
  origDist: Float64Array
): number {
  // Build union-find to track when pairs merge
  const parent = new Int32Array(2 * n - 1)
  for (let i = 0; i < 2 * n - 1; i++) parent[i] = i

  function find(x: number): number {
    while (parent[x] !== x) {
      parent[x] = parent[parent[x]!]!
      x = parent[x]!
    }
    return x
  }

  // Members of each active cluster node
  const members = new Map<number, number[]>()
  for (let i = 0; i < n; i++) members.set(i, [i])

  // Cophenetic distance matrix (flat upper triangle)
  const nPairs = n * (n - 1) / 2
  const cophDist = new Float64Array(nPairs)
  const origFlat = new Float64Array(nPairs)

  // Fill original flat distances
  let idx = 0
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      origFlat[idx] = origDist[i * n + j]!
      idx++
    }
  }

  // Process merges
  for (let s = 0; s < merges.length; s++) {
    const merge = merges[s]!
    const newCluster = n + s
    const rootA = find(merge.a)
    const rootB = find(merge.b)
    const membersA = members.get(rootA) ?? []
    const membersB = members.get(rootB) ?? []

    // Set cophenetic distance for all pairs (one from A, one from B)
    for (const a of membersA) {
      for (const b of membersB) {
        const i = Math.min(a, b)
        const j = Math.max(a, b)
        // Flat index: sum of (n-1) + (n-2) + ... + (n-i) + (j - i - 1)
        const flatIdx = i * n - i * (i + 1) / 2 + (j - i - 1)
        cophDist[flatIdx] = merge.height
      }
    }

    // Union
    parent[rootA] = newCluster
    parent[rootB] = newCluster
    members.set(newCluster, [...membersA, ...membersB])
    members.delete(rootA)
    members.delete(rootB)
  }

  // Pearson correlation between origFlat and cophDist
  return pearsonR(origFlat, cophDist, nPairs)
}

/** Pearson correlation between two Float64Arrays. */
function pearsonR(x: Float64Array, y: Float64Array, n: number): number {
  let sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0
  for (let i = 0; i < n; i++) {
    sx += x[i]!
    sy += y[i]!
    sxx += x[i]! * x[i]!
    syy += y[i]! * y[i]!
    sxy += x[i]! * y[i]!
  }
  const num = n * sxy - sx * sy
  const den = Math.sqrt((n * sxx - sx * sx) * (n * syy - sy * sy))
  return den > 0 ? num / den : 0
}
