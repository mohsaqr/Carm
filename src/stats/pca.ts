/**
 * PCA (Principal Component Analysis) via SVD.
 * Also provides varimax rotation and factor loading computation.
 */

import { Matrix } from '../core/matrix.js'
import type { PCAResult } from '../core/types.js'
import { preprocessData } from './preprocess.js'

// ─── PCA via SVD ──────────────────────────────────────────────────────────

/**
 * PCA via SVD on the standardized data matrix.
 * Equivalent to eigen-decomposition of the correlation matrix.
 *
 * Cross-validated with R:
 * > prcomp(data, scale. = TRUE)
 * > summary(pca)  # check proportion of variance explained
 */
export function runPCA(
  data: readonly (readonly number[])[],
  nComponents?: number,
  scale = true
): PCAResult {
  if (data.length < 2) throw new Error('runPCA: need at least 2 observations')
  const n = data.length
  const k = data[0]!.length
  if (k < 2) throw new Error('runPCA: need at least 2 variables')

  const pp = preprocessData(data, { method: scale ? 'standardize' : 'center' })
  const X = Matrix.fromArray(pp.data as number[][])

  // Compute correlation (or covariance) matrix R = X'X / (n-1)
  const Xt = X.transpose()
  const R = Xt.multiply(X).scale(1 / (n - 1))

  // Eigendecomposition of R — gives eigenvalues and eigenvectors directly
  const { values: eigVals, vectors: eigVecs } = R.eigen()

  const nc = nComponents ?? Math.min(n - 1, k)
  const eigenvalues = eigVals.slice(0, nc)
  // Use trace of R (sum of diagonal) for total variance — numerically exact for standardized data
  let totalVar = 0
  for (let i = 0; i < k; i++) totalVar += R.get(i, i)

  // Loadings = eigenvectors (columns are principal components, rows are variables)
  const loadings: number[][] = Array.from({ length: k }, (_, varIdx) =>
    Array.from({ length: nc }, (_, compIdx) => eigVecs.get(varIdx, compIdx))
  )

  // Scores = X * V (n × nc)
  const Vk = Matrix.fromArray(
    Array.from({ length: k }, (_, i) => Array.from({ length: nc }, (_, j) => eigVecs.get(i, j)))
  )
  const scoresM = X.multiply(Vk)  // n × nc
  const scores: number[][] = Array.from({ length: n }, (_, i) =>
    Array.from({ length: nc }, (_, j) => scoresM.get(i, j))
  )

  const varianceExplained = eigenvalues.map(e => totalVar > 0 ? e / totalVar : 0)
  const cumulativeVariance = varianceExplained.reduce<number[]>((acc, v, i) => {
    acc.push((acc[i - 1] ?? 0) + v)
    return acc
  }, [])

  return {
    loadings,
    scores,
    eigenvalues,
    varianceExplained,
    cumulativeVariance,
    nComponents: nc,
  }
}

// ─── Varimax rotation ─────────────────────────────────────────────────────

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
export function varimaxRotation(
  loadings: readonly (readonly number[])[],
  maxIter = 1000,
  tol = 1e-5,
  normalize = true
): { rotatedLoadings: number[][]; rotationMatrix: number[][] } {
  const p = loadings.length        // variables (rows)
  const nc = loadings[0]!.length   // components (columns)
  if (nc < 2) {
    return {
      rotatedLoadings: loadings.map(row => [...row]),
      rotationMatrix: [[1]],
    }
  }

  // Convert to mutable array
  let X: number[][] = loadings.map(row => [...row])

  // Kaiser normalization: divide each row by its communality (row norm)
  let sc: number[] | null = null
  if (normalize) {
    sc = X.map(row => Math.sqrt(row.reduce((s, v) => s + v * v, 0)))
    X = X.map((row, i) => {
      const h = sc![i]!
      return h > 0 ? row.map(v => v / h) : [...row]
    })
  }

  // SVD-based varimax iteration (matches R's stats::varimax exactly)
  let T: number[][] = Array.from({ length: nc }, (_, i) =>
    Array.from({ length: nc }, (_, j) => (i === j ? 1 : 0))
  )

  let d = 0
  for (let iter = 0; iter < maxIter; iter++) {
    // Z = X * T  (rotated loadings)
    const Z: number[][] = X.map(row => {
      const out = new Array(nc).fill(0) as number[]
      for (let j = 0; j < nc; j++) {
        let s = 0
        for (let l = 0; l < nc; l++) s += row[l]! * T[l]![j]!
        out[j] = s
      }
      return out
    })

    // B = X' * (Z^3 - Z * diag(colSums(Z^2) / p))
    // First compute colSums(Z^2) / p
    const colSumSq = new Array(nc).fill(0) as number[]
    for (let i = 0; i < p; i++) {
      for (let j = 0; j < nc; j++) colSumSq[j]! += Z[i]![j]! * Z[i]![j]!
    }
    for (let j = 0; j < nc; j++) colSumSq[j]! /= p

    // Compute the target: Z^3 - Z * diag(colSumSq)
    const target: number[][] = Z.map(row =>
      row.map((v, j) => v * v * v - v * colSumSq[j]!)
    )

    // B = X' * target (nc × nc)
    const B: number[][] = Array.from({ length: nc }, (_, i) => {
      const out = new Array(nc).fill(0) as number[]
      for (let j = 0; j < nc; j++) {
        let s = 0
        for (let r = 0; r < p; r++) s += X[r]![i]! * target[r]![j]!
        out[j] = s
      }
      return out
    })

    // Polar decomposition of B via eigendecomposition:
    // B = T · P where T is orthogonal, P = (B'B)^{1/2}
    // T = B · (B'B)^{-1/2}
    // Convergence criterion: d = sum(singular values) = trace(P)
    const Bm = Matrix.fromArray(B)
    const BtB = Bm.transpose().multiply(Bm)  // nc × nc, symmetric PSD
    const { values: eigVals, vectors: eigVecs } = BtB.eigen()

    // (B'B)^{-1/2} = Q · diag(1/sqrt(λ)) · Q'
    // Also compute sum(sqrt(λ)) for convergence (= sum of singular values)
    let dNew = 0
    const invSqrt: number[] = eigVals.map(lambda => {
      const sv = Math.sqrt(Math.max(0, lambda))
      dNew += sv
      return sv > 1e-15 ? 1 / sv : 0
    })

    // Compute (B'B)^{-1/2} = Q * diag(invSqrt) * Q'
    const invSqrtBtB: number[][] = Array.from({ length: nc }, (_, i) => {
      const row = new Array(nc).fill(0) as number[]
      for (let j = 0; j < nc; j++) {
        let s = 0
        for (let l = 0; l < nc; l++) s += eigVecs.get(i, l) * invSqrt[l]! * eigVecs.get(j, l)
        row[j] = s
      }
      return row
    })

    // T = B * (B'B)^{-1/2}
    const Tnew: number[][] = Array.from({ length: nc }, (_, i) => {
      const row = new Array(nc).fill(0) as number[]
      for (let j = 0; j < nc; j++) {
        let s = 0
        for (let l = 0; l < nc; l++) s += B[i]![l]! * invSqrtBtB[l]![j]!
        row[j] = s
      }
      return row
    })
    T = Tnew

    const dPast = d
    d = dNew
    if (d < dPast * (1 + tol)) break
  }

  // Final rotation: Z = X * T
  let Z: number[][] = X.map(row => {
    const out = new Array(nc).fill(0) as number[]
    for (let j = 0; j < nc; j++) {
      let s = 0
      for (let l = 0; l < nc; l++) s += row[l]! * T[l]![j]!
      out[j] = s
    }
    return out
  })

  // Kaiser de-normalization: multiply each row back by its communality
  if (normalize && sc) {
    Z = Z.map((row, i) => {
      const h = sc![i]!
      return row.map(v => v * h)
    })
  }

  return { rotatedLoadings: Z, rotationMatrix: T }
}

// ─── Scree data ───────────────────────────────────────────────────────────

export interface ScreeData {
  readonly components: readonly number[]
  readonly eigenvalues: readonly number[]
  readonly varianceExplained: readonly number[]
  readonly cumulativeVariance: readonly number[]
}

/** Extract scree plot data from a PCA result. */
export function screeData(pca: PCAResult): ScreeData {
  return {
    components: Array.from({ length: pca.nComponents }, (_, i) => i + 1),
    eigenvalues: pca.eigenvalues,
    varianceExplained: pca.varianceExplained,
    cumulativeVariance: pca.cumulativeVariance,
  }
}
