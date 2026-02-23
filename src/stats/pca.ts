/**
 * PCA (Principal Component Analysis) via SVD.
 * Also provides varimax rotation and factor loading computation.
 */

import { mean as _mean, sd as _sd, roundTo } from '../core/math.js'
import { Matrix } from '../core/matrix.js'
import type { PCAResult } from '../core/types.js'

// ─── Standardize matrix ───────────────────────────────────────────────────

/** Center and optionally scale each column. */
function standardize(data: readonly (readonly number[])[], scale = true): { X: Matrix; colMeans: number[]; colSDs: number[] } {
  const k = data[0]!.length
  const colMeans = Array.from({ length: k }, (_, j) => _mean(data.map(row => row[j] ?? 0)))
  const colSDs = Array.from({ length: k }, (_, j) => _sd(data.map(row => row[j] ?? 0)))

  const X = Matrix.fromArray(
    data.map(row =>
      row.map((v, j) => {
        const centered = v - (colMeans[j] ?? 0)
        const s = colSDs[j] ?? 1
        return scale && s !== 0 ? centered / s : centered
      })
    )
  )
  return { X, colMeans, colSDs }
}

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

  const { X } = standardize(data, scale)

  // Scale by 1/sqrt(n-1) so SVD gives principal components equivalent to eigen(cov)
  const Xs = X.scale(1 / Math.sqrt(n - 1))
  const { U: _U, S, V } = Xs.svd()

  const nc = nComponents ?? Math.min(n - 1, k)
  const eigenvalues = S.slice(0, nc).map(s => s * s)
  const totalVar = S.reduce((sum, s) => sum + s * s, 0)

  // Loadings = V (columns are principal components, rows are variables)
  const loadings: number[][] = Array.from({ length: k }, (_, varIdx) =>
    Array.from({ length: nc }, (_, compIdx) => V.get(varIdx, compIdx))
  )

  // Scores = X * V (n × nc)
  const Vk = Matrix.fromArray(
    Array.from({ length: k }, (_, i) => Array.from({ length: nc }, (_, j) => V.get(i, j)))
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
    eigenvalues: eigenvalues.map(e => roundTo(e, 6)),
    varianceExplained: varianceExplained.map(v => roundTo(v, 6)),
    cumulativeVariance: cumulativeVariance.map(v => roundTo(v, 6)),
    nComponents: nc,
  }
}

// ─── Varimax rotation ─────────────────────────────────────────────────────

/**
 * Varimax rotation of PCA loadings.
 * Maximizes the variance of squared loadings within each factor.
 * Reference: Kaiser (1958), Psychometrika 23:187-200
 *
 * Cross-validated with R:
 * > varimax(pca$rotation[, 1:3])
 */
export function varimaxRotation(
  loadings: readonly (readonly number[])[],
  maxIter = 1000,
  tol = 1e-6
): { rotatedLoadings: number[][]; rotationMatrix: number[][] } {
  const k = loadings.length        // variables
  const m = loadings[0]!.length   // components

  // Convert to 2D mutable array (k × m)
  let L: number[][] = loadings.map(row => [...row])

  // Rotation matrix starts as identity
  let T: number[][] = Array.from({ length: m }, (_, i) =>
    Array.from({ length: m }, (_, j) => (i === j ? 1 : 0))
  )

  for (let iter = 0; iter < maxIter; iter++) {
    let delta = 0
    for (let p = 0; p < m - 1; p++) {
      for (let q = p + 1; q < m; q++) {
        // Compute u and v for the rotation
        const u = L.map(row => (row[p] ?? 0) ** 2 - (row[q] ?? 0) ** 2)
        const v = L.map(row => 2 * (row[p] ?? 0) * (row[q] ?? 0))
        const A = u.reduce((s, ui) => s + ui, 0)
        const B = v.reduce((s, vi) => s + vi, 0)
        const C = u.reduce((s, ui, i) => s + ui ** 2 - (v[i] ?? 0) ** 2, 0)
        const D = u.reduce((s, ui, i) => s + ui * (v[i] ?? 0), 0) * 2

        const X_ = C - (A ** 2 - B ** 2) / k
        const Y_ = D - 2 * A * B / k
        const angle = Math.atan2(Y_, X_) / 4
        if (Math.abs(angle) < 1e-12) continue

        const cos = Math.cos(angle)
        const sin = Math.sin(angle)
        delta += Math.abs(angle)

        // Rotate columns p and q
        const newLp = L.map((row) => (row[p] ?? 0) * cos + (row[q] ?? 0) * sin)
        const newLq = L.map((row) => -(row[p] ?? 0) * sin + (row[q] ?? 0) * cos)
        L.forEach((row, i) => { row[p] = newLp[i]!; row[q] = newLq[i]! })

        // Update rotation matrix
        for (let r = 0; r < m; r++) {
          const tp = (T[r]?.[p] ?? 0) * cos + (T[r]?.[q] ?? 0) * sin
          const tq = -(T[r]?.[p] ?? 0) * sin + (T[r]?.[q] ?? 0) * cos
          T[r]![p] = tp
          T[r]![q] = tq
        }
      }
    }
    if (delta < tol) break
  }

  return { rotatedLoadings: L, rotationMatrix: T }
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
