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

import { mean as _mean, sd as _sd } from '../core/math.js'

// ─── Types ────────────────────────────────────────────────────────────────

export type PreprocessMethod = 'none' | 'center' | 'standardize' | 'log' | 'sqrt'

export interface PreprocessOptions {
  readonly method?: PreprocessMethod
}

export interface PreprocessResult {
  readonly data: readonly (readonly number[])[]
  readonly colMeans: readonly number[]
  readonly colSDs: readonly number[]
  readonly method: PreprocessMethod
  readonly centered: boolean
  readonly scaled: boolean
}

// ─── Preprocessing ────────────────────────────────────────────────────────

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
export function preprocessData(
  data: readonly (readonly number[])[],
  options?: PreprocessOptions
): PreprocessResult {
  const n = data.length
  if (n === 0) throw new Error('preprocessData: data cannot be empty')
  const d = data[0]!.length

  const method = options?.method ?? 'none'

  // For log/sqrt, transform first then compute means/SDs on transformed data
  if (method === 'log') {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < d; j++) {
        if (data[i]![j]! <= 0) {
          throw new Error(`preprocessData: log requires all values > 0, found ${data[i]![j]} at row ${i}, col ${j}`)
        }
      }
    }
    const transformed = data.map(row => row.map(v => Math.log(v)))
    const colMeans = computeColMeans(transformed, n, d)
    const colSDs = computeColSDs(transformed, colMeans, n, d)
    return { data: transformed, colMeans, colSDs, method, centered: false, scaled: false }
  }

  if (method === 'sqrt') {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < d; j++) {
        if (data[i]![j]! < 0) {
          throw new Error(`preprocessData: sqrt requires all values >= 0, found ${data[i]![j]} at row ${i}, col ${j}`)
        }
      }
    }
    const transformed = data.map(row => row.map(v => Math.sqrt(v)))
    const colMeans = computeColMeans(transformed, n, d)
    const colSDs = computeColSDs(transformed, colMeans, n, d)
    return { data: transformed, colMeans, colSDs, method, centered: false, scaled: false }
  }

  // Compute column means
  const colMeans = computeColMeans(data, n, d)

  // Compute column SDs (Bessel-corrected, matching R scale())
  const colSDs = computeColSDs(data, colMeans, n, d)

  if (method === 'none') {
    return { data, colMeans, colSDs, method, centered: false, scaled: false }
  }

  if (method === 'center') {
    const centered = data.map(row =>
      row.map((v, j) => v - colMeans[j]!)
    )
    return { data: centered, colMeans, colSDs, method, centered: true, scaled: false }
  }

  // method === 'standardize'
  // Zero-variance columns: SD = 1 (centered but not divided)
  const safeSDs = colSDs.map(s => s === 0 ? 1 : s)
  const standardized = data.map(row =>
    row.map((v, j) => (v - colMeans[j]!) / safeSDs[j]!)
  )
  return { data: standardized, colMeans, colSDs: safeSDs, method, centered: true, scaled: true }
}

// ─── Inverse Transform ───────────────────────────────────────────────────

/**
 * Inverse transform preprocessed data back to the original scale.
 *
 * @param data - N × D preprocessed matrix
 * @param params - the PreprocessResult containing transform parameters
 * @returns data in original scale
 */
export function inverseTransform(
  data: readonly (readonly number[])[],
  params: PreprocessResult
): readonly (readonly number[])[] {
  const { method, colMeans, colSDs } = params

  if (method === 'none') return data

  if (method === 'log') {
    return data.map(row => row.map(v => Math.exp(v)))
  }

  if (method === 'sqrt') {
    return data.map(row => row.map(v => v * v))
  }

  if (method === 'center') {
    return data.map(row =>
      row.map((v, j) => v + colMeans[j]!)
    )
  }

  // method === 'standardize'
  return data.map(row =>
    row.map((v, j) => v * colSDs[j]! + colMeans[j]!)
  )
}

// ─── Helpers ─────────────────────────────────────────────────────────────

function computeColMeans(
  data: readonly (readonly number[])[],
  n: number,
  d: number
): number[] {
  return Array.from({ length: d }, (_, j) => {
    let sum = 0
    for (let i = 0; i < n; i++) sum += data[i]![j]!
    return sum / n
  })
}

function computeColSDs(
  data: readonly (readonly number[])[],
  colMeans: readonly number[],
  n: number,
  d: number
): number[] {
  if (n < 2) return new Array(d).fill(0) as number[]
  return Array.from({ length: d }, (_, j) => {
    let ss = 0
    for (let i = 0; i < n; i++) {
      const diff = data[i]![j]! - colMeans[j]!
      ss += diff * diff
    }
    return Math.sqrt(ss / (n - 1))
  })
}
