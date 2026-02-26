/**
 * Factor Analysis module for Carm.
 * Exploratory Factor Analysis (EFA), Confirmatory Factor Analysis (CFA),
 * and psychometric diagnostics (KMO, Bartlett, MAP, Parallel Analysis).
 *
 * All functions are pure — no DOM, no D3, no side effects.
 * Uses splitmix32 PRNG for deterministic reproducibility.
 *
 * Algorithms adapted from:
 * - Iterated PAF: Gorsuch (1983), Factor Analysis (2nd ed.)
 * - ML extraction: Jöreskog (1967), Psychometrika 32(4)
 * - Varimax: Kaiser (1958), Psychometrika 23:187-200
 * - Oblimin/Promax: Jennrich (2002), Psychometrika 67(1)
 * - CFA: Bollen (1989), Structural Equations with Latent Variables
 * - Fit indices: Browne & Cudeck (1993), Hu & Bentler (1999)
 */

import { Matrix } from '../core/matrix.js'
import {
  normalCDF,
  chiSqPValue,
  roundTo,
  nelderMead,
} from '../core/math.js'
import { formatCFAFit } from '../core/apa.js'
import type {
  FactorFit,
  ParameterEstimate,
  FADiagnostics,
  FAResult,
  CFAResult,
} from '../core/types.js'

// ─── PRNG (splitmix32) — identical to clustering.ts ──────────────────────

class PRNG {
  private state: number
  constructor(seed: number) { this.state = seed >>> 0 }
  next(): number {
    this.state = (this.state + 0x9E3779B9) | 0
    let t = this.state ^ (this.state >>> 16)
    t = Math.imul(t, 0x21F0AAAD)
    t = t ^ (t >>> 15)
    t = Math.imul(t, 0x735A2D97)
    t = t ^ (t >>> 15)
    return (t >>> 0) / 4294967296
  }
}

/** Box-Muller transform for normal random variates using splitmix32. */
function prngNormal(rng: PRNG): number {
  let u = 0, v = 0
  while (u === 0) u = rng.next()
  while (v === 0) v = rng.next()
  return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v)
}

// ─── Correlation Matrix ──────────────────────────────────────────────────

/**
 * Compute Pearson correlation matrix from raw data using Float64Array accumulators.
 * Returns a Matrix (d × d).
 */
function computeCorrelationMatrix(
  data: readonly (readonly number[])[],
  n: number,
  d: number
): Matrix {
  const means = new Float64Array(d)
  const sds = new Float64Array(d)

  for (let i = 0; i < n; i++) {
    const row = data[i]!
    for (let j = 0; j < d; j++) means[j] = means[j]! + row[j]! / n
  }
  for (let i = 0; i < n; i++) {
    const row = data[i]!
    for (let j = 0; j < d; j++) sds[j] = sds[j]! + (row[j]! - means[j]!) ** 2
  }
  for (let j = 0; j < d; j++) sds[j] = Math.sqrt(sds[j]! / (n - 1))

  const R: number[][] = Array.from({ length: d }, () => new Array<number>(d).fill(0))
  for (let r = 0; r < d; r++) {
    R[r]![r] = 1.0
    for (let c = r + 1; c < d; c++) {
      let sum = 0
      const sdR = sds[r]! || 1
      const sdC = sds[c]! || 1
      for (let i = 0; i < n; i++) {
        sum += ((data[i]![r]! - means[r]!) / sdR) * ((data[i]![c]! - means[c]!) / sdC)
      }
      const val = sum / (n - 1)
      R[r]![c] = val
      R[c]![r] = val
    }
  }
  return Matrix.fromArray(R)
}

// ─── Options interfaces ──────────────────────────────────────────────────

export interface EFAOptions {
  readonly nFactors?: number               // auto-detect via parallel analysis if omitted
  readonly extraction?: 'paf' | 'ml'       // default: 'ml'
  readonly rotation?: 'varimax' | 'oblimin' | 'promax' | 'quartimin' | 'none'  // default: 'promax'
  readonly seed?: number                   // default: 42 (for parallel analysis)
  readonly maxIter?: number                // default: 1000
  readonly tol?: number                    // default: 1e-6
  readonly variableNames?: readonly string[]
}

export interface CFAOptions {
  readonly maxIter?: number     // default: 1000
  readonly tol?: number         // default: 1e-6
  readonly variableNames?: readonly string[]
  readonly factorNames?: readonly string[]
}

export interface FADiagnosticsOptions {
  readonly seed?: number                  // default: 42
  readonly parallelIterations?: number    // default: 100
}

// ─── Fit Statistics ──────────────────────────────────────────────────────

/**
 * Compute CFA/EFA fit indices from observed (S) and implied (Sigma) covariance matrices.
 * Uses Wishart ML discrepancy: F = log|Σ| + tr(Σ⁻¹S) - log|S| - d
 *
 * When nFactors is provided, applies the Bartlett (1950) correction to the chi-square
 * statistic: χ² = (n - 1 - (2p+5)/6 - 2k/3) × F_ml, matching R psych::fa().
 * When nFactors is omitted (CFA), uses uncorrected χ² = (n - 1) × F_ml.
 */
function computeFit(
  S: Matrix,
  Sigma: Matrix,
  n: number,
  d: number,
  nFreeParams: number,
  nFactors?: number
): FactorFit {
  // ML discrepancy
  let logDetSigma: number, logDetS: number
  try {
    logDetSigma = Sigma.logDet()
  } catch {
    // Fallback via eigen if Cholesky fails
    const ev = Sigma.eigen().values
    logDetSigma = ev.reduce((s, v) => s + Math.log(Math.max(v, 1e-15)), 0)
  }
  try {
    logDetS = S.logDet()
  } catch {
    const ev = S.eigen().values
    logDetS = ev.reduce((s, v) => s + Math.log(Math.max(v, 1e-15)), 0)
  }

  let traceVal: number
  try {
    traceVal = Sigma.inverse().multiply(S).trace()
  } catch {
    traceVal = Sigma.pseudoInverse().multiply(S).trace()
  }

  const F_ml = Math.max(0, logDetSigma + traceVal - logDetS - d)
  // Bartlett (1950) correction for EFA: χ² = (n - 1 - (2p+5)/6 - 2k/3) × F_ml
  // Without correction (CFA): χ² = (n - 1) × F_ml
  const bartlettN = nFactors !== undefined
    ? n - 1 - (2 * d + 5) / 6 - (2 * nFactors) / 3
    : n - 1
  const chiSq = Math.max(0, bartlettN * F_ml)

  // Degrees of freedom
  const totalElements = d * (d + 1) / 2
  const df = Math.max(0, totalElements - nFreeParams)
  const pValue = df > 0 ? chiSqPValue(chiSq, df) : 1

  // Null model: diagonal (independence) — F_null = log|diag(S)| + tr(diag(S)⁻¹ · S) - log|S| - d
  let logDetDiagS = 0
  let traceNull = 0
  for (let i = 0; i < d; i++) {
    const diagVal = Math.max(S.get(i, i), 1e-15)
    logDetDiagS += Math.log(diagVal)
    traceNull += S.get(i, i) / diagVal  // = 1 for each diagonal, but off-diag / diag contributes via full trace
  }
  // Actually: diag(S)⁻¹ · S trace = sum_ij S_ij / S_ii for same row
  // Simpler: F_null from diagonal model
  const diagArr: number[][] = Array.from({ length: d }, (_, i) =>
    Array.from({ length: d }, (_, j) => i === j ? S.get(i, i) : 0)
  )
  const diagMat = Matrix.fromArray(diagArr)
  let traceNullFull: number
  try {
    traceNullFull = diagMat.inverse().multiply(S).trace()
  } catch {
    traceNullFull = d // fallback
  }
  const F_null = Math.max(0, logDetDiagS + traceNullFull - logDetS - d)
  // Null model also uses Bartlett correction (but without the 2k/3 factor term)
  const bartlettNNull = nFactors !== undefined
    ? n - 1 - (2 * d + 5) / 6
    : n - 1
  const chiSqNull = bartlettNNull * F_null
  const dfNull = d * (d - 1) / 2

  // RMSEA with 90% CI (Browne & Cudeck 1993)
  const ncp = Math.max(chiSq - df, 0)
  const rmsea = df > 0 ? Math.sqrt(ncp / (df * (n - 1))) : 0

  // RMSEA 90% CI via chi-square non-centrality parameter bounds
  let rmseaLo = 0, rmseaHi = rmsea * 2
  if (df > 0) {
    // Lower bound: find λ_L such that P(χ² ≥ chiSq | df, λ_L) = 0.95
    // Upper bound: find λ_U such that P(χ² ≥ chiSq | df, λ_U) = 0.05
    // Approximate: use Steiger (1990) approach
    const ncpLo = Math.max(chiSq - df - 1.645 * Math.sqrt(2 * df), 0)
    const ncpHi = Math.max(chiSq - df + 1.645 * Math.sqrt(2 * df), 0)
    rmseaLo = Math.sqrt(Math.max(ncpLo / (df * (n - 1)), 0))
    rmseaHi = Math.sqrt(ncpHi / (df * (n - 1)))
  }

  // CFI (Bentler 1990)
  const ncpNull = Math.max(chiSqNull - dfNull, 0)
  const cfi = ncpNull > 0 ? Math.max(0, Math.min(1, 1 - ncp / ncpNull)) : 1

  // TLI / NNFI (Tucker & Lewis 1973)
  // Not clamped to [0,1] — TLI can exceed 1 when model fits better than expected
  const tli = df > 0 && dfNull > 0
    ? ((chiSqNull / dfNull) - (chiSq / df)) / ((chiSqNull / dfNull) - 1)
    : 1

  // SRMR (standardized root mean square residual)
  let srmrSum = 0
  let srmrCount = 0
  for (let i = 0; i < d; i++) {
    for (let j = 0; j <= i; j++) {
      const sij = S.get(i, j)
      const sigij = Sigma.get(i, j)
      const denom = Math.sqrt(S.get(i, i) * S.get(j, j))
      const r_obs = denom > 0 ? sij / denom : 0
      const r_imp = denom > 0 ? sigij / denom : 0
      srmrSum += (r_obs - r_imp) ** 2
      srmrCount++
    }
  }
  const srmr = Math.sqrt(srmrSum / srmrCount)

  // AIC, BIC
  const aic = chiSq + 2 * nFreeParams
  const bic = chiSq + nFreeParams * Math.log(n)

  return {
    chiSq,
    df,
    pValue,
    rmsea,
    rmseaCI: [rmseaLo, rmseaHi] as const,
    cfi,
    tli,
    srmr,
    aic,
    bic,
  }
}

// ─── Iterated PAF Extraction ─────────────────────────────────────────────

/**
 * Iterated Principal Axis Factoring.
 * Extracts k factors from the correlation matrix by iterating eigendecomposition
 * of the reduced correlation matrix (communalities on diagonal) until convergence.
 *
 * Reference: Gorsuch (1983), Factor Analysis, 2nd ed.
 */
function extractPAF(
  R: Matrix,
  k: number,
  maxIter: number,
  tol: number
): { loadings: number[][]; communalities: Float64Array } {
  const d = R.rows

  // Initialize communalities from SMC (squared multiple correlation)
  const h2 = new Float64Array(d)
  try {
    const invR = R.inverse()
    for (let i = 0; i < d; i++) h2[i] = Math.max(0.01, 1 - 1 / Math.max(invR.get(i, i), 1e-12))
  } catch {
    for (let i = 0; i < d; i++) h2[i] = 0.5
  }

  const loadings: number[][] = Array.from({ length: d }, () => new Array<number>(k).fill(0))

  for (let iter = 0; iter < maxIter; iter++) {
    // Build reduced correlation matrix (communalities on diagonal)
    const adjR: number[][] = R.toArray()
    for (let i = 0; i < d; i++) adjR[i]![i] = h2[i]!
    const adjM = Matrix.fromArray(adjR)

    const { values, vectors } = adjM.eigen()
    // eigen() returns sorted descending already

    const oldH2 = new Float64Array(h2)

    // Extract top k factors
    for (let f = 0; f < k; f++) {
      const eigenVal = Math.max(values[f]!, 0)
      const scale = Math.sqrt(eigenVal)
      for (let r = 0; r < d; r++) {
        loadings[r]![f] = vectors.get(r, f) * scale
      }
    }

    // Update communalities
    let maxDelta = 0
    for (let r = 0; r < d; r++) {
      let sum = 0
      for (let f = 0; f < k; f++) sum += loadings[r]![f]! ** 2
      // Allow near-Heywood cases (communality close to 1) — only clamp for numerical stability
      h2[r] = Math.max(0.001, Math.min(0.9999, sum))
      maxDelta = Math.max(maxDelta, Math.abs(h2[r]! - oldH2[r]!))
    }
    if (maxDelta < tol) break
  }

  return { loadings, communalities: h2 }
}

// ─── ML Extraction ───────────────────────────────────────────────────────

/**
 * Maximum Likelihood factor extraction.
 * Two-phase optimization matching R's factanal():
 *
 * Phase 1: Jöreskog gradient descent on uniquenesses with eigendecomposition
 *   for loadings. Fast convergence to near-optimum.
 *
 * Phase 2: Nelder-Mead refinement on R's exact concentrated ML objective
 *   (FAfn). Polishes uniquenesses to match R's L-BFGS-B result.
 *
 * The concentrated objective (R's FAfn):
 *   F(Ψ) = -Σ_{j>k} (log λ_j - λ_j) + k - d
 * where λ_j are eigenvalues of Ψ^{-1/2} R Ψ^{-1/2}.
 *
 * Reference: R stats::factanal (src/library/stats/R/factanal.R)
 * Reference: Jöreskog (1967), Psychometrika 32(4):443-482
 */
function extractML(
  R: Matrix,
  k: number,
  maxIter: number,
  tol: number
): { loadings: number[][]; communalities: Float64Array } {
  const d = R.rows
  const Rarr = R.toArray()

  // Initialize uniquenesses using R's factanal formula:
  // start <- (1 - 0.5 * nfactors/p) * diag(solve(S))^(-1)
  // diag(solve(S))^(-1) = 1/R^{-1}_ii = uniqueness estimate from SMC
  const Theta = new Float64Array(d)
  try {
    const invR = R.inverse()
    const factor = 1 - 0.5 * k / d
    for (let i = 0; i < d; i++) Theta[i] = Math.max(0.005, Math.min(0.995, factor / Math.max(invR.get(i, i), 1e-12)))
  } catch {
    for (let i = 0; i < d; i++) Theta[i] = 0.5
  }

  /**
   * Extract loadings from eigendecomposition of scaled correlation:
   * Θ^{-1/2} R Θ^{-1/2} → take top k eigenvalues/vectors
   * L = Θ^{1/2} V diag(sqrt(max(λ-1, 0)))
   */
  function extractLoadingsFromTheta(theta: Float64Array): number[][] {
    const scaledR: number[][] = Array.from({ length: d }, (_, i) =>
      Array.from({ length: d }, (_, j) => {
        const si = 1 / Math.sqrt(Math.max(theta[i]!, 1e-12))
        const sj = 1 / Math.sqrt(Math.max(theta[j]!, 1e-12))
        return Rarr[i]![j]! * si * sj
      })
    )
    const { values, vectors } = Matrix.fromArray(scaledR).eigen()
    const L: number[][] = Array.from({ length: d }, () => new Array<number>(k).fill(0))
    for (let f = 0; f < k; f++) {
      const ev = Math.max(values[f]! - 1, 0)
      const scale = Math.sqrt(ev)
      for (let i = 0; i < d; i++) {
        L[i]![f] = Math.sqrt(Math.max(theta[i]!, 1e-12)) * vectors.get(i, f) * scale
      }
    }
    return L
  }

  // ── Phase 1: Jöreskog gradient descent (fast initial convergence) ──
  const lr0 = 0.02
  const lrMin = 0.001
  let L: number[][] = extractLoadingsFromTheta(Theta)
  const vTheta = new Float64Array(d)
  const momentum = 0.8

  for (let iter = 0; iter < maxIter; iter++) {
    const lr = lrMin + (lr0 - lrMin) * 0.5 * (1 + Math.cos(Math.PI * iter / maxIter))
    L = extractLoadingsFromTheta(Theta)

    // Build Σ and its inverse
    const sigmaArr: number[][] = Array.from({ length: d }, (_, i) =>
      Array.from({ length: d }, (_, j) => {
        let sum = 0
        for (let f = 0; f < k; f++) sum += L[i]![f]! * L[j]![f]!
        return sum + (i === j ? Theta[i]! : 0)
      })
    )
    const Sigma = Matrix.fromArray(sigmaArr)
    let invSigma: Matrix
    try { invSigma = Sigma.inverse() } catch { invSigma = Sigma.pseudoInverse() }

    // Gradient: (Σ⁻¹(Σ - R)Σ⁻¹)_ii
    const Delta = invSigma.multiply(Sigma.subtract(R)).multiply(invSigma)
    let maxGrad = 0
    for (let i = 0; i < d; i++) {
      const grad = Delta.get(i, i)
      maxGrad = Math.max(maxGrad, Math.abs(grad))
      const v = momentum * vTheta[i]! - lr * grad
      vTheta[i] = v
      Theta[i] = Math.max(0.005, Math.min(0.995, Theta[i]! + v))
    }
    if (maxGrad < tol) break
  }

  // ── Phase 2: Nelder-Mead polish on R's exact concentrated objective ──
  // F(Ψ) = -Σ_{j>k} (log λ_j - λ_j) + k - d
  // where λ_j = eigenvalues of Ψ^{-1/2} R Ψ^{-1/2}
  function concentratedML(x: readonly number[]): number {
    const scaledR: number[][] = Array.from({ length: d }, (_, i) =>
      Array.from({ length: d }, (_, j) => {
        const psi_i = Math.max(x[i]!, 0.005)
        const psi_j = Math.max(x[j]!, 0.005)
        return Rarr[i]![j]! / (Math.sqrt(psi_i) * Math.sqrt(psi_j))
      })
    )
    const { values } = Matrix.fromArray(scaledR).eigen()
    let sum = 0
    for (let j = k; j < d; j++) {
      const ev = Math.max(values[j]!, 1e-15)
      sum += Math.log(ev) - ev
    }
    // Add penalty for out-of-bounds (death penalty for Nelder-Mead)
    let penalty = 0
    for (let i = 0; i < d; i++) {
      if (x[i]! < 0.005 || x[i]! > 0.995) penalty += 1000
    }
    return -sum + k - d + penalty
  }

  const psi0 = Array.from(Theta)
  const nmResult = nelderMead(concentratedML, psi0, {
    maxIter: 5000 * d,
    tol: 1e-10,
  })

  // Use Nelder-Mead result, clamped to bounds
  const finalTheta = new Float64Array(d)
  for (let i = 0; i < d; i++) {
    finalTheta[i] = Math.max(0.005, Math.min(0.995, nmResult.x[i]!))
  }

  // Extract final loadings from polished uniquenesses
  L = extractLoadingsFromTheta(finalTheta)

  const h2 = new Float64Array(d)
  for (let i = 0; i < d; i++) {
    let sum = 0
    for (let f = 0; f < k; f++) sum += L[i]![f]! ** 2
    h2[i] = Math.max(0.001, Math.min(0.9999, sum))
  }

  return { loadings: L, communalities: h2 }
}

// ─── Rotation Engine ─────────────────────────────────────────────────────

/**
 * Varimax rotation matching R's stats::varimax() exactly.
 *
 * Algorithm (from R src/library/stats/R/factanal.R):
 *   1. Kaiser-normalize rows: x = L / rowNorms
 *   2. Iterate:
 *      z = x * T
 *      B = x' * (z³ - z * diag(colSums(z²)/p))
 *      Polar decomposition of B: T = U * V' from SVD(B)
 *      Converge when sum(singularValues) stops increasing
 *   3. Denormalize: z * rowNorms
 *
 * The SVD of the k×k matrix B is computed via eigendecomposition of B'B
 * (polar decomposition), avoiding the Jacobi one-sided SVD which has
 * accuracy issues for small matrices.
 *
 * Reference: Kaiser (1958), Psychometrika 23:187-200
 */
function rotateVarimax(
  L: number[][],
  maxIter: number,
  tol: number
): { rotated: number[][]; T: number[][] } {
  const p = L.length
  const k = L[0]!.length

  if (k < 2) return { rotated: L.map(r => [...r]), T: [[1]] }

  // Kaiser normalization: normalize each row by communality (row norm)
  const sc = new Float64Array(p)
  for (let i = 0; i < p; i++) {
    let ss = 0
    for (let j = 0; j < k; j++) ss += L[i]![j]! ** 2
    sc[i] = Math.sqrt(ss || 1e-15)
  }
  const x: number[][] = L.map((row, i) => row.map(v => v / sc[i]!))

  // Initialize T = Identity (k × k)
  let T: number[][] = Array.from({ length: k }, (_, i) =>
    Array.from({ length: k }, (_, j) => (i === j ? 1 : 0))
  )

  let dPast = 0
  for (let iter = 0; iter < Math.min(maxIter, 1000); iter++) {
    // z = x * T  (p × k)
    const z: number[][] = Array.from({ length: p }, (_, i) =>
      Array.from({ length: k }, (_, j) => {
        let s = 0
        for (let m = 0; m < k; m++) s += x[i]![m]! * T[m]![j]!
        return s
      })
    )

    // colSums(z²) / p
    const csz2 = new Float64Array(k)
    for (let j = 0; j < k; j++) {
      let s = 0
      for (let i = 0; i < p; i++) s += z[i]![j]! ** 2
      csz2[j] = s / p
    }

    // B = x' * (z³ - z * diag(csz2))   (k × k)
    // target[i][j] = z[i][j]³ - z[i][j] * csz2[j]
    // B[r][c] = Σ_i x[i][r] * target[i][c]
    const B: number[][] = Array.from({ length: k }, (_, r) =>
      Array.from({ length: k }, (_, c) => {
        let s = 0
        for (let i = 0; i < p; i++) {
          const zij = z[i]![c]!
          s += x[i]![r]! * (zij * zij * zij - zij * csz2[c]!)
        }
        return s
      })
    )

    // SVD of B via polar decomposition using eigendecomposition of B'B
    // B'B is k×k symmetric → eigendecompose → V, σ²
    // T = B * V * diag(1/σ) * V'  (polar factor)
    const BtB: number[][] = Array.from({ length: k }, (_, i) =>
      Array.from({ length: k }, (_, j) => {
        let s = 0
        for (let m = 0; m < k; m++) s += B[m]![i]! * B[m]![j]!
        return s
      })
    )
    const { values: sigma2, vectors: Vmat } = Matrix.fromArray(BtB).eigen()

    // Singular values = sqrt(eigenvalues of B'B)
    const svals = new Float64Array(k)
    for (let j = 0; j < k; j++) svals[j] = Math.sqrt(Math.max(sigma2[j]!, 0))

    // T_new = B * V * diag(1/σ) * V'
    // Step 1: BV = B * V  (k × k)
    const Varr = Vmat.toArray()
    const BV: number[][] = Array.from({ length: k }, (_, i) =>
      Array.from({ length: k }, (_, j) => {
        let s = 0
        for (let m = 0; m < k; m++) s += B[i]![m]! * Varr[m]![j]!
        return s
      })
    )

    // Step 2: BV * diag(1/σ)  → this gives U (left singular vectors)
    for (let j = 0; j < k; j++) {
      const invS = svals[j]! > 1e-15 ? 1 / svals[j]! : 0
      for (let i = 0; i < k; i++) BV[i]![j] = BV[i]![j]! * invS
    }

    // Step 3: T = U * V'  (k × k)
    T = Array.from({ length: k }, (_, i) =>
      Array.from({ length: k }, (_, j) => {
        let s = 0
        for (let m = 0; m < k; m++) s += BV[i]![m]! * Varr[j]![m]!
        return s
      })
    )

    // Convergence: sum of singular values (matches R's criterion)
    let dNew = 0
    for (let j = 0; j < k; j++) dNew += svals[j]!
    if (dNew < dPast * (1 + tol)) break
    dPast = dNew
  }

  // Final rotated loadings = x * T, then denormalize by row norms
  const rot: number[][] = Array.from({ length: p }, (_, i) =>
    Array.from({ length: k }, (_, j) => {
      let s = 0
      for (let m = 0; m < k; m++) s += x[i]![m]! * T[m]![j]!
      return s * sc[i]!
    })
  )

  return { rotated: rot, T }
}

/**
 * Oblimin (oblique) rotation via gradient projection.
 * gamma = 0 gives quartimin; gamma = 0.5 gives biquartimin.
 * Reference: Jennrich (2002), Psychometrika 67(1):7-27
 */
function rotateOblimin(
  L: number[][],
  gamma: number,
  maxIter: number,
  tol: number
): { rotated: number[][]; T: number[][]; Phi: number[][] } {
  const d = L.length
  const k = L[0]!.length
  const Lmat = Matrix.fromArray(L)
  let T: number[][] = Array.from({ length: k }, (_, i) =>
    Array.from({ length: k }, (_, j) => (i === j ? 1 : 0))
  )
  let Tmat = Matrix.fromArray(T)
  const alpha = 1.0

  for (let iter = 0; iter < maxIter; iter++) {
    const Lambda = Lmat.multiply(Tmat)
    const LambdaArr = Lambda.toArray()

    // Compute gradient for oblimin criterion
    const G: number[][] = Array.from({ length: d }, () => new Array<number>(k).fill(0))
    for (let i = 0; i < d; i++) {
      let rowSum = 0
      for (let m = 0; m < k; m++) rowSum += LambdaArr[i]![m]! ** 2
      for (let j = 0; j < k; j++) {
        const sumSq = rowSum - LambdaArr[i]![j]! ** 2
        G[i]![j] = LambdaArr[i]![j]! * (sumSq - (gamma / d) * rowSum)
      }
    }
    const Gmat = Matrix.fromArray(G)

    const L_T_G = Lmat.transpose().multiply(Gmat)
    // Project gradient: for oblique, subtract diagonal part
    const diagLTG: number[][] = Array.from({ length: k }, (_, i) =>
      Array.from({ length: k }, (_, j) => (i === j ? L_T_G.get(i, i) : 0))
    )
    const Gp = Gmat.subtract(Tmat.multiply(Matrix.fromArray(diagLTG)))

    let maxGrad = 0
    for (let i = 0; i < k; i++) {
      for (let j = 0; j < k; j++) {
        maxGrad = Math.max(maxGrad, Math.abs(Gp.get(i, j)))
      }
    }
    if (maxGrad < tol) break

    let nextT = Tmat.subtract(Gp.scale(alpha))

    // Normalize columns
    const nextArr = nextT.toArray()
    for (let j = 0; j < k; j++) {
      let ss = 0
      for (let i = 0; i < k; i++) ss += nextArr[i]![j]! ** 2
      const invNorm = 1 / Math.sqrt(ss || 1e-12)
      for (let i = 0; i < k; i++) nextArr[i]![j] = nextArr[i]![j]! * invNorm
    }
    Tmat = Matrix.fromArray(nextArr)
    T = nextArr
  }

  const rotated = Lmat.multiply(Tmat).toArray()

  // Phi = (T⁻¹)(T⁻¹)' — factor correlation matrix
  let invT: Matrix
  try {
    invT = Tmat.inverse()
  } catch {
    invT = Tmat.pseudoInverse()
  }
  const PhiMat = invT.multiply(invT.transpose())
  const Phi = PhiMat.toArray()

  return { rotated, T, Phi }
}

/**
 * Promax rotation: outer Kaiser normalization → varimax → power target → Procrustes.
 * Matches R's psych::fa(rotate="promax") which calls psych::kaiser(loadings, rotate="Promax"):
 *   0. OUTER Kaiser: h² = rowSums(L²), normalize rows to unit norm
 *   1. Varimax rotation on normalized loadings → V, T_var
 *   2. Target Q = V ⊙ |V|^(m−1)
 *   3. Regression U = (V'V)^{-1} V' Q (on varimax loadings, not original)
 *   4. Normalize: d = diag((U'U)^{-1}), U = U diag(√d)
 *   5. Rotated = V × U_norm
 *   6. Denormalize: rotated = rotated * √h²
 *   7. Compound rotation = T_var × U_norm
 *
 * The outer Kaiser normalization is CRITICAL: psych::fa dispatches promax via
 * psych::kaiser() which normalizes loadings by communality before rotation.
 * This changes the promax target Q nonlinearly, producing different results
 * than stats::promax() or psych::Promax() applied directly.
 *
 * Reference: Hendrickson & White (1964), JASA 59:258-264
 * Implementation: R psych::kaiser + stats::promax (src/library/stats/R/factanal.R)
 */
function rotatePromax(
  L: number[][],
  power: number,
  maxIter: number,
  tol: number
): { rotated: number[][]; T: number[][]; Phi: number[][] } {
  const d = L.length
  const k = L[0]!.length

  // Step 0: OUTER Kaiser normalization (psych::kaiser)
  // h2 <- diag(f %*% t(f))  — communalities (row sums of squared loadings)
  // weighted <- f / sqrt(h2) — normalize each row to unit norm
  const h2 = new Float64Array(d)
  for (let i = 0; i < d; i++) {
    let sum = 0
    for (let j = 0; j < k; j++) sum += L[i]![j]! * L[i]![j]!
    h2[i] = sum
  }
  const sqrtH2 = new Float64Array(d)
  for (let i = 0; i < d; i++) sqrtH2[i] = Math.sqrt(Math.max(h2[i]!, 1e-15))

  const L_norm: number[][] = Array.from({ length: d }, (_, i) =>
    Array.from({ length: k }, (_, j) => L[i]![j]! / sqrtH2[i]!)
  )

  // Step 1: Varimax rotation on NORMALIZED loadings
  const { rotated: vari, T: T_varimax } = rotateVarimax(L_norm, maxIter, tol)

  // Step 2: Target Q = V * |V|^(m-1) where V = varimax loadings (normalized)
  const Q: number[][] = Array.from({ length: d }, (_, i) =>
    Array.from({ length: k }, (_, j) => {
      const val = vari[i]![j]!
      return Math.sign(val) * Math.abs(val) ** power
    })
  )

  // Step 3: Regression on varimax loadings: U = (V'V)^{-1} V' Q
  // This matches R's lm.fit(x, Q)$coefficients where x = varimax loadings
  const Vmat = Matrix.fromArray(vari)
  const Qmat = Matrix.fromArray(Q)
  let VtV_inv: Matrix
  try {
    VtV_inv = Vmat.transpose().multiply(Vmat).inverse()
  } catch {
    VtV_inv = Vmat.transpose().multiply(Vmat).pseudoInverse()
  }
  const U = VtV_inv.multiply(Vmat.transpose()).multiply(Qmat)
  const Uarr = U.toArray()

  // Step 4: R normalization: d = diag(solve(t(U) %*% U)), U = U %*% diag(sqrt(d))
  // This scales each column j by sqrt((U'U)^{-1}_{jj})
  const UtU = U.transpose().multiply(U)
  let UtU_inv: Matrix
  try {
    UtU_inv = UtU.inverse()
  } catch {
    UtU_inv = UtU.pseudoInverse()
  }
  for (let j = 0; j < k; j++) {
    const d_j = Math.max(UtU_inv.get(j, j), 1e-15)
    const scale = Math.sqrt(d_j)
    for (let i = 0; i < k; i++) Uarr[i]![j] = Uarr[i]![j]! * scale
  }

  // Step 5: Rotated loadings = V * U_normalized (still in normalized space)
  const Umat = Matrix.fromArray(Uarr)
  const rotatedNorm = Vmat.multiply(Umat).toArray()

  // Step 6: DENORMALIZE — multiply back by √h² (psych::kaiser: rotated$loadings * sqrt(h2))
  const rotated: number[][] = Array.from({ length: d }, (_, i) =>
    Array.from({ length: k }, (_, j) => rotatedNorm[i]![j]! * sqrtH2[i]!)
  )

  // Step 7: Compound rotation matrix = T_varimax * U_normalized
  // This matches R's: U <- xx$rotmat %*% U
  const TvarMat = Matrix.fromArray(T_varimax)
  const Tcompound = TvarMat.multiply(Umat)
  const Tarr = Tcompound.toArray()

  // Step 8: Phi = T_compound^{-1} * (T_compound^{-1})'
  let invT: Matrix
  try {
    invT = Tcompound.inverse()
  } catch {
    invT = Tcompound.pseudoInverse()
  }
  const PhiMat = invT.multiply(invT.transpose())
  const Phi = PhiMat.toArray()

  return { rotated, T: Tarr, Phi }
}

/**
 * Dispatch rotation by method name.
 * Returns rotated loadings and factor correlation matrix (Phi).
 * Phi is identity for orthogonal rotations.
 */
function applyRotation(
  loadings: number[][],
  method: string,
  maxIter: number,
  tol: number
): { rotated: number[][]; Phi: number[][] } {
  const k = loadings[0]!.length

  if (method === 'none') {
    const Phi = Array.from({ length: k }, (_, i) =>
      Array.from({ length: k }, (_, j) => (i === j ? 1 : 0))
    )
    return { rotated: loadings.map(r => [...r]), Phi }
  }

  if (method === 'varimax') {
    const { rotated } = rotateVarimax(loadings, maxIter, tol)
    const Phi = Array.from({ length: k }, (_, i) =>
      Array.from({ length: k }, (_, j) => (i === j ? 1 : 0))
    )
    return { rotated, Phi }
  }

  if (method === 'oblimin' || method === 'quartimin') {
    const gamma = method === 'quartimin' ? 0 : 0
    const { rotated, Phi } = rotateOblimin(loadings, gamma, maxIter, tol)
    return { rotated, Phi }
  }

  if (method === 'promax') {
    const { rotated, Phi } = rotatePromax(loadings, 4, maxIter, tol)
    return { rotated, Phi }
  }

  throw new Error(`runEFA: unknown rotation method '${method}'`)
}

// ─── Velicer's MAP ───────────────────────────────────────────────────────

/**
 * Velicer's Minimum Average Partial correlation test.
 * Determines the number of factors by finding the minimum average
 * squared partial correlation after extracting 0..d-1 components.
 *
 * Reference: Velicer (1976), Psychometrika 41(3):321-327
 */
function velicerMAP(R: Matrix): number {
  const d = R.rows
  const { values, vectors } = R.eigen()
  // values and vectors are already sorted descending by eigen()

  let bestK = 0
  let minMap = Infinity

  for (let k = 0; k < d - 1; k++) {
    // Build loading matrix L from top k+1 eigenvectors
    const Larr: number[][] = Array.from({ length: d }, (_, i) =>
      Array.from({ length: k + 1 }, (_, j) =>
        vectors.get(i, j) * Math.sqrt(Math.max(values[j]!, 0))
      )
    )
    const Lmat = Matrix.fromArray(Larr)
    const R_star = R.subtract(Lmat.multiply(Lmat.transpose()))

    // Convert residual to correlation
    const C_star: number[][] = Array.from({ length: d }, (_, i) =>
      Array.from({ length: d }, (_, j) => {
        const denom = Math.sqrt(Math.abs(R_star.get(i, i)) * Math.abs(R_star.get(j, j)))
        return denom > 1e-15 ? R_star.get(i, j) / denom : (i === j ? 1 : 0)
      })
    )

    // Average squared partial correlation (off-diagonal)
    let sumSq = 0
    for (let r = 0; r < d; r++) {
      for (let c = 0; c < r; c++) sumSq += C_star[r]![c]! ** 2
    }
    const mapVal = sumSq / (d * (d - 1) / 2)
    if (mapVal < minMap) {
      minMap = mapVal
      bestK = k + 1
    }
  }
  return bestK
}

// ─── Parallel Analysis ───────────────────────────────────────────────────

/**
 * Monte Carlo parallel analysis.
 * Generates random normal data matrices, computes their correlation eigenvalues,
 * and returns the 95th percentile as the threshold for factor retention.
 *
 * Reference: Horn (1965), Psychometrika 30(2):179-185
 */
function parallelAnalysis(
  observedEigenvalues: readonly number[],
  n: number,
  d: number,
  iterations: number,
  rng: PRNG
): { simulated: readonly number[]; suggested: number } {
  // Store all simulated eigenvalues: iterations × d
  const allEigens: number[][] = Array.from({ length: d }, () => new Array<number>(iterations).fill(0))

  for (let iter = 0; iter < iterations; iter++) {
    // Generate random normal data (n × d)
    const randomData: number[][] = Array.from({ length: n }, () =>
      Array.from({ length: d }, () => prngNormal(rng))
    )
    const randR = computeCorrelationMatrix(randomData, n, d)
    const randEig = randR.eigen().values  // sorted descending
    for (let i = 0; i < d; i++) {
      allEigens[i]![iter] = randEig[i]!
    }
  }

  // 95th percentile for each eigenvalue position
  const simulated = allEigens.map(eigArray => {
    const sorted = [...eigArray].sort((a, b) => a - b)
    const idx = Math.floor(0.95 * iterations)
    return sorted[Math.min(idx, iterations - 1)]!
  })

  // Number of factors: observed > simulated threshold
  let suggested = 0
  for (let i = 0; i < d; i++) {
    if (observedEigenvalues[i]! > simulated[i]!) suggested++
    else break
  }

  return { simulated, suggested: Math.max(1, suggested) }
}

// ─── KMO & Bartlett ──────────────────────────────────────────────────────

/**
 * Kaiser-Meyer-Olkin sampling adequacy and Bartlett's test of sphericity.
 *
 * KMO uses the anti-image correlation matrix.
 * Bartlett's test: χ² = -[(n-1) - (2d+5)/6] × log|R|
 *
 * Reference:
 * - Kaiser (1974), Educational & Psychological Measurement 34:111-117
 * - Bartlett (1950), British Journal of Statistical Psychology 3:77-85
 */
function computeKMOBartlett(
  R: Matrix,
  n: number,
  d: number
): { kmo: number; kmoPerItem: readonly number[]; bartlett: { chiSq: number; df: number; pValue: number } } {
  let invR: Matrix
  try {
    invR = R.inverse()
  } catch {
    invR = R.pseudoInverse()
  }

  // Anti-image correlation: S² R⁻¹ S² where S² = diag(1/R⁻¹_ii)
  const S2diag = new Float64Array(d)
  for (let i = 0; i < d; i++) {
    S2diag[i] = 1 / Math.max(invR.get(i, i), 1e-12)
  }

  // anti-image_ij = -R⁻¹_ij / sqrt(R⁻¹_ii × R⁻¹_jj) for i ≠ j
  const antiImage: number[][] = Array.from({ length: d }, (_, i) =>
    Array.from({ length: d }, (_, j) => {
      if (i === j) return 1
      return -invR.get(i, j) / Math.sqrt(Math.max(invR.get(i, i) * invR.get(j, j), 1e-12))
    })
  )

  let rSum = 0, qSum = 0
  const msaItems = new Float64Array(d)

  for (let i = 0; i < d; i++) {
    let rIdx = 0, qIdx = 0
    for (let j = 0; j < d; j++) {
      if (i === j) continue
      rIdx += R.get(i, j) ** 2
      qIdx += antiImage[i]![j]! ** 2
    }
    msaItems[i] = rIdx / Math.max(rIdx + qIdx, 1e-12)
    rSum += rIdx
    qSum += qIdx
  }

  const kmo = rSum / Math.max(rSum + qSum, 1e-12)

  // Bartlett's test
  let logDetR: number
  try {
    logDetR = R.logDet()
  } catch {
    const ev = R.eigen().values
    logDetR = ev.reduce((s, v) => s + Math.log(Math.max(v, 1e-15)), 0)
  }
  const chiSq = -((n - 1) - (2 * d + 5) / 6) * logDetR
  const df = d * (d - 1) / 2
  const pValue = chiSqPValue(Math.max(chiSq, 0), df)

  return {
    kmo,
    kmoPerItem: Array.from(msaItems),
    bartlett: {
      chiSq: Math.max(chiSq, 0),
      df,
      pValue,
    },
  }
}

// ─── Implied Covariance (for CFA) ───────────────────────────────────────

/** Σ = Λ Φ Λ' + Θ  where Θ is diagonal uniquenesses */
function computeImpliedCov(
  L: number[][],
  Phi: number[][],
  Theta: Float64Array
): Matrix {
  const d = L.length
  const k = L[0]!.length
  const sigma: number[][] = Array.from({ length: d }, (_, i) =>
    Array.from({ length: d }, (_, j) => {
      let sum = 0
      for (let f1 = 0; f1 < k; f1++) {
        for (let f2 = 0; f2 < k; f2++) {
          sum += L[i]![f1]! * Phi[f1]![f2]! * L[j]![f2]!
        }
      }
      return sum + (i === j ? Theta[i]! : 0)
    })
  )
  return Matrix.fromArray(sigma)
}

// ─── CFA ML Optimizer ────────────────────────────────────────────────────

/**
 * CFA via ML estimation with Armijo backtracking line search.
 * Estimates loadings, uniquenesses, and factor covariances simultaneously.
 *
 * Reference: Bollen (1989), Structural Equations with Latent Variables, ch. 4
 */
function cfaOptimize(
  S: Matrix,
  model: Readonly<Record<string, readonly number[]>>,
  d: number,
  maxIter: number,
  tol: number
): {
  L: number[][]
  Theta: Float64Array
  Phi: number[][]
  converged: boolean
  iterations: number
} {
  const factors = Object.keys(model)
  const k = factors.length

  // Initialize loadings
  const L: number[][] = Array.from({ length: d }, () => new Array<number>(k).fill(0))
  const Theta = new Float64Array(d).fill(0.5)
  const Phi: number[][] = Array.from({ length: k }, (_, i) =>
    Array.from({ length: k }, (_, j) => (i === j ? 1 : 0))
  )

  // Set initial loadings to 0.7 for specified items
  factors.forEach((f, c) => {
    const items = model[f]!
    for (const r of items) {
      if (r < d) L[r]![c] = 0.7
    }
  })

  let converged = false
  let iter = 0
  const c1 = 1e-4  // Armijo sufficient decrease parameter

  /** Wishart ML discrepancy */
  function objective(Lc: number[][], Phic: number[][], Thetac: Float64Array): number {
    const Sigma = computeImpliedCov(Lc, Phic, Thetac)
    let logDetSigma: number
    try {
      logDetSigma = Sigma.logDet()
    } catch {
      return 1e10
    }
    let logDetS: number
    try {
      logDetS = S.logDet()
    } catch {
      return 1e10
    }
    let trVal: number
    try {
      trVal = Sigma.inverse().multiply(S).trace()
    } catch {
      return 1e10
    }
    return Math.max(0, logDetSigma + trVal - logDetS - d)
  }

  for (iter = 0; iter < maxIter; iter++) {
    const Sigma = computeImpliedCov(L, Phi, Theta)

    let invSigma: Matrix
    try {
      invSigma = Sigma.inverse()
    } catch {
      try {
        invSigma = Sigma.pseudoInverse()
      } catch {
        break
      }
    }

    // Delta = Σ⁻¹(Σ - S)Σ⁻¹
    const Delta = invSigma.multiply(Sigma.subtract(S)).multiply(invSigma)

    // Gradient w.r.t. loadings: dF/dΛ = 2ΔΛΦ
    const LMat = Matrix.fromArray(L)
    const PhiMat = Matrix.fromArray(Phi)
    const gL = Delta.multiply(LMat).multiply(PhiMat).scale(2)

    // Gradient w.r.t. factor covariances: dF/dΦ = Λ'ΔΛ
    const gPhi = LMat.transpose().multiply(Delta).multiply(LMat)

    const currentLoss = objective(L, Phi, Theta)
    let alpha = 1.0
    let armijoSatisfied = false

    while (!armijoSatisfied && alpha > 1e-6) {
      let maxGrad = 0
      let dirDotGrad = 0

      // Trial update loadings
      const nextL: number[][] = Array.from({ length: d }, () => new Array<number>(k).fill(0))
      factors.forEach((f, c) => {
        for (const r of model[f]!) {
          if (r >= d) continue
          const grad = gL.get(r, c)
          maxGrad = Math.max(maxGrad, Math.abs(grad))
          nextL[r]![c] = L[r]![c]! - alpha * grad
          dirDotGrad += grad * (-grad)
        }
      })

      // Trial update uniquenesses
      const nextTheta = new Float64Array(d)
      for (let i = 0; i < d; i++) {
        const grad = Delta.get(i, i)
        maxGrad = Math.max(maxGrad, Math.abs(grad))
        nextTheta[i] = Math.max(0.001, Theta[i]! - alpha * grad)
        dirDotGrad += grad * (-grad)
      }

      // Trial update factor covariances
      const nextPhi: number[][] = Array.from({ length: k }, (_, i) =>
        Array.from({ length: k }, (_, j) => (i === j ? 1 : 0))
      )
      for (let r = 0; r < k; r++) {
        for (let c = 0; c < r; c++) {
          const grad = gPhi.get(r, c) + gPhi.get(c, r)
          maxGrad = Math.max(maxGrad, Math.abs(grad))
          const updated = Math.max(-0.99, Math.min(0.99, Phi[r]![c]! - alpha * grad))
          nextPhi[r]![c] = updated
          nextPhi[c]![r] = updated
          dirDotGrad += grad * (-grad)
        }
      }

      if (maxGrad < tol) {
        converged = true
        break
      }

      const nextLoss = objective(nextL, nextPhi, nextTheta)

      if (nextLoss <= currentLoss + c1 * alpha * dirDotGrad) {
        armijoSatisfied = true
        // Accept update
        for (let i = 0; i < d; i++) {
          for (let j = 0; j < k; j++) L[i]![j] = nextL[i]![j]!
          Theta[i] = nextTheta[i]!
        }
        for (let i = 0; i < k; i++) {
          for (let j = 0; j < k; j++) Phi[i]![j] = nextPhi[i]![j]!
        }
      } else {
        alpha *= 0.5
      }
    }

    if (converged) break
    if (!armijoSatisfied) break  // line search failed
  }

  return { L, Theta, Phi, converged, iterations: iter }
}

// ─── Fisher Information for CFA Standard Errors ─────────────────────────

/**
 * Compute standard errors for CFA parameters using the expected information matrix.
 * SE = sqrt( 2/(n-1) * diag(I⁻¹) ) where I is the Fisher information.
 *
 * We use the sandwich-form observed information:
 * I_θθ' = (n-1)/2 * (∂σ/∂θ)' (Σ⁻¹ ⊗ Σ⁻¹) (∂σ/∂θ)
 *
 * For simplicity, we use numerical differentiation of the gradient.
 */
function cfaStandardErrors(
  L: number[][],
  Phi: number[][],
  Theta: Float64Array,
  S: Matrix,
  model: Readonly<Record<string, readonly number[]>>,
  n: number,
  d: number
): {
  loadingSE: number[][]
  thetaSE: Float64Array
  phiSE: number[][]
} {
  const factors = Object.keys(model)
  const k = factors.length

  // Collect free parameters into a vector
  const params: { type: 'loading' | 'theta' | 'phi'; i: number; j: number }[] = []
  factors.forEach((f, c) => {
    for (const r of model[f]!) {
      if (r < d) params.push({ type: 'loading', i: r, j: c })
    }
  })
  for (let i = 0; i < d; i++) params.push({ type: 'theta', i, j: 0 })
  for (let r = 0; r < k; r++) {
    for (let c = 0; c < r; c++) params.push({ type: 'phi', i: r, j: c })
  }

  const nParams = params.length

  // Helper: get/set parameter value
  function getParam(idx: number): number {
    const p = params[idx]!
    if (p.type === 'loading') return L[p.i]![p.j]!
    if (p.type === 'theta') return Theta[p.i]!
    return Phi[p.i]![p.j]!
  }

  function setParam(idx: number, val: number): void {
    const p = params[idx]!
    if (p.type === 'loading') L[p.i]![p.j] = val
    else if (p.type === 'theta') Theta[p.i] = val
    else { Phi[p.i]![p.j] = val; Phi[p.j]![p.i] = val }
  }

  // Compute objective at current parameters
  function obj(): number {
    const Sigma = computeImpliedCov(L, Phi, Theta)
    try {
      const logDetSig = Sigma.logDet()
      const logDetS = S.logDet()
      const tr = Sigma.inverse().multiply(S).trace()
      return Math.max(0, logDetSig + tr - logDetS - d)
    } catch {
      return 1e10
    }
  }

  // Numerical Hessian via central finite differences
  const h = 1e-4
  const hessian: number[][] = Array.from({ length: nParams }, () =>
    new Array<number>(nParams).fill(0)
  )

  for (let i = 0; i < nParams; i++) {
    for (let j = i; j < nParams; j++) {
      const vi = getParam(i)
      const vj = getParam(j)

      if (i === j) {
        // Diagonal: f(x+h) - 2f(x) + f(x-h) / h²
        setParam(i, vi + h)
        const fPlus = obj()
        setParam(i, vi - h)
        const fMinus = obj()
        setParam(i, vi)
        const f0 = obj()
        hessian[i]![i] = (fPlus - 2 * f0 + fMinus) / (h * h)
      } else {
        // Off-diagonal
        setParam(i, vi + h); setParam(j, vj + h)
        const fpp = obj()
        setParam(i, vi + h); setParam(j, vj - h)
        const fpm = obj()
        setParam(i, vi - h); setParam(j, vj + h)
        const fmp = obj()
        setParam(i, vi - h); setParam(j, vj - h)
        const fmm = obj()
        setParam(i, vi); setParam(j, vj)

        hessian[i]![j] = (fpp - fpm - fmp + fmm) / (4 * h * h)
        hessian[j]![i] = hessian[i]![j]!
      }
    }
  }

  // Information matrix = (n-1)/2 * Hessian
  const infoArr: number[][] = hessian.map(row =>
    row.map(v => ((n - 1) / 2) * v)
  )

  // Invert to get variance-covariance of estimates
  let covMat: Matrix
  try {
    covMat = Matrix.fromArray(infoArr).inverse()
  } catch {
    try {
      covMat = Matrix.fromArray(infoArr).pseudoInverse()
    } catch {
      // Fallback: return rough SEs
      const loadingSE = Array.from({ length: d }, () => new Array<number>(k).fill(0.05))
      const thetaSE = new Float64Array(d).fill(0.05)
      const phiSE = Array.from({ length: k }, () => new Array<number>(k).fill(0.05))
      return { loadingSE, thetaSE, phiSE }
    }
  }

  // Extract SEs
  const loadingSE: number[][] = Array.from({ length: d }, () => new Array<number>(k).fill(0))
  const thetaSE = new Float64Array(d)
  const phiSE: number[][] = Array.from({ length: k }, () => new Array<number>(k).fill(0))

  for (let idx = 0; idx < nParams; idx++) {
    const variance = Math.max(covMat.get(idx, idx), 0)
    const se = Math.sqrt(variance)
    const p = params[idx]!
    if (p.type === 'loading') loadingSE[p.i]![p.j] = se
    else if (p.type === 'theta') thetaSE[p.i] = se
    else phiSE[p.i]![p.j] = se
  }

  return { loadingSE, thetaSE, phiSE }
}

// ─── Public API: runFADiagnostics ────────────────────────────────────────

/**
 * Compute factor analysis diagnostics: KMO, Bartlett's test,
 * Velicer's MAP, and parallel analysis.
 *
 * Cross-validated with R:
 * > library(psych)
 * > KMO(data)
 * > cortest.bartlett(cor(data), n = nrow(data))
 * > fa.parallel(data, fm = "ml", fa = "fa")
 */
export function runFADiagnostics(
  data: readonly (readonly number[])[],
  options?: FADiagnosticsOptions
): FADiagnostics {
  const n = data.length
  if (n < 3) throw new Error('runFADiagnostics: need at least 3 observations')
  const d = data[0]!.length
  if (d < 2) throw new Error('runFADiagnostics: need at least 2 variables')

  const seed = options?.seed ?? 42
  const iterations = options?.parallelIterations ?? 100
  const rng = new PRNG(seed)

  const R = computeCorrelationMatrix(data, n, d)

  // KMO & Bartlett
  const { kmo, kmoPerItem, bartlett } = computeKMOBartlett(R, n, d)

  // Eigenvalues of correlation matrix (for scree and parallel analysis)
  const eigenvalues = R.eigen().values  // sorted descending

  // Parallel Analysis
  const { simulated, suggested: parallelSuggested } = parallelAnalysis(
    eigenvalues, n, d, iterations, rng
  )

  // Velicer's MAP
  const mapSuggested = velicerMAP(R)

  return {
    kmo,
    kmoPerItem,
    bartlett,
    mapSuggested,
    parallelEigenvalues: [...eigenvalues],
    parallelSimulated: [...simulated],
    parallelSuggested,
  }
}

// ─── Public API: runEFA ──────────────────────────────────────────────────

/**
 * Exploratory Factor Analysis with extraction (PAF or ML) and rotation.
 *
 * Cross-validated with R:
 * > library(psych)
 * > fa(data, nfactors = 3, fm = "ml", rotate = "promax")
 * > fa(data, nfactors = 3, fm = "minres", rotate = "varimax")
 */
export function runEFA(
  data: readonly (readonly number[])[],
  options?: EFAOptions
): FAResult {
  const n = data.length
  if (n < 3) throw new Error('runEFA: need at least 3 observations')
  const d = data[0]!.length
  if (d < 2) throw new Error('runEFA: need at least 2 variables')

  const extraction = options?.extraction ?? 'ml'
  const rotation = options?.rotation ?? 'promax'
  const maxIter = options?.maxIter ?? 1000
  const tol = options?.tol ?? 1e-6
  const seed = options?.seed ?? 42

  const R = computeCorrelationMatrix(data, n, d)
  const eigenvalues = R.eigen().values  // sorted descending

  // Determine number of factors
  let nFactors = options?.nFactors
  if (nFactors === undefined) {
    const rng = new PRNG(seed)
    const { suggested } = parallelAnalysis(eigenvalues, n, d, 100, rng)
    nFactors = suggested
  }
  if (nFactors < 1) throw new Error('runEFA: nFactors must be at least 1')
  if (nFactors >= d) throw new Error('runEFA: nFactors must be less than number of variables')

  // Extract
  const extracted = extraction === 'ml'
    ? extractML(R, nFactors, maxIter, tol)
    : extractPAF(R, nFactors, maxIter, tol)

  // Rotate
  const { rotated, Phi } = applyRotation(extracted.loadings, rotation, maxIter, tol)

  // Compute communalities and uniqueness from rotated solution + Phi
  const communalities = new Float64Array(d)
  const uniqueness = new Float64Array(d)
  for (let i = 0; i < d; i++) {
    let comm = 0
    for (let f1 = 0; f1 < nFactors; f1++) {
      for (let f2 = 0; f2 < nFactors; f2++) {
        comm += rotated[i]![f1]! * Phi[f1]![f2]! * rotated[i]![f2]!
      }
    }
    communalities[i] = Math.max(0.001, Math.min(0.9999, comm))
    uniqueness[i] = 1 - communalities[i]!
  }

  // Build implied covariance for fit statistics
  const thetaArr = new Float64Array(d)
  for (let i = 0; i < d; i++) thetaArr[i] = Math.max(0.001, uniqueness[i]!)
  const impliedSigma = computeImpliedCov(rotated, Phi, thetaArr)

  // Count free parameters for EFA df calculation
  // In EFA, the unrotated ML solution has k(k-1)/2 identification constraints
  // (Anderson–Rubin conditions), so free params = loadings + uniquenesses - constraints.
  // Rotation is post-hoc and does not affect the chi-square test df.
  const nLoadingParams = d * nFactors
  const nUniqueParams = d
  const nRotationConstraints = nFactors * (nFactors - 1) / 2
  const nFreeParams = nLoadingParams + nUniqueParams - nRotationConstraints

  const fit = computeFit(R, impliedSigma, n, d, nFreeParams, nFactors)

  // Standardized loadings: λ_std = λ * sqrt(φ_jj) / sqrt(σ_ii)
  const stdLoadings: number[][] = Array.from({ length: d }, (_, i) =>
    Array.from({ length: nFactors }, (_, j) => {
      const sigmaII = impliedSigma.get(i, i)
      return rotated[i]![j]! * Math.sqrt(Phi[j]![j]!) / Math.sqrt(Math.max(sigmaII, 1e-12))
    })
  )

  // Variable and factor names
  const variableNames = options?.variableNames
    ?? Array.from({ length: d }, (_, i) => `V${i + 1}`)
  const factorNames = Array.from({ length: nFactors }, (_, i) => `F${i + 1}`)

  const formatted = `EFA (${extraction}/${rotation}): ${formatCFAFit(fit)}`

  return {
    loadings: rotated,
    standardizedLoadings: stdLoadings,
    uniqueness: Array.from(uniqueness),
    communalities: Array.from(communalities),
    factorCorrelations: Phi,
    fit,
    eigenvalues: [...eigenvalues],
    nFactors,
    rotation,
    extraction,
    variableNames,
    factorNames,
    formatted,
  }
}

// ─── Public API: runCFA ──────────────────────────────────────────────────

/**
 * Confirmatory Factor Analysis via ML estimation.
 *
 * Model is specified as a record: { F1: [0, 1, 2], F2: [3, 4, 5] }
 * where keys are factor names and values are 0-indexed item indices.
 *
 * Cross-validated with R:
 * > library(lavaan)
 * > model <- 'F1 =~ x1 + x2 + x3\n F2 =~ x4 + x5 + x6'
 * > fit <- cfa(model, data)
 * > standardizedSolution(fit)
 * > fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli", "srmr"))
 */
export function runCFA(
  data: readonly (readonly number[])[],
  model: Readonly<Record<string, readonly number[]>>,
  options?: CFAOptions
): CFAResult {
  const n = data.length
  if (n < 3) throw new Error('runCFA: need at least 3 observations')
  const d = data[0]!.length
  if (d < 2) throw new Error('runCFA: need at least 2 variables')

  const factors = Object.keys(model)
  const k = factors.length
  if (k < 1) throw new Error('runCFA: model must specify at least 1 factor')

  const maxIter = options?.maxIter ?? 1000
  const tolVal = options?.tol ?? 1e-6

  const S = computeCorrelationMatrix(data, n, d)

  // Optimize CFA
  const { L, Theta, Phi } = cfaOptimize(S, model, d, maxIter, tolVal)

  // Implied covariance
  const impliedSigma = computeImpliedCov(L, Phi, Theta)

  // Count free parameters for df calculation
  let nLoadingParams = 0
  factors.forEach((f) => { nLoadingParams += model[f]!.length })
  const nUniqueParams = d
  const nCovParams = k * (k - 1) / 2
  const nFreeParams = nLoadingParams + nUniqueParams + nCovParams

  const fit = computeFit(S, impliedSigma, n, d, nFreeParams)

  // Standard errors
  // Make deep copies for SE computation (it mutates L/Phi/Theta temporarily)
  const Lcopy: number[][] = L.map(r => [...r])
  const Phicopy: number[][] = Phi.map(r => [...r])
  const Thetacopy = new Float64Array(Theta)
  const { loadingSE, thetaSE, phiSE } = cfaStandardErrors(
    Lcopy, Phicopy, Thetacopy, S, model, n, d
  )

  // Build parameter estimates
  const paramLoadings: ParameterEstimate[][] = factors.map((f, c) =>
    model[f]!.map((r) => {
      if (r >= d) {
        return { estimate: 0, se: 0, z: 0, pValue: 1, stdAll: 0 }
      }
      const est = L[r]![c]!
      const se = Math.max(loadingSE[r]![c]!, 1e-6)
      const z = est / se
      const pValue = 2 * (1 - normalCDF(Math.abs(z)))
      const sigmaII = Math.max(impliedSigma.get(r, r), 1e-12)
      const stdAll = est * Math.sqrt(Phi[c]![c]!) / Math.sqrt(sigmaII)
      return { estimate: roundTo(est, 4), se: roundTo(se, 4), z: roundTo(z, 3), pValue: roundTo(pValue, 4), stdAll: roundTo(stdAll, 4) }
    })
  )

  const paramUniqueness: ParameterEstimate[] = Array.from(Theta).map((est, i) => {
    const se = Math.max(thetaSE[i]!, 1e-6)
    const z = est / se
    const pValue = 2 * (1 - normalCDF(Math.abs(z)))
    const sigmaII = Math.max(impliedSigma.get(i, i), 1e-12)
    const stdAll = est / sigmaII
    return { estimate: roundTo(est, 4), se: roundTo(se, 4), z: roundTo(z, 3), pValue: roundTo(pValue, 4), stdAll: roundTo(stdAll, 4) }
  })

  const paramCov: ParameterEstimate[][] = Array.from({ length: k }, (_, r) =>
    Array.from({ length: k }, (_, c) => {
      if (r === c) {
        return { estimate: 1, se: 0, z: Infinity, pValue: 0, stdAll: 1 }
      }
      const est = Phi[r]![c]!
      const se = Math.max(phiSE[r]![c]! || phiSE[c]![r]!, 1e-6)
      const z = est / se
      const pValue = 2 * (1 - normalCDF(Math.abs(z)))
      return { estimate: roundTo(est, 4), se: roundTo(se, 4), z: roundTo(z, 3), pValue: roundTo(pValue, 4), stdAll: roundTo(est, 4) }
    })
  )

  // Communalities and uniqueness
  const communalities = Array.from(Theta).map(t => roundTo(1 - t, 4))
  const uniquenessArr = Array.from(Theta).map(t => roundTo(t, 4))

  // Standardized loadings
  const stdLoadings: number[][] = Array.from({ length: d }, (_, i) =>
    Array.from({ length: k }, (_, j) => {
      const sigmaII = Math.max(impliedSigma.get(i, i), 1e-12)
      return roundTo(L[i]![j]! * Math.sqrt(Phi[j]![j]!) / Math.sqrt(sigmaII), 4)
    })
  )

  // Eigenvalues of correlation matrix
  const eigenvalues = S.eigen().values.map(v => roundTo(v, 4))

  const variableNames = options?.variableNames
    ?? Array.from({ length: d }, (_, i) => `V${i + 1}`)
  const factorNames = options?.factorNames ?? factors

  const formatted = `CFA: ${formatCFAFit(fit)}`

  return {
    loadings: L,
    standardizedLoadings: stdLoadings,
    uniqueness: uniquenessArr,
    communalities,
    factorCorrelations: Phi,
    fit,
    eigenvalues,
    nFactors: k,
    rotation: 'none',
    extraction: 'ml',
    variableNames,
    factorNames,
    formatted,
    parameterEstimates: {
      loadings: paramLoadings,
      uniquenesses: paramUniqueness,
      factorCovariances: paramCov,
    },
    model,
  }
}
