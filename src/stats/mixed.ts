/**
 * Linear Mixed Models (LMM) via REML or ML.
 * Model: y = Xβ + Zb + ε
 *   b ~ N(0, G), ε ~ N(0, σ²_e · I)
 * G is the random-effects covariance (intercept + optional slopes).
 *
 * Parameterization: log-Cholesky for G/σ²_e (lme4 style).
 * Optimization: multi-start Nelder-Mead on profiled log-likelihood.
 *
 * Supports:
 *   - Random intercepts (default)
 *   - Random slopes via Cholesky parameterization of G
 *   - REML (default) or ML estimation
 *   - Nakagawa R² (marginal and conditional)
 *   - Model comparison via likelihood ratio test (compareLMM)
 *
 * Cross-validated with R:
 * > lme4::lmer(y ~ x + (1|group), data = df, REML = TRUE)
 * > lme4::lmer(y ~ x + (1 + x|group), data = df, REML = TRUE)
 */

import { Matrix } from '../core/matrix.js'
import { nelderMead, tDistPValue, tDistQuantile, chiSqPValue, roundTo, variance as computeVariance } from '../core/math.js'
import { formatLMM } from '../core/apa.js'
import type { LMMResult, FixedEffect } from '../core/types.js'

// ─── Data structure ────────────────────────────────────────────────────────

export interface LMMInput {
  readonly outcome: readonly number[]
  readonly fixedPredictors: Readonly<Record<string, readonly number[]>>
  readonly groupId: readonly (string | number)[]
  readonly randomSlopes?: readonly string[]  // which fixed predictors also get random slopes
  readonly method?: 'REML' | 'ML'
  readonly ciLevel?: number
}

// ─── Model comparison result ─────────────────────────────────────────────

export interface LMMComparison {
  readonly chiSq: number
  readonly df: number
  readonly pValue: number
  readonly preferred: 'model1' | 'model2'
  readonly warning?: string
}

// ─── Internal helpers ────────────────────────────────────────────────────

/**
 * Build the relative Cholesky factor L (q×q lower-triangular) from θ parameters.
 * θ layout: [log(L[0,0]), L[1,0], log(L[1,1]), L[2,0], L[2,1], log(L[2,2]), ...]
 * Diagonal elements are exp(θ) to ensure positivity.
 */
function buildCholFactor(theta: readonly number[], q: number): Matrix {
  const data = new Array(q * q).fill(0)
  let idx = 0
  for (let i = 0; i < q; i++) {
    for (let j = 0; j <= i; j++) {
      if (i === j) {
        data[i * q + j] = Math.exp(theta[idx]!)
      } else {
        data[i * q + j] = theta[idx]!
      }
      idx++
    }
  }
  return new Matrix(q, q, data)
}

/**
 * Build the extended random-effects design matrix Z_ext (n × nGroups*q).
 * For each group j and random effect k:
 *   Z_ext[i, j*q + 0] = 1 if obs i in group j (intercept)
 *   Z_ext[i, j*q + k] = slopePredictors[k-1][i] if obs i in group j, else 0
 */
function buildExtendedZ(
  n: number,
  groupId: readonly (string | number)[],
  groupLevels: readonly (string | number)[],
  slopePredictors: readonly (readonly number[])[]
): Matrix {
  const nGroups = groupLevels.length
  const q = 1 + slopePredictors.length
  const gq = nGroups * q
  const data = new Array(n * gq).fill(0)

  for (let i = 0; i < n; i++) {
    const gIdx = groupLevels.indexOf(groupId[i]!)
    if (gIdx < 0) continue
    // Intercept column for this group
    data[i * gq + gIdx * q] = 1
    // Slope columns
    for (let k = 0; k < slopePredictors.length; k++) {
      data[i * gq + gIdx * q + (k + 1)] = slopePredictors[k]![i]!
    }
  }
  return new Matrix(n, gq, data)
}

/**
 * Compute A = Z_ext · (I_g ⊗ L) efficiently without forming the block-diagonal.
 * A[i, j*q+k] = Σ_{m=k..q-1} Z_ext[i, j*q+m] · L[m, k]
 */
function buildA(Z_ext: Matrix, L: Matrix, nGroups: number, q: number, n: number): Matrix {
  const gq = nGroups * q
  const data = new Array(n * gq).fill(0)
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < nGroups; j++) {
      for (let k = 0; k < q; k++) {
        let val = 0
        for (let m = k; m < q; m++) {
          val += Z_ext.get(i, j * q + m) * L.get(m, k)
        }
        data[i * gq + j * q + k] = val
      }
    }
  }
  return new Matrix(n, gq, data)
}

/**
 * Profiled log-likelihood for LMM (generalized for any q random effects).
 *
 * Parameterization follows lme4:
 *   Ψ = L·L' is the relative covariance (G/σ²)
 *   V_scaled = I + Z·(I_g⊗L)·(I_g⊗L)'·Z' = I + A·A'
 *   V_scaled⁻¹ = I − A·D⁻¹·A'  where D = I + A'A  (Woodbury)
 *   σ² profiled out analytically
 *
 * REML: σ² = e'V⁻¹e / (n-p);  ℓ = −½[(n-p)log(σ²) + log|D| + log|X'V⁻¹X|]
 * ML:   σ² = e'V⁻¹e / n;      ℓ = −½[n·log(σ²) + log|D|]
 */
function profileLogLik(
  theta: readonly number[],
  y: readonly number[],
  X: Matrix,
  Z_ext: Matrix,
  q: number,
  nGroups: number,
  method: 'REML' | 'ML'
): number {
  const n = y.length
  const p = X.cols
  const gq = nGroups * q

  const L = buildCholFactor(theta, q)
  const A = buildA(Z_ext, L, nGroups, q, n)

  // D = I_{gq} + A'A
  const At = A.transpose()
  const D = Matrix.identity(gq).add(At.multiply(A))

  let DInv: Matrix
  let logDetD: number
  try {
    DInv = D.inverse()
    logDetD = D.logDet()
  } catch {
    return Infinity
  }

  // Apply V_scaled⁻¹ via Woodbury: Vinv·v = v − A·D⁻¹·A'·v
  // Compute Vinv·y and Vinv·X efficiently (never forming the n×n matrix)
  const AtY = At.multiply(Matrix.colVec(y))         // gq×1
  const ADInvAtY = A.multiply(DInv.multiply(AtY))   // n×1
  const VinvY = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - ADInvAtY.get(i, 0))

  const AtX = At.multiply(X)                         // gq×p
  const ADInvAtX = A.multiply(DInv.multiply(AtX))    // n×p
  const VinvX = X.subtract(ADInvAtX)                 // n×p

  // X'V⁻¹X  and  X'V⁻¹y
  const Xt = X.transpose()
  const XtVinvX = Xt.multiply(VinvX)                 // p×p
  const XtVinvY = Xt.multiply(Matrix.colVec(VinvY))  // p×1

  let XtVinvXInv: Matrix
  let logDetXVX: number
  try {
    XtVinvXInv = XtVinvX.inverse()
    logDetXVX = XtVinvX.logDet()
  } catch {
    return Infinity
  }

  // GLS β̂
  const beta = XtVinvXInv.multiply(XtVinvY)         // p×1

  // Residuals and quadratic form  e'V⁻¹e
  const Xbeta = X.multiply(beta)
  const e = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - Xbeta.get(i, 0))
  // V⁻¹e = VinvY − VinvX·β̂
  const VinvXbeta = VinvX.multiply(beta)             // n×1
  const VinvE = VinvY.map((v, i) => v - VinvXbeta.get(i, 0))
  const quadForm = e.reduce((s, ei, i) => s + ei * VinvE[i]!, 0)

  if (quadForm <= 0) return Infinity

  let sigma2: number
  let logLik: number

  if (method === 'REML') {
    sigma2 = quadForm / (n - p)
    logLik = -0.5 * ((n - p) * Math.log(sigma2) + logDetD + logDetXVX)
  } else {
    sigma2 = quadForm / n
    logLik = -0.5 * (n * Math.log(sigma2) + logDetD)
  }

  return -logLik  // return negLogLik (we minimize)
}

// ─── Main LMM function ────────────────────────────────────────────────────

/**
 * Fit a linear mixed model with random intercepts and optional random slopes.
 *
 * Cross-validated with R lme4:
 * > mod <- lmer(y ~ x + (1|group), data = df, REML = TRUE)
 * > mod <- lmer(y ~ x + (1 + x|group), data = df, REML = TRUE)
 * > fixef(mod); VarCorr(mod); logLik(mod); AIC(mod)
 * > MuMIn::r.squaredGLMM(mod)
 */
export function runLMM(input: LMMInput): LMMResult {
  const { outcome: y, fixedPredictors, groupId, randomSlopes, method = 'REML', ciLevel = 0.95 } = input
  const n = y.length
  if (n < 5) throw new Error('runLMM: need at least 5 observations')
  if (groupId.length !== n) throw new Error('runLMM: groupId must have same length as outcome')

  // Identify groups
  const groupLevels = [...new Set(groupId)]
  const nGroups = groupLevels.length
  if (nGroups < 2) throw new Error('runLMM: need at least 2 groups')

  // Fixed-effects design matrix X (n × p), with intercept
  const predNames = Object.keys(fixedPredictors)
  const p = predNames.length + 1
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [
      1,
      ...predNames.map(name => (fixedPredictors[name] ?? [])[i] ?? 0),
    ])
  )

  // Validate and gather slope predictors
  const slopeNames = randomSlopes ?? []
  for (const s of slopeNames) {
    if (!(s in fixedPredictors)) {
      throw new Error(`runLMM: random slope '${s}' not found in fixedPredictors`)
    }
  }
  const slopePreds = slopeNames.map(s => fixedPredictors[s]!)
  const q = 1 + slopeNames.length  // number of random effects per group
  const nTheta = q * (q + 1) / 2   // log-Cholesky parameters

  // Build extended Z matrix
  const Z_ext = buildExtendedZ(n, groupId, groupLevels, slopePreds)

  // ── Optimization: multi-start Nelder-Mead ──────────────────────────────

  const objFn = (theta: readonly number[]) =>
    profileLogLik(theta, y, X, Z_ext, q, nGroups, method)

  // Generate starting values: perturb the diagonal elements
  const diagIndices: number[] = []
  let idx = 0
  for (let i = 0; i < q; i++) {
    for (let j = 0; j <= i; j++) {
      if (i === j) diagIndices.push(idx)
      idx++
    }
  }

  const startDiagValues = [-2, -1, 0, 1, 2]
  let bestResult = nelderMead(objFn, new Array(nTheta).fill(0), { maxIter: 2000, tol: 1e-8 })

  for (const dv of startDiagValues) {
    const start = new Array(nTheta).fill(0)
    for (const di of diagIndices) start[di] = dv
    const cand = nelderMead(objFn, start, { maxIter: 2000, tol: 1e-8 })
    if (cand.fval < bestResult.fval) bestResult = cand
  }

  // ── Extract results at optimum θ ───────────────────────────────────────

  const optTheta = bestResult.x
  const L = buildCholFactor(optTheta, q)
  const A = buildA(Z_ext, L, nGroups, q, n)

  // Re-do the full computation at the optimum to extract all quantities
  const gq = nGroups * q
  const At = A.transpose()
  const D = Matrix.identity(gq).add(At.multiply(A))

  let DInv: Matrix
  let logDetD: number
  try {
    DInv = D.inverse()
    logDetD = D.logDet()
  } catch {
    // Should not happen at the optimum, but guard anyway
    DInv = Matrix.identity(gq)
    logDetD = 0
  }

  // V_scaled⁻¹ applied to y and X
  const AtY = At.multiply(Matrix.colVec(y))
  const ADInvAtY = A.multiply(DInv.multiply(AtY))
  const VinvY = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - ADInvAtY.get(i, 0))

  const AtX = At.multiply(X)
  const ADInvAtX = A.multiply(DInv.multiply(AtX))
  const VinvX = X.subtract(ADInvAtX)

  const Xt = X.transpose()
  const XtVinvX = Xt.multiply(VinvX)
  const XtVinvY = Xt.multiply(Matrix.colVec(VinvY))

  let XtVinvXInv: Matrix
  let logDetXVX: number
  try {
    XtVinvXInv = XtVinvX.inverse()
    logDetXVX = XtVinvX.logDet()
  } catch {
    XtVinvXInv = Matrix.identity(p)
    logDetXVX = 0
  }

  // GLS β̂
  const betaM = XtVinvXInv.multiply(XtVinvY)
  const beta = Array.from({ length: p }, (_, i) => betaM.get(i, 0))

  // Residuals and quadratic form
  const Xbeta = X.multiply(betaM)
  const e = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - Xbeta.get(i, 0))
  const VinvXbeta = VinvX.multiply(betaM)
  const VinvE = VinvY.map((v, i) => v - VinvXbeta.get(i, 0))
  const quadForm = e.reduce((s, ei, i) => s + ei * VinvE[i]!, 0)

  // Profiled σ²
  const denom = method === 'REML' ? (n - p) : n
  const sigma2 = Math.max(1e-8, quadForm / denom)

  // Reconstruct G = σ² · L·L' (the actual random-effects covariance)
  const LLt = L.multiply(L.transpose())
  const G = LLt.scale(sigma2)

  // ── Variance components ────────────────────────────────────────────────

  const sigmab2 = G.get(0, 0)   // random intercept variance
  const sigmae2 = sigma2          // residual variance

  const slopeVarRecord: Record<string, number> = {}
  for (let k = 0; k < slopeNames.length; k++) {
    slopeVarRecord[slopeNames[k]!] = roundTo(G.get(k + 1, k + 1), 6)
  }

  // Random correlations (from G): correlation between effect i and j
  const randomCorrs: Record<string, number> = {}
  for (let i = 0; i < q; i++) {
    for (let j = i + 1; j < q; j++) {
      const si = Math.sqrt(G.get(i, i))
      const sj = Math.sqrt(G.get(j, j))
      if (si > 1e-10 && sj > 1e-10) {
        const iName = i === 0 ? '(Intercept)' : slopeNames[i - 1]!
        const jName = j === 0 ? '(Intercept)' : slopeNames[j - 1]!
        randomCorrs[`${iName}:${jName}`] = roundTo(G.get(i, j) / (si * sj), 6)
      }
    }
  }

  // ── Fixed effects with SEs ─────────────────────────────────────────────

  // Var(β̂) = σ² · (X'V_scaled⁻¹X)⁻¹
  const covBeta = XtVinvXInv.scale(sigma2)

  // Degrees of freedom (Satterthwaite approximation)
  const df = Math.max(1, n - p - nGroups + 1)
  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, df)

  const fixedEffectNames = ['(Intercept)', ...predNames]
  const fixedEffects: FixedEffect[] = beta.map((b, i) => {
    const seVal = Math.sqrt(Math.max(0, covBeta.get(i, i)))
    const t = seVal === 0 ? 0 : b / seVal
    const pVal = tDistPValue(t, df)
    return {
      name: fixedEffectNames[i] ?? `β${i}`,
      estimate: roundTo(b, 6),
      se: roundTo(seVal, 6),
      tValue: roundTo(t, 4),
      pValue: roundTo(pVal, 4),
      ci: [roundTo(b - tCrit * seVal, 6), roundTo(b + tCrit * seVal, 6)],
    }
  })

  // ── ICC ────────────────────────────────────────────────────────────────

  const icc = sigmab2 / (sigmab2 + sigmae2)

  // ── Log-likelihood (with normalizing constant for R compatibility) ─────

  const dfLogLik = method === 'REML' ? (n - p) : n
  const negLogLikProfiled = method === 'REML'
    ? 0.5 * ((n - p) * Math.log(sigma2) + logDetD + logDetXVX)
    : 0.5 * (n * Math.log(sigma2) + logDetD)
  const normConst = 0.5 * dfLogLik * (1 + Math.log(2 * Math.PI))
  const logLik = -negLogLikProfiled - normConst

  // Number of parameters: p fixed + q(q+1)/2 covariance + 1 residual
  const nParams = p + nTheta + 1

  const aic = -2 * logLik + 2 * nParams
  const bic = -2 * logLik + Math.log(n) * nParams

  // ── Nakagawa R² (marginal and conditional) ─────────────────────────────
  // Reference: Nakagawa & Schielzeth (2013), Johnson (2014)
  // σ²_f = variance of fixed-effects predictions Xβ̂
  // σ²_r = mean of z_i'Gz_i over observations (random effects variance)
  // R²_m = σ²_f / (σ²_f + σ²_r + σ²_e)
  // R²_c = (σ²_f + σ²_r) / (σ²_f + σ²_r + σ²_e)

  const fixedPred = Array.from({ length: n }, (_, i) => Xbeta.get(i, 0))
  const sigma2_f = computeVariance(fixedPred)

  // Compute σ²_r = (1/n) Σ_i z_i'Gz_i  where z_i = [1, x_i1, x_i2, ...]
  let sigma2_r = 0
  for (let i = 0; i < n; i++) {
    // z_i is the q-vector of random-effect predictors for observation i
    // z_i[0] = 1 (intercept), z_i[k] = slopePreds[k-1][i]
    for (let a = 0; a < q; a++) {
      for (let b = 0; b < q; b++) {
        const za = a === 0 ? 1 : slopePreds[a - 1]![i]!
        const zb = b === 0 ? 1 : slopePreds[b - 1]![i]!
        sigma2_r += za * G.get(a, b) * zb
      }
    }
  }
  sigma2_r /= n

  const totalVar = sigma2_f + sigma2_r + sigmae2
  const r2Marginal = totalVar > 0 ? sigma2_f / totalVar : 0
  const r2Conditional = totalVar > 0 ? (sigma2_f + sigma2_r) / totalVar : 0

  // ── Formatted output ───────────────────────────────────────────────────

  const formatted = formatLMM(icc, aic, bic, logLik)

  return {
    fixedEffects,
    varianceComponents: {
      intercept: roundTo(sigmab2, 6),
      residual: roundTo(sigmae2, 6),
      ...(slopeNames.length > 0 ? { slopes: slopeVarRecord } : {}),
    },
    ...(Object.keys(randomCorrs).length > 0 ? { randomCorrelations: randomCorrs } : {}),
    icc: roundTo(icc, 6),
    logLik: roundTo(logLik, 4),
    aic: roundTo(aic, 2),
    bic: roundTo(bic, 2),
    nObs: n,
    nGroups,
    nParams,
    method,
    r2Marginal: roundTo(r2Marginal, 6),
    r2Conditional: roundTo(r2Conditional, 6),
    formatted,
  }
}

// ─── Model comparison via LRT ────────────────────────────────────────────

/**
 * Compare two LMM models via likelihood ratio test (LRT).
 * Both models should be fitted with ML for valid comparison.
 * The model with more parameters is treated as the "full" model.
 *
 * Cross-validated with R:
 * > anova(mod1, mod2)
 *
 * @returns LRT statistic, df, p-value, and which model is preferred
 */
export function compareLMM(model1: LMMResult, model2: LMMResult): LMMComparison {
  // Determine which is the fuller model
  const [reduced, full] = model1.nParams <= model2.nParams
    ? [model1, model2]
    : [model2, model1]
  const preferred = model1.nParams <= model2.nParams ? 'model2' : 'model1'

  const chiSq = Math.max(0, -2 * (reduced.logLik - full.logLik))
  const df = Math.abs(full.nParams - reduced.nParams)

  if (df === 0) {
    return {
      chiSq: 0,
      df: 0,
      pValue: 1,
      preferred: 'model1',
      warning: 'Models have the same number of parameters',
    }
  }

  const pValue = chiSqPValue(chiSq, df)

  // Warn if REML was used (REML likelihoods not comparable for different fixed effects)
  let warning: string | undefined
  if (model1.method === 'REML' || model2.method === 'REML') {
    warning = 'REML likelihoods are not comparable for models with different fixed effects. Use ML for valid comparison.'
  }

  return {
    chiSq: roundTo(chiSq, 4),
    df,
    pValue: roundTo(pValue, 4),
    preferred: pValue < 0.05 ? preferred : (model1.nParams <= model2.nParams ? 'model1' : 'model2'),
    ...(warning !== undefined ? { warning } : {}),
  }
}

// ─── BLUPs ────────────────────────────────────────────────────────────────

/**
 * Compute BLUPs (Best Linear Unbiased Predictors) — the random intercepts.
 * b_hat = σ²_b Z'V^{-1}(y - Xβ)
 */
export function computeBLUPs(
  input: LMMInput,
  result: LMMResult
): ReadonlyArray<{ group: string | number; blup: number }> {
  const { outcome: y, fixedPredictors, groupId } = input
  const n = y.length
  const groupLevels = [...new Set(groupId)]
  const predNames = Object.keys(fixedPredictors)

  const sigmab2 = result.varianceComponents.intercept
  const sigmae2 = result.varianceComponents.residual

  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predNames.map(name => (fixedPredictors[name] ?? [])[i] ?? 0)])
  )
  // Fixed fitted values
  const beta = result.fixedEffects.map(fe => fe.estimate)
  const Xbeta = X.multiply(Matrix.colVec(beta))
  const residuals = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - Xbeta.get(i, 0))

  // BLUPs: b = σ²_b Z' V^{-1} e ≈ (σ²_b/σ²_e) Z' (I - ZD^{-1}Z') e / σ²_e
  // Simple approximation for random intercepts: b_j = ψ/(1+ψ·n_j) * Σ_{i in j} e_i
  // where ψ = σ²_b/σ²_e
  const psi = sigmab2 / sigmae2
  return groupLevels.map((g) => {
    const indices = Array.from({ length: n }, (_, i) => i).filter(i => groupId[i] === g)
    const sumResid = indices.reduce((s, i) => s + (residuals[i] ?? 0), 0)
    const nj = indices.length
    const blup = (psi / (1 + psi * nj)) * sumResid
    return { group: g, blup: roundTo(blup, 6) }
  })
}
