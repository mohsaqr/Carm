/**
 * Logistic Generalized Linear Mixed Model (GLMM) via Laplace approximation.
 * Model: logit(μ_i) = x_i'β + z_i'b, b ~ N(0, G), y_i | b ~ Bernoulli(μ_i)
 *
 * Algorithm: Cholesky-based PLS (Penalized Least Squares) matching lme4's
 * internal PIRLS. Uses the u-parameterization (u = Λ⁻¹b, u ~ N(0,I))
 * with Cholesky forward/back substitution for the Schur complement.
 * Outer loop: multi-start Nelder-Mead over log-Cholesky θ.
 *
 * Key difference from a naive augmented system solve: lme4's PLS computes
 * β via the Schur complement using the Cholesky factor of Z*'WZ*+I,
 * which gives a specific β-b decomposition that all R optimizers agree on.
 * A direct 12×12 system solve gives a numerically different decomposition
 * along the β₀-b ridge.
 *
 * Cross-validated with R:
 * > lme4::glmer(y ~ x + (1|group), data = df, family = binomial)
 */

import { Matrix } from '../core/matrix.js'
import { nelderMead, normalCDF, normalQuantile, roundTo } from '../core/math.js'
import { formatGLMM } from '../core/apa.js'
import { buildCholFactor, buildExtendedZ } from './mixed.js'
import type { GLMMResult, GLMMFixedEffect } from '../core/types.js'

// ─── Data structure ────────────────────────────────────────────────────────

export interface GLMMInput {
  readonly outcome: readonly number[]
  readonly fixedPredictors: Readonly<Record<string, readonly number[]>>
  readonly groupId: readonly (string | number)[]
  readonly randomSlopes?: readonly string[]
  readonly ciLevel?: number
  readonly maxIter?: number
  readonly tol?: number
}

// ─── Logistic helpers ──────────────────────────────────────────────────────

function logistic(eta: number): number {
  if (eta > 20) return 1 - 1e-9
  if (eta < -20) return 1e-9
  return 1 / (1 + Math.exp(-eta))
}

function binomialDeviance(y: readonly number[], mu: readonly number[]): number {
  let dev = 0
  for (let i = 0; i < y.length; i++) {
    const yi = y[i]!
    const mi = Math.max(1e-10, Math.min(1 - 1e-10, mu[i]!))
    dev += -2 * (yi * Math.log(mi) + (1 - yi) * Math.log(1 - mi))
  }
  return dev
}

// ─── Triangular solves ─────────────────────────────────────────────────────

/** Forward solve: L · x = B where L is lower-triangular. Returns x. */
function forwardSolve(L: Matrix, B: Matrix): Matrix {
  const n = L.rows
  const m = B.cols
  const data = new Array(n * m).fill(0)
  for (let col = 0; col < m; col++) {
    for (let i = 0; i < n; i++) {
      let val = B.get(i, col)
      for (let j = 0; j < i; j++) val -= L.get(i, j) * data[j * m + col]!
      data[i * m + col] = val / L.get(i, i)
    }
  }
  return new Matrix(n, m, data)
}

/** Back solve: L' · x = B where L is lower-triangular. Returns x. */
function backSolve(L: Matrix, B: Matrix): Matrix {
  const n = L.rows
  const m = B.cols
  const data = new Array(n * m).fill(0)
  for (let col = 0; col < m; col++) {
    for (let i = n - 1; i >= 0; i--) {
      let val = B.get(i, col)
      for (let j = i + 1; j < n; j++) val -= L.get(j, i) * data[j * m + col]!
      data[i * m + col] = val / L.get(i, i)
    }
  }
  return new Matrix(n, m, data)
}

// ─── Build Z* = Z · blkdiag(Λ,...,Λ) ──────────────────────────────────────

function buildZstar(Z_ext: Matrix, Lambda: Matrix, n: number, nGroups: number, q: number): Matrix {
  const gq = nGroups * q
  const data = new Array(n * gq).fill(0)
  for (let i = 0; i < n; i++) {
    for (let g = 0; g < nGroups; g++) {
      for (let k = 0; k < q; k++) {
        let val = 0
        for (let m = 0; m < q; m++) val += Z_ext.get(i, g * q + m) * Lambda.get(m, k)
        data[i * gq + g * q + k] = val
      }
    }
  }
  return new Matrix(n, gq, data)
}

// ─── Fixed-effects logistic regression init ──────────────────────────────

function logisticInit(y: readonly number[], X: Matrix): number[] {
  const n = y.length
  const p = X.cols
  const yMean = y.reduce((s, v) => s + v, 0) / n
  const betaArr = new Array(p).fill(0)
  betaArr[0] = Math.log(Math.max(1e-4, yMean) / Math.max(1e-4, 1 - yMean))
  for (let iter = 0; iter < 15; iter++) {
    const eta = Array.from({ length: n }, (_, i) => {
      let v = 0
      for (let j = 0; j < p; j++) v += X.get(i, j) * betaArr[j]!
      return v
    })
    const mu = eta.map(logistic)
    const w = mu.map(m => Math.max(1e-10, m * (1 - m)))
    const z = eta.map((e, i) => e + (y[i]! - mu[i]!) / w[i]!)
    const WXd = new Array(n * p).fill(0)
    const Wzd = new Array(n).fill(0)
    for (let i = 0; i < n; i++) {
      Wzd[i] = w[i]! * z[i]!
      for (let j = 0; j < p; j++) WXd[i * p + j] = w[i]! * X.get(i, j)
    }
    const WX = new Matrix(n, p, WXd)
    const Xt = X.transpose()
    try {
      const newB = Xt.multiply(WX).inverse().multiply(Xt.multiply(Matrix.colVec(Wzd)))
      for (let j = 0; j < p; j++) betaArr[j] = newB.get(j, 0)
    } catch { break }
  }
  return betaArr
}

// ─── PIRLS inner loop (PLS formulation, matching lme4) ────────────────────

interface PIRLSResult {
  readonly beta: readonly number[]
  readonly b: readonly number[]
  readonly mu: readonly number[]
  readonly penalizedDeviance: number
  readonly converged: boolean
}

/**
 * Cholesky-based PLS (lme4's PIRLS formulation).
 *
 * Works with u = Λ⁻¹b (standardized random effects, u ~ N(0,I)).
 * At each step:
 *   1. L_z = chol(Z*'WZ* + I)
 *   2. RZX = L_z⁻¹ Z*'WX  (forward solve)
 *   3. β = (X'WX - RZX'RZX)⁻¹ (X'Wz - RZX' L_z⁻¹ Z*'Wz)
 *   4. u = L_z⁻ᵀ (L_z⁻¹ Z*'Wz_resid)
 *
 * This matches R's lme4 exactly: same Schur complement, same
 * Cholesky-based forward/back substitution.
 */
function pirlsInner(
  theta: readonly number[],
  y: readonly number[],
  X: Matrix,
  Z_ext: Matrix,
  q: number,
  nGroups: number,
  maxIter: number,
  tol: number,
  initBeta?: readonly number[],
  initU?: readonly number[],
): PIRLSResult {
  const n = y.length
  const p = X.cols
  const gq = nGroups * q

  const Lambda = buildCholFactor(theta, q)
  const Zstar = buildZstar(Z_ext, Lambda, n, nGroups, q)

  // Initialize
  let betaArr = initBeta ? [...initBeta] : (() => {
    const yMean = y.reduce((s, v) => s + v, 0) / n
    const b = new Array(p).fill(0)
    b[0] = Math.log(Math.max(1e-4, yMean) / Math.max(1e-4, 1 - yMean))
    return b
  })()
  let uArr = initU ? [...initU] : new Array(gq).fill(0)

  let mu = new Array(n).fill(0.5)
  let prevPenDev = Infinity
  let converged = false

  for (let iter = 0; iter < maxIter; iter++) {
    // η = Xβ + Z*u
    const betaM = Matrix.colVec(betaArr)
    const uM = Matrix.colVec(uArr)
    const Xbeta = X.multiply(betaM)
    const Zstaru = Zstar.multiply(uM)
    const eta = Array.from({ length: n }, (_, i) => Xbeta.get(i, 0) + Zstaru.get(i, 0))

    mu = eta.map(logistic)
    const w = mu.map(m => Math.max(1e-10, m * (1 - m)))
    const zWork = eta.map((e, i) => e + (y[i]! - mu[i]!) / w[i]!)

    // Weighted matrices
    const WZstarData = new Array(n * gq).fill(0)
    const WXdata = new Array(n * p).fill(0)
    const WzData = new Array(n).fill(0)
    for (let i = 0; i < n; i++) {
      const wi = w[i]!
      WzData[i] = wi * zWork[i]!
      for (let j = 0; j < gq; j++) WZstarData[i * gq + j] = wi * Zstar.get(i, j)
      for (let j = 0; j < p; j++) WXdata[i * p + j] = wi * X.get(i, j)
    }
    const WZstar = new Matrix(n, gq, WZstarData)
    const WX = new Matrix(n, p, WXdata)
    const Wz = Matrix.colVec(WzData)
    const Zstart = Zstar.transpose()
    const Xt = X.transpose()

    // 1. H = Z*'WZ* + I
    const ZtWZ = Zstart.multiply(WZstar)
    const H = ZtWZ.add(Matrix.identity(gq))

    // 2. Cholesky: L_z L_z' = H
    let Lz: Matrix
    try {
      Lz = H.cholesky()
    } catch {
      return { beta: betaArr, b: new Array(gq).fill(0), mu, penalizedDeviance: Infinity, converged: false }
    }

    // 3. Forward solve: RZX = L_z⁻¹ Z*'WX, cu = L_z⁻¹ Z*'Wz
    const ZtWX = Zstart.multiply(WX)
    const ZtWz = Zstart.multiply(Wz)
    const RZX = forwardSolve(Lz, ZtWX)
    const cu = forwardSolve(Lz, ZtWz)

    // 4. β via Schur complement: (X'WX - RZX'RZX)⁻¹ (X'Wz - RZX'cu)
    const XtWX = Xt.multiply(WX)
    const XtWz = Xt.multiply(Wz)
    const RZXt = RZX.transpose()
    const schurLHS = XtWX.subtract(RZXt.multiply(RZX))
    const schurRHS = XtWz.subtract(RZXt.multiply(cu))

    let newBetaM: Matrix
    try {
      newBetaM = schurLHS.inverse().multiply(schurRHS)
    } catch {
      return { beta: betaArr, b: new Array(gq).fill(0), mu, penalizedDeviance: Infinity, converged: false }
    }
    const newBeta = Array.from({ length: p }, (_, i) => newBetaM.get(i, 0))

    // 5. u via back-solve: L_z' u = cu - RZX β
    const cuMinusRZXbeta = cu.subtract(RZX.multiply(newBetaM))
    const newUM = backSolve(Lz, cuMinusRZXbeta)
    const newU = Array.from({ length: gq }, (_, i) => newUM.get(i, 0))

    // Penalized deviance: dev(y,μ) + ||u||²
    const newEta = Array.from({ length: n }, (_, i) => {
      let val = 0
      for (let j = 0; j < p; j++) val += X.get(i, j) * newBeta[j]!
      for (let j = 0; j < gq; j++) val += Zstar.get(i, j) * newU[j]!
      return val
    })
    const newMu = newEta.map(logistic)
    const dev = binomialDeviance(y, newMu)
    const uNormSq = newU.reduce((s, v) => s + v * v, 0)
    let penDev = dev + uNormSq

    const oldBeta = [...betaArr]

    // Step-halving
    let accepted = false
    if (penDev <= prevPenDev + 1e-4) {
      betaArr = newBeta
      uArr = newU
      mu = newMu
      accepted = true
    } else {
      let stepSize = 0.5
      const oldU = [...uArr]
      for (let attempt = 0; attempt < 10; attempt++) {
        const hBeta = oldBeta.map((b, i) => b + stepSize * (newBeta[i]! - b))
        const hU = oldU.map((u, i) => u + stepSize * (newU[i]! - u))
        const hEta = Array.from({ length: n }, (_, i) => {
          let val = 0
          for (let j = 0; j < p; j++) val += X.get(i, j) * hBeta[j]!
          for (let j = 0; j < gq; j++) val += Zstar.get(i, j) * hU[j]!
          return val
        })
        const hMu = hEta.map(logistic)
        const hDev = binomialDeviance(y, hMu)
        const hUNorm = hU.reduce((s, v) => s + v * v, 0)
        const hPenDev = hDev + hUNorm
        if (hPenDev <= prevPenDev + 1e-4) {
          betaArr = hBeta; uArr = hU; mu = hMu; penDev = hPenDev; accepted = true; break
        }
        stepSize *= 0.5
      }
      if (!accepted) penDev = prevPenDev
    }

    // Convergence
    const devConv = Math.abs(prevPenDev - penDev) < tol * (Math.abs(penDev) + 0.1)
    let maxBC = 0
    for (let i = 0; i < p; i++) {
      const rd = Math.abs(betaArr[i]! - oldBeta[i]!) / (1 + Math.abs(oldBeta[i]!))
      if (rd > maxBC) maxBC = rd
    }
    if (devConv || (accepted && maxBC < tol)) {
      converged = true; prevPenDev = penDev; break
    }
    prevPenDev = penDev
  }

  // Convert u to b: b_g = Λ · u_g
  const finalB = new Array(gq).fill(0)
  for (let g = 0; g < nGroups; g++) {
    for (let k = 0; k < q; k++) {
      let val = 0
      for (let m = 0; m < q; m++) val += Lambda.get(k, m) * uArr[g * q + m]!
      finalB[g * q + k] = val
    }
  }

  return { beta: betaArr, b: finalB, mu, penalizedDeviance: prevPenDev, converged }
}

// ─── Log-determinant (lme4 parameterization) ──────────────────────────────

function computeLogDetLme4(
  Z_ext: Matrix, L: Matrix, mu: readonly number[],
  n: number, gq: number, nGroups: number, q: number,
): number {
  const Zstar = buildZstar(Z_ext, L, n, nGroups, q)
  const wHat = mu.map(m => m * (1 - m))
  const WZdata = new Array(n * gq).fill(0)
  for (let i = 0; i < n; i++) {
    const wi = wHat[i]!
    for (let j = 0; j < gq; j++) WZdata[i * gq + j] = wi * Zstar.get(i, j)
  }
  const WZ = new Matrix(n, gq, WZdata)
  const H = Zstar.transpose().multiply(WZ).add(Matrix.identity(gq))
  return H.logDet()
}

// ─── Main GLMM function ───────────────────────────────────────────────────

export function runGLMM(input: GLMMInput): GLMMResult {
  const {
    outcome: y, fixedPredictors, groupId,
    randomSlopes, ciLevel = 0.95, maxIter = 25, tol = 1e-8,
  } = input
  const n = y.length
  if (n < 5) throw new Error('runGLMM: need at least 5 observations')
  if (groupId.length !== n) throw new Error('runGLMM: groupId must have same length as outcome')
  for (let i = 0; i < n; i++) {
    if (y[i] !== 0 && y[i] !== 1) throw new Error(`runGLMM: outcome must be 0/1 binary, got ${y[i]} at index ${i}`)
  }

  const groupLevels = [...new Set(groupId)]
  const nGroups = groupLevels.length
  if (nGroups < 2) throw new Error('runGLMM: need at least 2 groups')

  const predNames = Object.keys(fixedPredictors)
  const p = predNames.length + 1
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predNames.map(name => (fixedPredictors[name] ?? [])[i] ?? 0)])
  )

  const slopeNames = randomSlopes ?? []
  for (const s of slopeNames) {
    if (!(s in fixedPredictors)) throw new Error(`runGLMM: random slope '${s}' not found in fixedPredictors`)
  }
  const slopePreds = slopeNames.map(s => fixedPredictors[s]!)
  const q = 1 + slopeNames.length
  const nTheta = q * (q + 1) / 2
  const gq = nGroups * q

  const Z_ext = buildExtendedZ(n, groupId, groupLevels, slopePreds)
  const initBeta = logisticInit(y, X)

  // ── Optimization ─────────────────────────────────────────────────────────

  let cachedBeta = [...initBeta]
  let cachedU: number[] = new Array(gq).fill(0)

  const objFn = (theta: readonly number[]): number => {
    const pirls = pirlsInner(theta, y, X, Z_ext, q, nGroups, maxIter, tol, cachedBeta, cachedU)
    if (!pirls.converged && pirls.penalizedDeviance === Infinity) return Infinity
    cachedBeta = [...pirls.beta]
    // Convert b back to u for warm-start
    const Lambda = buildCholFactor(theta, q)
    let LambdaInv: Matrix
    try { LambdaInv = Lambda.inverse() } catch { return Infinity }
    cachedU = new Array(gq).fill(0)
    for (let g = 0; g < nGroups; g++) {
      for (let k = 0; k < q; k++) {
        let val = 0
        for (let m = 0; m < q; m++) val += LambdaInv.get(k, m) * pirls.b[g * q + m]!
        cachedU[g * q + k] = val
      }
    }
    try {
      return pirls.penalizedDeviance + computeLogDetLme4(Z_ext, Lambda, pirls.mu, n, gq, nGroups, q)
    } catch { return Infinity }
  }

  const diagIndices: number[] = []
  let idx = 0
  for (let i = 0; i < q; i++) {
    for (let j = 0; j <= i; j++) { if (i === j) diagIndices.push(idx); idx++ }
  }

  // Primary run: start from θ=0 (variance=1, matching R's default start)
  // with continuous warm-starting to track the solution path like lme4
  cachedBeta = [...initBeta]; cachedU = new Array(gq).fill(0)
  let bestResult = nelderMead(objFn, new Array(nTheta).fill(0), { maxIter: 5000, tol: 1e-10 })

  // Multi-start safety net: only try alternative starts if the primary
  // run found a clearly suboptimal solution. Reset warm-start for each
  // alternative to avoid poisoning, but prefer the primary run's decomposition.
  const altStarts = [-4, -2, 1, 2]
  for (const dv of altStarts) {
    cachedBeta = [...initBeta]; cachedU = new Array(gq).fill(0)
    const start = new Array(nTheta).fill(0)
    for (const di of diagIndices) start[di] = dv
    const cand = nelderMead(objFn, start, { maxIter: 3000, tol: 1e-10 })
    if (cand.fval < bestResult.fval - 0.1) bestResult = cand  // only switch if substantially better
  }

  // ── Extract results ────────────────────────────────────────────────────

  const optTheta = bestResult.x
  const finalPirls = pirlsInner(optTheta, y, X, Z_ext, q, nGroups, 100, 1e-12, initBeta)
  const beta = finalPirls.beta
  const bOpt = finalPirls.b
  const muOpt = finalPirls.mu

  const L = buildCholFactor(optTheta, q)
  const G = L.multiply(L.transpose())
  const sigmab2 = G.get(0, 0)

  const slopeVarRecord: Record<string, number> = {}
  for (let k = 0; k < slopeNames.length; k++) {
    slopeVarRecord[slopeNames[k]!] = roundTo(G.get(k + 1, k + 1), 6)
  }

  const randomCorrs: Record<string, number> = {}
  for (let i = 0; i < q; i++) {
    for (let j = i + 1; j < q; j++) {
      const si = Math.sqrt(G.get(i, i)), sj = Math.sqrt(G.get(j, j))
      if (si > 1e-10 && sj > 1e-10) {
        const iN = i === 0 ? '(Intercept)' : slopeNames[i - 1]!
        const jN = j === 0 ? '(Intercept)' : slopeNames[j - 1]!
        randomCorrs[`${iN}:${jN}`] = roundTo(G.get(i, j) / (si * sj), 6)
      }
    }
  }

  // ── SEs via Schur complement ──────────────────────────────────────────

  const wHat = muOpt.map(m => m * (1 - m))
  const Zstar = buildZstar(Z_ext, L, n, nGroups, q)
  const WZsData = new Array(n * gq).fill(0)
  const WXdata = new Array(n * p).fill(0)
  for (let i = 0; i < n; i++) {
    const wi = wHat[i]!
    for (let j = 0; j < gq; j++) WZsData[i * gq + j] = wi * Zstar.get(i, j)
    for (let j = 0; j < p; j++) WXdata[i * p + j] = wi * X.get(i, j)
  }
  const WZs = new Matrix(n, gq, WZsData)
  const WX = new Matrix(n, p, WXdata)
  const Xt = X.transpose()
  const Zst = Zstar.transpose()
  const XtWX = Xt.multiply(WX)
  const ZtWZ = Zst.multiply(WZs)
  const H = ZtWZ.add(Matrix.identity(gq))
  let Lz: Matrix
  try { Lz = H.cholesky() } catch { Lz = Matrix.identity(gq) }
  const RZX = forwardSolve(Lz, Zst.multiply(WX))
  const schurLHS = XtWX.subtract(RZX.transpose().multiply(RZX))
  let covBeta: Matrix
  try { covBeta = schurLHS.inverse() } catch { covBeta = Matrix.identity(p) }

  const zCrit = normalQuantile(1 - (1 - ciLevel) / 2)
  const fixedEffectNames = ['(Intercept)', ...predNames]

  const fixedEffects: GLMMFixedEffect[] = beta.map((b, i) => {
    const seVal = Math.sqrt(Math.max(0, covBeta.get(i, i)))
    const z = seVal === 0 ? 0 : b / seVal
    const pVal = 2 * (1 - normalCDF(Math.abs(z)))
    const ciLo = b - zCrit * seVal, ciHi = b + zCrit * seVal
    return {
      name: fixedEffectNames[i] ?? `β${i}`,
      estimate: roundTo(b, 6), se: roundTo(seVal, 6),
      zValue: roundTo(z, 4), pValue: roundTo(Math.max(0, pVal), 4),
      ci: [roundTo(ciLo, 6), roundTo(ciHi, 6)] as readonly [number, number],
      or: roundTo(Math.exp(b), 4),
      orCI: [roundTo(Math.exp(ciLo), 4), roundTo(Math.exp(ciHi), 4)] as readonly [number, number],
    }
  })

  const PI_SQ_OVER_3 = Math.PI * Math.PI / 3
  const icc = sigmab2 / (sigmab2 + PI_SQ_OVER_3)

  // Log-likelihood (Laplace)
  // Convert b to u for ||u||²
  let LambdaInv: Matrix
  try { LambdaInv = L.inverse() } catch { LambdaInv = Matrix.identity(q) }
  const uFinal = new Array(gq).fill(0)
  for (let g = 0; g < nGroups; g++) {
    for (let k = 0; k < q; k++) {
      let val = 0
      for (let m = 0; m < q; m++) val += LambdaInv.get(k, m) * bOpt[g * q + m]!
      uFinal[g * q + k] = val
    }
  }
  const uNormSq = uFinal.reduce((s, v) => s + v * v, 0)
  const dev = binomialDeviance(y, muOpt)
  let logDetH: number
  try { logDetH = computeLogDetLme4(Z_ext, L, muOpt, n, gq, nGroups, q) } catch { logDetH = 0 }
  const neg2LL = dev + uNormSq + logDetH
  const logLik = -0.5 * neg2LL

  const nParams = p + nTheta
  const aic = -2 * logLik + 2 * nParams
  const bic = -2 * logLik + Math.log(n) * nParams
  const deviance = -2 * logLik
  const formatted = formatGLMM(icc, aic, deviance)

  return {
    fixedEffects,
    varianceComponents: { intercept: roundTo(sigmab2, 6), ...(slopeNames.length > 0 ? { slopes: slopeVarRecord } : {}) },
    ...(Object.keys(randomCorrs).length > 0 ? { randomCorrelations: randomCorrs } : {}),
    icc: roundTo(icc, 6), logLik: roundTo(logLik, 4), deviance: roundTo(deviance, 4),
    aic: roundTo(aic, 2), bic: roundTo(bic, 2),
    nObs: n, nGroups, nParams, family: 'binomial', link: 'logit', formatted,
  }
}
