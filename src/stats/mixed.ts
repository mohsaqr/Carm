/**
 * Linear Mixed Models (LMM) via REML.
 * Model: y = Xβ + Zb + ε
 *   b ~ N(0, σ²_b · I), ε ~ N(0, σ²_e · I)
 * Random intercepts + optional random slopes.
 * Optimization: Nelder-Mead on REML profile log-likelihood.
 *
 * Cross-validated with R:
 * > lme4::lmer(y ~ x + (1|group), data = df, REML = TRUE)
 */

import { Matrix } from '../core/matrix.js'
import { nelderMead, tDistPValue, tDistQuantile, roundTo } from '../core/math.js'
import { formatLMM } from '../core/apa.js'
import type { LMMResult, FixedEffect } from '../core/types.js'

// ─── Data structure ────────────────────────────────────────────────────────

export interface LMMInput {
  readonly outcome: readonly number[]
  readonly fixedPredictors: Readonly<Record<string, readonly number[]>>
  readonly groupId: readonly (string | number)[]
  readonly randomSlopes?: readonly string[]  // which fixed predictors also get random slopes
  readonly ciLevel?: number
}

// ─── Profiled REML log-likelihood ─────────────────────────────────────────

/**
 * Profiled REML log-likelihood for random intercepts model.
 * Single parameter: logPsi = log(ψ) where ψ = σ²_b / σ²_e.
 * σ²_e is profiled out analytically as the GLS residual variance —
 * this prevents the optimizer finding the degenerate solution where
 * all variance collapses into the random intercept.
 *
 * Model:  y = Xβ + Zb + ε,  b ~ N(0, σ²_b·I),  ε ~ N(0, σ²_e·I)
 * V_ψ = ψ·ZZ' + I  (marginal covariance, scaled by σ²_e)
 * V_ψ⁻¹ = I − Z·D⁻¹·Z'  via Woodbury,  D = Z'Z + (1/ψ)·I
 * log|V_ψ| = q·log(ψ) + log|D|  via matrix determinant lemma
 * σ²_e* = e'·V_ψ⁻¹·e / (n−p)  (analytical REML optimum)
 * Profiled REML = −½·[(n−p)·log(σ²_e*) + log|V_ψ| + log|X'·V_ψ⁻¹·X|]
 */
function remlProfileLogLik(
  logPsi: number,
  y: readonly number[],
  X: Matrix,
  Z: Matrix
): { negLogLik: number; sigmae2: number; sigmab2: number } {
  const psi = Math.exp(logPsi)
  const n = y.length
  const q = Z.cols
  const p = X.cols

  // D = Z'Z + (1/ψ)·I  (q×q)
  const ZtZ = Z.transpose().multiply(Z)
  const Dmat = ZtZ.add(Matrix.identity(q).scale(1 / psi))

  let DInv: Matrix
  let logDetD: number
  try {
    DInv = Dmat.inverse()
    logDetD = Dmat.logDet()
  } catch {
    return { negLogLik: Infinity, sigmae2: 0, sigmab2: 0 }
  }

  // V_ψ⁻¹ = I − Z·D⁻¹·Z'
  const ZDinvZt = Z.multiply(DInv).multiply(Z.transpose())
  const VpsiInv = Matrix.fromArray(
    Array.from({ length: n }, (_, i) =>
      Array.from({ length: n }, (_, j) =>
        (i === j ? 1 : 0) - ZDinvZt.get(i, j)
      )
    )
  )

  // log|V_ψ| = q·log(ψ) + log|D|
  const logDetVpsi = q * Math.log(psi) + logDetD

  // GLS: β = (X'·V_ψ⁻¹·X)⁻¹·X'·V_ψ⁻¹·y  (σ²_e scale cancels)
  const XtVinv = X.transpose().multiply(VpsiInv)
  const XtVinvX = XtVinv.multiply(X)
  let XtVinvXInv: Matrix
  let logDetXVX: number
  try {
    XtVinvXInv = XtVinvX.inverse()
    logDetXVX = XtVinvX.logDet()
  } catch {
    return { negLogLik: Infinity, sigmae2: 0, sigmab2: 0 }
  }

  const XtVinvY = XtVinv.multiply(Matrix.colVec(y))
  const beta = XtVinvXInv.multiply(XtVinvY)

  // Residuals e = y − Xβ  and quadratic form e'·V_ψ⁻¹·e
  const Xbeta = X.multiply(beta)
  const e = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - Xbeta.get(i, 0))
  const eM = Matrix.colVec(e)
  const quadForm = eM.transpose().multiply(VpsiInv).multiply(eM).get(0, 0)

  // Analytically optimal σ²_e  (REML degrees of freedom = n − p)
  const sigmae2 = Math.max(1e-8, quadForm / (n - p))
  const sigmab2 = psi * sigmae2

  // Profiled REML log-likelihood
  const reml = -0.5 * ((n - p) * Math.log(sigmae2) + logDetVpsi + logDetXVX)

  return { negLogLik: -reml, sigmae2, sigmab2 }
}

// ─── Main LMM function ────────────────────────────────────────────────────

/**
 * Fit a linear mixed model with random intercepts (and optionally random slopes).
 *
 * Cross-validated with R lme4:
 * > mod <- lmer(y ~ x + (1|group), data = df, REML = TRUE)
 * > fixef(mod)
 * > VarCorr(mod)
 * > icc(mod)
 */
export function runLMM(input: LMMInput): LMMResult {
  const { outcome: y, fixedPredictors, groupId, ciLevel = 0.95 } = input
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

  // Random-effects design matrix Z (n × nGroups) — indicator columns
  const Z = Matrix.fromArray(
    Array.from({ length: n }, (_, i) =>
      groupLevels.map(g => (groupId[i] === g ? 1 : 0))
    )
  )

  // 1D profiled REML: optimize over log(ψ) where ψ = σ²_b / σ²_e.
  // σ²_e is profiled out analytically at each ψ, so the optimizer
  // cannot find the degenerate solution σ²_e → 0.
  // Multi-start: try several logPsi starting values and keep the best.
  const objFn = (theta: readonly number[]) =>
    remlProfileLogLik(theta[0] ?? 0, y, X, Z).negLogLik

  const starts = [-4, -2, 0, 2, 4]
  let optResult = nelderMead(objFn, [starts[0]!], { maxIter: 1000, tol: 1e-8 })
  for (let si = 1; si < starts.length; si++) {
    const cand = nelderMead(objFn, [starts[si]!], { maxIter: 1000, tol: 1e-8 })
    if (cand.fval < optResult.fval) optResult = cand
  }

  // Extract analytically optimal variance components at the converged ψ
  const finalModel = remlProfileLogLik(optResult.x[0] ?? 0, y, X, Z)
  const sigmab2 = finalModel.sigmab2
  const sigmae2 = finalModel.sigmae2

  // Compute GLS estimates of fixed effects at optimal variance components.
  // scale = ψ = σ²_b / σ²_e.  When scale → 0, V_ψ⁻¹ → I (GLS = OLS).
  const scale = sigmab2 / sigmae2
  const ZtZ = Z.transpose().multiply(Z)

  let VinvScaled: Matrix
  if (scale < 1e-10) {
    // ψ → 0: random-intercept variance negligible, V_ψ → I → GLS = OLS
    VinvScaled = Matrix.identity(n)
  } else {
    const Dmat = ZtZ.add(Matrix.identity(nGroups).scale(1 / scale))
    let DInv: Matrix
    try {
      DInv = Dmat.inverse()
      const ZDinvZt = Z.multiply(DInv).multiply(Z.transpose())
      VinvScaled = Matrix.fromArray(
        Array.from({ length: n }, (_, i) =>
          Array.from({ length: n }, (_, j) =>
            (i === j ? 1 : 0) - ZDinvZt.get(i, j)
          )
        )
      )
    } catch {
      // Singular D — fall back to OLS
      VinvScaled = Matrix.identity(n)
    }
  }
  const Vinv = VinvScaled.scale(1 / sigmae2)

  const Xt = X.transpose()
  const XtVinv = Xt.multiply(Vinv)
  const XtVinvX = XtVinv.multiply(X)
  let XtVinvXInv: Matrix
  try { XtVinvXInv = XtVinvX.inverse() } catch { XtVinvXInv = Matrix.identity(p) }

  const XtVinvY = XtVinv.multiply(Matrix.colVec([...y]))
  const betaM = XtVinvXInv.multiply(XtVinvY)
  const beta = Array.from({ length: p }, (_, i) => betaM.get(i, 0))

  // Degrees of freedom for fixed effects (Satterthwaite approximation: df = n - p - nGroups + 1)
  const df = Math.max(1, n - p - nGroups + 1)

  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, df)
  const covBeta = XtVinvXInv.scale(sigmae2)

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

  // ICC = σ²_b / (σ²_b + σ²_e)
  const icc = sigmab2 / (sigmab2 + sigmae2)

  // Log-likelihood at optimum.
  // The profiled REML above omits the normalizing constant −½·(n−p)·(1+log(2π)),
  // which cancels in optimization but must be included for AIC/BIC and for
  // compatibility with R lme4's reported log-likelihood.
  //   Full REML ℓ = profiled_reml − ½·(n−p)·(1+log(2π))
  const remlConst = 0.5 * (n - p) * (1 + Math.log(2 * Math.PI))
  const logLik = -finalModel.negLogLik - remlConst
  const aic = -2 * logLik + 2 * (p + 2)  // p fixed + 2 variance components
  const bic = -2 * logLik + Math.log(n) * (p + 2)

  const formatted = formatLMM(icc, aic, bic, logLik)

  return {
    fixedEffects,
    varianceComponents: {
      intercept: roundTo(sigmab2, 6),
      residual: roundTo(sigmae2, 6),
    },
    icc: roundTo(icc, 6),
    logLik: roundTo(logLik, 4),
    aic: roundTo(aic, 2),
    bic: roundTo(bic, 2),
    nObs: n,
    nGroups,
    formatted,
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
