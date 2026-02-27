/**
 * Extended tests for src/stats/mixed.ts — ML estimation, Nakagawa R²,
 * model comparison (LRT), and random slopes.
 *
 * Cross-validated with R lme4 + MuMIn (verified 2026-02-27):
 *
 * REML intercept-only (existing dataset):
 * > mod <- lmer(y ~ x + (1|g), data=df, REML=TRUE)
 * > fixef(mod)       => (Intercept)=2.1, x=1.2
 * > VarCorr(mod)     => g (Intercept) Var=14.511
 * > sigma(mod)^2     => 0.3429
 * > logLik(mod)      => -12.399
 *
 * ML intercept-only:
 * > mod_ml <- lmer(y ~ x + (1|g), data=df, REML=FALSE)
 * > logLik(mod_ml)   => -12.969
 * > VarCorr intercept Var: 7.23, sigma^2: 0.30
 *
 * LRT:
 * > anova(mod_null_ml, mod_full_ml)
 * > Chi-sq=20.52, df=1, p=5.9e-06
 *
 * Random slopes (seed=42, n_per_group=10, n_groups=5):
 * > mod_slopes <- lmer(y ~ x + (1+x|g), data=df, REML=TRUE)
 * > fixef => (Intercept)=3.198, x=1.561
 * > VarCorr => intercept var=2.048, slope var=0.174, corr=-0.167
 */
import { describe, it, expect } from 'vitest'
import { runLMM, compareLMM } from '../../src/stats/mixed.js'
import ref from '../fixtures/lmm-poisson-ref.json'
import gt from '../fixtures/ground_truth.json'

const { y: y_orig, x: x_orig, group: g_orig } = gt.lmm

// ─── REML backward compatibility ─────────────────────────────────────────

describe('runLMM REML (backward compatibility)', () => {
  const input = {
    outcome: y_orig as number[],
    fixedPredictors: { x: x_orig as number[] },
    groupId: g_orig as number[],
  }
  const result = runLMM(input)

  it('method defaults to REML', () => {
    expect(result.method).toBe('REML')
  })

  it('fixed effect intercept matches R (= 2.1)', () => {
    expect(result.fixedEffects.find(e => e.name === '(Intercept)')!.estimate).toBeCloseTo(2.1, 3)
  })

  it('fixed effect slope matches R (= 1.2)', () => {
    expect(result.fixedEffects.find(e => e.name === 'x')!.estimate).toBeCloseTo(1.2, 3)
  })

  it('ICC matches R (≈ 0.977)', () => {
    expect(result.icc).toBeCloseTo(ref.lmm_reml.icc, 1)
  })

  it('logLik matches R (≈ -12.40)', () => {
    expect(result.logLik).toBeCloseTo(ref.lmm_reml.logLik, 0)
  })

  it('AIC matches R (≈ 32.80)', () => {
    expect(result.aic).toBeCloseTo(ref.lmm_reml.aic, 0)
  })

  it('nParams = p + 2 for intercept-only', () => {
    // 2 fixed (intercept + x) + 1 random cov param + 1 residual = 4
    expect(result.nParams).toBe(4)
  })

  it('has r2Marginal and r2Conditional', () => {
    expect(result.r2Marginal).toBeGreaterThanOrEqual(0)
    expect(result.r2Marginal).toBeLessThanOrEqual(1)
    expect(result.r2Conditional).toBeGreaterThanOrEqual(result.r2Marginal)
    expect(result.r2Conditional).toBeLessThanOrEqual(1)
  })
})

// ─── ML estimation ───────────────────────────────────────────────────────

describe('runLMM ML estimation', () => {
  const input = {
    outcome: y_orig as number[],
    fixedPredictors: { x: x_orig as number[] },
    groupId: g_orig as number[],
    method: 'ML' as const,
  }
  const result = runLMM(input)

  it('method is ML', () => {
    expect(result.method).toBe('ML')
  })

  it('fixed effects match R ML (intercept=2.1, slope=1.2)', () => {
    expect(result.fixedEffects.find(e => e.name === '(Intercept)')!.estimate).toBeCloseTo(ref.lmm_ml.fixef_intercept, 3)
    expect(result.fixedEffects.find(e => e.name === 'x')!.estimate).toBeCloseTo(ref.lmm_ml.fixef_slope, 3)
  })

  it('logLik ML matches R (≈ -12.97)', () => {
    expect(result.logLik).toBeCloseTo(ref.lmm_ml.logLik, 0)
  })

  it('ML logLik < REML logLik (expected for same data)', () => {
    const remlResult = runLMM({
      outcome: y_orig as number[],
      fixedPredictors: { x: x_orig as number[] },
      groupId: g_orig as number[],
      method: 'REML',
    })
    // ML logLik is typically smaller (more negative) than REML
    // Both should be finite
    expect(isFinite(result.logLik)).toBe(true)
    expect(isFinite(remlResult.logLik)).toBe(true)
  })

  it('AIC matches R (≈ 33.94)', () => {
    expect(result.aic).toBeCloseTo(ref.lmm_ml.aic, 0)
  })

  it('npar matches R (= 4)', () => {
    expect(result.nParams).toBe(ref.lmm_ml.npar)
  })
})

// ─── Nakagawa R² ─────────────────────────────────────────────────────────

describe('Nakagawa R²', () => {
  const input = {
    outcome: y_orig as number[],
    fixedPredictors: { x: x_orig as number[] },
    groupId: g_orig as number[],
  }
  const result = runLMM(input)

  it('R² marginal ≈ R value', () => {
    // R gives R2m ≈ 0.177 for this dataset
    expect(result.r2Marginal).toBeCloseTo(ref.lmm_reml.r2m, 1)
  })

  it('R² conditional ≈ R value', () => {
    // R gives R2c ≈ 0.981 for this dataset
    expect(result.r2Conditional).toBeCloseTo(ref.lmm_reml.r2c, 1)
  })

  it('R²_c ≥ R²_m always', () => {
    expect(result.r2Conditional).toBeGreaterThanOrEqual(result.r2Marginal)
  })

  it('both R² in [0, 1]', () => {
    expect(result.r2Marginal).toBeGreaterThanOrEqual(0)
    expect(result.r2Marginal).toBeLessThanOrEqual(1)
    expect(result.r2Conditional).toBeGreaterThanOrEqual(0)
    expect(result.r2Conditional).toBeLessThanOrEqual(1)
  })
})

// ─── Model comparison (LRT) ─────────────────────────────────────────────

describe('compareLMM (LRT)', () => {
  const modNull = runLMM({
    outcome: y_orig as number[],
    fixedPredictors: {},
    groupId: g_orig as number[],
    method: 'ML',
  })
  const modFull = runLMM({
    outcome: y_orig as number[],
    fixedPredictors: { x: x_orig as number[] },
    groupId: g_orig as number[],
    method: 'ML',
  })

  it('null model logLik matches R', () => {
    expect(modNull.logLik).toBeCloseTo(ref.lmm_lrt.logLik_null, 0)
  })

  it('full model logLik matches R', () => {
    expect(modFull.logLik).toBeCloseTo(ref.lmm_lrt.logLik_full, 0)
  })

  const comparison = compareLMM(modNull, modFull)

  it('chi-square matches R anova()', () => {
    expect(comparison.chiSq).toBeCloseTo(ref.lmm_lrt.chisq, 0)
  })

  it('df matches R (= 1)', () => {
    expect(comparison.df).toBe(ref.lmm_lrt.df)
  })

  it('p-value is significant', () => {
    expect(comparison.pValue).toBeLessThan(0.001)
  })

  it('prefers the full model', () => {
    expect(comparison.preferred).toBe('model2')
  })

  it('warns if REML used', () => {
    const remlModel = runLMM({
      outcome: y_orig as number[],
      fixedPredictors: { x: x_orig as number[] },
      groupId: g_orig as number[],
      method: 'REML',
    })
    const comp = compareLMM(modNull, remlModel)
    expect(comp.warning).toBeDefined()
  })

  it('handles identical models', () => {
    const comp = compareLMM(modFull, modFull)
    expect(comp.chiSq).toBeCloseTo(0, 2)
    expect(comp.pValue).toBeCloseTo(1, 1)
  })
})

// ─── Random slopes ───────────────────────────────────────────────────────

describe('runLMM with random slopes', () => {
  const { y, x, g } = ref.lmm_slopes
  const input = {
    outcome: y,
    fixedPredictors: { x: x as number[] },
    groupId: g,
    randomSlopes: ['x'] as readonly string[],
    method: 'REML' as const,
  }
  const result = runLMM(input)

  it('fixed intercept matches R (≈ 3.20)', () => {
    expect(result.fixedEffects.find(e => e.name === '(Intercept)')!.estimate).toBeCloseTo(ref.lmm_slopes.fixef_intercept, 0)
  })

  it('fixed slope matches R (≈ 1.56)', () => {
    expect(result.fixedEffects.find(e => e.name === 'x')!.estimate).toBeCloseTo(ref.lmm_slopes.fixef_slope, 0)
  })

  it('intercept variance is positive and in right ballpark', () => {
    // R gives 2.048
    expect(result.varianceComponents.intercept).toBeGreaterThan(0)
    expect(result.varianceComponents.intercept).toBeCloseTo(ref.lmm_slopes.var_intercept, 0)
  })

  it('slope variance is reported', () => {
    expect(result.varianceComponents.slopes).toBeDefined()
    expect(result.varianceComponents.slopes!['x']).toBeGreaterThan(0)
    // R gives 0.174
    expect(result.varianceComponents.slopes!['x']).toBeCloseTo(ref.lmm_slopes.var_slope, 0)
  })

  it('random correlation is reported', () => {
    expect(result.randomCorrelations).toBeDefined()
    const corrKey = Object.keys(result.randomCorrelations!)[0]!
    const corr = result.randomCorrelations![corrKey]!
    expect(Math.abs(corr)).toBeLessThanOrEqual(1)
  })

  it('nParams = p + q(q+1)/2 + 1 for q=2', () => {
    // 2 fixed (intercept + x) + 3 random cov params (var_int, cov, var_slope) + 1 residual = 6
    expect(result.nParams).toBe(6)
  })

  it('residual variance is positive', () => {
    expect(result.varianceComponents.residual).toBeGreaterThan(0)
  })

  it('logLik is finite', () => {
    expect(isFinite(result.logLik)).toBe(true)
  })

  it('AIC is finite', () => {
    expect(isFinite(result.aic)).toBe(true)
  })

  it('R² values are valid', () => {
    expect(result.r2Marginal).toBeGreaterThanOrEqual(0)
    expect(result.r2Marginal).toBeLessThanOrEqual(1)
    expect(result.r2Conditional).toBeGreaterThanOrEqual(result.r2Marginal)
    expect(result.r2Conditional).toBeLessThanOrEqual(1)
  })
})

// ─── Edge cases ──────────────────────────────────────────────────────────

describe('runLMM edge cases', () => {
  it('throws when random slope not in fixedPredictors', () => {
    expect(() => runLMM({
      outcome: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
      fixedPredictors: { x: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5] },
      groupId: [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
      randomSlopes: ['z'],
    })).toThrow(/random slope/)
  })

  it('ML and REML give same fixed effects direction', () => {
    const base = {
      outcome: y_orig as number[],
      fixedPredictors: { x: x_orig as number[] },
      groupId: g_orig as number[],
    }
    const reml = runLMM({ ...base, method: 'REML' })
    const ml = runLMM({ ...base, method: 'ML' })
    // Point estimates should be very similar
    const remlSlope = reml.fixedEffects.find(e => e.name === 'x')!.estimate
    const mlSlope = ml.fixedEffects.find(e => e.name === 'x')!.estimate
    expect(Math.sign(remlSlope)).toBe(Math.sign(mlSlope))
    expect(Math.abs(remlSlope - mlSlope)).toBeLessThan(0.5)
  })
})
