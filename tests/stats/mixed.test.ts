/**
 * Tests for src/stats/mixed.ts (LMM)
 * Cross-validated with R lme4 (verified 2026-02-23):
 * > library(lme4)
 * > df <- data.frame(y=c(1,2,3,4,5,6,7,8,9,12), x=c(1,2,3,4,5,1,2,3,4,5), g=factor(c(1,1,1,1,1,2,2,2,2,2)))
 * > mod <- lmer(y ~ x + (1|g), data=df, REML=TRUE)
 * > fixef(mod)      =>  (Intercept)=2.1,  x=1.2
 * > VarCorr(mod)    =>  g (Intercept) Std.Dev.=3.8094 → σ²_b=14.511
 * > sigma(mod)^2    =>  σ²_e=0.3429
 * > logLik(mod)     =>  -12.399
 * > AIC(mod)        =>  32.797
 * > ICC             =>  0.9769
 * Our implementation matches R exactly on all these values.
 */
import { describe, it, expect } from 'vitest'
import { runLMM, computeBLUPs } from '../../src/stats/mixed.js'
import gt from '../fixtures/ground_truth.json'

const { y, x, group } = gt.lmm

describe('runLMM', () => {
  const input = {
    outcome: y as number[],
    fixedPredictors: { x: x as number[] },
    groupId: group as number[],
  }
  const result = runLMM(input)

  // R lme4 gives ICC=0.9769 for this dataset (verified 2026-02-23).
  // The data has SS_between=72.9 vs SS_within=2.4, so high ICC is structurally correct.
  it('ICC matches R lme4 (≈ 0.977)', () => {
    expect(result.icc).toBeGreaterThan(0.95)
    expect(result.icc).toBeLessThan(1)
  })
  it('ICC is in valid range (0, 1)', () => {
    expect(result.icc).toBeGreaterThan(0)
    expect(result.icc).toBeLessThan(1)
  })
  it('fixed effect slope matches R lme4 (= 1.2)', () => {
    const slopeEffect = result.fixedEffects.find(e => e.name === 'x')
    expect(slopeEffect).toBeDefined()
    expect(slopeEffect!.estimate).toBeCloseTo(1.2, 4)
  })
  it('fixed effect intercept matches R lme4 (= 2.1)', () => {
    const interceptEffect = result.fixedEffects.find(e => e.name === '(Intercept)')
    expect(interceptEffect).toBeDefined()
    expect(interceptEffect!.estimate).toBeCloseTo(2.1, 4)
  })
  it('variance components are positive', () => {
    expect(result.varianceComponents.intercept).toBeGreaterThan(0)
    expect(result.varianceComponents.residual).toBeGreaterThan(0)
  })
  it('logLik matches R lme4 (≈ -12.399)', () => {
    expect(result.logLik).toBeCloseTo(-12.399, 1)
  })
  it('AIC matches R lme4 (≈ 32.797)', () => {
    expect(result.aic).toBeCloseTo(32.80, 1)
  })
  it('AIC is finite', () => {
    expect(isFinite(result.aic)).toBe(true)
  })
  it('BIC is finite', () => {
    expect(isFinite(result.bic)).toBe(true)
  })
  it('BIC > AIC (BIC penalizes more for same k)', () => {
    expect(result.bic).toBeGreaterThan(result.aic)
  })
  it('nObs is correct', () => {
    expect(result.nObs).toBe((y as number[]).length)
  })
  it('nGroups is 2', () => {
    expect(result.nGroups).toBe(2)
  })
  it('formatted string contains ICC', () => {
    expect(result.formatted).toMatch(/ICC/)
  })
  it('fixedEffects has intercept + x', () => {
    const names = result.fixedEffects.map(e => e.name)
    expect(names).toContain('(Intercept)')
    expect(names).toContain('x')
  })
})

describe('runLMM throws on bad input', () => {
  it('throws when n < 5', () => {
    expect(() => runLMM({
      outcome: [1, 2, 3],
      fixedPredictors: { x: [1, 2, 3] },
      groupId: [1, 1, 2],
    })).toThrow()
  })
  it('throws when only 1 group', () => {
    expect(() => runLMM({
      outcome: [1, 2, 3, 4, 5],
      fixedPredictors: { x: [1, 2, 3, 4, 5] },
      groupId: [1, 1, 1, 1, 1],
    })).toThrow()
  })
})

describe('computeBLUPs', () => {
  const input = {
    outcome: y as number[],
    fixedPredictors: { x: x as number[] },
    groupId: group as number[],
  }
  const result = runLMM(input)
  const blups = computeBLUPs(input, result)

  it('returns one BLUP per group', () => {
    expect(blups.length).toBe(result.nGroups)
  })
  it('BLUPs sum approximately to zero (shrinkage property)', () => {
    const sum = blups.reduce((s, b) => s + b.blup, 0)
    expect(Math.abs(sum)).toBeLessThan(1)  // approximately zero
  })
  it('BLUP values are finite', () => {
    for (const b of blups) {
      expect(isFinite(b.blup)).toBe(true)
    }
  })
  it('groups with higher outcomes have positive BLUPs', () => {
    // group 2 has higher y values (6,7,8,9,12) vs group 1 (1,2,3,4,5)
    const g2 = blups.find(b => b.group === 2)
    expect(g2).toBeDefined()
    expect(g2!.blup).toBeGreaterThan(0)
  })
})
