/**
 * 100-DATASET CROSS-VALIDATION — glmm.ts vs R lme4 + glmmTMB
 *
 * Fits our runGLMM on 100 synthetic datasets with varying:
 *   - n: 50–600 observations
 *   - nGroups: 5–20
 *   - β₀: -2 to +2
 *   - β₁: -1.5 to +1.5
 *   - σ_b: 0.3 to 2.5
 *
 * Compares against two R packages:
 *   - lme4::glmer (BOBYQA, Laplace)
 *   - glmmTMB::glmmTMB (TMB automatic differentiation, Laplace)
 *
 * Reports median and 95th percentile absolute errors across all 100 datasets.
 *
 * R script: validation/r-reference/glmm-100-datasets.R
 * Fixture:  tests/fixtures/glmm-100-ref.json
 */

import { describe, it, expect } from 'vitest'
import { runGLMM } from '../../src/stats/glmm.js'
import fixture from '../fixtures/glmm-100-ref.json'

interface DatasetRef {
  id: number
  params: { n: number; nGroups: number; nPerGroup: number; trueBeta0: number; trueBeta1: number; trueSigmaB: number }
  data: { y: number[]; x: number[]; group: number[] }
  lme4: { intercept: number; x: number; interceptSE: number; xSE: number; variance: number; logLik: number; aic: number }
  glmmTMB: { intercept: number; x: number; interceptSE: number; xSE: number; variance: number; logLik: number; aic: number }
}

const datasets = fixture as DatasetRef[]

// Pre-compute all results once (expensive)
const allResults = datasets.map(ds => {
  try {
    return runGLMM({
      outcome: ds.data.y,
      fixedPredictors: { x: ds.data.x },
      groupId: ds.data.group,
    })
  } catch {
    return null
  }
})

// Collect diffs
const diffs = {
  vsLme4: {
    logLik: [] as number[],
    aic: [] as number[],
    variance: [] as number[],
    intercept: [] as number[],
    xCoef: [] as number[],
    xSE: [] as number[],
  },
  vsTMB: {
    logLik: [] as number[],
    aic: [] as number[],
    variance: [] as number[],
    intercept: [] as number[],
    xCoef: [] as number[],
    xSE: [] as number[],
  },
  lme4VsTMB: {
    logLik: [] as number[],
    intercept: [] as number[],
    xCoef: [] as number[],
    variance: [] as number[],
  },
}

let nFitted = 0
for (let i = 0; i < datasets.length; i++) {
  const r = allResults[i]
  if (!r) continue
  const ds = datasets[i]!
  nFitted++

  // vs lme4
  diffs.vsLme4.logLik.push(Math.abs(r.logLik - ds.lme4.logLik))
  diffs.vsLme4.aic.push(Math.abs(r.aic - ds.lme4.aic))
  diffs.vsLme4.variance.push(Math.abs(r.varianceComponents.intercept - ds.lme4.variance))
  diffs.vsLme4.intercept.push(Math.abs(r.fixedEffects[0]!.estimate - ds.lme4.intercept))
  diffs.vsLme4.xCoef.push(Math.abs(r.fixedEffects[1]!.estimate - ds.lme4.x))
  diffs.vsLme4.xSE.push(Math.abs(r.fixedEffects[1]!.se - ds.lme4.xSE))

  // vs glmmTMB
  diffs.vsTMB.logLik.push(Math.abs(r.logLik - ds.glmmTMB.logLik))
  diffs.vsTMB.aic.push(Math.abs(r.aic - ds.glmmTMB.aic))
  diffs.vsTMB.variance.push(Math.abs(r.varianceComponents.intercept - ds.glmmTMB.variance))
  diffs.vsTMB.intercept.push(Math.abs(r.fixedEffects[0]!.estimate - ds.glmmTMB.intercept))
  diffs.vsTMB.xCoef.push(Math.abs(r.fixedEffects[1]!.estimate - ds.glmmTMB.x))
  diffs.vsTMB.xSE.push(Math.abs(r.fixedEffects[1]!.se - ds.glmmTMB.xSE))

  // lme4 vs glmmTMB (baseline: how much do the two R packages disagree?)
  diffs.lme4VsTMB.logLik.push(Math.abs(ds.lme4.logLik - ds.glmmTMB.logLik))
  diffs.lme4VsTMB.intercept.push(Math.abs(ds.lme4.intercept - ds.glmmTMB.intercept))
  diffs.lme4VsTMB.xCoef.push(Math.abs(ds.lme4.x - ds.glmmTMB.x))
  diffs.lme4VsTMB.variance.push(Math.abs(ds.lme4.variance - ds.glmmTMB.variance))
}

function median(arr: number[]): number {
  const s = [...arr].sort((a, b) => a - b)
  const mid = Math.floor(s.length / 2)
  return s.length % 2 === 0 ? (s[mid - 1]! + s[mid]!) / 2 : s[mid]!
}
function p95(arr: number[]): number {
  const s = [...arr].sort((a, b) => a - b)
  return s[Math.floor(s.length * 0.95)]!
}

describe('100-dataset GLMM cross-validation', () => {

  it('should fit at least 95 of 100 datasets', () => {
    expect(nFitted).toBeGreaterThanOrEqual(95)
  })

  // ── vs lme4 ──────────────────────────────────────────────────────────────

  describe('carm vs lme4 (across 100 datasets)', () => {
    it('median logLik diff < 0.5', () => {
      const m = median(diffs.vsLme4.logLik)
      console.log(`  logLik — median: ${m.toFixed(4)}, p95: ${p95(diffs.vsLme4.logLik).toFixed(4)}`)
      expect(m).toBeLessThan(0.5)
    })

    it('median AIC diff < 1.0', () => {
      const m = median(diffs.vsLme4.aic)
      console.log(`  AIC    — median: ${m.toFixed(4)}, p95: ${p95(diffs.vsLme4.aic).toFixed(4)}`)
      expect(m).toBeLessThan(1.0)
    })

    it('median variance diff < 0.1', () => {
      const m = median(diffs.vsLme4.variance)
      console.log(`  Var    — median: ${m.toFixed(4)}, p95: ${p95(diffs.vsLme4.variance).toFixed(4)}`)
      expect(m).toBeLessThan(0.1)
    })

    it('median x coefficient diff < 0.15', () => {
      const m = median(diffs.vsLme4.xCoef)
      console.log(`  β_x   — median: ${m.toFixed(4)}, p95: ${p95(diffs.vsLme4.xCoef).toFixed(4)}`)
      expect(m).toBeLessThan(0.15)
    })

    it('median x SE diff < 0.05', () => {
      const m = median(diffs.vsLme4.xSE)
      console.log(`  SE_x  — median: ${m.toFixed(4)}, p95: ${p95(diffs.vsLme4.xSE).toFixed(4)}`)
      expect(m).toBeLessThan(0.05)
    })
  })

  // ── vs glmmTMB ───────────────────────────────────────────────────────────

  describe('carm vs glmmTMB (across 100 datasets)', () => {
    it('median logLik diff < 0.5', () => {
      const m = median(diffs.vsTMB.logLik)
      console.log(`  logLik — median: ${m.toFixed(4)}, p95: ${p95(diffs.vsTMB.logLik).toFixed(4)}`)
      expect(m).toBeLessThan(0.5)
    })

    it('median AIC diff < 1.0', () => {
      const m = median(diffs.vsTMB.aic)
      console.log(`  AIC    — median: ${m.toFixed(4)}, p95: ${p95(diffs.vsTMB.aic).toFixed(4)}`)
      expect(m).toBeLessThan(1.0)
    })

    it('median variance diff < 0.1', () => {
      const m = median(diffs.vsTMB.variance)
      console.log(`  Var    — median: ${m.toFixed(4)}, p95: ${p95(diffs.vsTMB.variance).toFixed(4)}`)
      expect(m).toBeLessThan(0.1)
    })

    it('median x coefficient diff < 0.15', () => {
      const m = median(diffs.vsTMB.xCoef)
      console.log(`  β_x   — median: ${m.toFixed(4)}, p95: ${p95(diffs.vsTMB.xCoef).toFixed(4)}`)
      expect(m).toBeLessThan(0.15)
    })
  })

  // ── lme4 vs glmmTMB baseline ─────────────────────────────────────────────

  describe('baseline: lme4 vs glmmTMB disagreement', () => {
    it('reports how much the two R packages disagree', () => {
      console.log(`\n  === R packages vs each other (baseline) ===`)
      console.log(`  logLik — median: ${median(diffs.lme4VsTMB.logLik).toFixed(4)}, p95: ${p95(diffs.lme4VsTMB.logLik).toFixed(4)}`)
      console.log(`  β₀    — median: ${median(diffs.lme4VsTMB.intercept).toFixed(4)}, p95: ${p95(diffs.lme4VsTMB.intercept).toFixed(4)}`)
      console.log(`  β_x   — median: ${median(diffs.lme4VsTMB.xCoef).toFixed(4)}, p95: ${p95(diffs.lme4VsTMB.xCoef).toFixed(4)}`)
      console.log(`  Var    — median: ${median(diffs.lme4VsTMB.variance).toFixed(4)}, p95: ${p95(diffs.lme4VsTMB.variance).toFixed(4)}`)
      expect(true).toBe(true)
    })
  })

  // ── Summary ──────────────────────────────────────────────────────────────

  describe('summary', () => {
    it('prints full comparison table', () => {
      const pad = (s: string, n: number) => s.padStart(n)
      console.log(`\n  ╔═══════════╤═══════════════════════════╤═══════════════════════════╤═══════════════════════════╗`)
      console.log(`  ║           │    carm vs lme4            │    carm vs glmmTMB         │    lme4 vs glmmTMB         ║`)
      console.log(`  ║  Metric   │  median       p95          │  median       p95          │  median       p95          ║`)
      console.log(`  ╠═══════════╪═══════════════════════════╪═══════════════════════════╪═══════════════════════════╣`)

      const rows = [
        ['logLik', diffs.vsLme4.logLik, diffs.vsTMB.logLik, diffs.lme4VsTMB.logLik],
        ['β_x', diffs.vsLme4.xCoef, diffs.vsTMB.xCoef, diffs.lme4VsTMB.xCoef],
        ['Var', diffs.vsLme4.variance, diffs.vsTMB.variance, diffs.lme4VsTMB.variance],
        ['SE_x', diffs.vsLme4.xSE, diffs.vsTMB.xSE, [] as number[]],
      ] as const

      for (const [name, a, b, c] of rows) {
        const ma = a.length ? pad(median(a).toFixed(4), 8) : pad('—', 8)
        const pa = a.length ? pad(p95(a).toFixed(4), 8) : pad('—', 8)
        const mb = b.length ? pad(median(b).toFixed(4), 8) : pad('—', 8)
        const pb = b.length ? pad(p95(b).toFixed(4), 8) : pad('—', 8)
        const mc = c.length ? pad(median(c).toFixed(4), 8) : pad('—', 8)
        const pc = c.length ? pad(p95(c).toFixed(4), 8) : pad('—', 8)
        console.log(`  ║  ${(name as string).padEnd(9)}│${ma}    ${pa}    │${mb}    ${pb}    │${mc}    ${pc}    ║`)
      }
      console.log(`  ╚═══════════╧═══════════════════════════╧═══════════════════════════╧═══════════════════════════╝`)
      console.log(`  n_fitted = ${nFitted} / ${datasets.length}`)
      expect(true).toBe(true)
    })
  })
})
