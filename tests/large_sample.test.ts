/**
 * Large-sample cross-validation tests for Carm vs R.
 *
 * All datasets are fully deterministic (no random numbers).
 * Expected values generated with R 4.3.2 — see inline R code comments.
 *
 * Larger n surfaces numerical issues that small fixtures hide:
 *  - n=200 arithmetic sequences expose quantile / skewness / kurtosis accuracy
 *  - n=100/group t-test checks Welch-Satterthwaite df at scale
 *  - n=50/group ANOVA with true overlap checks F accuracy
 *  - n=200 regression checks R², adj-R², slope accuracy
 *  - n=300 LMM (10 groups × 30 obs) checks ICC / variance decomposition
 */

import { describe, it, expect } from 'vitest'
import { describe as describeStats, skewness, kurtosis, shapiroWilk } from '../src/stats/descriptive.js'
import {
  tTestIndependent,
  oneWayANOVA,
  mannWhitneyU,
  kruskalWallis,
} from '../src/stats/comparison.js'
import { pearsonCorrelation, spearmanCorrelation } from '../src/stats/correlation.js'
import { linearRegression } from '../src/stats/regression.js'
import { runLMM } from '../src/stats/mixed.js'

// ─── Data generators (deterministic) ─────────────────────────────────────────

function range(start: number, end: number): number[] {
  return Array.from({ length: end - start + 1 }, (_, i) => start + i)
}

// ─── DESCRIPTIVE (n=200, x = 1:200) ──────────────────────────────────────────
// R:
//   x200 <- 1:200
//   mean(x200)                        → 100.5
//   sd(x200)                          → 57.87918
//   var(x200)                         → 3350
//   median(x200)                      → 100.5
//   quantile(x200, 0.25, type=7)      → 50.75
//   quantile(x200, 0.75, type=7)      → 150.25
//   e1071::skewness(x200, type=2)     → 0
//   e1071::kurtosis(x200, type=2)     → -1.2

describe('Descriptive — n=200 arithmetic sequence', () => {
  const x = range(1, 200)
  const r = describeStats(x)

  it('mean = 100.5', () => {
    expect(r.mean).toBeCloseTo(100.5, 6)
  })
  it('sd = 57.879 (R ground truth)', () => {
    expect(r.sd).toBeCloseTo(57.87918, 3)
  })
  it('variance = 3350 (R ground truth)', () => {
    expect(r.variance).toBeCloseTo(3350, 2)
  })
  it('median = 100.5', () => {
    expect(r.median).toBeCloseTo(100.5, 6)
  })
  it('Q1 = 50.75 (R type-7)', () => {
    expect(r.q1).toBeCloseTo(50.75, 4)
  })
  it('Q3 = 150.25 (R type-7)', () => {
    expect(r.q3).toBeCloseTo(150.25, 4)
  })
  it('skewness = 0 (symmetric sequence)', () => {
    expect(Math.abs(skewness(x))).toBeLessThan(0.01)
  })
  it('kurtosis = -1.2 (R e1071 type=2, exact for uniform seq)', () => {
    expect(kurtosis(x)).toBeCloseTo(-1.2, 2)
  })
})

// ─── SHAPIRO-WILK (n=50 subset) ──────────────────────────────────────────────
// R:
//   shapiro.test(1:50)   → W=0.9556, p=0.0581

describe('Shapiro-Wilk — n=50 arithmetic sequence', () => {
  const x50 = range(1, 50)
  const sw = shapiroWilk(x50)

  it('W ≈ 0.9556 (R ground truth, tol=0.01)', () => {
    expect(sw.statistic).toBeCloseTo(0.9556, 2)
  })
  it('p ≈ 0.058 (R ground truth, tol=0.02)', () => {
    expect(sw.pValue).toBeCloseTo(0.0581, 2)
  })
  it('W is in [0,1]', () => {
    expect(sw.statistic).toBeGreaterThan(0)
    expect(sw.statistic).toBeLessThanOrEqual(1)
  })
})

// ─── T-TEST (n=100/group) ─────────────────────────────────────────────────────
// R:
//   g1 <- 1:100; g2 <- 11:110    # shift = 10
//   t.test(g1, g2, var.equal=FALSE)
//     t = -2.437, df = 198, p = 0.01568
//     95% CI: [-18.091, -1.909]
//   effsize::cohen.d(g1, g2)$estimate = -0.3447

describe('t-test independent — n=100/group', () => {
  const g1 = range(1, 100)
  const g2 = range(11, 110)
  const r = tTestIndependent(g1, g2)

  it('t = -2.437 (R ground truth, tol=0.01)', () => {
    expect(r.statistic).toBeCloseTo(-2.4373, 2)
  })
  it('df = 198 (R Welch-Satterthwaite)', () => {
    expect(r.df as number).toBeCloseTo(198, 0)
  })
  it('p = 0.0157 (R ground truth, tol=0.001)', () => {
    expect(r.pValue).toBeCloseTo(0.01568, 3)
  })
  it('95% CI lower ≈ -18.09 (tol=0.1)', () => {
    expect(r.ci[0]).toBeCloseTo(-18.091, 1)
  })
  it('95% CI upper ≈ -1.91 (tol=0.1)', () => {
    expect(r.ci[1]).toBeCloseTo(-1.909, 1)
  })
  it('Cohen d ≈ -0.345 (R effsize, tol=0.01)', () => {
    expect(r.effectSize.value).toBeCloseTo(-0.3447, 2)
  })
})

// ─── ONE-WAY ANOVA (n=50/group, 3 groups) ────────────────────────────────────
// R:
//   gA=1:50; gB=26:75; gC=51:100
//   aov(vals ~ groups)
//     F = 147.06, df = (2, 147), p ≈ 8.4e-36
//     eta² = 0.6668

describe('One-way ANOVA — n=50/group × 3 groups', () => {
  const gA = range(1, 50)
  const gB = range(26, 75)
  const gC = range(51, 100)
  const r = oneWayANOVA([
    { values: gA, label: 'A' },
    { values: gB, label: 'B' },
    { values: gC, label: 'C' },
  ])

  it('F = 147.06 (R ground truth, tol=1)', () => {
    expect(r.statistic).toBeCloseTo(147.06, 0)
  })
  it('df between = 2', () => {
    expect((r.df as readonly [number, number])[0]).toBe(2)
  })
  it('df within = 147', () => {
    expect((r.df as readonly [number, number])[1]).toBe(147)
  })
  it('p < 1e-20 (highly significant)', () => {
    expect(r.pValue).toBeLessThan(1e-20)
  })
  it('omega² ≈ 0.661 (Carm returns ω², not η²; R ω²=0.6607, tol=0.01)', () => {
    // R: eta²=0.6668; omega²=(SS_B - df_B*MS_W)/(SS_T+MS_W)=0.6607
    expect(r.effectSize.value).toBeCloseTo(0.6607, 2)
  })
})

// ─── MANN-WHITNEY (n=100/group) ───────────────────────────────────────────────
// R:
//   wilcox.test(1:100, 11:110)
//     W = 4050, p = 0.02034

describe('Mann-Whitney U — n=100/group', () => {
  const g1 = range(1, 100)
  const g2 = range(11, 110)
  const r = mannWhitneyU(g1, g2)

  it('W = 4050 (R ground truth, tol=1)', () => {
    expect(r.statistic).toBeCloseTo(4050, 0)
  })
  it('p ≈ 0.0203 (R ground truth, tol=0.002)', () => {
    expect(r.pValue).toBeCloseTo(0.02034, 3)
  })
})

// ─── KRUSKAL-WALLIS (n=50/group × 3) ─────────────────────────────────────────
// R:
//   kruskal.test(list(1:50, 26:75, 51:100))
//     H = 101.42, df = 2, p ≈ 9.5e-23

describe('Kruskal-Wallis — n=50/group × 3 groups', () => {
  const gA = range(1, 50)
  const gB = range(26, 75)
  const gC = range(51, 100)
  const r = kruskalWallis([
    { values: gA, label: 'A' },
    { values: gB, label: 'B' },
    { values: gC, label: 'C' },
  ])

  it('H ≈ 101.4 (R ground truth, tol=1)', () => {
    expect(r.statistic).toBeCloseTo(101.42, 0)
  })
  it('df = 2', () => {
    expect(r.df as number).toBe(2)
  })
  it('p < 1e-15 (highly significant)', () => {
    expect(r.pValue).toBeLessThan(1e-15)
  })
})

// ─── PEARSON (n=200) ──────────────────────────────────────────────────────────
// R:
//   x <- 1:200; y <- 2*x + 5 + rep(c(-3,3), 100)
//   cor.test(x, y)
//     r = 0.9996628, t = 541.74, p ≈ 6.36e-316
//     95% CI: [0.99955, 0.99975]

describe('Pearson correlation — n=200', () => {
  const x = range(1, 200)
  const y = x.map((v, i) => 2 * v + 5 + (i % 2 === 0 ? -3 : 3))
  const r = pearsonCorrelation(x, y)

  it('r ≈ 0.9997 (R ground truth, tol=0.0001)', () => {
    expect(r.statistic).toBeCloseTo(0.9996628, 4)
  })
  it('p < 1e-100 (effectively zero)', () => {
    expect(r.pValue).toBeLessThan(1e-100)
  })
  it('95% CI lower > 0.999', () => {
    expect(r.ci[0]).toBeGreaterThan(0.999)
  })
  it('95% CI upper < 1.0', () => {
    expect(r.ci[1]).toBeLessThan(1.0)
  })
})

// ─── SPEARMAN (n=200) ─────────────────────────────────────────────────────────
// R:
//   cor.test(x, y, method="spearman")   (same data)
//     rho = 0.9996677 (ties → normal approx p)

describe('Spearman correlation — n=200', () => {
  const x = range(1, 200)
  const y = x.map((v, i) => 2 * v + 5 + (i % 2 === 0 ? -3 : 3))
  const r = spearmanCorrelation(x, y)

  it('rho ≈ 0.9997 (R ground truth, tol=0.001)', () => {
    expect(r.statistic).toBeCloseTo(0.9997, 3)
  })
  it('p < 1e-100', () => {
    expect(r.pValue).toBeLessThan(1e-100)
  })
})

// ─── LINEAR REGRESSION (n=200) ────────────────────────────────────────────────
// R:
//   lm(y ~ x)  (same data: y = 2x+5 ±3)
//     intercept = 4.9548, slope = 2.0005
//     R² = 0.99933, adj-R² = 0.99932
//     F = 293480, F_p ≈ 6.36e-316

describe('Linear regression — n=200', () => {
  const x = range(1, 200)
  const y = x.map((v, i) => 2 * v + 5 + (i % 2 === 0 ? -3 : 3))
  const r = linearRegression(x, y)

  it('slope ≈ 2.0005 (R ground truth, tol=0.001)', () => {
    const slope = r.coefficients.find(c => c.name === 'x')
    expect(slope?.estimate).toBeCloseTo(2.0005, 3)
  })
  it('intercept ≈ 4.955 (R ground truth, tol=0.01)', () => {
    const intercept = r.coefficients.find(c => c.name === '(Intercept)')
    expect(intercept?.estimate).toBeCloseTo(4.9548, 2)
  })
  it('R² ≈ 0.9993 (R ground truth, tol=0.001)', () => {
    expect(r.r2).toBeCloseTo(0.9993258, 3)
  })
  it('adj-R² ≈ 0.9993 (R ground truth, tol=0.001)', () => {
    expect(r.adjR2).toBeCloseTo(0.9993224, 3)
  })
  it('F p-value < 1e-100', () => {
    expect(r.fPValue).toBeLessThan(1e-100)
  })
  it('residuals sum to zero (OLS property)', () => {
    const sumRes = r.residuals.reduce((s, v) => s + v, 0)
    expect(Math.abs(sumRes)).toBeLessThan(1e-8)
  })
})

// ─── LMM (n=300, 10 groups × 30 obs) ─────────────────────────────────────────
// R:
//   library(lme4); library(performance)
//   offsets <- c(10,5,0,-5,-10,8,-8,3,-3,0)
//   x_lmm <- rep(1:30, 10); g_lmm <- rep(1:10, each=30)
//   noise <- rep(c(2,-2), 150)
//   y_lmm <- 2*x_lmm + rep(offsets, each=30) + noise
//   mod <- lmer(y ~ x + (1|g), data=df_lmm, REML=TRUE)
//   fixef(mod)      → intercept=0.2069, slope=1.9867
//   VarCorr(mod)    → sigma_b²=43.862, sigma_e²=4.138
//   icc(mod)        → ICC_adjusted=0.9138
//   AIC(mod)        → 1344.55
//   logLik(mod)     → -668.274

describe('LMM — n=300 (10 groups × 30 obs)', () => {
  const nGroups = 10
  const nPer = 30
  const offsets = [10, 5, 0, -5, -10, 8, -8, 3, -3, 0]
  const xArr: number[] = []
  const yArr: number[] = []
  const gArr: number[] = []

  for (let gi = 0; gi < nGroups; gi++) {
    for (let i = 1; i <= nPer; i++) {
      const noise = i % 2 === 1 ? 2 : -2
      xArr.push(i)
      yArr.push(2 * i + (offsets[gi] ?? 0) + noise)
      gArr.push(gi + 1)
    }
  }

  const r = runLMM({ outcome: yArr, fixedPredictors: { x: xArr }, groupId: gArr })

  it('slope ≈ 1.987 (R lme4, tol=0.1)', () => {
    const slope = r.fixedEffects.find(e => e.name === 'x')
    expect(slope?.estimate).toBeCloseTo(1.9867, 1)
  })
  it('ICC ≈ 0.914 (R lme4, tol=0.05)', () => {
    expect(r.icc).toBeCloseTo(0.9138, 2)
  })
  it('sigma_b² ≈ 43.86 (R lme4, tol=2)', () => {
    expect(r.varianceComponents.intercept).toBeCloseTo(43.86, 0)
  })
  it('sigma_e² ≈ 4.138 (R lme4, tol=0.5)', () => {
    expect(r.varianceComponents.residual).toBeCloseTo(4.138, 1)
  })
  it('AIC ≈ 1344.5 (R lme4, tol=5)', () => {
    // Full REML ℓ = profiled_reml − ½(n-p)(1+log(2π)); AIC = -2ℓ + 2k
    expect(r.aic).toBeCloseTo(1344.55, 0)
  })
  it('logLik ≈ -668.3 (R lme4, tol=2)', () => {
    expect(r.logLik).toBeCloseTo(-668.274, 0)
  })
  it('nObs = 300', () => {
    expect(r.nObs).toBe(300)
  })
  it('nGroups = 10', () => {
    expect(r.nGroups).toBe(10)
  })
})
