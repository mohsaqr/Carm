/**
 * LMM diagnostic: deterministic datasets, no random numbers.
 * Run after: npm run build
 * Usage: node tmp/lmm-verify.mjs
 *
 * R cross-validation code at bottom of file.
 */

import { runLMM } from '../dist/index.js'

// ─── Dataset helpers ──────────────────────────────────────────────────────────

function makeDataset(groups, xVals, slopes, offsets, withinNoise, label) {
  const y = [], x = [], g = []
  for (let gi = 0; gi < groups; gi++) {
    for (let i = 0; i < xVals.length; i++) {
      const xi = xVals[i]
      const noise = withinNoise[i % withinNoise.length]
      y.push(slopes[gi] * xi + offsets[gi] + noise)
      x.push(xi)
      g.push(gi + 1)
    }
  }
  return { y, x, g, label }
}

// ─── Dataset definitions (fully deterministic) ───────────────────────────────

/**
 * DS1: Low ICC
 * 3 groups × 10 obs, slope=1, offsets=[0,0,0], within-noise=±2 alternating
 * True σ²_b = 0, σ²_e = 4. True ICC = 0.
 * Expected: ICC ≈ 0, slope ≈ 1
 *
 * R:
 * library(lme4)
 * x <- rep(1:10, 3); g <- rep(1:3, each=10)
 * noise <- rep(c(2,-2), 15)
 * y <- x + noise
 * mod <- lmer(y ~ x + (1|g), REML=TRUE)
 * fixef(mod); icc(mod)  # slope~1, ICC~0
 */
const ds1 = makeDataset(
  3, [1,2,3,4,5,6,7,8,9,10],
  [1,1,1], [0,0,0],
  [2,-2,2,-2,2,-2,2,-2,2,-2],
  'DS1: ICC≈0 (no between-group var)'
)

/**
 * DS2: High ICC
 * 3 groups × 10 obs, slope=1, offsets=[0,10,-10], within-noise=±0.5
 * True σ²_b ≈ 100 (offsets [0,10,-10] → variance = 200/3 ≈ 66.7), σ²_e = 0.25
 * True ICC ≈ 66.7/(66.7+0.25) ≈ 0.996
 * Expected: ICC > 0.9, slope ≈ 1
 *
 * R:
 * offsets <- c(0,10,-10); x <- rep(1:10,3); g <- rep(1:3, each=10)
 * noise <- rep(c(0.5,-0.5), 15)
 * y <- x + rep(offsets, each=10) + noise
 * mod <- lmer(y ~ x + (1|g), REML=TRUE)
 */
const ds2 = makeDataset(
  3, [1,2,3,4,5,6,7,8,9,10],
  [1,1,1], [0,10,-10],
  [0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5],
  'DS2: ICC≈0.99 (large between-group var)'
)

/**
 * DS3: Medium ICC
 * 4 groups × 8 obs, slope=2, offsets=[0,3,-3,0], within-noise=±3
 * σ²_b (from offsets [0,3,-3,0]) = mean=0, var=(0+9+9+0)/4 = 4.5
 * σ²_e = 9
 * True ICC ≈ 4.5/(4.5+9) ≈ 0.33
 * Expected: ICC 0.2–0.5, slope ≈ 2
 *
 * R:
 * offsets <- c(0,3,-3,0); x <- rep(1:8,4); g <- rep(1:4, each=8)
 * noise <- rep(c(3,-3,3,-3,3,-3,3,-3), 4)
 * y <- 2*x + rep(offsets, each=8) + noise
 * mod <- lmer(y ~ x + (1|g), REML=TRUE)
 */
const ds3 = makeDataset(
  4, [1,2,3,4,5,6,7,8],
  [2,2,2,2], [0,3,-3,0],
  [3,-3,3,-3,3,-3,3,-3],
  'DS3: ICC≈0.33 (medium between-group var, slope=2)'
)

/**
 * DS4: Slope recovery (slope=3)
 * 2 groups × 12 obs, slope=3, offsets=[0,5], within-noise=±1
 * σ²_b (from offsets [0,5]) ≈ 12.5, σ²_e = 1
 * True ICC ≈ 12.5/13.5 ≈ 0.93
 * Expected: slope ≈ 3, ICC > 0.8
 *
 * R:
 * x <- rep(1:12,2); g <- rep(1:2, each=12)
 * noise <- rep(c(1,-1),12)
 * y <- 3*x + rep(c(0,5),each=12) + noise
 * mod <- lmer(y ~ x + (1|g), REML=TRUE)
 */
const ds4 = makeDataset(
  2, [1,2,3,4,5,6,7,8,9,10,11,12],
  [3,3], [0,5],
  [1,-1,1,-1,1,-1,1,-1,1,-1,1,-1],
  'DS4: slope=3, ICC≈0.93'
)

/**
 * DS5: OLS only (no random intercept effect)
 * Same groups but identical offsets → should give same result as OLS
 */
const ds5 = makeDataset(
  4, [1,2,3,4,5,6,7,8,9,10],
  [1.5,1.5,1.5,1.5], [0,0,0,0],
  [2,-2,2,-2,2,-2,2,-2,2,-2],
  'DS5: ICC=0 (identical group offsets, slope=1.5)'
)

// ─── OLS reference (simple linear regression) ────────────────────────────────
function ols(y, x) {
  const n = y.length
  const xm = x.reduce((s,v)=>s+v,0)/n
  const ym = y.reduce((s,v)=>s+v,0)/n
  const sxx = x.reduce((s,v)=>s+(v-xm)**2,0)
  const sxy = x.reduce((s,v,i)=>s+(v-xm)*(y[i]-ym),0)
  const slope = sxy/sxx
  const int = ym - slope*xm
  return { slope, int }
}

// ─── Run all datasets ─────────────────────────────────────────────────────────
const datasets = [ds1, ds2, ds3, ds4, ds5]

console.log('='.repeat(72))
console.log('LMM VERIFICATION — Deterministic datasets')
console.log('='.repeat(72))

for (const ds of datasets) {
  console.log(`\n${ds.label}`)
  console.log('-'.repeat(60))

  const olsEst = ols(ds.y, ds.x)
  console.log(`  OLS slope: ${olsEst.slope.toFixed(4)}, OLS intercept: ${olsEst.int.toFixed(4)}`)

  try {
    const r = runLMM({
      outcome: ds.y,
      fixedPredictors: { x: ds.x },
      groupId: ds.g,
    })
    const slopeEff = r.fixedEffects.find(e => e.name === 'x')
    const intEff   = r.fixedEffects.find(e => e.name === '(Intercept)')
    console.log(`  LMM slope: ${slopeEff?.estimate.toFixed(4)},  intercept: ${intEff?.estimate.toFixed(4)}`)
    console.log(`  ICC: ${r.icc.toFixed(4)},  σ²_b: ${r.varianceComponents.intercept.toFixed(4)},  σ²_e: ${r.varianceComponents.residual.toFixed(4)}`)
    console.log(`  AIC: ${r.aic.toFixed(2)},  BIC: ${r.bic.toFixed(2)},  logLik: ${r.logLik.toFixed(4)}`)
    console.log(`  nObs: ${r.nObs}, nGroups: ${r.nGroups}`)
  } catch(e) {
    console.log(`  ERROR: ${e.message}`)
  }
}

console.log('\n' + '='.repeat(72))
console.log('Expected summary:')
console.log('  DS1: ICC≈0, slope≈1 (OLS)')
console.log('  DS2: ICC>0.9, slope≈1')
console.log('  DS3: ICC≈0.3, slope≈2')
console.log('  DS4: ICC>0.8, slope≈3')
console.log('  DS5: ICC≈0, slope≈1.5 (OLS)')
console.log('='.repeat(72))
