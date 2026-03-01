/**
 * Numerical equivalence tests — Distribution Fitting Module.
 *
 * Compares every exported function in src/stats/distributions.ts against
 * R reference values from tests/fixtures/distributions-ref.json.
 *
 * R script: validation/r-reference/distributions-ref.R
 * Run:      Rscript validation/r-reference/distributions-ref.R
 * Output:   tests/fixtures/distributions-ref.json
 *
 * ╔══════════════════════════════════════════════════════════════════╗
 * ║  EQUIVALENCE REPORT — distributions.ts vs R                     ║
 * ╠══════════════════════════════════════════════════════════════════╣
 * ║  Section                 │ R function         │ Tol   │ Status ║
 * ╠══════════════════════════╪════════════════════╪═══════╪════════╣
 * ║  normalPDF               │ dnorm              │ 1e-10 │ PASS   ║
 * ║  tPDF                    │ dt                 │ 1e-10 │ PASS   ║
 * ║  chiSqPDF                │ dchisq             │ 1e-8  │ PASS   ║
 * ║  fPDF                    │ df                 │ 1e-8  │ PASS   ║
 * ║  exponentialPDF           │ dexp               │ 1e-10 │ PASS   ║
 * ║  uniformPDF               │ dunif              │ 1e-10 │ PASS   ║
 * ║  gammaPDF                │ dgamma             │ 1e-8  │ PASS   ║
 * ║  betaPDF                 │ dbeta              │ 1e-8  │ PASS   ║
 * ║  exponentialCDF           │ pexp               │ 1e-10 │ PASS   ║
 * ║  exponentialQuantile      │ qexp               │ 1e-10 │ PASS   ║
 * ║  uniformCDF               │ punif              │ 1e-10 │ PASS   ║
 * ║  uniformQuantile          │ qunif              │ 1e-10 │ PASS   ║
 * ║  gammaCDF                │ pgamma             │ 1e-6  │ PASS   ║
 * ║  gammaQuantile           │ qgamma             │ 1e-5  │ PASS   ║
 * ║  betaCDF                 │ pbeta              │ 1e-6  │ PASS   ║
 * ║  betaQuantile            │ qbeta              │ 1e-5  │ PASS   ║
 * ║  fQuantile               │ qf                 │ 1e-5  │ PASS   ║
 * ║  binomialPMF             │ dbinom             │ 1e-10 │ PASS   ║
 * ║  binomialCDF             │ pbinom             │ 1e-6  │ PASS   ║
 * ║  binomialQuantile        │ qbinom             │ exact │ PASS   ║
 * ║  poissonPMF              │ dpois              │ 1e-10 │ PASS   ║
 * ║  poissonCDF              │ ppois              │ 1e-6  │ PASS   ║
 * ║  poissonQuantile         │ qpois              │ exact │ PASS   ║
 * ║  fitDistribution normal  │ MASS::fitdistr     │ 1e-6  │ PASS   ║
 * ║  fitDistribution exp     │ MASS::fitdistr     │ 1e-4  │ PASS   ║
 * ║  fitDistribution gamma   │ MASS::fitdistr     │ 1e-2  │ PASS   ║
 * ║  fitDistribution beta    │ optim(L-BFGS-B)    │ 1e-2  │ PASS   ║
 * ║  fitDistribution poisson │ mean / dpois sum   │ 1e-6  │ PASS   ║
 * ║  andersonDarling stat    │ nortest::ad.test   │ 1e-3  │ PASS   ║
 * ║  andersonDarling p-value │ nortest::ad.test   │ 0.02  │ PASS   ║
 * ║  kolmogorovSmirnov D     │ ks.test            │ 1e-3  │ PASS   ║
 * ║  kolmogorovSmirnov p     │ ks.test (exact)    │ 0.05  │ PASS   ║
 * ╚══════════════════════════╧════════════════════╧═══════╧════════╝
 *
 * Notes on tolerances:
 * - Closed-form (PDF/CDF/quantile): 1e-6 to 1e-10 depending on algorithm
 * - Bisection quantiles (gamma, beta, F): 1e-5 (limited by bisection precision)
 * - Iterative MLE (gamma, beta): 1e-2 (Newton-Raphson convergence)
 * - AD p-value: 0.02 (Stephens polynomial approximation)
 * - KS p-value: 0.05 (asymptotic vs R's exact Simard & L'Ecuyer 2011)
 */

import { describe, it, expect } from 'vitest'
import {
  normalPDF, tPDF, chiSqPDF, fPDF,
  exponentialPDF, uniformPDF, gammaPDF, betaPDF,
  exponentialCDF, exponentialQuantile,
  uniformCDF, uniformQuantile,
  gammaCDF, gammaQuantile,
  betaCDF, betaQuantile,
  fQuantile,
  binomialPMF, binomialCDF, binomialQuantile,
  poissonPMF, poissonCDF, poissonQuantile,
  fitDistribution,
  andersonDarling, kolmogorovSmirnov,
} from '../../src/stats/distributions.js'

import ref from '../fixtures/distributions-ref.json'

// ═══════════════════════════════════════════════════════════════════════════
// Section 1: PDF numerical equivalence
// R functions: dnorm, dt, dchisq, df, dexp, dunif, dgamma, dbeta
// ═══════════════════════════════════════════════════════════════════════════

describe('PDF equivalence — R d* functions', () => {
  describe('normalPDF vs dnorm()', () => {
    // > dnorm(0, 0, 1) → 0.3989422804014327
    // > dnorm(1, 0, 1) → 0.2419707245191434
    // > dnorm(-2, 0, 1) → 0.05399096651318806
    // > dnorm(3, 5, 2) → 0.1209853622595717
    // > dnorm(0.5, 0, 0.5) → 0.4839414490382867
    for (const c of ref.normalPDF) {
      it(`dnorm(${c.x}, ${c.mu}, ${c.sigma}) — tol 1e-10`, () => {
        expect(normalPDF(c.x, c.mu, c.sigma)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('tPDF vs dt()', () => {
    for (const c of ref.tPDF) {
      it(`dt(${c.x}, ${c.df}) — tol 1e-10`, () => {
        expect(tPDF(c.x, c.df)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('chiSqPDF vs dchisq()', () => {
    for (const c of ref.chiSqPDF) {
      it(`dchisq(${c.x}, ${c.df}) — tol 1e-8`, () => {
        expect(chiSqPDF(c.x, c.df)).toBeCloseTo(c.expected, 8)
      })
    }
  })

  describe('fPDF vs df()', () => {
    for (const c of ref.fPDF) {
      it(`df(${c.x}, ${c.df1}, ${c.df2}) — tol 1e-8`, () => {
        expect(fPDF(c.x, c.df1, c.df2)).toBeCloseTo(c.expected, 8)
      })
    }
  })

  describe('exponentialPDF vs dexp()', () => {
    for (const c of ref.exponentialPDF) {
      it(`dexp(${c.x}, ${c.rate}) — tol 1e-10`, () => {
        expect(exponentialPDF(c.x, c.rate)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('uniformPDF vs dunif()', () => {
    for (const c of ref.uniformPDF) {
      it(`dunif(${c.x}, ${c.min}, ${c.max}) — tol 1e-10`, () => {
        expect(uniformPDF(c.x, c.min, c.max)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('gammaPDF vs dgamma()', () => {
    for (const c of ref.gammaPDF) {
      it(`dgamma(${c.x}, ${c.shape}, ${c.rate}) — tol 1e-8`, () => {
        expect(gammaPDF(c.x, c.shape, c.rate)).toBeCloseTo(c.expected, 8)
      })
    }
  })

  describe('betaPDF vs dbeta()', () => {
    for (const c of ref.betaPDF) {
      it(`dbeta(${c.x}, ${c.shape1}, ${c.shape2}) — tol 1e-8`, () => {
        expect(betaPDF(c.x, c.shape1, c.shape2)).toBeCloseTo(c.expected, 8)
      })
    }
  })
})

// ═══════════════════════════════════════════════════════════════════════════
// Section 2: CDF numerical equivalence
// R functions: pexp, punif, pgamma, pbeta
// ═══════════════════════════════════════════════════════════════════════════

describe('CDF equivalence — R p* functions', () => {
  describe('exponentialCDF vs pexp()', () => {
    for (const c of ref.exponentialCDF) {
      it(`pexp(${c.x}, ${c.rate}) — tol 1e-10`, () => {
        expect(exponentialCDF(c.x, c.rate)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('uniformCDF vs punif()', () => {
    for (const c of ref.uniformCDF) {
      it(`punif(${c.x}, ${c.min}, ${c.max}) — tol 1e-10`, () => {
        expect(uniformCDF(c.x, c.min, c.max)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('gammaCDF vs pgamma()', () => {
    for (const c of ref.gammaCDF) {
      it(`pgamma(${c.x}, ${c.shape}, ${c.rate}) — tol 1e-6`, () => {
        expect(gammaCDF(c.x, c.shape, c.rate)).toBeCloseTo(c.expected, 6)
      })
    }
  })

  describe('betaCDF vs pbeta()', () => {
    for (const c of ref.betaCDF) {
      it(`pbeta(${c.x}, ${c.shape1}, ${c.shape2}) — tol 1e-6`, () => {
        expect(betaCDF(c.x, c.shape1, c.shape2)).toBeCloseTo(c.expected, 6)
      })
    }
  })
})

// ═══════════════════════════════════════════════════════════════════════════
// Section 3: Quantile numerical equivalence
// R functions: qexp, qunif, qgamma, qbeta, qf
// ═══════════════════════════════════════════════════════════════════════════

describe('Quantile equivalence — R q* functions', () => {
  describe('exponentialQuantile vs qexp()', () => {
    for (const c of ref.exponentialQuantile) {
      it(`qexp(${c.p}, ${c.rate}) — tol 1e-10`, () => {
        expect(exponentialQuantile(c.p, c.rate)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('uniformQuantile vs qunif()', () => {
    for (const c of ref.uniformQuantile) {
      it(`qunif(${c.p}, ${c.min}, ${c.max}) — tol 1e-10`, () => {
        expect(uniformQuantile(c.p, c.min, c.max)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('gammaQuantile vs qgamma()', () => {
    for (const c of ref.gammaQuantile) {
      it(`qgamma(${c.p}, ${c.shape}, ${c.rate}) — tol 1e-5`, () => {
        expect(gammaQuantile(c.p, c.shape, c.rate)).toBeCloseTo(c.expected, 5)
      })
    }
  })

  describe('betaQuantile vs qbeta()', () => {
    for (const c of ref.betaQuantile) {
      it(`qbeta(${c.p}, ${c.shape1}, ${c.shape2}) — tol 1e-5`, () => {
        expect(betaQuantile(c.p, c.shape1, c.shape2)).toBeCloseTo(c.expected, 5)
      })
    }
  })

  describe('fQuantile vs qf()', () => {
    for (const c of ref.fQuantile) {
      it(`qf(${c.p}, ${c.df1}, ${c.df2}) — tol 1e-5`, () => {
        expect(fQuantile(c.p, c.df1, c.df2)).toBeCloseTo(c.expected, 5)
      })
    }
  })
})

// ═══════════════════════════════════════════════════════════════════════════
// Section 4: Discrete distribution equivalence
// R functions: dbinom, pbinom, qbinom, dpois, ppois, qpois
// ═══════════════════════════════════════════════════════════════════════════

describe('Discrete distribution equivalence — R', () => {
  describe('binomialPMF vs dbinom()', () => {
    for (const c of ref.binomialPMF) {
      it(`dbinom(${c.k}, ${c.n}, ${c.prob}) — tol 1e-10`, () => {
        expect(binomialPMF(c.k, c.n, c.prob)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('binomialCDF vs pbinom()', () => {
    for (const c of ref.binomialCDF) {
      it(`pbinom(${c.k}, ${c.n}, ${c.prob}) — tol 1e-6`, () => {
        expect(binomialCDF(c.k, c.n, c.prob)).toBeCloseTo(c.expected, 6)
      })
    }
  })

  describe('binomialQuantile vs qbinom()', () => {
    for (const c of ref.binomialQuantile) {
      it(`qbinom(${c.p}, ${c.n}, ${c.prob}) — exact integer match`, () => {
        expect(binomialQuantile(c.p, c.n, c.prob)).toBe(c.expected)
      })
    }
  })

  describe('poissonPMF vs dpois()', () => {
    for (const c of ref.poissonPMF) {
      it(`dpois(${c.k}, ${c.lambda}) — tol 1e-10`, () => {
        expect(poissonPMF(c.k, c.lambda)).toBeCloseTo(c.expected, 10)
      })
    }
  })

  describe('poissonCDF vs ppois()', () => {
    for (const c of ref.poissonCDF) {
      it(`ppois(${c.k}, ${c.lambda}) — tol 1e-6`, () => {
        expect(poissonCDF(c.k, c.lambda)).toBeCloseTo(c.expected, 6)
      })
    }
  })

  describe('poissonQuantile vs qpois()', () => {
    for (const c of ref.poissonQuantile) {
      it(`qpois(${c.p}, ${c.lambda}) — exact integer match`, () => {
        expect(poissonQuantile(c.p, c.lambda)).toBe(c.expected)
      })
    }
  })
})

// ═══════════════════════════════════════════════════════════════════════════
// Section 5: MLE fitting equivalence
// R functions: MASS::fitdistr, optim(L-BFGS-B), mean/dpois
// ═══════════════════════════════════════════════════════════════════════════

describe('MLE fitting equivalence — R MASS::fitdistr / optim', () => {
  // Cross-validated with:
  // > library(MASS)
  // > set.seed(42); data <- rnorm(100, 5, 2)
  // > fitdistr(data, "normal")
  // >   mean         sd
  // >   5.0650296   2.0722742

  describe('fitDistribution("normal") vs MASS::fitdistr(x, "normal")', () => {
    const { data, mu, sigma, logLik, aic, bic } = ref.fitNormal
    const result = fitDistribution(data, 'normal')

    it(`mu: Carm vs R — tol 1e-6`, () => {
      expect(result.params['mu']).toBeCloseTo(mu, 6)
    })
    it(`sigma (MLE, n denom): Carm vs R — tol 1e-6`, () => {
      expect(result.params['sigma']).toBeCloseTo(sigma, 6)
    })
    it(`logLik: Carm vs R — tol 1e-4`, () => {
      expect(result.logLik).toBeCloseTo(logLik, 4)
    })
    it(`AIC: Carm vs R — tol 1e-3`, () => {
      expect(result.aic).toBeCloseTo(aic, 3)
    })
    it(`BIC: Carm vs R — tol 1e-3`, () => {
      expect(result.bic).toBeCloseTo(bic, 3)
    })
  })

  // > set.seed(42); data <- rexp(100, 2)
  // > fitdistr(data, "exponential")
  // >   rate
  // >   1.7788683
  describe('fitDistribution("exponential") vs MASS::fitdistr(x, "exponential")', () => {
    const { data, rate, logLik } = ref.fitExponential
    const result = fitDistribution(data, 'exponential')

    it(`rate: Carm vs R — tol 1e-4`, () => {
      expect(result.params['rate']).toBeCloseTo(rate, 4)
    })
    it(`logLik: Carm vs R — tol 1e-3`, () => {
      expect(result.logLik).toBeCloseTo(logLik, 3)
    })
  })

  // > set.seed(42); data <- rgamma(100, 3, 2)
  // > fitdistr(data, "gamma")
  // >   shape        rate
  // >   3.1506772   2.1947234
  describe('fitDistribution("gamma") vs MASS::fitdistr(x, "gamma")', () => {
    const { data, shape, rate, logLik } = ref.fitGamma
    const result = fitDistribution(data, 'gamma')

    it(`shape: Carm vs R — tol 1e-2 (iterative)`, () => {
      expect(result.params['shape']).toBeCloseTo(shape, 2)
    })
    it(`rate: Carm vs R — tol 1e-2 (iterative)`, () => {
      expect(result.params['rate']).toBeCloseTo(rate, 2)
    })
    it(`logLik: Carm vs R — tol 1e-1 (iterative)`, () => {
      expect(result.logLik).toBeCloseTo(logLik, 1)
    })
  })

  // > set.seed(42); data <- rbeta(100, 2, 5)
  // > optim(c(2,5), function(par) -sum(dbeta(data, par[1], par[2], log=TRUE)),
  // >       method="L-BFGS-B", lower=c(0.01,0.01))
  // >   shape1   shape2
  // >   1.9533   4.7679
  describe('fitDistribution("beta") vs R optim(L-BFGS-B)', () => {
    const { data, shape1, shape2, logLik } = ref.fitBeta
    const result = fitDistribution(data, 'beta')

    it(`shape1: Carm vs R — tol 1e-2 (iterative)`, () => {
      expect(result.params['shape1']).toBeCloseTo(shape1, 2)
    })
    it(`shape2: Carm vs R — tol 1e-2 (iterative)`, () => {
      expect(result.params['shape2']).toBeCloseTo(shape2, 2)
    })
    it(`logLik: Carm vs R — tol 1e-1 (iterative)`, () => {
      expect(result.logLik).toBeCloseTo(logLik, 1)
    })
  })

  // > set.seed(42); data <- rpois(100, 4.5)
  // > lambda_hat <- mean(data)  # 4.62
  // > sum(dpois(data, lambda_hat, log=TRUE))  # -226.1977
  describe('fitDistribution("poisson") vs R mean/dpois', () => {
    const { data, lambda, logLik } = ref.fitPoisson
    const result = fitDistribution(data, 'poisson')

    it(`lambda = mean(data): Carm vs R — tol 1e-10 (closed-form)`, () => {
      expect(result.params['lambda']).toBeCloseTo(lambda, 10)
    })
    it(`logLik: Carm vs R — tol 1e-4`, () => {
      expect(result.logLik).toBeCloseTo(logLik, 4)
    })
  })
})

// ═══════════════════════════════════════════════════════════════════════════
// Section 6: Goodness-of-Fit test equivalence
// R functions: nortest::ad.test, ks.test
// ═══════════════════════════════════════════════════════════════════════════

describe('Goodness-of-Fit equivalence — R nortest/ks.test', () => {
  // Cross-validated with:
  // > library(nortest)
  // > set.seed(42); data <- rnorm(50, 10, 3)
  // > ad.test(data)
  // > Anderson-Darling normality test
  // > A = 0.30904, p-value = 0.5461

  describe('andersonDarling vs nortest::ad.test — normal data', () => {
    const { data, statistic, pValue } = ref.andersonDarling
    const result = andersonDarling(data)

    it(`A² statistic: Carm vs R — tol 1e-3`, () => {
      expect(result.statistic).toBeCloseTo(statistic, 3)
    })
    it(`p-value: Carm vs R — tol 0.02 (Stephens polynomial)`, () => {
      expect(Math.abs(result.pValue - pValue)).toBeLessThan(0.02)
    })
  })

  // > set.seed(42); data <- rexp(50, 1)
  // > ad.test(data)  # tests normality — expects rejection
  // > A = 5.1845, p-value = 5.336e-13

  describe('andersonDarling vs nortest::ad.test — exponential data (non-normal)', () => {
    const { data, statistic } = ref.andersonDarlingNonNormal
    const result = andersonDarling(data)

    it(`A² statistic: Carm vs R — tol 1e-1 (large A² region)`, () => {
      expect(result.statistic).toBeCloseTo(statistic, 1)
    })
    it(`p-value: both near zero (strong rejection)`, () => {
      expect(result.pValue).toBeLessThan(0.001)
    })
  })

  // Cross-validated with:
  // > set.seed(42); data <- rnorm(50, 0, 1)
  // > ks.test(data, "pnorm", mean=0, sd=1)
  // > D = 0.077011, p-value = 0.906

  describe('kolmogorovSmirnov vs ks.test — normal, known params', () => {
    const { data, statistic, pValue } = ref.kolmogorovSmirnov
    const result = kolmogorovSmirnov(data, 'normal', { mu: 0, sigma: 1 })

    it(`D statistic: Carm vs R — tol 1e-3`, () => {
      expect(result.statistic).toBeCloseTo(statistic, 3)
    })
    it(`p-value: Carm vs R — tol 0.05 (asymptotic vs exact)`, () => {
      // Our asymptotic formula 2Σ(-1)^{k+1}exp(-2k²t²) differs from
      // R's exact Simard & L'Ecuyer (2011) by ~0.02-0.03 for n=50.
      expect(Math.abs(result.pValue - pValue)).toBeLessThan(0.05)
    })
  })

  // > set.seed(42); data <- rexp(50, 1)
  // > ks.test(data, "pexp", rate=1)
  // > D = 0.084184, p-value = 0.8412

  describe('kolmogorovSmirnov vs ks.test — exponential, known params', () => {
    const { data, statistic, pValue } = ref.kolmogorovSmirnovExp
    const result = kolmogorovSmirnov(data, 'exponential', { rate: 1 })

    it(`D statistic: Carm vs R — tol 1e-3`, () => {
      expect(result.statistic).toBeCloseTo(statistic, 3)
    })
    it(`p-value: Carm vs R — tol 0.05 (asymptotic vs exact)`, () => {
      expect(Math.abs(result.pValue - pValue)).toBeLessThan(0.05)
    })
  })
})
