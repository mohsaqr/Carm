/**
 * Tests for src/stats/distributions.ts — Distribution fitting module.
 *
 * Cross-validated with R:
 * > dnorm(0, 0, 1)  # normalPDF
 * > dt(0, 1)        # tPDF
 * > dchisq(1, 1)    # chiSqPDF
 * > df(1, 5, 10)    # fPDF
 * > dexp(0, 1)      # exponentialPDF
 * > dunif(0.5, 0, 1) # uniformPDF
 * > dgamma(1, 1, 1) # gammaPDF
 * > dbeta(0.5, 1, 1) # betaPDF
 * > pexp(1, 1), qexp(0.5, 1) # exponentialCDF/Quantile
 * > pgamma(1, 1, 1), qgamma(0.5, 1, 1) # gammaCDF/Quantile
 * > pbeta(0.5, 2, 5), qbeta(0.5, 2, 5) # betaCDF/Quantile
 * > qf(0.5, 5, 10)  # fQuantile
 * > dbinom(3, 10, 0.5), pbinom(3, 10, 0.5), qbinom(0.5, 10, 0.5) # binomial
 * > dpois(3, 5), ppois(3, 5), qpois(0.5, 5) # poisson
 * > MASS::fitdistr(...) # MLE fitting
 * > nortest::ad.test(...) # Anderson-Darling
 * > ks.test(...) # Kolmogorov-Smirnov
 *
 * R script: validation/r-reference/distributions-ref.R
 * Fixture:  tests/fixtures/distributions-ref.json
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

// ─── Helper types for fixture entries ────────────────────────────────────────

interface PDFEntry { x: number; mu?: number; sigma?: number; df?: number; df1?: number; df2?: number; rate?: number; min?: number; max?: number; shape?: number; shape1?: number; shape2?: number; expected: number }
interface CDFEntry { x: number; rate?: number; min?: number; max?: number; shape?: number; shape1?: number; shape2?: number; expected: number }
interface QuantileEntry { p: number; rate?: number; min?: number; max?: number; shape?: number; shape1?: number; shape2?: number; df1?: number; df2?: number; expected: number }
interface DiscreteEntry { k: number; n?: number; prob?: number; lambda?: number; expected: number }

// ─── Section A: Continuous PDFs ──────────────────────────────────────────────

describe('normalPDF', () => {
  const cases = ref.normalPDF as PDFEntry[]
  cases.forEach(({ x, mu, sigma, expected }) => {
    it(`normalPDF(${x}, ${mu}, ${sigma}) ≈ ${expected}`, () => {
      expect(normalPDF(x, mu!, sigma!)).toBeCloseTo(expected, 10)
    })
  })
  it('throws for sigma ≤ 0', () => {
    expect(() => normalPDF(0, 0, 0)).toThrow()
    expect(() => normalPDF(0, 0, -1)).toThrow()
  })
})

describe('tPDF', () => {
  const cases = ref.tPDF as PDFEntry[]
  cases.forEach(({ x, df, expected }) => {
    it(`tPDF(${x}, ${df}) ≈ ${expected}`, () => {
      expect(tPDF(x, df!)).toBeCloseTo(expected, 10)
    })
  })
})

describe('chiSqPDF', () => {
  const cases = ref.chiSqPDF as PDFEntry[]
  cases.forEach(({ x, df, expected }) => {
    it(`chiSqPDF(${x}, ${df}) ≈ ${expected}`, () => {
      expect(chiSqPDF(x, df!)).toBeCloseTo(expected, 8)
    })
  })
  it('returns 0 for x < 0', () => {
    expect(chiSqPDF(-1, 5)).toBe(0)
  })
})

describe('fPDF', () => {
  const cases = ref.fPDF as PDFEntry[]
  cases.forEach(({ x, df1, df2, expected }) => {
    it(`fPDF(${x}, ${df1}, ${df2}) ≈ ${expected}`, () => {
      expect(fPDF(x, df1!, df2!)).toBeCloseTo(expected, 8)
    })
  })
  it('returns 0 for x < 0', () => {
    expect(fPDF(-1, 5, 10)).toBe(0)
  })
})

describe('exponentialPDF', () => {
  const cases = ref.exponentialPDF as PDFEntry[]
  cases.forEach(({ x, rate, expected }) => {
    it(`exponentialPDF(${x}, ${rate}) ≈ ${expected}`, () => {
      expect(exponentialPDF(x, rate!)).toBeCloseTo(expected, 10)
    })
  })
  it('returns 0 for x < 0', () => {
    expect(exponentialPDF(-1, 1)).toBe(0)
  })
})

describe('uniformPDF', () => {
  const cases = ref.uniformPDF as PDFEntry[]
  cases.forEach(({ x, min, max, expected }) => {
    it(`uniformPDF(${x}, ${min}, ${max}) ≈ ${expected}`, () => {
      expect(uniformPDF(x, min!, max!)).toBeCloseTo(expected, 10)
    })
  })
})

describe('gammaPDF', () => {
  const cases = ref.gammaPDF as PDFEntry[]
  cases.forEach(({ x, shape, rate, expected }) => {
    it(`gammaPDF(${x}, ${shape}, ${rate}) ≈ ${expected}`, () => {
      expect(gammaPDF(x, shape!, rate!)).toBeCloseTo(expected, 8)
    })
  })
  it('returns 0 for x < 0', () => {
    expect(gammaPDF(-1, 2, 1)).toBe(0)
  })
})

describe('betaPDF', () => {
  const cases = ref.betaPDF as PDFEntry[]
  cases.forEach(({ x, shape1, shape2, expected }) => {
    it(`betaPDF(${x}, ${shape1}, ${shape2}) ≈ ${expected}`, () => {
      expect(betaPDF(x, shape1!, shape2!)).toBeCloseTo(expected, 8)
    })
  })
  it('returns 0 outside [0,1]', () => {
    expect(betaPDF(-0.1, 2, 5)).toBe(0)
    expect(betaPDF(1.1, 2, 5)).toBe(0)
  })
})

// ─── Section B: CDFs and Quantiles ───────────────────────────────────────────

describe('exponentialCDF', () => {
  const cases = ref.exponentialCDF as CDFEntry[]
  cases.forEach(({ x, rate, expected }) => {
    it(`exponentialCDF(${x}, ${rate}) ≈ ${expected}`, () => {
      expect(exponentialCDF(x, rate!)).toBeCloseTo(expected, 10)
    })
  })
})

describe('exponentialQuantile', () => {
  const cases = ref.exponentialQuantile as QuantileEntry[]
  cases.forEach(({ p, rate, expected }) => {
    it(`exponentialQuantile(${p}, ${rate}) ≈ ${expected}`, () => {
      expect(exponentialQuantile(p, rate!)).toBeCloseTo(expected, 10)
    })
  })
  it('returns Infinity for p=1', () => {
    expect(exponentialQuantile(1, 1)).toBe(Infinity)
  })
})

describe('uniformCDF', () => {
  const cases = ref.uniformCDF as CDFEntry[]
  cases.forEach(({ x, min, max, expected }) => {
    it(`uniformCDF(${x}, ${min}, ${max}) ≈ ${expected}`, () => {
      expect(uniformCDF(x, min!, max!)).toBeCloseTo(expected, 10)
    })
  })
})

describe('uniformQuantile', () => {
  const cases = ref.uniformQuantile as QuantileEntry[]
  cases.forEach(({ p, min, max, expected }) => {
    it(`uniformQuantile(${p}, ${min}, ${max}) ≈ ${expected}`, () => {
      expect(uniformQuantile(p, min!, max!)).toBeCloseTo(expected, 10)
    })
  })
})

describe('gammaCDF', () => {
  const cases = ref.gammaCDF as CDFEntry[]
  cases.forEach(({ x, shape, rate, expected }) => {
    it(`gammaCDF(${x}, ${shape}, ${rate}) ≈ ${expected}`, () => {
      expect(gammaCDF(x, shape!, rate!)).toBeCloseTo(expected, 6)
    })
  })
})

describe('gammaQuantile', () => {
  const cases = ref.gammaQuantile as QuantileEntry[]
  cases.forEach(({ p, shape, rate, expected }) => {
    it(`gammaQuantile(${p}, ${shape}, ${rate}) ≈ ${expected}`, () => {
      expect(gammaQuantile(p, shape!, rate!)).toBeCloseTo(expected, 5)
    })
  })
})

describe('betaCDF', () => {
  const cases = ref.betaCDF as CDFEntry[]
  cases.forEach(({ x, shape1, shape2, expected }) => {
    it(`betaCDF(${x}, ${shape1}, ${shape2}) ≈ ${expected}`, () => {
      expect(betaCDF(x, shape1!, shape2!)).toBeCloseTo(expected, 6)
    })
  })
})

describe('betaQuantile', () => {
  const cases = ref.betaQuantile as QuantileEntry[]
  cases.forEach(({ p, shape1, shape2, expected }) => {
    it(`betaQuantile(${p}, ${shape1}, ${shape2}) ≈ ${expected}`, () => {
      expect(betaQuantile(p, shape1!, shape2!)).toBeCloseTo(expected, 5)
    })
  })
})

describe('fQuantile', () => {
  const cases = ref.fQuantile as QuantileEntry[]
  cases.forEach(({ p, df1, df2, expected }) => {
    it(`fQuantile(${p}, ${df1}, ${df2}) ≈ ${expected}`, () => {
      expect(fQuantile(p, df1!, df2!)).toBeCloseTo(expected, 5)
    })
  })
})

// ─── Section C: Discrete distributions ───────────────────────────────────────

describe('binomialPMF', () => {
  const cases = ref.binomialPMF as DiscreteEntry[]
  cases.forEach(({ k, n, prob, expected }) => {
    it(`binomialPMF(${k}, ${n}, ${prob}) ≈ ${expected}`, () => {
      expect(binomialPMF(k, n!, prob!)).toBeCloseTo(expected, 10)
    })
  })
  it('returns 0 for k out of range', () => {
    expect(binomialPMF(-1, 10, 0.5)).toBe(0)
    expect(binomialPMF(11, 10, 0.5)).toBe(0)
  })
  it('handles edge cases: prob=0, prob=1', () => {
    expect(binomialPMF(0, 5, 0)).toBe(1)
    expect(binomialPMF(1, 5, 0)).toBe(0)
    expect(binomialPMF(5, 5, 1)).toBe(1)
    expect(binomialPMF(4, 5, 1)).toBe(0)
  })
})

describe('binomialCDF', () => {
  const cases = ref.binomialCDF as DiscreteEntry[]
  cases.forEach(({ k, n, prob, expected }) => {
    it(`binomialCDF(${k}, ${n}, ${prob}) ≈ ${expected}`, () => {
      expect(binomialCDF(k, n!, prob!)).toBeCloseTo(expected, 6)
    })
  })
})

describe('binomialQuantile', () => {
  const cases = ref.binomialQuantile as DiscreteEntry[]
  cases.forEach(({ p, n, prob, expected }) => {
    it(`binomialQuantile(${p}, ${n}, ${prob}) = ${expected}`, () => {
      expect(binomialQuantile(p as unknown as number, n!, prob!)).toBe(expected)
    })
  })
})

describe('poissonPMF', () => {
  const cases = ref.poissonPMF as DiscreteEntry[]
  cases.forEach(({ k, lambda, expected }) => {
    it(`poissonPMF(${k}, ${lambda}) ≈ ${expected}`, () => {
      expect(poissonPMF(k, lambda!)).toBeCloseTo(expected, 10)
    })
  })
  it('returns 0 for k < 0', () => {
    expect(poissonPMF(-1, 5)).toBe(0)
  })
  it('handles lambda=0', () => {
    expect(poissonPMF(0, 0)).toBe(1)
    expect(poissonPMF(1, 0)).toBe(0)
  })
})

describe('poissonCDF', () => {
  const cases = ref.poissonCDF as DiscreteEntry[]
  cases.forEach(({ k, lambda, expected }) => {
    it(`poissonCDF(${k}, ${lambda}) ≈ ${expected}`, () => {
      expect(poissonCDF(k, lambda!)).toBeCloseTo(expected, 6)
    })
  })
})

describe('poissonQuantile', () => {
  const cases = ref.poissonQuantile as DiscreteEntry[]
  cases.forEach(({ p, lambda, expected }) => {
    it(`poissonQuantile(${p}, ${lambda}) = ${expected}`, () => {
      expect(poissonQuantile(p as unknown as number, lambda!)).toBe(expected)
    })
  })
})

// ─── Section D: MLE Fitting ─────────────────────────────────────────────────

describe('fitDistribution — normal', () => {
  const { data, mu, sigma, logLik, aic, bic } = ref.fitNormal
  const result = fitDistribution(data, 'normal')

  it('mu matches R MASS::fitdistr', () => {
    expect(result.params['mu']).toBeCloseTo(mu, 6)
  })
  it('sigma matches R (MLE denominator n)', () => {
    expect(result.params['sigma']).toBeCloseTo(sigma, 6)
  })
  it('logLik matches R', () => {
    expect(result.logLik).toBeCloseTo(logLik, 4)
  })
  it('AIC matches R', () => {
    expect(result.aic).toBeCloseTo(aic, 3)
  })
  it('BIC matches R', () => {
    expect(result.bic).toBeCloseTo(bic, 3)
  })
  it('formatted contains distribution name', () => {
    expect(result.formatted).toContain('Normal')
  })
})

describe('fitDistribution — exponential', () => {
  const { data, rate, logLik } = ref.fitExponential
  const result = fitDistribution(data, 'exponential')

  it('rate matches R', () => {
    expect(result.params['rate']).toBeCloseTo(rate, 4)
  })
  it('logLik matches R', () => {
    expect(result.logLik).toBeCloseTo(logLik, 3)
  })
})

describe('fitDistribution — gamma', () => {
  const { data, shape, rate, logLik } = ref.fitGamma
  const result = fitDistribution(data, 'gamma')

  it('shape matches R MASS::fitdistr', () => {
    // Iterative — allow 1e-3 tolerance
    expect(result.params['shape']).toBeCloseTo(shape, 2)
  })
  it('rate matches R', () => {
    expect(result.params['rate']).toBeCloseTo(rate, 2)
  })
  it('logLik matches R', () => {
    expect(result.logLik).toBeCloseTo(logLik, 1)
  })
})

describe('fitDistribution — beta', () => {
  const { data, shape1, shape2, logLik } = ref.fitBeta
  const result = fitDistribution(data, 'beta')

  it('shape1 matches R optim', () => {
    expect(result.params['shape1']).toBeCloseTo(shape1, 2)
  })
  it('shape2 matches R optim', () => {
    expect(result.params['shape2']).toBeCloseTo(shape2, 2)
  })
  it('logLik matches R', () => {
    expect(result.logLik).toBeCloseTo(logLik, 1)
  })
})

describe('fitDistribution — poisson', () => {
  const { data, lambda, logLik } = ref.fitPoisson
  const result = fitDistribution(data, 'poisson')

  it('lambda = mean (closed form)', () => {
    expect(result.params['lambda']).toBeCloseTo(lambda, 10)
  })
  it('logLik matches R', () => {
    expect(result.logLik).toBeCloseTo(logLik, 4)
  })
})

describe('fitDistribution — uniform', () => {
  const data = [0.1, 0.3, 0.5, 0.7, 0.9]
  const result = fitDistribution(data, 'uniform')

  it('min = min(data)', () => {
    expect(result.params['min']).toBe(0.1)
  })
  it('max = max(data)', () => {
    expect(result.params['max']).toBe(0.9)
  })
  it('logLik = -n * log(max - min)', () => {
    expect(result.logLik).toBeCloseTo(-5 * Math.log(0.8), 10)
  })
})

describe('fitDistribution — edge cases', () => {
  it('throws on empty data', () => {
    expect(() => fitDistribution([], 'normal')).toThrow()
  })
  it('throws on negative data for exponential', () => {
    expect(() => fitDistribution([-1, 2, 3], 'exponential')).toThrow()
  })
  it('throws on non-integer data for poisson', () => {
    expect(() => fitDistribution([1.5, 2, 3], 'poisson')).toThrow()
  })
  it('throws on data outside (0,1) for beta', () => {
    expect(() => fitDistribution([0, 0.5, 1], 'beta')).toThrow()
  })
})

// ─── Section E: Goodness-of-Fit Tests ────────────────────────────────────────

describe('andersonDarling — normal data', () => {
  const { data, statistic, pValue } = ref.andersonDarling
  const result = andersonDarling(data)

  it('test statistic matches R nortest::ad.test', () => {
    expect(result.statistic).toBeCloseTo(statistic, 3)
  })
  it('p-value matches R nortest::ad.test', () => {
    expect(result.pValue).toBeCloseTo(pValue, 2)
  })
  it('high p-value (data is normal)', () => {
    expect(result.pValue).toBeGreaterThan(0.05)
  })
  it('testName is Anderson-Darling', () => {
    expect(result.testName).toBe('Anderson-Darling')
  })
  it('formatted string contains A²', () => {
    expect(result.formatted).toContain('A²')
  })
})

describe('andersonDarling — non-normal data', () => {
  const { data, statistic, pValue } = ref.andersonDarlingNonNormal
  const result = andersonDarling(data)

  it('test statistic matches R', () => {
    expect(result.statistic).toBeCloseTo(statistic, 2)
  })
  it('very low p-value (data is exponential, tested as normal)', () => {
    expect(result.pValue).toBeLessThan(0.01)
  })
})

describe('andersonDarling — edge cases', () => {
  it('throws for n < 3', () => {
    expect(() => andersonDarling([1, 2])).toThrow()
  })
})

describe('kolmogorovSmirnov — normal data with known params', () => {
  const { data, statistic, pValue } = ref.kolmogorovSmirnov
  const result = kolmogorovSmirnov(data, 'normal', { mu: 0, sigma: 1 })

  it('D statistic matches R ks.test', () => {
    expect(result.statistic).toBeCloseTo(statistic, 3)
  })
  it('p-value matches R ks.test (asymptotic vs exact, tol=0.05)', () => {
    // Our asymptotic formula differs from R's exact Simard & L'Ecuyer (2011).
    // For n=50 with small D, the difference can be ~0.02-0.03.
    expect(Math.abs(result.pValue - pValue)).toBeLessThan(0.05)
  })
  it('high p-value (data is from the tested distribution)', () => {
    expect(result.pValue).toBeGreaterThan(0.05)
  })
})

describe('kolmogorovSmirnov — exponential data with known params', () => {
  const { data, statistic, pValue } = ref.kolmogorovSmirnovExp
  const result = kolmogorovSmirnov(data, 'exponential', { rate: 1 })

  it('D statistic matches R ks.test', () => {
    expect(result.statistic).toBeCloseTo(statistic, 3)
  })
  it('p-value matches R ks.test (asymptotic vs exact, tol=0.05)', () => {
    expect(Math.abs(result.pValue - pValue)).toBeLessThan(0.05)
  })
})

describe('kolmogorovSmirnov — defaults to normal with MLE params', () => {
  const data = ref.andersonDarling.data
  const result = kolmogorovSmirnov(data)

  it('returns a valid StatResult', () => {
    expect(result.testName).toBe('Kolmogorov-Smirnov')
    expect(result.statistic).toBeGreaterThan(0)
    expect(result.pValue).toBeGreaterThanOrEqual(0)
    expect(result.pValue).toBeLessThanOrEqual(1)
  })
})

// ─── Inverse consistency tests ──────────────────────────────────────────────

describe('CDF-Quantile inverse consistency', () => {
  it('exponential: CDF(Q(p)) ≈ p', () => {
    for (const p of [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]) {
      expect(exponentialCDF(exponentialQuantile(p, 2), 2)).toBeCloseTo(p, 10)
    }
  })

  it('uniform: CDF(Q(p)) ≈ p', () => {
    for (const p of [0.01, 0.1, 0.5, 0.9, 0.99]) {
      expect(uniformCDF(uniformQuantile(p, 3, 7), 3, 7)).toBeCloseTo(p, 10)
    }
  })

  it('gamma: CDF(Q(p)) ≈ p', () => {
    for (const p of [0.1, 0.5, 0.9]) {
      expect(gammaCDF(gammaQuantile(p, 3, 2), 3, 2)).toBeCloseTo(p, 5)
    }
  })

  it('beta: CDF(Q(p)) ≈ p', () => {
    for (const p of [0.1, 0.5, 0.9]) {
      expect(betaCDF(betaQuantile(p, 2, 5), 2, 5)).toBeCloseTo(p, 5)
    }
  })

  it('binomial: CDF(Q(p)) ≥ p and CDF(Q(p)-1) < p', () => {
    const n = 20, prob = 0.3
    for (const p of [0.1, 0.5, 0.9]) {
      const k = binomialQuantile(p, n, prob)
      expect(binomialCDF(k, n, prob)).toBeGreaterThanOrEqual(p - 1e-10)
      if (k > 0) expect(binomialCDF(k - 1, n, prob)).toBeLessThan(p + 1e-10)
    }
  })

  it('poisson: CDF(Q(p)) ≥ p', () => {
    const lambda = 5
    for (const p of [0.1, 0.5, 0.9]) {
      const k = poissonQuantile(p, lambda)
      expect(poissonCDF(k, lambda)).toBeGreaterThanOrEqual(p - 1e-10)
    }
  })
})

// ─── PDF integrates to 1 (numerical check) ──────────────────────────────────

describe('PDF integrates to ≈ 1', () => {
  /** Simple trapezoidal rule integration. */
  function integrate(f: (x: number) => number, lo: number, hi: number, n = 10000): number {
    const dx = (hi - lo) / n
    let sum = (f(lo) + f(hi)) / 2
    for (let i = 1; i < n; i++) sum += f(lo + i * dx)
    return sum * dx
  }

  it('normalPDF integrates to 1', () => {
    expect(integrate(x => normalPDF(x, 0, 1), -8, 8)).toBeCloseTo(1, 4)
  })
  it('tPDF(df=5) integrates to 1', () => {
    expect(integrate(x => tPDF(x, 5), -20, 20)).toBeCloseTo(1, 4)
  })
  it('chiSqPDF(df=3) integrates to 1', () => {
    expect(integrate(x => chiSqPDF(x, 3), 0.001, 30, 20000)).toBeCloseTo(1, 3)
  })
  it('fPDF(5,10) integrates to 1', () => {
    expect(integrate(x => fPDF(x, 5, 10), 0.001, 20, 20000)).toBeCloseTo(1, 3)
  })
  it('exponentialPDF(rate=2) integrates to 1', () => {
    expect(integrate(x => exponentialPDF(x, 2), 0, 20, 20000)).toBeCloseTo(1, 4)
  })
  it('gammaPDF(3, 2) integrates to 1', () => {
    expect(integrate(x => gammaPDF(x, 3, 2), 0.001, 15, 20000)).toBeCloseTo(1, 3)
  })
  it('betaPDF(2, 5) integrates to 1', () => {
    expect(integrate(x => betaPDF(x, 2, 5), 0.001, 0.999, 20000)).toBeCloseTo(1, 3)
  })
})

// ─── Binomial PMF sums to 1 ─────────────────────────────────────────────────

describe('binomial PMF sums to 1', () => {
  it('n=10, p=0.5', () => {
    let sum = 0
    for (let k = 0; k <= 10; k++) sum += binomialPMF(k, 10, 0.5)
    expect(sum).toBeCloseTo(1, 10)
  })
  it('n=20, p=0.25', () => {
    let sum = 0
    for (let k = 0; k <= 20; k++) sum += binomialPMF(k, 20, 0.25)
    expect(sum).toBeCloseTo(1, 10)
  })
})

// ─── Poisson PMF sums to ≈ 1 ────────────────────────────────────────────────

describe('poisson PMF sums to ≈ 1', () => {
  it('lambda=5', () => {
    let sum = 0
    for (let k = 0; k <= 50; k++) sum += poissonPMF(k, 5)
    expect(sum).toBeCloseTo(1, 10)
  })
  it('lambda=0.5', () => {
    let sum = 0
    for (let k = 0; k <= 20; k++) sum += poissonPMF(k, 0.5)
    expect(sum).toBeCloseTo(1, 10)
  })
})
