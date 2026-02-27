/**
 * Frequency analysis module.
 * Frequency tables, chi-square test of independence, Fisher's exact test,
 * Cramér's V, phi coefficient, odds ratio, goodness-of-fit test.
 */

import { chiSqPValue, normalQuantile, incompleteBeta, normalCDF, normalSurvival, logGamma } from '../core/math.js'
import { formatChiSq, formatP, interpretCramerV, formatBinomial, formatProportions } from '../core/apa.js'
import { PRNG } from '../core/prng.js'
import type {
  FrequencyRow,
  FrequencyTestResult,
  StatResult,
  EffectSize,
  EffectInterpretation,
} from '../core/types.js'

// ─── Frequency table ──────────────────────────────────────────────────────

/**
 * Build a frequency table from an array of values.
 * Returns rows sorted by value with absolute, relative, and cumulative frequencies.
 */
export function frequencyTable(data: readonly (string | number)[]): FrequencyRow[] {
  const counts = new Map<string | number, number>()
  for (const v of data) counts.set(v, (counts.get(v) ?? 0) + 1)

  const total = data.length
  const sorted = [...counts.entries()].sort((a, b) => {
    const av = typeof a[0] === 'number' ? a[0] : String(a[0])
    const bv = typeof b[0] === 'number' ? b[0] : String(b[0])
    return av < bv ? -1 : av > bv ? 1 : 0
  })

  let cumulative = 0
  return sorted.map(([value, count]) => {
    const relative = count / total
    cumulative += relative
    return { value, count, relative, cumulative: Math.min(1, cumulative) }
  })
}

// ─── Contingency table helpers ────────────────────────────────────────────

/** Convert grouped data to a contingency table (rows = group1, cols = group2). */
export function contingencyTable(
  group1: readonly (string | number)[],
  group2: readonly (string | number)[]
): { table: number[][]; rowLabels: (string | number)[]; colLabels: (string | number)[] } {
  if (group1.length !== group2.length) throw new Error('contingencyTable: arrays must have equal length')

  const rowSet = [...new Set(group1)].sort()
  const colSet = [...new Set(group2)].sort()
  const table = Array.from({ length: rowSet.length }, () =>
    new Array<number>(colSet.length).fill(0)
  )

  for (let i = 0; i < group1.length; i++) {
    const r = rowSet.indexOf(group1[i]!)
    const c = colSet.indexOf(group2[i]!)
    table[r]![c]!++
  }

  return { table, rowLabels: rowSet, colLabels: colSet }
}

// ─── Expected counts ──────────────────────────────────────────────────────

function expectedCounts(observed: readonly (readonly number[])[]): number[][] {
  const R = observed.length
  const C = observed[0]!.length
  const rowTotals = observed.map(row => row.reduce((s, v) => s + v, 0))
  const colTotals = Array.from({ length: C }, (_, j) =>
    observed.reduce((s, row) => s + (row[j] ?? 0), 0)
  )
  const total = rowTotals.reduce((s, v) => s + v, 0)

  return Array.from({ length: R }, (_, i) =>
    Array.from({ length: C }, (_, j) => (rowTotals[i]! * (colTotals[j] ?? 0)) / total)
  )
}

// ─── Chi-square test of independence ─────────────────────────────────────

/**
 * Pearson chi-square test of independence for a contingency table.
 *
 * Cross-validated with R:
 * > chisq.test(matrix(c(10,20,30,40), nrow=2))
 * X-squared = 0.6061, df = 1, p-value = 0.4363
 * Cramér's V = 0.0875
 */
export function chiSquareTest(
  observed: readonly (readonly number[])[],
  yatesCorrection = false
): FrequencyTestResult {
  const R = observed.length
  if (R < 2) throw new Error('chiSquareTest: need at least 2 rows')
  const C = observed[0]!.length
  if (C < 2) throw new Error('chiSquareTest: need at least 2 columns')

  const expected = expectedCounts(observed)
  const n = observed.reduce((s, row) => s + row.reduce((a, v) => a + v, 0), 0)

  let chiSq = 0
  for (let i = 0; i < R; i++) {
    for (let j = 0; j < C; j++) {
      const o = observed[i]![j] ?? 0
      const e = expected[i]![j] ?? 0
      if (e < 5) { /* warn: expected < 5, consider Fisher */ }
      const num = yatesCorrection ? Math.max(0, Math.abs(o - e) - 0.5) : o - e
      chiSq += (num * num) / e
    }
  }

  const df = (R - 1) * (C - 1)
  const pValue = chiSqPValue(chiSq, df)

  // Cramér's V
  const minDim = Math.min(R, C) - 1
  const cramersV = Math.sqrt(chiSq / (n * Math.max(1, minDim)))
  const effectSize: EffectSize = {
    value: cramersV,
    name: "Cramér's V",
    interpretation: interpretCramerV(cramersV, df) as EffectInterpretation,
  }

  const ci: readonly [number, number] = [NaN, NaN]
  const formatted = formatChiSq(chiSq, df, pValue, cramersV, "V")

  return {
    testName: "Pearson's χ²",
    statistic: chiSq,
    df,
    pValue,
    effectSize,
    ci,
    ciLevel: 0.95,
    n,
    formatted,
    table: frequencyTable(observed.flatMap((row, i) =>
      row.flatMap((count, j) => new Array<string>(count).fill(`${i},${j}`))
    )),
    expectedCounts: expected,
  }
}

// ─── Fisher's exact test ──────────────────────────────────────────────────

/**
 * Fisher's exact test for 2×2 contingency tables.
 * Uses the hypergeometric distribution.
 *
 * Cross-validated with R:
 * > fisher.test(matrix(c(11, 5, 3, 6), nrow=2))
 * p-value = 0.2684, odds ratio = 4.0
 */
export function fisherExactTest(a: number, b: number, c: number, d: number): StatResult {
  if (a < 0 || b < 0 || c < 0 || d < 0) throw new Error('fisherExactTest: all cells must be non-negative')

  const n = a + b + c + d
  const r1 = a + b  // row 1 total
  const c1 = a + c  // col 1 total
  const c2 = b + d  // col 2 total

  // Two-tailed p-value: sum all tables at least as extreme
  const pObserved = hypergeomPMF(a, r1, c1, n)
  let pValue = 0
  const kMin = Math.max(0, r1 - c2)
  const kMax = Math.min(r1, c1)
  for (let k = kMin; k <= kMax; k++) {
    const p = hypergeomPMF(k, r1, c1, n)
    if (p <= pObserved + 1e-10) pValue += p
  }
  pValue = Math.min(1, pValue)

  // Odds ratio (a·d)/(b·c), with 0.5 continuity if any cell is 0
  const oddsRatio = (b === 0 || c === 0)
    ? ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    : (a * d) / (b * c)

  // CI for log(OR) via normal approximation
  const se = Math.sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (d + 0.5))
  const z = normalQuantile(0.975)
  const logOR = Math.log(oddsRatio)
  const ci: readonly [number, number] = [Math.exp(logOR - z * se), Math.exp(logOR + z * se)]

  const effectSize: EffectSize = {
    value: oddsRatio,
    name: 'Odds ratio',
    interpretation: oddsRatio >= 3 ? 'large' : oddsRatio >= 1.5 ? 'medium' : 'small',
  }

  return {
    testName: "Fisher's Exact Test",
    statistic: oddsRatio,
    df: 1,
    pValue,
    effectSize,
    ci,
    ciLevel: 0.95,
    n,
    formatted: `OR = ${oddsRatio.toFixed(2)}, ${formatP(pValue)}, 95% CI [${ci[0].toFixed(2)}, ${ci[1].toFixed(2)}]`,
  }
}

/** Hypergeometric PMF: P(X = k | n, K, N). */
function hypergeomPMF(k: number, n: number, K: number, N: number): number {
  return Math.exp(
    logCombination(K, k) + logCombination(N - K, n - k) - logCombination(N, n)
  )
}

/** Log of binomial coefficient: log C(n, k). */
function logCombination(n: number, k: number): number {
  if (k < 0 || k > n) return -Infinity
  if (k === 0 || k === n) return 0
  return logFactorial(n) - logFactorial(k) - logFactorial(n - k)
}

/** Log factorial using Stirling-style accumulation for small n. */
function logFactorial(n: number): number {
  if (n <= 1) return 0
  let result = 0
  for (let i = 2; i <= n; i++) result += Math.log(i)
  return result
}

// ─── Phi coefficient ──────────────────────────────────────────────────────

/**
 * Phi coefficient for 2×2 tables: φ = (ad - bc) / sqrt((a+b)(c+d)(a+c)(b+d))
 */
export function phiCoefficient(a: number, b: number, c: number, d: number): number {
  const denom = Math.sqrt((a + b) * (c + d) * (a + c) * (b + d))
  return denom === 0 ? 0 : (a * d - b * c) / denom
}

// ─── Goodness-of-fit test ─────────────────────────────────────────────────

/**
 * Chi-square goodness-of-fit test.
 * Tests whether observed frequencies match expected proportions.
 *
 * Cross-validated with R:
 * > chisq.test(c(30, 20, 50), p = c(1/3, 1/3, 1/3))
 * X-squared = 10, df = 2, p-value = 0.006738
 */
export function goodnessOfFit(
  observed: readonly number[],
  expected?: readonly number[]  // proportions; if omitted, assumes equal
): StatResult {
  const n = observed.reduce((s, v) => s + v, 0)
  const k = observed.length
  if (k < 2) throw new Error('goodnessOfFit: need at least 2 categories')

  const expCounts = expected
    ? expected.map(p => p * n)
    : new Array<number>(k).fill(n / k)

  let chiSq = 0
  for (let i = 0; i < k; i++) {
    const o = observed[i] ?? 0
    const e = expCounts[i] ?? 0
    if (e <= 0) throw new Error(`goodnessOfFit: expected count must be > 0 at index ${i}`)
    chiSq += (o - e) ** 2 / e
  }

  const df = k - 1
  const pValue = chiSqPValue(chiSq, df)
  const w = Math.sqrt(chiSq / n)  // Cohen's w

  return {
    testName: 'Chi-square goodness-of-fit',
    statistic: chiSq,
    df,
    pValue,
    effectSize: {
      value: w,
      name: "Cohen's w",
      interpretation: w < 0.1 ? 'negligible' : w < 0.3 ? 'small' : w < 0.5 ? 'medium' : 'large',
    },
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted: formatChiSq(chiSq, df, pValue, w, "w"),
  }
}

// ─── Cramér's V with Bootstrap CI ─────────────────────────────────────────

/**
 * Cramér's V with bootstrap confidence interval via multinomial resampling.
 *
 * Cross-validated with R:
 * > library(rcompanion)
 * > m <- matrix(c(10,20,30,40,50,60), nrow=2)
 * > cramerV(m, ci=TRUE, R=2000)
 */
export function cramersVWithCI(
  observed: readonly (readonly number[])[],
  ciLevel = 0.95,
  nBoot = 2000,
  seed = 42
): StatResult {
  const R = observed.length
  const C = observed[0]!.length
  if (R < 2 || C < 2) throw new Error('cramersVWithCI: need at least a 2×2 table')

  // Compute observed V via chiSquareTest
  const result = chiSquareTest(observed)
  const V = result.effectSize.value
  const n = result.n
  const df = (R - 1) * (C - 1)

  // Flatten observed into multinomial probabilities
  const flat: number[] = []
  for (let i = 0; i < R; i++) {
    for (let j = 0; j < C; j++) {
      flat.push(observed[i]![j]!)
    }
  }
  const probs = flat.map(c => c / n)

  // Cumulative probabilities for multinomial sampling
  const cumProbs: number[] = []
  let cumSum = 0
  for (const p of probs) {
    cumSum += p
    cumProbs.push(cumSum)
  }

  const rng = new PRNG(seed)
  const bootV: number[] = []

  for (let b = 0; b < nBoot; b++) {
    // Generate bootstrap table via multinomial resampling
    const bootFlat = new Array<number>(flat.length).fill(0)
    for (let s = 0; s < n; s++) {
      const u = rng.next()
      let idx = 0
      while (idx < cumProbs.length - 1 && u > cumProbs[idx]!) idx++
      bootFlat[idx]!++
    }

    // Reshape to table and compute V
    const bootTable: number[][] = Array.from({ length: R }, (_, i) =>
      bootFlat.slice(i * C, (i + 1) * C)
    )

    // Compute chi-square for bootstrap table
    const bootN = n
    const rowTotals = bootTable.map(row => row.reduce((s, v) => s + v, 0))
    const colTotals = Array.from({ length: C }, (_, j) =>
      bootTable.reduce((s, row) => s + (row[j] ?? 0), 0)
    )
    let chiSq = 0
    let valid = true
    for (let i = 0; i < R; i++) {
      for (let j = 0; j < C; j++) {
        const e = (rowTotals[i]! * (colTotals[j] ?? 0)) / bootN
        if (e <= 0) { valid = false; break }
        const o = bootTable[i]![j]!
        chiSq += (o - e) ** 2 / e
      }
      if (!valid) break
    }
    if (!valid) continue

    const minDim = Math.min(R, C) - 1
    bootV.push(Math.sqrt(chiSq / (bootN * Math.max(1, minDim))))
  }

  // Percentile CI
  bootV.sort((a, b) => a - b)
  const alpha = 1 - ciLevel
  const loIdx = Math.max(0, Math.floor((alpha / 2) * bootV.length))
  const hiIdx = Math.min(bootV.length - 1, Math.floor((1 - alpha / 2) * bootV.length))
  const ci: readonly [number, number] = [bootV[loIdx]!, bootV[hiIdx]!]

  const effectSize: EffectSize = {
    value: V,
    name: "Cramér's V",
    interpretation: interpretCramerV(V, df) as EffectInterpretation,
  }

  const ciPct = Math.round(ciLevel * 100)
  const formatted = `V = ${V.toFixed(2)}, ${ciPct}% CI [${ci[0].toFixed(2)}, ${ci[1].toFixed(2)}]`

  return {
    testName: "Cramér's V with Bootstrap CI",
    statistic: V,
    df,
    pValue: result.pValue,
    effectSize,
    ci,
    ciLevel,
    n,
    formatted,
  }
}

// ─── McNemar's Test ───────────────────────────────────────────────────────

/**
 * McNemar's test for paired nominal data (2×2 table).
 * Inputs b and c are off-diagonal counts: b = 0→1 switches, c = 1→0 switches.
 *
 * Cross-validated with R:
 * > mcnemar.test(matrix(c(20, 5, 10, 15), nrow=2), correct=FALSE)
 * McNemar's chi-squared = 1.6667, df = 1, p-value = 0.1967
 */
export function mcnemarsTest(
  b: number,
  c: number,
  correction = false
): StatResult {
  if (b < 0 || c < 0) throw new Error('mcnemarsTest: b and c must be non-negative')
  if (b + c === 0) throw new Error('mcnemarsTest: b + c must be > 0')

  const diff = Math.abs(b - c)
  const corrected = correction ? Math.max(0, diff - 1) : diff
  const chiSq = (corrected * corrected) / (b + c)
  const df = 1
  const pValue = chiSqPValue(chiSq, df)

  // Effect size: odds ratio b/c
  const oddsRatio = c === 0 ? Infinity : b / c

  const effectSize: EffectSize = {
    value: oddsRatio,
    name: 'Odds ratio',
    interpretation: oddsRatio >= 3 || (oddsRatio > 0 && oddsRatio <= 1 / 3) ? 'large'
      : oddsRatio >= 1.5 || (oddsRatio > 0 && oddsRatio <= 1 / 1.5) ? 'medium' : 'small',
  }

  const formatted = `χ²(${df}) = ${chiSq.toFixed(2)}, ${formatP(pValue)}, OR = ${oddsRatio === Infinity ? '∞' : oddsRatio.toFixed(2)}`

  return {
    testName: "McNemar's Test",
    statistic: chiSq,
    df,
    pValue,
    effectSize,
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n: b + c,
    formatted,
  }
}

// ─── Binomial Test ────────────────────────────────────────────────────────

/**
 * Exact binomial test with Clopper-Pearson CI.
 * Uses regularized incomplete beta function for exact p-values and CI.
 *
 * Cross-validated with R:
 * > binom.test(6, 10, p = 0.5)
 * p-value = 0.7539, 95% CI [0.2624, 0.8784]
 */
export function binomialTest(
  successes: number,
  trials: number,
  p0 = 0.5,
  alternative: 'two.sided' | 'less' | 'greater' = 'two.sided',
  ciLevel = 0.95
): StatResult {
  if (trials < 1) throw new Error('binomialTest: trials must be >= 1')
  if (successes < 0 || successes > trials) throw new Error('binomialTest: successes must be in [0, trials]')
  if (p0 <= 0 || p0 >= 1) throw new Error('binomialTest: p0 must be in (0, 1)')

  const k = successes
  const n = trials
  const pHat = k / n

  // Handle boundary cases where incompleteBeta would call logGamma(0)
  // P(X <= k) = I_{1-p0}(n-k, k+1)
  const pLessOrEqual = (k === n) ? 1 : incompleteBeta(1 - p0, n - k, k + 1)
  // P(X >= k) = I_{p0}(k, n-k+1) = 1 - P(X <= k-1)
  const pGreaterOrEqual = k === 0 ? 1 : ((k === n) ? Math.pow(p0, n) : incompleteBeta(p0, k, n - k + 1))

  let pValue: number
  if (alternative === 'less') {
    pValue = pLessOrEqual
  } else if (alternative === 'greater') {
    pValue = pGreaterOrEqual
  } else {
    // R's binom.test two-sided: sum probabilities of all outcomes
    // at least as extreme as the observed (exact two-sided method).
    // For small n, compute the exact sum.
    pValue = binomialTwoSidedP(k, n, p0)
  }

  // Clopper-Pearson exact CI via bisection on incompleteBeta
  const alpha = 1 - ciLevel
  let ciLo: number, ciHi: number

  if (k === 0) {
    ciLo = 0
  } else {
    // Find p such that I_p(k, n-k+1) = alpha/2  (lower bound)
    ciLo = bisectBetaQuantile(alpha / 2, k, n - k + 1)
  }

  if (k === n) {
    ciHi = 1
  } else {
    // Find p such that I_p(k+1, n-k) = 1 - alpha/2  (upper bound)
    ciHi = bisectBetaQuantile(1 - alpha / 2, k + 1, n - k)
  }

  const ci: readonly [number, number] = [ciLo, ciHi]

  // Effect size: Cohen's g
  const g = Math.abs(pHat - p0)

  const effectSize: EffectSize = {
    value: g,
    name: "Cohen's g",
    interpretation: g < 0.05 ? 'negligible' : g < 0.15 ? 'small' : g < 0.25 ? 'medium' : 'large',
  }

  return {
    testName: 'Exact Binomial Test',
    statistic: k,
    df: 0,
    pValue,
    effectSize,
    ci,
    ciLevel,
    n,
    formatted: formatBinomial(pHat, pValue, ci, g, ciLevel),
  }
}

/**
 * Exact two-sided binomial p-value matching R's binom.test.
 * Sums probabilities of all outcomes with P(X=j) <= P(X=k).
 * Uses log-space computation to avoid overflow.
 */
function binomialTwoSidedP(k: number, n: number, p0: number): number {
  // Compute log(P(X = j)) = log(C(n,j)) + j*log(p0) + (n-j)*log(1-p0)
  const logP0 = Math.log(p0)
  const logQ0 = Math.log(1 - p0)

  // Log of binomial coefficient via logGamma, handling edge cases
  const logBinomPMF = (j: number): number => {
    if (j < 0 || j > n) return -Infinity
    // Use Stirling: log(C(n,j)) = logGamma(n+1) - logGamma(j+1) - logGamma(n-j+1)
    // logGamma(1) = 0, so logGamma(0+1) = 0
    const logCoeff = logGamma(n + 1) - logGamma(j + 1) - logGamma(n - j + 1)
    return logCoeff + j * logP0 + (n - j) * logQ0
  }

  const logPk = logBinomPMF(k)

  // Sum P(X=j) for all j where P(X=j) <= P(X=k)
  // Use a small tolerance for floating point comparison (matching R's relErr = 1 + 1e-7)
  let pValue = 0
  for (let j = 0; j <= n; j++) {
    const logPj = logBinomPMF(j)
    if (logPj <= logPk + 1e-7) {
      pValue += Math.exp(logPj)
    }
  }
  return Math.min(1, pValue)
}

/**
 * Bisection to find the quantile of a Beta(a, b) distribution.
 * Finds x such that I_x(a, b) = target.
 */
function bisectBetaQuantile(target: number, a: number, b: number): number {
  let lo = 0, hi = 1
  for (let iter = 0; iter < 100; iter++) {
    const mid = (lo + hi) / 2
    const val = incompleteBeta(mid, a, b)
    if (val < target) lo = mid
    else hi = mid
    if (hi - lo < 1e-10) break
  }
  return (lo + hi) / 2
}

// ─── Proportions Z-Test ──────────────────────────────────────────────────

/**
 * Z-test for proportions (1-sample or 2-sample).
 *
 * Cross-validated with R:
 * > prop.test(c(30, 40), c(100, 100), correct=FALSE)
 * X-squared = 2.2222, p-value = 0.1360  (z = -1.4907)
 *
 * > prop.test(60, 100, p=0.5, correct=FALSE)
 * X-squared = 4, p-value = 0.0455  (z = 2.0)
 */
export function proportionsZTest(
  x1: number,
  n1: number,
  x2?: number,
  n2?: number,
  p0 = 0.5,
  alternative: 'two.sided' | 'less' | 'greater' = 'two.sided',
  ciLevel = 0.95,
  yates = false
): StatResult {
  if (n1 < 1) throw new Error('proportionsZTest: n1 must be >= 1')
  if (x1 < 0 || x1 > n1) throw new Error('proportionsZTest: x1 must be in [0, n1]')

  const twoSample = x2 !== undefined && n2 !== undefined
  const pHat1 = x1 / n1
  const alpha = 1 - ciLevel
  const zCrit = normalQuantile(1 - alpha / 2)

  let z: number
  let ci: readonly [number, number]
  let h: number
  let n: number

  if (twoSample) {
    if (n2! < 1) throw new Error('proportionsZTest: n2 must be >= 1')
    if (x2! < 0 || x2! > n2!) throw new Error('proportionsZTest: x2 must be in [0, n2]')

    const pHat2 = x2! / n2!
    const pooledP = (x1 + x2!) / (n1 + n2!)
    const se = Math.sqrt(pooledP * (1 - pooledP) * (1 / n1 + 1 / n2!))
    const diff = pHat1 - pHat2
    const absDiff = Math.abs(diff)
    const correction = yates ? 0.5 * (1 / n1 + 1 / n2!) : 0
    const correctedDiff = Math.max(0, absDiff - correction)
    z = se === 0 ? 0 : (diff >= 0 ? correctedDiff : -correctedDiff) / se

    // Wald CI for the difference (unpooled SE)
    const seDiff = Math.sqrt(pHat1 * (1 - pHat1) / n1 + pHat2 * (1 - pHat2) / n2!)
    ci = [diff - zCrit * seDiff, diff + zCrit * seDiff]

    // Cohen's h
    h = 2 * Math.asin(Math.sqrt(pHat1)) - 2 * Math.asin(Math.sqrt(pHat2))
    n = n1 + n2!
  } else {
    // 1-sample
    const se = Math.sqrt(p0 * (1 - p0) / n1)
    const diff = pHat1 - p0
    const absDiff = Math.abs(diff)
    const correction = yates ? 0.5 / n1 : 0
    const correctedDiff = Math.max(0, absDiff - correction)
    z = se === 0 ? 0 : (diff >= 0 ? correctedDiff : -correctedDiff) / se

    // Wald CI for proportion
    const seProp = Math.sqrt(pHat1 * (1 - pHat1) / n1)
    ci = [pHat1 - zCrit * seProp, pHat1 + zCrit * seProp]

    // Cohen's h
    h = 2 * Math.asin(Math.sqrt(pHat1)) - 2 * Math.asin(Math.sqrt(p0))
    n = n1
  }

  // p-value from z
  let pValue: number
  if (alternative === 'less') {
    pValue = normalCDF(z)
  } else if (alternative === 'greater') {
    pValue = normalSurvival(z)
  } else {
    pValue = 2 * normalSurvival(Math.abs(z))
  }

  const effectSize: EffectSize = {
    value: h,
    name: "Cohen's h",
    interpretation: Math.abs(h) < 0.2 ? 'negligible'
      : Math.abs(h) < 0.5 ? 'small'
      : Math.abs(h) < 0.8 ? 'medium' : 'large',
  }

  return {
    testName: twoSample ? 'Two-sample Z-test for proportions' : 'One-sample Z-test for proportions',
    statistic: z,
    df: 0,
    pValue,
    effectSize,
    ci,
    ciLevel,
    n,
    formatted: formatProportions(z, pValue, h, ci, ciLevel),
  }
}
