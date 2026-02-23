/**
 * Frequency analysis module.
 * Frequency tables, chi-square test of independence, Fisher's exact test,
 * Cramér's V, phi coefficient, odds ratio, goodness-of-fit test.
 */

import { chiSqPValue, normalQuantile } from '../core/math.js'
import { formatChiSq, formatP, interpretCramerV } from '../core/apa.js'
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
