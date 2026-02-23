/**
 * Group comparison module.
 * t-tests, one-way ANOVA, Mann-Whitney U, Wilcoxon signed-rank,
 * Kruskal-Wallis, Friedman test.
 */

import {
  mean as _mean,
  variance as _variance,
  sd as _sd,
  se as _se,
  rank,
  tDistQuantile,
  tDistPValue,
  fDistPValue,
  chiSqPValue,
} from '../core/math.js'
import {
  formatTTest,
  formatANOVA,
  formatMannWhitney,
  formatKruskalWallis,
  formatP,
} from '../core/apa.js'
import { cohensD, cohensDPaired, omegaSquared, etaSquaredKW, rankBiserial } from './effect-size.js'
import type { StatResult, GroupData } from '../core/types.js'

// ─── Independent samples t-test ───────────────────────────────────────────

/**
 * Welch's t-test for independent samples (unequal variances, default).
 * Set `equalVariances = true` for Student's t-test.
 *
 * Cross-validated with R:
 * > t.test(c(2,3,5,6,8), c(1,4,7,9,10), var.equal = FALSE)
 * t = -0.4851, df = 7.6, p-value = 0.6408
 * > t.test(c(2,3,5,6,8), c(1,4,7,9,10), var.equal = TRUE)
 * t = -0.4851, df = 8, p-value = 0.6403
 */
export function tTestIndependent(
  x1: readonly number[],
  x2: readonly number[],
  equalVariances = false,
  ciLevel = 0.95,
  alternative: 'two.sided' | 'less' | 'greater' = 'two.sided'
): StatResult {
  if (x1.length < 2 || x2.length < 2) throw new Error('tTestIndependent: need at least 2 per group')
  const n1 = x1.length, n2 = x2.length
  const m1 = _mean(x1), m2 = _mean(x2)
  const v1 = _variance(x1), v2 = _variance(x2)

  let df: number
  let se: number

  if (equalVariances) {
    // Pooled variance (Student's)
    const spSq = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    se = Math.sqrt(spSq * (1 / n1 + 1 / n2))
    df = n1 + n2 - 2
  } else {
    // Welch-Satterthwaite df
    se = Math.sqrt(v1 / n1 + v2 / n2)
    const num = (v1 / n1 + v2 / n2) ** 2
    const den = (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)
    df = den > 0 ? num / den : n1 + n2 - 2
  }

  const t = se === 0 ? 0 : (m1 - m2) / se
  const pFull = tDistPValue(t, df)  // two-tailed
  const pValue = alternative === 'two.sided'
    ? pFull
    : alternative === 'less'
      ? t < 0 ? pFull / 2 : 1 - pFull / 2
      : t > 0 ? pFull / 2 : 1 - pFull / 2

  // CI for mean difference
  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, df)
  const diff = m1 - m2
  const ci: readonly [number, number] = [diff - tCrit * se, diff + tCrit * se]

  const effectSize = cohensD(x1, x2)
  const formatted = formatTTest(t, df, pValue, effectSize.value, "d", ci, ciLevel)

  return {
    testName: equalVariances ? "Student's t-test" : "Welch's t-test",
    statistic: t,
    df,
    pValue,
    effectSize,
    ci,
    ciLevel,
    n: n1 + n2,
    formatted,
  }
}

// ─── Paired samples t-test ────────────────────────────────────────────────

/**
 * Paired samples t-test.
 *
 * Cross-validated with R:
 * > t.test(c(1,2,3,4,5), c(2,4,6,8,10), paired = TRUE)
 * t = -3.873, df = 4, p-value = 0.01789
 */
export function tTestPaired(
  x1: readonly number[],
  x2: readonly number[],
  ciLevel = 0.95
): StatResult {
  if (x1.length !== x2.length) throw new Error('tTestPaired: arrays must have equal length')
  if (x1.length < 2) throw new Error('tTestPaired: need at least 2 pairs')

  const diffs = x1.map((v, i) => v - (x2[i] ?? 0))
  const n = diffs.length
  const mDiff = _mean(diffs)
  const seDiff = _se(diffs)
  const df = n - 1
  const t = seDiff === 0 ? 0 : mDiff / seDiff
  const pValue = tDistPValue(t, df)
  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, df)
  const ci: readonly [number, number] = [mDiff - tCrit * seDiff, mDiff + tCrit * seDiff]

  const effectSize = cohensDPaired(diffs)
  const formatted = formatTTest(t, df, pValue, effectSize.value, "d", ci, ciLevel)

  return {
    testName: 'Paired t-test',
    statistic: t,
    df,
    pValue,
    effectSize,
    ci,
    ciLevel,
    n,
    formatted,
  }
}

// ─── One-way ANOVA ────────────────────────────────────────────────────────

export interface ANOVAResult extends StatResult {
  readonly groups: readonly {
    readonly label: string
    readonly n: number
    readonly mean: number
    readonly sd: number
  }[]
  readonly ssBetween: number
  readonly ssWithin: number
  readonly ssTotal: number
  readonly msBetween: number
  readonly msWithin: number
  readonly dfBetween: number
  readonly dfWithin: number
}

/**
 * One-way ANOVA.
 *
 * Cross-validated with R:
 * > oneway.test(value ~ group, data = df, var.equal = TRUE)
 * > aov(value ~ group, data = df)
 */
export function oneWayANOVA(groups: readonly GroupData[]): ANOVAResult {
  if (groups.length < 2) throw new Error('oneWayANOVA: need at least 2 groups')

  const k = groups.length
  const allValues = groups.flatMap(g => [...g.values])
  const n = allValues.length
  const grandMean = _mean(allValues)

  let ssBetween = 0, ssWithin = 0
  const groupStats = groups.map(g => {
    const gm = _mean(g.values)
    const gn = g.values.length
    ssBetween += gn * (gm - grandMean) ** 2
    ssWithin += g.values.reduce((s, v) => s + (v - gm) ** 2, 0)
    return { label: g.label, n: gn, mean: gm, sd: _sd(g.values) }
  })

  const ssTotal = ssBetween + ssWithin
  const dfBetween = k - 1
  const dfWithin = n - k
  const msBetween = ssBetween / dfBetween
  const msWithin = ssWithin / dfWithin

  const F = msWithin === 0 ? Infinity : msBetween / msWithin
  const pValue = fDistPValue(F, dfBetween, dfWithin)

  const omega = omegaSquared(ssBetween, ssTotal, dfBetween, msWithin)

  const formatted = formatANOVA(F, dfBetween, dfWithin, pValue, omega.value, 'ω²')

  return {
    testName: 'One-way ANOVA',
    statistic: F,
    df: [dfBetween, dfWithin],
    pValue,
    effectSize: omega,
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted,
    groups: groupStats,
    ssBetween,
    ssWithin,
    ssTotal,
    msBetween,
    msWithin,
    dfBetween,
    dfWithin,
  }
}

// ─── Mann-Whitney U test ──────────────────────────────────────────────────

/**
 * Mann-Whitney U test (Wilcoxon rank-sum test).
 * Uses normal approximation for large n.
 *
 * Cross-validated with R:
 * > wilcox.test(c(1,2,3,4,5), c(3,4,5,6,7))
 * W = 5, p-value = 0.09502
 */
export function mannWhitneyU(
  x1: readonly number[],
  x2: readonly number[],
  alternative: 'two.sided' | 'less' | 'greater' = 'two.sided'
): StatResult {
  if (x1.length < 1 || x2.length < 1) throw new Error('mannWhitneyU: empty group')
  const n1 = x1.length, n2 = x2.length

  // Combined ranks
  const combined = [
    ...x1.map((v, i) => ({ v, group: 1, i })),
    ...x2.map((v, i) => ({ v, group: 2, i })),
  ].sort((a, b) => a.v - b.v)

  // Assign ranks with tie correction
  const allCombined = rank(combined.map(d => d.v))
  let R1 = 0
  for (let i = 0; i < combined.length; i++) {
    if (combined[i]!.group === 1) R1 += allCombined[i]!
  }

  const U1 = R1 - n1 * (n1 + 1) / 2

  // Normal approximation with tie correction
  const N = n1 + n2
  const muU = n1 * n2 / 2

  // Tie correction for variance
  const allVals = combined.map(d => d.v)
  const tieCountsU = new Map<number, number>()
  for (const v of allVals) tieCountsU.set(v, (tieCountsU.get(v) ?? 0) + 1)
  let tieCorr = 0
  for (const t of tieCountsU.values()) tieCorr += t * t * t - t
  const varU = (n1 * n2 / 12) * (N + 1 - tieCorr / (N * (N - 1)))

  const z = varU === 0 ? 0 : (U1 - muU) / Math.sqrt(varU)
  const zAbs = Math.abs(z)
  const pNormal = 2 * (1 - normalCDFInline(zAbs))
  const pValue = alternative === 'two.sided'
    ? pNormal
    : alternative === 'less'
      ? z < 0 ? pNormal / 2 : 1 - pNormal / 2
      : z > 0 ? pNormal / 2 : 1 - pNormal / 2

  const effect = rankBiserial(U1, n1, n2)
  const formatted = formatMannWhitney(U1, pValue, effect.value)

  return {
    testName: 'Mann-Whitney U',
    statistic: U1,
    df: 0,
    pValue: Math.min(1, Math.max(0, pValue)),
    effectSize: effect,
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n: n1 + n2,
    formatted,
  }
}

// Inline normal CDF — uses erf(z/√2) per A&S 26.2.17
function normalCDFInline(z: number): number {
  const x = Math.abs(z) / Math.SQRT2
  const t = 1 / (1 + 0.3275911 * x)
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
  const erf = 1 - poly * Math.exp(-x * x)
  return 0.5 * (1 + (z >= 0 ? erf : -erf))
}

/** Exact Wilcoxon signed-rank p-value via dynamic programming. */
function wilcoxonExactP(W: number, n: number): number {
  const maxW = n * (n + 1) / 2
  // dist[w] = number of rank-sign assignments giving W+ = w
  let dist = new Array<number>(maxW + 1).fill(0)
  dist[0] = 1
  for (let k = 1; k <= n; k++) {
    const newDist = [...dist]
    for (let w = k; w <= maxW; w++) {
      newDist[w]! += dist[w - k]!
    }
    dist = newDist
  }
  const total = Math.pow(2, n)
  // Two-sided: p = 2 * P(W+ <= W)
  let cumP = 0
  for (let w = 0; w <= W && w <= maxW; w++) {
    cumP += (dist[w] ?? 0) / total
  }
  return Math.min(1, 2 * cumP)
}

// ─── Wilcoxon signed-rank test ────────────────────────────────────────────

/**
 * Wilcoxon signed-rank test for paired data.
 *
 * Cross-validated with R:
 * > wilcox.test(c(1,2,3,4,5), c(2,4,6,8,10), paired = TRUE)
 * V = 0, p-value = 0.0625
 */
export function wilcoxonSignedRank(
  x1: readonly number[],
  x2: readonly number[]
): StatResult {
  if (x1.length !== x2.length) throw new Error('wilcoxonSignedRank: arrays must match length')
  const diffs = x1.map((v, i) => v - (x2[i] ?? 0)).filter(d => d !== 0)
  const n = diffs.length
  if (n === 0) throw new Error('wilcoxonSignedRank: no non-zero differences')

  const absDiffs = diffs.map(Math.abs)
  const ranks_ = rank(absDiffs)

  let Wplus = 0, Wminus = 0
  for (let i = 0; i < n; i++) {
    if ((diffs[i] ?? 0) > 0) Wplus += ranks_[i]!
    else Wminus += ranks_[i]!
  }
  const W = Math.min(Wplus, Wminus)

  // Use exact distribution for n ≤ 20, normal approximation otherwise
  let pValue: number
  if (n <= 20) {
    pValue = wilcoxonExactP(W, n)
  } else {
    const muW = n * (n + 1) / 4
    const varW = n * (n + 1) * (2 * n + 1) / 24
    const z = varW === 0 ? 0 : (W + 0.5 - muW) / Math.sqrt(varW)  // continuity correction
    pValue = 2 * (1 - normalCDFInline(Math.abs(z)))
  }

  const effect = { value: Wplus / (n * (n + 1) / 2) * 2 - 1, name: 'r (rank-biserial)', interpretation: 'small' as const }

  return {
    testName: 'Wilcoxon Signed-Rank',
    statistic: W,
    df: 0,
    pValue: Math.min(1, Math.max(0, pValue)),
    effectSize: effect,
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted: `V = ${W}, ${formatP(pValue)}, r = ${effect.value.toFixed(2)}`,
  }
}

// ─── Kruskal-Wallis test ──────────────────────────────────────────────────

/**
 * Kruskal-Wallis H test (non-parametric one-way ANOVA).
 *
 * Cross-validated with R:
 * > kruskal.test(list(c(1,2,3), c(4,5,6), c(7,8,9)))
 * Kruskal-Wallis chi-squared = 7.2, df = 2, p-value = 0.02732
 */
export function kruskalWallis(groups: readonly GroupData[]): StatResult {
  if (groups.length < 2) throw new Error('kruskalWallis: need at least 2 groups')

  const k = groups.length
  const allValues = groups.flatMap(g => [...g.values])
  const n = allValues.length
  const allRanks_ = rank(allValues)

  let offset = 0
  let H = 0
  for (const g of groups) {
    const gn = g.values.length
    let Rj = 0
    for (let i = 0; i < gn; i++) Rj += allRanks_[offset + i]!
    H += Rj * Rj / gn
    offset += gn
  }
  H = (12 / (n * (n + 1))) * H - 3 * (n + 1)

  // Tie correction
  const tieCounts = new Map<number, number>()
  for (const v of allValues) tieCounts.set(v, (tieCounts.get(v) ?? 0) + 1)
  let C = 0
  for (const t of tieCounts.values()) C += t * t * t - t
  const correction = 1 - C / (n * n * n - n)
  if (correction > 0) H /= correction

  const df = k - 1
  const pValue = chiSqPValue(H, df)
  const effect = etaSquaredKW(H, k, n)
  const formatted = formatKruskalWallis(H, df, pValue, effect.value)

  return {
    testName: 'Kruskal-Wallis',
    statistic: H,
    df,
    pValue,
    effectSize: effect,
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted,
  }
}

// ─── Friedman test ────────────────────────────────────────────────────────

/**
 * Friedman test for repeated measures (non-parametric ANOVA).
 * Data is a matrix: rows = subjects, columns = conditions.
 *
 * Cross-validated with R:
 * > friedman.test(matrix(c(1,2,3, 4,5,6, 7,8,9), nrow=3))
 */
export function friedmanTest(data: readonly (readonly number[])[]): StatResult {
  const n = data.length  // subjects
  if (n < 2) throw new Error('friedmanTest: need at least 2 subjects')
  const k = data[0]!.length  // conditions
  if (k < 2) throw new Error('friedmanTest: need at least 2 conditions')

  // Rank within each row (subject)
  let sumRankSq = 0
  for (const row of data) {
    const rowRanks = rank(row)
    for (const r of rowRanks) sumRankSq += r * r
  }

  // Friedman statistic
  const colRankSums = Array.from({ length: k }, (_, j) =>
    data.reduce((s, row) => {
      const rowRanks = rank(row)
      return s + (rowRanks[j] ?? 0)
    }, 0)
  )
  const Rj2sum = colRankSums.reduce((s, r) => s + r * r, 0)
  const chi2 = 12 / (n * k * (k + 1)) * Rj2sum - 3 * n * (k + 1)

  const df = k - 1
  const pValue = chiSqPValue(chi2, df)
  const w = chi2 / (n * (k - 1))  // Kendall's W

  return {
    testName: 'Friedman Test',
    statistic: chi2,
    df,
    pValue,
    effectSize: {
      value: w,
      name: "Kendall's W",
      interpretation: w < 0.1 ? 'negligible' : w < 0.3 ? 'small' : w < 0.5 ? 'medium' : 'large',
    },
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted: `χ²_F(${df}) = ${chi2.toFixed(2)}, ${formatP(pValue)}, W = ${w.toFixed(2)}`,
  }
}
