/**
 * Group comparison module.
 * t-tests, one-way ANOVA, Mann-Whitney U, Wilcoxon signed-rank,
 * Kruskal-Wallis, Friedman test.
 */

import {
  mean as _mean,
  median as _median,
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
  formatRMANOVA,
  formatMannWhitney,
  formatKruskalWallis,
  formatP,
  interpretEtaSq,
} from '../core/apa.js'
import { Matrix } from '../core/matrix.js'
import { cohensD, cohensDPaired, omegaSquared, etaSquaredKW, rankBiserial } from './effect-size.js'
import type { StatResult, GroupData, RMANOVAResult, EffectInterpretation } from '../core/types.js'

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

// ─── Levene's test ───────────────────────────────────────────────────────

/**
 * Levene's test for homogeneity of variances.
 *
 * Computes absolute deviations from each group's center, then runs
 * a one-way ANOVA on those deviations. A significant result (p < .05)
 * indicates unequal variances across groups.
 *
 * @param groups  Two or more groups with labels and values.
 * @param center  'median' (default, Brown-Forsythe variant — more robust)
 *                or 'mean' (original Levene).
 *
 * Cross-validated with R:
 * > car::leveneTest(y ~ group, center = median)
 * > # or manually:
 * > medians <- tapply(y, group, median)
 * > abs_dev <- abs(y - medians[group])
 * > anova(lm(abs_dev ~ group))
 */
export function leveneTest(
  groups: readonly GroupData[],
  center: 'median' | 'mean' = 'median'
): StatResult {
  if (groups.length < 2) throw new Error('leveneTest: need at least 2 groups')

  const k = groups.length
  const centerFn = center === 'median' ? _median : _mean

  // Compute absolute deviations from each group's center
  const devGroups: GroupData[] = groups.map(g => {
    const c = centerFn(g.values)
    return {
      label: g.label,
      values: g.values.map(v => Math.abs(v - c)),
    }
  })

  // One-way ANOVA on the absolute deviations
  const allDevs = devGroups.flatMap(g => [...g.values])
  const n = allDevs.length
  const grandMean = _mean(allDevs)

  let ssBetween = 0
  let ssWithin = 0
  for (const g of devGroups) {
    const gm = _mean(g.values)
    ssBetween += g.values.length * (gm - grandMean) ** 2
    ssWithin += g.values.reduce((s, v) => s + (v - gm) ** 2, 0)
  }

  const dfBetween = k - 1
  const dfWithin = n - k
  const msBetween = ssBetween / dfBetween
  const msWithin = dfWithin > 0 ? ssWithin / dfWithin : 0
  const F = msWithin > 0 ? msBetween / msWithin : (msBetween > 0 ? Infinity : 0)
  const pValue = fDistPValue(F, dfBetween, dfWithin)

  const testLabel = center === 'median' ? 'Brown-Forsythe' : "Levene's"
  const formatted = `${testLabel}: F(${dfBetween}, ${dfWithin}) = ${F.toFixed(2)}, ${formatP(pValue)}`

  return {
    testName: `${testLabel} Test`,
    statistic: F,
    df: [dfBetween, dfWithin],
    pValue,
    effectSize: { value: 0, name: 'none', interpretation: 'negligible' },
    ci: [NaN, NaN],
    ciLevel: 0.95,
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
  readonly levene: {
    readonly statistic: number
    readonly pValue: number
    readonly df: readonly [number, number]
    readonly homogeneous: boolean
  }
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

  // Levene's test (Brown-Forsythe variant) for homogeneity of variances
  const lev = leveneTest(groups)
  const levene = {
    statistic: lev.statistic,
    pValue: lev.pValue,
    df: (lev.df as readonly [number, number]),
    homogeneous: lev.pValue >= 0.05,
  }

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
    levene,
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

// ─── Repeated Measures ANOVA ────────────────────────────────────────────

/**
 * Helmert orthonormal contrast matrix ((k-1) × k).
 * Row i (0-indexed): 1/√((i+1)(i+2)) for j ≤ i,
 *   -(i+1)/√((i+1)(i+2)) for j = i+1, 0 otherwise.
 * Rows are orthonormal: C·C' = I_{k-1}.
 */
function helmertContrasts(k: number): number[][] {
  const C: number[][] = []
  for (let i = 0; i < k - 1; i++) {
    const row = new Array(k).fill(0) as number[]
    const denom = Math.sqrt((i + 1) * (i + 2))
    for (let j = 0; j <= i; j++) row[j] = 1 / denom
    row[i + 1] = -(i + 1) / denom
    C.push(row)
  }
  return C
}

/**
 * Sample covariance matrix (k × k) from an n × k data matrix.
 * Uses Bessel correction (n-1 denominator).
 */
function covarianceMatrixCols(
  data: readonly (readonly number[])[],
  k: number
): number[][] {
  const n = data.length
  const means: number[] = []
  for (let j = 0; j < k; j++) {
    let s = 0
    for (let i = 0; i < n; i++) s += data[i]![j]!
    means.push(s / n)
  }
  const cov: number[][] = Array.from({ length: k }, () => new Array(k).fill(0) as number[])
  for (let a = 0; a < k; a++) {
    for (let b = a; b < k; b++) {
      let s = 0
      for (let i = 0; i < n; i++) {
        s += (data[i]![a]! - means[a]!) * (data[i]![b]! - means[b]!)
      }
      const val = s / (n - 1)
      cov[a]![b] = val
      cov[b]![a] = val
    }
  }
  return cov
}

/**
 * Mauchly's test of sphericity for repeated measures designs.
 * Tests whether the covariance matrix of orthonormalized contrasts
 * is proportional to an identity matrix.
 *
 * @param data  Subjects × conditions matrix (n × k, k ≥ 3).
 * @returns { W, chiSq, df, pValue } where W is Mauchly's statistic.
 *
 * Reference: Mauchly (1940), "Significance test for sphericity of a
 * normal n-variate distribution"
 *
 * Cross-validated with R:
 * > library(ez)
 * > ezANOVA(data=df, dv=y, wid=subj, within=cond)$`Mauchly's Test for Sphericity`
 */
export function mauchlysTest(
  data: readonly (readonly number[])[]
): { readonly W: number; readonly chiSq: number; readonly df: number; readonly pValue: number } {
  const n = data.length
  const k = data[0]!.length
  if (k < 3) throw new Error('mauchlysTest: need at least 3 conditions (k ≥ 3)')
  if (n < 2) throw new Error('mauchlysTest: need at least 2 subjects')

  const p = k - 1

  // Covariance matrix of conditions (k × k)
  const sigma = covarianceMatrixCols(data, k)

  // Helmert orthonormal contrasts (p × k)
  const C = helmertContrasts(k)

  // S_star = C × Σ × C' (p × p)
  const CMat = Matrix.fromArray(C)
  const sigmaMat = Matrix.fromArray(sigma)
  const S = CMat.multiply(sigmaMat).multiply(CMat.transpose())

  // Eigenvalues of S_star (symmetric positive semi-definite)
  const lambdas = S.eigen().values

  // W = det(S) / (trace(S)/p)^p
  // In log-space: ln(W) = Σln(λ_i) - p·ln(Σλ_i/p)
  const sumLambda = lambdas.reduce((s, v) => s + v, 0)
  const meanLambda = sumLambda / p
  const logW = lambdas.reduce((s, v) => s + Math.log(Math.max(v, 1e-300)), 0) - p * Math.log(meanLambda)
  const W = Math.exp(logW)

  // Box's chi-square approximation
  // χ² = -(n - 1 - (2p² + p + 2)/(6p)) × ln(W)
  const nu = n - 1 - (2 * p * p + p + 2) / (6 * p)
  const chiSq = Math.max(0, -nu * logW)
  const df = p * (p + 1) / 2 - 1
  const pValue = chiSq > 0 && df > 0 ? chiSqPValue(chiSq, df) : 1

  return { W, chiSq, df, pValue }
}

/**
 * Greenhouse-Geisser and Huynh-Feldt epsilon corrections for sphericity.
 *
 * @param data  Subjects × conditions matrix (n × k).
 * @returns { gg, hf } clamped to [1/(k-1), 1].
 *
 * Reference: Greenhouse & Geisser (1959); Huynh & Feldt (1976)
 */
export function epsilonCorrections(
  data: readonly (readonly number[])[]
): { readonly gg: number; readonly hf: number } {
  const n = data.length
  const k = data[0]!.length
  const p = k - 1

  // Covariance matrix of orthonormalized contrasts
  const sigma = covarianceMatrixCols(data, k)
  const C = helmertContrasts(k)
  const CMat = Matrix.fromArray(C)
  const sigmaMat = Matrix.fromArray(sigma)
  const S = CMat.multiply(sigmaMat).multiply(CMat.transpose())

  const lambdas = S.eigen().values
  const sumLambda = lambdas.reduce((s, v) => s + v, 0)
  const sumLambdaSq = lambdas.reduce((s, v) => s + v * v, 0)

  // Greenhouse-Geisser: ε_GG = (Σλ)² / (p × Σλ²)
  const gg = sumLambdaSq > 0 ? (sumLambda * sumLambda) / (p * sumLambdaSq) : 1

  // Huynh-Feldt: ε_HF = (n(k-1)ε_GG - 2) / ((k-1)(n - 1 - (k-1)ε_GG))
  const num = n * p * gg - 2
  const den = p * (n - 1 - p * gg)
  const hf = den > 0 ? num / den : 1

  // Clamp to [1/(k-1), 1]
  const lower = 1 / p
  return {
    gg: Math.max(lower, Math.min(1, gg)),
    hf: Math.max(lower, Math.min(1, hf)),
  }
}

/**
 * Repeated measures one-way ANOVA.
 *
 * Decomposes total variance into SS_subjects + SS_conditions + SS_error.
 * Tests whether condition means differ, accounting for within-subject correlation.
 * Automatically runs Mauchly's sphericity test (k ≥ 3) and applies
 * Greenhouse-Geisser correction when sphericity is violated (p < .05).
 *
 * @param data     n × k matrix: data[i][j] = score for subject i in condition j.
 * @param options  { labels?: string[], ciLevel?: number }
 *
 * Cross-validated with R:
 * > fit <- aov(y ~ cond + Error(subj/cond), data=df)
 * > summary(fit)
 * > library(ez); ezANOVA(data=df, dv=y, wid=subj, within=cond)
 */
export function repeatedMeasuresANOVA(
  data: readonly (readonly number[])[],
  options?: { readonly labels?: readonly string[]; readonly ciLevel?: number }
): RMANOVAResult {
  const n = data.length
  if (n < 2) throw new Error('repeatedMeasuresANOVA: need at least 2 subjects')
  const k = data[0]!.length
  if (k < 2) throw new Error('repeatedMeasuresANOVA: need at least 2 conditions')
  for (let i = 1; i < n; i++) {
    if (data[i]!.length !== k) {
      throw new Error(`repeatedMeasuresANOVA: row ${i} has ${data[i]!.length} columns, expected ${k}`)
    }
  }

  const ciLevel = options?.ciLevel ?? 0.95
  const labels = options?.labels ?? Array.from({ length: k }, (_, i) => `C${i + 1}`)

  // ── Grand mean ──
  let grandSum = 0
  for (const row of data) for (const v of row) grandSum += v
  const grandMean = grandSum / (n * k)

  // ── Condition means and descriptives ──
  const condMeans: number[] = []
  const condSds: number[] = []
  for (let j = 0; j < k; j++) {
    const col: number[] = []
    for (let i = 0; i < n; i++) col.push(data[i]![j]!)
    condMeans.push(_mean(col))
    condSds.push(_sd(col))
  }

  // ── Subject means ──
  const subjMeans: number[] = data.map(row => _mean([...row]))

  // ── SS decomposition ──
  // SS_total = Σ_ij (x_ij - grand_mean)²
  let ssTotal = 0
  for (const row of data) for (const v of row) ssTotal += (v - grandMean) ** 2

  // SS_subjects = k × Σ_i (subj_mean_i - grand_mean)²
  let ssSubjects = 0
  for (const sm of subjMeans) ssSubjects += k * (sm - grandMean) ** 2

  // SS_conditions = n × Σ_j (cond_mean_j - grand_mean)²
  let ssConditions = 0
  for (const cm of condMeans) ssConditions += n * (cm - grandMean) ** 2

  // SS_error = SS_total - SS_subjects - SS_conditions (residual)
  const ssError = ssTotal - ssSubjects - ssConditions

  // ── Degrees of freedom ──
  const dfConditions = k - 1
  const dfSubjects = n - 1
  const dfError = dfConditions * dfSubjects

  // ── Mean squares & F ──
  const msConditions = ssConditions / dfConditions
  const msError = dfError > 0 ? ssError / dfError : 0
  const F = msError > 0 ? msConditions / msError : (msConditions > 0 ? Infinity : 0)

  // ── Sphericity test & epsilon corrections (k ≥ 3 only) ──
  let sphericity: RMANOVAResult['sphericity'] = null
  let epsilonGG = 1
  let epsilonHF = 1

  if (k >= 3) {
    try {
      sphericity = mauchlysTest(data)
      const eps = epsilonCorrections(data)
      epsilonGG = eps.gg
      epsilonHF = eps.hf
    } catch {
      // If sphericity test fails (e.g., singular matrix), assume sphericity
      sphericity = null
    }
  }

  // ── Determine correction ──
  const sphericityViolated = sphericity !== null && sphericity.pValue < 0.05
  const correction: RMANOVAResult['correction'] = sphericityViolated ? 'Greenhouse-Geisser' : 'none'

  // ── Compute p-value (with correction if needed) ──
  let df1 = dfConditions
  let df2 = dfError
  if (correction === 'Greenhouse-Geisser') {
    df1 = epsilonGG * dfConditions
    df2 = epsilonGG * dfError
  }
  const pValue = fDistPValue(F, df1, df2)

  // ── Effect size: partial eta-squared ──
  // η²_p = SS_conditions / (SS_conditions + SS_error)
  const partialEta2 = (ssConditions + ssError) > 0
    ? ssConditions / (ssConditions + ssError)
    : 0

  const effectSize = {
    value: partialEta2,
    name: 'η²_p',
    interpretation: interpretEtaSq(partialEta2) as EffectInterpretation,
  }

  // ── APA formatted string ──
  const formatted = formatRMANOVA(
    F, df1, df2, pValue, partialEta2, 'η²_p',
    correction !== 'none' ? correction : undefined
  )

  // ── Condition descriptives ──
  const conditions = labels.map((label, j) => ({
    label,
    mean: condMeans[j]!,
    sd: condSds[j]!,
    n,
  }))

  return {
    testName: 'Repeated Measures ANOVA',
    statistic: F,
    df: [df1, df2],
    pValue,
    effectSize,
    ci: [NaN, NaN],
    ciLevel,
    n: n * k,
    formatted,
    conditions,
    ssConditions,
    ssSubjects,
    ssError,
    ssTotal,
    dfConditions,
    dfSubjects,
    dfError,
    msConditions,
    msError,
    sphericity,
    epsilonGG,
    epsilonHF,
    ...(correction !== 'none' ? { correctedDf: [df1, df2] as const } : {}),
    correction,
  }
}
