/**
 * Post-hoc tests for group comparisons.
 * Tukey HSD, Games-Howell, Dunn's test.
 */

import {
  mean as _mean,
  variance as _variance,
  tDistQuantile,
  rank,
  adjustPValues,
  pValueStudentizedRangeApprox,
  normalSurvival,
} from '../core/math.js'
import type { PairwiseResult, GroupData, PAdjMethod } from '../core/types.js'

// ─── Tukey HSD ───────────────────────────────────────────────────────────

/**
 * Tukey's Honest Significant Difference test.
 * Assumes equal variances (uses pooled MSE from ANOVA).
 * Uses Studentized range distribution approximation.
 *
 * Cross-validated with R:
 * > TukeyHSD(aov(value ~ group, data = df))
 */
export function tukeyHSD(
  groups: readonly GroupData[],
  msWithin: number,
  dfWithin: number,
  ciLevel = 0.95
): PairwiseResult[] {
  const results: PairwiseResult[] = []
  const alpha = 1 - ciLevel
  const k = groups.length

  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const g1 = groups[i]!
      const g2 = groups[j]!
      const n1 = g1.values.length
      const n2 = g2.values.length
      const m1 = _mean(g1.values)
      const m2 = _mean(g2.values)
      const diff = m1 - m2

      // SE = sqrt(MSE/2 * (1/n1 + 1/n2))
      const se = Math.sqrt(msWithin / 2 * (1 / n1 + 1 / n2))
      const q = se === 0 ? 0 : Math.abs(diff) / se  // studentized range statistic

      // p-value via proper studentized range distribution
      const pValue = pValueStudentizedRangeApprox(q, k, dfWithin)

      // CI uses studentized range critical value approximation
      // Use Bonferroni-corrected t as approximation for qtukey critical value
      const tCrit = tDistQuantile(1 - alpha / (k * (k - 1)), dfWithin)
      const ciHalf = tCrit * se
      const ci: readonly [number, number] = [diff - ciHalf, diff + ciHalf]

      results.push({
        group1: g1.label,
        group2: g2.label,
        meanDiff: diff,
        se,
        statistic: q,
        pValue,
        pValueAdj: pValue,  // already family-wise corrected
        ci,
        significant: pValue < alpha,
      })
    }
  }
  return results
}

// ─── Games-Howell ────────────────────────────────────────────────────────

/**
 * Games-Howell test — does not assume equal variances.
 * Use when Levene's test is significant or group sizes differ substantially.
 *
 * Cross-validated with R:
 * > rstatix::games_howell_test(df, value ~ group)
 */
export function gamesHowell(
  groups: readonly GroupData[],
  ciLevel = 0.95
): PairwiseResult[] {
  const results: PairwiseResult[] = []
  const alpha = 1 - ciLevel
  const k = groups.length

  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const g1 = groups[i]!
      const g2 = groups[j]!
      const n1 = g1.values.length
      const n2 = g2.values.length
      const m1 = _mean(g1.values)
      const m2 = _mean(g2.values)
      const v1 = _variance(g1.values)
      const v2 = _variance(g2.values)
      const diff = m1 - m2

      // Welch SE (with /2 to match studentized range convention used by rstatix)
      const se = Math.sqrt((v1 / n1 + v2 / n2) / 2)
      const q = se === 0 ? 0 : Math.abs(diff) / se

      // Welch-Satterthwaite df
      const dfNum = (v1 / n1 + v2 / n2) ** 2
      const dfDen = (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)
      const df = dfDen > 0 ? dfNum / dfDen : n1 + n2 - 2

      const pValue = pValueStudentizedRangeApprox(q, k, df)
      const tCrit = tDistQuantile(1 - alpha / (k * (k - 1)), df)
      const ciHalf = tCrit * se

      results.push({
        group1: g1.label,
        group2: g2.label,
        meanDiff: diff,
        se,
        statistic: q,
        pValue,
        pValueAdj: pValue,
        ci: [diff - ciHalf, diff + ciHalf],
        significant: pValue < alpha,
      })
    }
  }
  return results
}

// ─── Dunn's test ──────────────────────────────────────────────────────────

/**
 * Dunn's post-hoc test following Kruskal-Wallis.
 * Uses rank sums and z-scores with adjustable p-value correction.
 *
 * Cross-validated with R:
 * > dunn.test::dunn.test(values, groups, method = "bonferroni")
 */
export function dunnTest(
  groups: readonly GroupData[],
  method: PAdjMethod = 'bonferroni'
): PairwiseResult[] {
  const k = groups.length
  const allValues = groups.flatMap(g => [...g.values])
  const n = allValues.length
  const allRanks = rank(allValues)

  // Rank sums and averages per group
  let offset = 0
  const groupRankMeans: number[] = []
  const groupNs: number[] = []

  for (const g of groups) {
    const gn = g.values.length
    let Rj = 0
    for (let i = 0; i < gn; i++) Rj += allRanks[offset + i]!
    groupRankMeans.push(Rj / gn)
    groupNs.push(gn)
    offset += gn
  }

  // Tie correction
  const tieCounts = new Map<number, number>()
  for (const v of allValues) tieCounts.set(v, (tieCounts.get(v) ?? 0) + 1)
  let C = 0
  for (const t of tieCounts.values()) C += t * t * t - t
  const tieAdj = C / (12 * (n - 1))

  // Variance of rank sums (with tie correction matching R's dunn.test)
  const baseVar = n * (n + 1) / 12 - tieAdj

  // Raw p-values
  const rawPValues: number[] = []
  const rawZs: number[] = []
  const rawSEs: number[] = []
  const pairs: Array<{ i: number; j: number }> = []

  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const diff = (groupRankMeans[i] ?? 0) - (groupRankMeans[j] ?? 0)
      const se = Math.sqrt(
        baseVar * (1 / (groupNs[i] ?? 1) + 1 / (groupNs[j] ?? 1))
      )
      const z = se === 0 ? 0 : diff / se
      // One-tailed p-value matching R's dunn.test default (altp=FALSE)
      // Uses normalSurvival for numerical stability in the tails
      const p = normalSurvival(Math.abs(z))
      rawPValues.push(p)
      rawZs.push(z)
      rawSEs.push(se)
      pairs.push({ i, j })
    }
  }

  // R's dunn.test implements Holm/BH without enforcing monotonicity (no cummax/cummin).
  // Use adjustPValues for bonferroni (identical), but match dunn.test behavior for holm/BH.
  const adjPValues = method === 'holm'
    ? dunnStyleHolm(rawPValues)
    : method === 'BH'
    ? dunnStyleBH(rawPValues)
    : adjustPValues(rawPValues, method)

  return pairs.map(({ i, j }, idx) => ({
    group1: groups[i]!.label,
    group2: groups[j]!.label,
    meanDiff: (groupRankMeans[i] ?? 0) - (groupRankMeans[j] ?? 0),
    se: rawSEs[idx]!,
    statistic: rawZs[idx]!,
    pValue: rawPValues[idx]!,
    pValueAdj: adjPValues[idx]!,
    ci: [NaN, NaN] as readonly [number, number],
    significant: (adjPValues[idx] ?? 1) < 0.05,
  }))
}

// ─── Dunn-style Holm (no cummax) ─────────────────────────────────────────

/**
 * Holm adjustment matching R's dunn.test behavior.
 * Unlike standard Holm (p.adjust), dunn.test does NOT enforce
 * monotonicity via cumulative maximum. It simply multiplies
 * each sorted p-value by (m+1-i) and caps at 1.
 */
function dunnStyleHolm(pValues: readonly number[]): number[] {
  const m = pValues.length
  if (m === 0) return []
  const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => a.p - b.p)
  const adj = new Array<number>(m)
  for (let k = 0; k < m; k++) {
    const { p, i } = order[k]!
    adj[i] = Math.min(1, p * (m - k))
  }
  return adj
}

/**
 * BH (Benjamini-Hochberg) adjustment matching R's dunn.test behavior.
 * Unlike standard BH (p.adjust), dunn.test does NOT enforce
 * monotonicity via cumulative minimum.
 */
function dunnStyleBH(pValues: readonly number[]): number[] {
  const m = pValues.length
  if (m === 0) return []
  const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => b.p - a.p) // descending
  const adj = new Array<number>(m)
  for (let k = 0; k < m; k++) {
    const { p, i } = order[k]!
    const rank = m - k // rank from smallest (1-indexed)
    adj[i] = Math.min(1, p * m / rank)
  }
  return adj
}
