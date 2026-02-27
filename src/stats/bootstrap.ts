/**
 * Bootstrap confidence interval estimation.
 * Percentile and BCa methods with seeded PRNG for reproducibility.
 */

import { quantile, normalQuantile, normalCDF, mean as _mean } from '../core/math.js'
import { PRNG } from '../core/prng.js'
import { formatStat, formatCI } from '../core/apa.js'
import type { BootstrapCIResult } from '../core/types.js'

/**
 * Bootstrap CI for a single-sample statistic.
 * @param data       Input data array
 * @param statFn     Function that computes the statistic from a data array
 * @param options    { nBoot, ciLevel, method, seed }
 *
 * Cross-validated with R:
 * > library(boot)
 * > b <- boot(data, function(d,i) mean(d[i]), R=2000)
 * > boot.ci(b, type="perc")
 * > boot.ci(b, type="bca")
 */
export function bootstrapCI(
  data: readonly number[],
  statFn: (d: readonly number[]) => number,
  options?: {
    readonly nBoot?: number
    readonly ciLevel?: number
    readonly method?: 'percentile' | 'bca'
    readonly seed?: number
  }
): BootstrapCIResult {
  const nBoot = options?.nBoot ?? 2000
  const ciLevel = options?.ciLevel ?? 0.95
  const method = options?.method ?? 'percentile'
  const seed = options?.seed ?? 42
  const n = data.length
  if (n === 0) throw new Error('bootstrapCI: empty data')

  const rng = new PRNG(seed)
  const estimate = statFn(data)

  // Generate bootstrap replicates
  const bootStats: number[] = []
  const sample = new Array<number>(n)
  for (let b = 0; b < nBoot; b++) {
    for (let i = 0; i < n; i++) {
      sample[i] = data[Math.floor(rng.next() * n)]!
    }
    bootStats.push(statFn(sample))
  }
  bootStats.sort((a, b) => a - b)

  const bootMean = _mean(bootStats)
  const se = Math.sqrt(bootStats.reduce((s, v) => s + (v - bootMean) ** 2, 0) / (nBoot - 1))
  const alpha = 1 - ciLevel

  let ci: readonly [number, number]
  if (method === 'bca') {
    // BCa: bias-correction and acceleration
    // z0 = Phi^-1(proportion of boot stats < estimate)
    const nBelow = bootStats.filter(v => v < estimate).length
    const z0 = normalQuantile(Math.max(0.0001, Math.min(0.9999, nBelow / nBoot)))

    // acceleration: jackknife
    const jackStats: number[] = []
    for (let i = 0; i < n; i++) {
      const jk = [...data.slice(0, i), ...data.slice(i + 1)]
      jackStats.push(statFn(jk))
    }
    const jkMean = _mean(jackStats)
    let num3 = 0, denom2 = 0
    for (const j of jackStats) {
      const diff = jkMean - j
      num3 += diff ** 3
      denom2 += diff ** 2
    }
    const acc = denom2 > 0 ? num3 / (6 * Math.pow(denom2, 1.5)) : 0

    // Adjusted percentiles
    const zAlpha = normalQuantile(alpha / 2)
    const zUpper = normalQuantile(1 - alpha / 2)
    const a1 = normalCDF(z0 + (z0 + zAlpha) / (1 - acc * (z0 + zAlpha)))
    const a2 = normalCDF(z0 + (z0 + zUpper) / (1 - acc * (z0 + zUpper)))

    const lo = quantile(bootStats, Math.max(0, Math.min(1, a1)))
    const hi = quantile(bootStats, Math.max(0, Math.min(1, a2)))
    ci = [lo, hi]
  } else {
    // Percentile method
    ci = [
      quantile(bootStats, alpha / 2),
      quantile(bootStats, 1 - alpha / 2),
    ]
  }

  const ciPct = Math.round(ciLevel * 100)
  const formatted = `estimate = ${formatStat(estimate)}, ${ciPct}% CI ${formatCI(ci)}, SE_boot = ${formatStat(se)}`

  return { estimate, ci, se, ciLevel, nBoot, method, formatted }
}

/**
 * Bootstrap CI for a two-sample statistic.
 * @param x1       First sample
 * @param x2       Second sample
 * @param statFn   Function that computes the statistic from two data arrays
 * @param options  { nBoot, ciLevel, seed }
 *
 * Cross-validated with R:
 * > library(boot)
 * > b <- boot(cbind(x1,x2), function(d,i) mean(d[i,1])-mean(d[i,2]), R=2000)
 * > boot.ci(b, type="perc")
 */
export function bootstrapCITwoSample(
  x1: readonly number[],
  x2: readonly number[],
  statFn: (a: readonly number[], b: readonly number[]) => number,
  options?: {
    readonly nBoot?: number
    readonly ciLevel?: number
    readonly seed?: number
  }
): BootstrapCIResult {
  const nBoot = options?.nBoot ?? 2000
  const ciLevel = options?.ciLevel ?? 0.95
  const seed = options?.seed ?? 42
  const n1 = x1.length, n2 = x2.length
  if (n1 === 0 || n2 === 0) throw new Error('bootstrapCITwoSample: empty sample')

  const rng = new PRNG(seed)
  const estimate = statFn(x1, x2)

  const bootStats: number[] = []
  const s1 = new Array<number>(n1)
  const s2 = new Array<number>(n2)
  for (let b = 0; b < nBoot; b++) {
    for (let i = 0; i < n1; i++) s1[i] = x1[Math.floor(rng.next() * n1)]!
    for (let i = 0; i < n2; i++) s2[i] = x2[Math.floor(rng.next() * n2)]!
    bootStats.push(statFn(s1, s2))
  }
  bootStats.sort((a, b) => a - b)

  const bootMean = _mean(bootStats)
  const se = Math.sqrt(bootStats.reduce((s, v) => s + (v - bootMean) ** 2, 0) / (nBoot - 1))
  const alpha = 1 - ciLevel
  const ci: readonly [number, number] = [
    quantile(bootStats, alpha / 2),
    quantile(bootStats, 1 - alpha / 2),
  ]

  const ciPct = Math.round(ciLevel * 100)
  const formatted = `estimate = ${formatStat(estimate)}, ${ciPct}% CI ${formatCI(ci)}, SE_boot = ${formatStat(se)}`

  return { estimate, ci, se, ciLevel, nBoot, method: 'percentile', formatted }
}
