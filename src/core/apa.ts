/**
 * APA 7th edition string formatting for statistical results.
 * Every stat function uses these formatters to produce the `formatted` string
 * that gets embedded directly into plots as subtitles/captions.
 */

import { roundTo } from './math.js'

/**
 * Format a p-value in APA style.
 * < .001 → "p < .001"
 * otherwise → "p = .025" (leading zero dropped, 3 decimal places)
 */
export function formatP(p: number): string {
  if (p < 0.001) return 'p < .001'
  const rounded = roundTo(p, 3).toFixed(3).replace('0.', '.')
  return `p = ${rounded}`
}

/**
 * Format a statistic value: 2 decimal places, no leading zero for absolute values.
 * Used for t, F, z, χ², etc.
 */
export function formatStat(v: number, decimals = 2): string {
  return roundTo(v, decimals).toFixed(decimals)
}

/**
 * Format a confidence interval: [lo, hi] → "[0.08, 1.22]"
 */
export function formatCI(ci: readonly [number, number], decimals = 2): string {
  return `[${roundTo(ci[0], decimals).toFixed(decimals)}, ${roundTo(ci[1], decimals).toFixed(decimals)}]`
}

/**
 * Format degrees of freedom — integer if whole number, else 2 decimal places.
 */
export function formatDF(df: number | readonly [number, number]): string {
  if (typeof df === 'number') {
    return Number.isInteger(df) ? String(df) : df.toFixed(2)
  }
  const [df1, df2] = df
  const f1 = Number.isInteger(df1) ? String(df1) : df1.toFixed(2)
  const f2 = Number.isInteger(df2) ? String(df2) : df2.toFixed(2)
  return `${f1}, ${f2}`
}

/**
 * Full APA string for a t-test result.
 * e.g. "t(48) = 2.31, p = .025, d = 0.65, 95% CI [0.08, 1.22]"
 */
export function formatTTest(
  t: number,
  df: number,
  pValue: number,
  effectSize: number,
  effectName: string,
  ci: readonly [number, number],
  ciLevel = 0.95
): string {
  const ciPct = Math.round(ciLevel * 100)
  return `t(${formatDF(df)}) = ${formatStat(t)}, ${formatP(pValue)}, ${effectName} = ${formatStat(effectSize)}, ${ciPct}% CI ${formatCI(ci)}`
}

/**
 * Full APA string for ANOVA.
 * e.g. "F(2, 87) = 4.21, p = .018, η² = 0.09"
 */
export function formatANOVA(
  F: number,
  df1: number,
  df2: number,
  pValue: number,
  effectSize: number,
  effectName = 'η²'
): string {
  return `F(${formatDF([df1, df2])}) = ${formatStat(F)}, ${formatP(pValue)}, ${effectName} = ${formatStat(effectSize)}`
}

/**
 * APA string for repeated measures ANOVA.
 * Sphericity assumed: "F(2, 18) = 4.21, p = .018, η²_p = 0.32"
 * GG corrected: "F_GG(1.23, 11.07) = 4.21, p = .042, η²_p = 0.32"
 */
export function formatRMANOVA(
  F: number,
  df1: number,
  df2: number,
  pValue: number,
  effectSize: number,
  effectName = 'η²_p',
  correction?: 'Greenhouse-Geisser' | 'Huynh-Feldt'
): string {
  const prefix = correction === 'Greenhouse-Geisser' ? 'F_GG' : correction === 'Huynh-Feldt' ? 'F_HF' : 'F'
  return `${prefix}(${formatDF([df1, df2])}) = ${formatStat(F)}, ${formatP(pValue)}, ${effectName} = ${formatStat(effectSize)}`
}

/**
 * Full APA string for chi-square test.
 * e.g. "χ²(3) = 8.46, p = .037, V = 0.29"
 */
export function formatChiSq(
  chiSq: number,
  df: number,
  pValue: number,
  effectSize: number,
  effectName = "V"
): string {
  return `χ²(${formatDF(df)}) = ${formatStat(chiSq)}, ${formatP(pValue)}, ${effectName} = ${formatStat(effectSize)}`
}

/**
 * APA string for Pearson/Spearman/Kendall correlation.
 * e.g. "r(48) = 0.54, p = .003, 95% CI [0.28, 0.73]"
 */
export function formatCorrelation(
  r: number,
  df: number,
  pValue: number,
  ci: readonly [number, number],
  method = 'r',
  ciLevel = 0.95
): string {
  const ciPct = Math.round(ciLevel * 100)
  return `${method}(${formatDF(df)}) = ${formatStat(r)}, ${formatP(pValue)}, ${ciPct}% CI ${formatCI(ci)}`
}

/**
 * APA string for linear regression.
 * e.g. "R² = 0.42, F(3, 96) = 23.2, p < .001"
 */
export function formatRegression(
  r2: number,
  adjR2: number,
  F: number,
  df1: number,
  df2: number,
  pValue: number
): string {
  return `R² = ${formatStat(r2)}, adj. R² = ${formatStat(adjR2)}, F(${formatDF([df1, df2])}) = ${formatStat(F)}, ${formatP(pValue)}`
}

/**
 * APA string for Mann-Whitney / Wilcoxon.
 * e.g. "W = 423, p = .031, r = 0.27"
 */
export function formatMannWhitney(
  W: number,
  pValue: number,
  r: number
): string {
  return `W = ${formatStat(W, 0)}, ${formatP(pValue)}, r = ${formatStat(r)}`
}

/**
 * APA string for Kruskal-Wallis.
 * e.g. "H(2) = 12.3, p = .002, η²_H = 0.18"
 */
export function formatKruskalWallis(
  H: number,
  df: number,
  pValue: number,
  etaSq: number
): string {
  return `H(${formatDF(df)}) = ${formatStat(H)}, ${formatP(pValue)}, η²_H = ${formatStat(etaSq)}`
}

/**
 * APA string for LMM.
 * e.g. "ICC = 0.42, AIC = 1234.5"
 */
export function formatLMM(icc: number, aic: number, bic: number, logLik: number): string {
  return `ICC = ${formatStat(icc)}, AIC = ${formatStat(aic, 1)}, BIC = ${formatStat(bic, 1)}, logLik = ${formatStat(logLik, 1)}`
}

/**
 * Interpret Cohen's d effect size.
 */
export function interpretCohensD(d: number): string {
  const abs = Math.abs(d)
  if (abs < 0.2) return 'negligible'
  if (abs < 0.5) return 'small'
  if (abs < 0.8) return 'medium'
  return 'large'
}

/**
 * Interpret eta-squared / omega-squared.
 */
export function interpretEtaSq(eta: number): string {
  if (eta < 0.01) return 'negligible'
  if (eta < 0.06) return 'small'
  if (eta < 0.14) return 'medium'
  return 'large'
}

/**
 * Interpret Pearson r.
 */
export function interpretR(r: number): string {
  const abs = Math.abs(r)
  if (abs < 0.1) return 'negligible'
  if (abs < 0.3) return 'small'
  if (abs < 0.5) return 'medium'
  if (abs < 0.7) return 'large'
  return 'very large'
}

/**
 * Interpret Cramér's V.
 */
export function interpretCramerV(v: number, df: number): string {
  // Benchmarks depend on df (Cohen 1988 adjusted for contingency tables)
  const small = df === 1 ? 0.1 : df === 2 ? 0.07 : 0.06
  const medium = df === 1 ? 0.3 : df === 2 ? 0.21 : 0.17
  if (v < small) return 'negligible'
  if (v < medium) return 'small'
  if (v < medium * 1.5) return 'medium'
  return 'large'
}

import type { EffectInterpretation, FactorFit } from './types.js'

/**
 * APA string for CFA/EFA fit indices.
 * e.g. "χ²(24) = 28.42, p = .241, RMSEA = 0.042 [0.000, 0.085], CFI = 0.987, TLI = 0.982, SRMR = 0.038"
 */
export function formatCFAFit(fit: FactorFit): string {
  return `χ²(${formatDF(fit.df)}) = ${formatStat(fit.chiSq)}, ${formatP(fit.pValue)}, RMSEA = ${formatStat(fit.rmsea, 3)} ${formatCI(fit.rmseaCI, 3)}, CFI = ${formatStat(fit.cfi, 3)}, TLI = ${formatStat(fit.tli, 3)}, SRMR = ${formatStat(fit.srmr, 3)}`
}

/** Generic effect size interpretation helper. */
export function interpretEffect(value: number, thresholds: readonly [number, number, number]): EffectInterpretation {
  const abs = Math.abs(value)
  if (abs < thresholds[0]) return 'negligible'
  if (abs < thresholds[1]) return 'small'
  if (abs < thresholds[2]) return 'medium'
  return 'large'
}
