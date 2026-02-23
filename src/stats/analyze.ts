/**
 * Field-based analysis dispatch.
 * Pass named fields with declared types; analyze() selects and runs the
 * right statistical test automatically.
 */

import type {
  StatResult,
  FrequencyTestResult,
  DescriptiveResult,
  PairwiseResult,
  GroupData,
  FieldType,
  Field,
  NumericField,
  GroupField,
  AnalyzeOptions,
  AnalysisResult,
} from '../core/types.js'

import {
  tTestIndependent,
  tTestPaired,
  oneWayANOVA,
  kruskalWallis,
  mannWhitneyU,
  wilcoxonSignedRank,
} from './comparison.js'
import type { ANOVAResult } from './comparison.js'
import { tukeyHSD, dunnTest } from './post-hoc.js'
import { chiSquareTest, fisherExactTest, contingencyTable } from './frequency.js'
import { describe as computeDescribe, shapiroWilk } from './descriptive.js'

// ─── Internal resolved options ────────────────────────────────────────────

interface ResolvedOptions {
  ciLevel: number
  paired: boolean
  pAdjMethod: NonNullable<AnalyzeOptions['pAdjMethod']>
  forceTest: string | undefined
  equalVariances: boolean
  normalityAlpha: number
}

const DEFAULTS: ResolvedOptions = {
  ciLevel: 0.95,
  paired: false,
  pAdjMethod: 'holm',
  forceTest: undefined,
  equalVariances: false,
  normalityAlpha: 0.05,
}

// ─── detectFieldType ──────────────────────────────────────────────────────

/**
 * Infer the FieldType of an array of raw values.
 *
 * Rules (in order):
 *   1. If every value is a finite number and exactly 2 unique values both
 *      in {0, 1} → 'binary'
 *   2. If every value is a finite number → 'numeric'
 *   3. If there are exactly 2 unique string/mixed values → 'binary'
 *   4. Otherwise → 'categorical'
 *
 * @param values - Raw column values (string | number), length ≥ 1
 * @returns Inferred FieldType
 *
 * @example
 * detectFieldType([0, 1, 0, 1])      // → 'binary'
 * detectFieldType([1.2, 3.4, 5.6])   // → 'numeric'
 * detectFieldType(['A','B','A','C'])  // → 'categorical'
 * detectFieldType(['yes','no'])       // → 'binary'
 */
export function detectFieldType(values: readonly (string | number)[]): FieldType {
  if (values.length === 0) return 'categorical'

  const allFiniteNumber = values.every(v => typeof v === 'number' && isFinite(v as number))

  if (allFiniteNumber) {
    const nums = values as readonly number[]
    const unique = new Set(nums)
    if (unique.size === 2) {
      const sorted = [...unique].sort((a, b) => a - b)
      if (sorted[0] === 0 && sorted[1] === 1) return 'binary'
    }
    return 'numeric'
  }

  const unique = new Set(values)
  return unique.size === 2 ? 'binary' : 'categorical'
}

// ─── splitGroups ──────────────────────────────────────────────────────────

/**
 * Split a numeric outcome array into per-group vectors using a parallel
 * label array. Returns GroupData[] sorted by label (string sort).
 *
 * @param outcome - Numeric observations, length n
 * @param labels  - Group label for each observation, length n (parallel)
 * @returns Array of { label, values }, one per unique label, sorted
 * @throws Error if outcome.length !== labels.length
 */
function splitGroups(
  outcome: readonly number[],
  labels: readonly (string | number)[]
): GroupData[] {
  if (outcome.length !== labels.length) {
    throw new Error(
      `analyze: outcome length (${outcome.length}) !== predictor length (${labels.length})`
    )
  }

  const map = new Map<string, number[]>()
  for (let i = 0; i < labels.length; i++) {
    const key = String(labels[i])
    if (!map.has(key)) map.set(key, [])
    map.get(key)!.push(outcome[i]!)
  }

  return [...map.entries()]
    .sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0))
    .map(([label, values]) => ({ label, values }))
}

// ─── checkNormality ───────────────────────────────────────────────────────

/**
 * Run Shapiro-Wilk on each group. Returns per-group W and p-value,
 * plus a boolean indicating whether ALL groups pass normality.
 *
 * Groups with n < 3 or n > 50 are treated as normal (SW undefined /
 * overly sensitive at large n; CLT applies).
 *
 * @param groups - GroupData[] from splitGroups()
 * @param alpha  - Significance threshold (default 0.05)
 */
function checkNormality(
  groups: readonly GroupData[],
  alpha = 0.05
): { allNormal: boolean; results: ReadonlyArray<{ group: string; W: number; p: number }> } {
  const results = groups.map(g => {
    const n = g.values.length
    if (n < 3 || n > 50) {
      // Treat as normal: SW undefined for n<3, overly sensitive for large n
      return { group: g.label, W: 1, p: 1 }
    }
    const sw = shapiroWilk(g.values)
    return { group: g.label, W: sw.statistic, p: sw.pValue }
  })

  const allNormal = results.every(r => r.p >= alpha)
  return { allNormal, results }
}

// ─── selectTest ───────────────────────────────────────────────────────────

/**
 * Choose the appropriate test name given field types, groups, and options.
 *
 * Dispatch table:
 *   numeric + binary, paired=false → 't-test-independent' (or 'mann-whitney')
 *   numeric + binary, paired=true  → 't-test-paired'      (or 'wilcoxon')
 *   numeric + categorical (3+)     → 'one-way-anova'      (or 'kruskal-wallis')
 *   binary/categorical + any       → 'chi-square'
 *   numeric + none                 → 'describe-only'
 */
function selectTest(
  outcomeFt: FieldType,
  predictorFt: FieldType | undefined,
  groups: readonly GroupData[],
  opts: ResolvedOptions
): string {
  // forceTest overrides all auto-selection
  if (opts.forceTest) return opts.forceTest

  if (outcomeFt === 'numeric') {
    if (predictorFt === undefined) return 'describe-only'

    if (predictorFt === 'binary') {
      const { allNormal } = checkNormality(groups, opts.normalityAlpha)
      if (opts.paired) {
        return allNormal ? 't-test-paired' : 'wilcoxon'
      }
      return allNormal ? 't-test-independent' : 'mann-whitney'
    }

    if (predictorFt === 'categorical') {
      const { allNormal } = checkNormality(groups, opts.normalityAlpha)
      return allNormal ? 'one-way-anova' : 'kruskal-wallis'
    }

    throw new Error(
      `analyze: unsupported predictor type '${predictorFt}' for numeric outcome. ` +
        `Use 'binary' or 'categorical'.`
    )
  }

  // binary or categorical outcome → frequency test
  return 'chi-square'
}

// ─── analyze ──────────────────────────────────────────────────────────────

/**
 * High-level statistical dispatch: pass fields, get the right test result.
 *
 * Automatically:
 *   - Detects field types from the declared .type property
 *   - Splits numeric outcome by group labels
 *   - Selects parametric vs non-parametric via Shapiro-Wilk
 *   - Runs the selected test with remaining options forwarded
 *   - Computes descriptive statistics for numeric outcomes
 *   - Runs post-hoc tests when 3+ groups are present
 *
 * @param outcome   - The outcome/dependent variable field
 * @param predictor - The grouping/independent variable field (optional)
 * @param opts      - Tuning options (all optional)
 *
 * @returns AnalysisResult with test name, StatResult, optional descriptives,
 *          optional posthoc, and the normality check used for routing.
 *
 * @throws Error if outcome and predictor have different lengths
 * @throws Error if paired=true but group sizes are unequal
 * @throws Error if forceTest names an unknown test
 *
 * @example — independent t-test (auto-detected)
 * analyze(
 *   { type: 'numeric', name: 'score', values: [72, 85, 90, 68, 77] },
 *   { type: 'binary',  name: 'group', values: ['A','B','A','B','A'] }
 * )
 *
 * @example — force Kruskal-Wallis
 * analyze(
 *   { type: 'numeric',     name: 'rt',   values: [...] },
 *   { type: 'categorical', name: 'cond', values: [...] },
 *   { forceTest: 'kruskal-wallis' }
 * )
 *
 * @example — paired t-test
 * analyze(
 *   { type: 'numeric', name: 'post', values: [80, 85, 90] },
 *   { type: 'binary',  name: 'time', values: ['pre','post','pre'] },
 *   { paired: true }
 * )
 */
export function analyze(
  outcome: Field,
  predictor?: Field,
  opts?: AnalyzeOptions
): AnalysisResult {
  const options: ResolvedOptions = { ...DEFAULTS, ...opts }

  // ── Categorical / binary outcome → frequency test ──────────────────────
  if (outcome.type === 'binary' || outcome.type === 'categorical') {
    if (predictor === undefined) {
      throw new Error(
        `analyze: categorical/binary outcome requires a predictor for frequency tests`
      )
    }

    const groupOutcome = outcome as GroupField
    const groupPredictor = predictor as GroupField | NumericField

    const { table } = contingencyTable(
      groupOutcome.values,
      groupPredictor.values as readonly (string | number)[]
    )

    const testName = options.forceTest ?? 'chi-square'
    let result: StatResult | FrequencyTestResult

    if (testName === 'fisher') {
      if (table.length !== 2 || table[0]!.length !== 2) {
        throw new Error(`analyze: Fisher's exact test requires a 2×2 table`)
      }
      result = fisherExactTest(table[0]![0]!, table[0]![1]!, table[1]![0]!, table[1]![1]!)
    } else {
      // chi-square; auto-fallback to Fisher for 2×2 with expected counts < 5
      if (!options.forceTest && table.length === 2 && table[0]!.length === 2) {
        const n = table[0]![0]! + table[0]![1]! + table[1]![0]! + table[1]![1]!
        const r0 = table[0]![0]! + table[0]![1]!
        const r1 = n - r0
        const c0 = table[0]![0]! + table[1]![0]!
        const c1 = n - c0
        const anyLowExpected =
          (r0 * c0) / n < 5 ||
          (r0 * c1) / n < 5 ||
          (r1 * c0) / n < 5 ||
          (r1 * c1) / n < 5

        if (anyLowExpected) {
          result = fisherExactTest(table[0]![0]!, table[0]![1]!, table[1]![0]!, table[1]![1]!)
          return {
            test: result.testName,
            outcome: outcome.name,
            predictor: predictor.name,
            result,
          }
        }
      }
      result = chiSquareTest(table)
    }

    return {
      test: result.testName,
      outcome: outcome.name,
      predictor: predictor.name,
      result,
    }
  }

  // ── Numeric outcome ────────────────────────────────────────────────────
  const numericOutcome = outcome as NumericField
  const numericValues = numericOutcome.values

  // No predictor → descriptive statistics only
  if (predictor === undefined) {
    const desc = computeDescribe(numericValues, options.ciLevel)
    const dummyResult: StatResult = {
      testName: 'Descriptive statistics',
      statistic: NaN,
      df: 0,
      pValue: NaN,
      effectSize: { value: NaN, name: 'none', interpretation: 'negligible' },
      ci: [NaN, NaN],
      ciLevel: options.ciLevel,
      n: numericValues.length,
      formatted: desc.formatted,
    }
    return {
      test: 'descriptive',
      outcome: outcome.name,
      result: dummyResult,
      descriptives: [desc],
    }
  }

  // Split outcome values by predictor labels
  const groups = splitGroups(numericValues, predictor.values as readonly (string | number)[])

  // Normality check (used for routing; always computed for the result)
  const normalityCheck = checkNormality(groups, options.normalityAlpha)

  // Select test
  const selectedTest = selectTest(outcome.type, predictor.type, groups, options)

  // Run selected test
  let result: StatResult | FrequencyTestResult
  let posthoc: PairwiseResult[] | undefined

  if (selectedTest === 't-test-independent' || selectedTest === 'mann-whitney') {
    if (groups.length !== 2) {
      throw new Error(
        `analyze: '${selectedTest}' requires exactly 2 groups, got ${groups.length}`
      )
    }
    const [g1, g2] = groups as [GroupData, GroupData]

    if (selectedTest === 't-test-independent') {
      result = tTestIndependent(g1.values, g2.values, options.equalVariances, options.ciLevel)
    } else {
      result = mannWhitneyU(g1.values, g2.values)
    }
  } else if (selectedTest === 't-test-paired' || selectedTest === 'wilcoxon') {
    if (groups.length !== 2) {
      throw new Error(
        `analyze: '${selectedTest}' requires exactly 2 groups, got ${groups.length}`
      )
    }
    const [g1, g2] = groups as [GroupData, GroupData]

    if (g1.values.length !== g2.values.length) {
      throw new Error(
        `analyze: paired=true but group sizes are unequal ` +
          `(${g1.label}: ${g1.values.length}, ${g2.label}: ${g2.values.length})`
      )
    }

    if (selectedTest === 't-test-paired') {
      result = tTestPaired(g1.values, g2.values, options.ciLevel)
    } else {
      result = wilcoxonSignedRank(g1.values, g2.values)
    }
  } else if (selectedTest === 'one-way-anova') {
    const anovaResult = oneWayANOVA(groups) as ANOVAResult
    result = anovaResult
    posthoc = tukeyHSD(groups, anovaResult.msWithin, anovaResult.dfWithin, options.ciLevel)
  } else if (selectedTest === 'kruskal-wallis') {
    result = kruskalWallis(groups)
    posthoc = dunnTest(groups, options.pAdjMethod)
  } else if (selectedTest === 'chi-square' || selectedTest === 'fisher') {
    throw new Error(
      `analyze: '${selectedTest}' is not applicable to numeric outcomes`
    )
  } else {
    throw new Error(`analyze: unknown test '${selectedTest}'`)
  }

  // Descriptive statistics for each group
  const descriptives: DescriptiveResult[] = groups.map(g =>
    computeDescribe(g.values, options.ciLevel)
  )

  return {
    test: result.testName,
    outcome: outcome.name,
    predictor: predictor.name,
    result,
    descriptives,
    // Only include optional fields when they have a value (exactOptionalPropertyTypes)
    ...(posthoc !== undefined && { posthoc }),
    normality: normalityCheck.results,
  }
}
