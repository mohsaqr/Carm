/**
 * Field-based analysis dispatch.
 * Pass named fields with declared types; analyze() selects and runs the
 * right statistical test automatically.
 */
import type { FieldType, Field, AnalyzeOptions, AnalysisResult } from '../core/types.js';
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
export declare function detectFieldType(values: readonly (string | number)[]): FieldType;
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
export declare function analyze(outcome: Field, predictor?: Field, opts?: AnalyzeOptions): AnalysisResult;
//# sourceMappingURL=analyze.d.ts.map