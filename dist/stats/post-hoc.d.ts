/**
 * Post-hoc tests for group comparisons.
 * Tukey HSD, Games-Howell, Dunn's test.
 */
import type { PairwiseResult, GroupData, PAdjMethod } from '../core/types.js';
/**
 * Tukey's Honest Significant Difference test.
 * Assumes equal variances (uses pooled MSE from ANOVA).
 * Uses Studentized range distribution approximation.
 *
 * Cross-validated with R:
 * > TukeyHSD(aov(value ~ group, data = df))
 */
export declare function tukeyHSD(groups: readonly GroupData[], msWithin: number, dfWithin: number, ciLevel?: number): PairwiseResult[];
/**
 * Games-Howell test — does not assume equal variances.
 * Use when Levene's test is significant or group sizes differ substantially.
 *
 * Cross-validated with R:
 * > rstatix::games_howell_test(df, value ~ group)
 */
export declare function gamesHowell(groups: readonly GroupData[], ciLevel?: number): PairwiseResult[];
/**
 * Dunn's post-hoc test following Kruskal-Wallis.
 * Uses rank sums and z-scores with adjustable p-value correction.
 *
 * Cross-validated with R:
 * > dunn.test::dunn.test(values, groups, method = "bonferroni")
 */
export declare function dunnTest(groups: readonly GroupData[], method?: PAdjMethod): PairwiseResult[];
//# sourceMappingURL=post-hoc.d.ts.map