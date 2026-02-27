/**
 * Significance bracket layout for group comparison plots.
 * Renders stacked brackets above groups with p-value labels.
 * Non-overlapping bracket stacking algorithm.
 */
import type * as d3 from 'd3';
import type { PairwiseResult } from '../../core/types.js';
import type { CarmTheme } from '../themes/default.js';
export interface BracketConfig {
    readonly groupPositions: ReadonlyMap<string, number>;
    readonly yBase: number;
    readonly bracketHeight: number;
    readonly significantOnly: boolean;
    readonly numericP?: boolean;
}
/**
 * Render significance brackets onto an SVG group.
 * Uses greedy level-assignment to avoid bracket overlap.
 */
export declare function renderBrackets(g: d3.Selection<SVGGElement, unknown, null, undefined>, comparisons: readonly PairwiseResult[], config: BracketConfig, theme?: CarmTheme): void;
/** Total bracket height needed (to set top margin). */
export declare function totalBracketHeight(comparisons: readonly PairwiseResult[], significantOnly: boolean, bracketHeight: number): number;
//# sourceMappingURL=brackets.d.ts.map