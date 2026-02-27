/**
 * Smart axis rendering for Carm plots.
 * Auto tick formatting, minimal grid, clean style.
 */
import type * as d3 from 'd3';
import type { CarmTheme } from '../themes/default.js';
type AnyScale = d3.ScaleLinear<number, number> | d3.ScaleBand<string>;
/** Render a styled x-axis with label. */
export declare function renderXAxis(g: d3.Selection<SVGGElement, unknown, null, undefined>, _scale: AnyScale, height: number, label: string, width: number, theme?: CarmTheme): void;
/** Render a styled y-axis with label. */
export declare function renderYAxis(g: d3.Selection<SVGGElement, unknown, null, undefined>, _scale: d3.ScaleLinear<number, number>, height: number, label: string, theme?: CarmTheme): void;
/** Render horizontal grid lines. */
export declare function renderGridLines(g: d3.Selection<SVGGElement, unknown, null, undefined>, scale: d3.ScaleLinear<number, number>, width: number, theme?: CarmTheme): void;
export {};
//# sourceMappingURL=axis.d.ts.map