/**
 * Parallel coordinates plot — one vertical axis per variable, evenly spaced
 * horizontally. Each data row is a polyline through its per-axis values.
 * Axes are individually min-max scaled. Rows are optionally coloured by group.
 * Lines are semi-transparent to show density.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface ParallelCoordsConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface ParallelCoordsData {
    /** rows[obsIndex][varIndex] */
    readonly rows: readonly (readonly number[])[];
    readonly axes: readonly string[];
    /** optional group index per row (0-based) for colour coding */
    readonly groups?: readonly number[];
    readonly testResult?: StatResult;
}
/**
 * Render a parallel coordinates plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - observations (rows × vars), axis names, optional groups + stat result
 * @param config - visual configuration
 */
export declare function renderParallelCoords(container: HTMLElement, data: ParallelCoordsData, config?: ParallelCoordsConfig): void;
//# sourceMappingURL=parallel-coords.d.ts.map