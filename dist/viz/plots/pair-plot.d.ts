/**
 * Scatter matrix (pairs plot) — n×n grid of mini-plots.
 * Off-diagonal cells: scatter plot of variable i vs variable j.
 * Diagonal cells: histogram of variable i.
 * Each cell is a self-contained SVG <g> fitting in the available space.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface PairPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface PairPlotData {
    /** data[varIndex][obsIndex] — each inner array is one variable's observations */
    readonly data: readonly (readonly number[])[];
    readonly labels: readonly string[];
    readonly testResult?: StatResult;
}
/**
 * Render a scatter matrix (pairs plot).
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - variable arrays + labels + optional stat result
 * @param config - visual configuration
 */
export declare function renderPairPlot(container: HTMLElement, data: PairPlotData, config?: PairPlotConfig): void;
//# sourceMappingURL=pair-plot.d.ts.map