/**
 * Forest plot: per-study CI lines with square points + pooled diamond.
 * Vertical reference line at 0. Numeric estimates shown on the right.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface ForestPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface ForestStudy {
    readonly label: string;
    readonly estimate: number;
    readonly ciLow: number;
    readonly ciHigh: number;
    readonly weight?: number;
}
export interface ForestPooled {
    readonly estimate: number;
    readonly ciLow: number;
    readonly ciHigh: number;
}
export interface ForestPlotData {
    readonly studies: readonly ForestStudy[];
    readonly pooled?: ForestPooled;
    readonly testResult?: StatResult;
}
export declare function renderForestPlot(container: HTMLElement, data: ForestPlotData, config?: ForestPlotConfig): void;
//# sourceMappingURL=forest-plot.d.ts.map