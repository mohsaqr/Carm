/**
 * Cleveland dot plot: horizontal layout with labels on y-axis.
 * Supports single dots or paired dots connected by a line (group comparison).
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface DotPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface DotPlotData {
    readonly labels: readonly string[];
    readonly values: readonly number[];
    readonly group2?: readonly number[];
    readonly group1Label?: string;
    readonly group2Label?: string;
    readonly testResult?: StatResult;
}
export declare function renderDotPlot(container: HTMLElement, data: DotPlotData, config?: DotPlotConfig): void;
//# sourceMappingURL=dot-plot.d.ts.map