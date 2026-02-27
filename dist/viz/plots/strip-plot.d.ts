/**
 * Strip plot: pure jittered dot plot per group.
 * Shows individual data points with a horizontal mean line per group.
 * No violin, no box — raw distribution only.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface StripPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showN?: boolean;
}
export interface StripGroup {
    readonly label: string;
    readonly values: readonly number[];
}
export interface StripPlotData {
    readonly groups: readonly StripGroup[];
    readonly testResult?: StatResult;
}
export declare function renderStripPlot(container: HTMLElement, data: StripPlotData, config?: StripPlotConfig): void;
//# sourceMappingURL=strip-plot.d.ts.map