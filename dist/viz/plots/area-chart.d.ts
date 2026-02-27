/**
 * Area chart with optional stacking.
 * Non-stacked: overlapping semi-transparent areas + lines per series.
 * Stacked: cumulative areas per series.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface AreaChartConfig {
    readonly stacked?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface AreaSeries {
    readonly label: string;
    readonly x: readonly number[];
    readonly y: readonly number[];
}
export interface AreaChartData {
    readonly series: readonly AreaSeries[];
    readonly testResult?: StatResult;
}
export declare function renderAreaChart(container: HTMLElement, data: AreaChartData, config?: AreaChartConfig): void;
//# sourceMappingURL=area-chart.d.ts.map