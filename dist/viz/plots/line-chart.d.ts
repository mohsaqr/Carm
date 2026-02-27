/**
 * Multi-line chart with optional area fill.
 * Each series has its own line + color. Points at each data location.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface LineChartConfig {
    readonly showArea?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface LineChartSeries {
    readonly label: string;
    readonly x: readonly number[];
    readonly y: readonly number[];
}
export interface LineChartData {
    readonly series: readonly LineChartSeries[];
    readonly testResult?: StatResult;
}
export declare function renderLineChart(container: HTMLElement, data: LineChartData, config?: LineChartConfig): void;
//# sourceMappingURL=line-chart.d.ts.map