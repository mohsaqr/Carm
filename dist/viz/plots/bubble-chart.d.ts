/**
 * Bubble chart: scatter plot where bubble size is proportional to r value.
 * Color by group if group field is present. Tooltip on hover.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface BubbleChartConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface BubblePoint {
    readonly x: number;
    readonly y: number;
    readonly r: number;
    readonly label?: string;
    readonly group?: string;
}
export interface BubbleChartData {
    readonly points: readonly BubblePoint[];
    readonly testResult?: StatResult;
}
export declare function renderBubbleChart(container: HTMLElement, data: BubbleChartData, config?: BubbleChartConfig): void;
//# sourceMappingURL=bubble-chart.d.ts.map