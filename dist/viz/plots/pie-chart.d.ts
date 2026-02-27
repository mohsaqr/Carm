/**
 * Pie / donut chart with polyline labels and optional test result subtitle.
 * Donut mode uses innerRadius = outerRadius * 0.55.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface PieChartConfig {
    readonly donut?: boolean;
    readonly showPercentages?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface PieSlice {
    readonly label: string;
    readonly value: number;
}
export interface PieChartData {
    readonly slices: readonly PieSlice[];
    readonly testResult?: StatResult;
}
export declare function renderPieChart(container: HTMLElement, data: PieChartData, config?: PieChartConfig): void;
//# sourceMappingURL=pie-chart.d.ts.map