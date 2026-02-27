/**
 * Waffle chart — 10×10 grid of small squares where each square represents 1%
 * of the total. Squares are coloured by slice proportionally (floor-rounded,
 * remainder assigned to largest slice). Legend below shows colour, label, %.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface WaffleChartConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface WaffleChartData {
    readonly slices: readonly {
        readonly label: string;
        readonly value: number;
    }[];
    readonly testResult?: StatResult;
}
/**
 * Render a waffle chart (10×10 proportional grid).
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - named slices with values + optional stat result
 * @param config - visual configuration
 */
export declare function renderWaffleChart(container: HTMLElement, data: WaffleChartData, config?: WaffleChartConfig): void;
//# sourceMappingURL=waffle-chart.d.ts.map