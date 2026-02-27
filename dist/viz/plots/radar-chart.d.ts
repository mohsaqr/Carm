/**
 * Radar / spider chart — one polygon per series, one spoke per axis.
 * Values are normalized per axis to [0,1] using min/max across all series.
 * Spokes are uniformly spaced angularly. Filled polygons with opacity.
 * Axis labels at spoke tips. Legend for series.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface RadarChartConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly levels?: number;
}
export interface RadarChartData {
    readonly series: readonly {
        readonly label: string;
        readonly values: readonly number[];
    }[];
    readonly axes: readonly string[];
    readonly testResult?: StatResult;
}
/**
 * Render a radar / spider chart.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - series (each with values per axis) + axis names + optional stat result
 * @param config - visual configuration
 */
export declare function renderRadarChart(container: HTMLElement, data: RadarChartData, config?: RadarChartConfig): void;
//# sourceMappingURL=radar-chart.d.ts.map