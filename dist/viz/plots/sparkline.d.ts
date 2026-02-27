/**
 * Sparkline — minimal inline trend line for dashboard embedding.
 * Tiny margins, no axes, no labels, no grid. Just the line, an optional
 * shaded area, and a highlighted end-point dot.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface SparklineConfig {
    readonly showArea?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface SparklineData {
    readonly values: readonly number[];
    readonly testResult?: StatResult;
}
/**
 * Render a sparkline (minimal inline trend line).
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - numeric time series values + optional stat result
 * @param config - visual configuration (very minimal by design)
 */
export declare function renderSparkline(container: HTMLElement, data: SparklineData, config?: SparklineConfig): void;
//# sourceMappingURL=sparkline.d.ts.map