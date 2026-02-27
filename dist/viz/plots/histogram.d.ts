/**
 * Histogram with density overlay and optional normality curve.
 */
import type { CarmTheme } from '../themes/default.js';
import type { DescriptiveResult } from '../../core/types.js';
export interface HistogramConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly bins?: number;
    readonly showDensity?: boolean;
    readonly showNormalCurve?: boolean;
    readonly color?: string;
}
export interface HistogramData {
    readonly values: readonly number[];
    readonly descriptives?: DescriptiveResult;
}
export declare function renderHistogram(container: HTMLElement, data: HistogramData, config?: HistogramConfig): void;
//# sourceMappingURL=histogram.d.ts.map