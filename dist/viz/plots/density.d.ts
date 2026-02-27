/**
 * Density plot: KDE curves for one or more series with optional rug ticks.
 * Uses Silverman's rule-of-thumb bandwidth: 1.06 * σ * n^(-1/5)
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface DensityConfig {
    readonly bandwidth?: number;
    readonly showRug?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showLegend?: boolean;
}
export interface DensitySeries {
    readonly label: string;
    readonly values: readonly number[];
}
export interface DensityData {
    readonly series: readonly DensitySeries[];
    readonly testResult?: StatResult;
}
export declare function renderDensity(container: HTMLElement, data: DensityData, config?: DensityConfig): void;
//# sourceMappingURL=density.d.ts.map