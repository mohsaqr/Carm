/**
 * Scatter plot with regression line, CI band, and marginal distributions.
 * ggscatterstats style: stat result in subtitle, regression equation on plot.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult, RegressionResult } from '../../core/types.js';
export interface ScatterStatsConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showMarginals?: boolean;
    readonly showCI?: boolean;
    readonly pointSize?: number;
    readonly showEquation?: boolean;
}
export interface ScatterStatsData {
    readonly x: readonly number[];
    readonly y: readonly number[];
    readonly labels?: readonly string[];
    readonly correlationResult?: StatResult;
    readonly regressionResult?: RegressionResult;
}
export declare function renderScatterStats(container: HTMLElement, data: ScatterStatsData, config?: ScatterStatsConfig): void;
//# sourceMappingURL=scatter-stats.d.ts.map