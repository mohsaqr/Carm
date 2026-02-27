/**
 * QQ plot (quantile-quantile plot) with confidence band.
 * Tests normality assumption visually.
 */
import type { CarmTheme } from '../themes/default.js';
export interface QQPlotConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showCI?: boolean;
}
export declare function renderQQPlot(container: HTMLElement, values: readonly number[], config?: QQPlotConfig): void;
//# sourceMappingURL=qq-plot.d.ts.map