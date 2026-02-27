/**
 * Coefficient dot-and-whisker plot (ggcoefstats style).
 * Shows regression/LMM coefficients with CI bars and p-value annotations.
 */
import type { CarmTheme } from '../themes/default.js';
import type { RegressionCoef } from '../../core/types.js';
export interface CoefPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showZeroLine?: boolean;
    readonly excludeIntercept?: boolean;
}
export declare function renderCoefPlot(container: HTMLElement, coefficients: readonly RegressionCoef[], config?: CoefPlotConfig): void;
//# sourceMappingURL=coef-plot.d.ts.map