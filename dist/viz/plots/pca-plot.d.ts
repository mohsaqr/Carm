/**
 * PCA visualizations: biplot, scree plot, loadings heatmap.
 */
import type { CarmTheme } from '../themes/default.js';
import type { PCAResult } from '../../core/types.js';
export interface PCAPlotConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly variableLabels?: readonly string[];
    readonly observationLabels?: readonly string[];
    readonly type?: 'biplot' | 'scree' | 'loadings';
}
export declare function renderPCAPlot(container: HTMLElement, pca: PCAResult, config?: PCAPlotConfig): void;
//# sourceMappingURL=pca-plot.d.ts.map