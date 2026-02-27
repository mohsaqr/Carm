/**
 * Correlation heatmap (correlogram) with significance stars and clustering.
 */
import type { CarmTheme } from '../themes/default.js';
import type { CorrelationMatrix } from '../../stats/correlation.js';
export interface CorrelogramConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showValues?: boolean;
    readonly showSignificance?: boolean;
    readonly showLegend?: boolean;
}
export declare function renderCorrelogram(container: HTMLElement, data: CorrelationMatrix, config?: CorrelogramConfig): void;
//# sourceMappingURL=correlogram.d.ts.map