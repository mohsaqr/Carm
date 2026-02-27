/**
 * Standalone box plot: Q1/Q3/median/whiskers (1.5 IQR)/outliers per group.
 * No violin. Uses quantile() from core/math.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface BoxplotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showN?: boolean;
    readonly showMean?: boolean;
    readonly showOutliers?: boolean;
    readonly showMedian?: boolean;
}
export interface BoxplotGroup {
    readonly label: string;
    readonly values: readonly number[];
}
export interface BoxplotData {
    readonly groups: readonly BoxplotGroup[];
    readonly testResult?: StatResult;
}
export declare function renderBoxplot(container: HTMLElement, data: BoxplotData, config?: BoxplotConfig): void;
//# sourceMappingURL=boxplot.d.ts.map