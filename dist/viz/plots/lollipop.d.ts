/**
 * Lollipop chart: stem + circle for each category.
 * Horizontal layout (categories on y-axis, values on x-axis).
 * Optional descending sort.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface LollipopConfig {
    readonly sorted?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface LollipopData {
    readonly labels: readonly string[];
    readonly values: readonly number[];
    readonly testResult?: StatResult;
}
export declare function renderLollipop(container: HTMLElement, data: LollipopData, config?: LollipopConfig): void;
//# sourceMappingURL=lollipop.d.ts.map