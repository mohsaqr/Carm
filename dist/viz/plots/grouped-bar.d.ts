/**
 * Grouped or stacked bar chart: multiple series per category.
 * type='grouped' renders bars side-by-side; type='stacked' stacks them.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface GroupedBarConfig {
    readonly type?: 'grouped' | 'stacked';
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface GroupedBarSeries {
    readonly label: string;
    readonly values: readonly number[];
}
export interface GroupedBarData {
    readonly categories: readonly string[];
    readonly series: readonly GroupedBarSeries[];
    readonly testResult?: StatResult;
}
export declare function renderGroupedBar(container: HTMLElement, data: GroupedBarData, config?: GroupedBarConfig): void;
//# sourceMappingURL=grouped-bar.d.ts.map