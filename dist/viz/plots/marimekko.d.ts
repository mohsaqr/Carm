/**
 * Marimekko (Mekko) chart — a 100% stacked bar chart where column WIDTHS are
 * proportional to each column's total value, so both axes encode data.
 *
 * - Column width  = that category's total as a fraction of the grand total.
 * - Within each column, stacked segments sum to 100%, coloured by series.
 * - Cells receive percentage (or raw value) labels when large enough.
 * - A "spine" bar below the x-axis shows relative column widths at a glance.
 *
 * Usage:
 *   renderMarimekko(container, data, config)
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface MarimekkoData {
    /** Labels for each column (X variable). */
    readonly categories: readonly string[];
    /** One series per stacked segment. values[j] corresponds to categories[j]. */
    readonly series: readonly {
        readonly label: string;
        readonly values: readonly number[];
    }[];
    readonly testResult?: StatResult;
}
export interface MarimekkoConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    /** Show percentage labels inside cells when the cell is large enough. Default: true. */
    readonly showPercentLabels?: boolean;
    /** Show raw values instead of percentages in cell labels. Default: false. */
    readonly showValueLabels?: boolean;
}
export declare function renderMarimekko(container: HTMLElement, data: MarimekkoData, config?: MarimekkoConfig): void;
//# sourceMappingURL=marimekko.d.ts.map