/**
 * Annotated bar chart with counts, percentages, and significance brackets.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult, FrequencyRow } from '../../core/types.js';
export interface BarStatsConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showPercentages?: boolean;
    readonly showCounts?: boolean;
    readonly orientation?: 'vertical' | 'horizontal';
}
export interface BarStatsData {
    readonly rows: readonly FrequencyRow[];
    readonly testResult?: StatResult;
}
export declare function renderBarStats(container: HTMLElement, data: BarStatsData, config?: BarStatsConfig): void;
//# sourceMappingURL=bar-stats.d.ts.map