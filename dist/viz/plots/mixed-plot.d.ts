/**
 * LMM visualization: caterpillar plot of BLUPs + variance component chart.
 */
import type { CarmTheme } from '../themes/default.js';
import type { LMMResult } from '../../core/types.js';
export interface MixedPlotConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export declare function renderMixedPlot(container: HTMLElement, result: LMMResult, blups: ReadonlyArray<{
    group: string | number;
    blup: number;
}>, config?: MixedPlotConfig): void;
//# sourceMappingURL=mixed-plot.d.ts.map