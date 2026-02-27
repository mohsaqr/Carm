/**
 * Beeswarm plot — individual points per group arranged to avoid overlap.
 * Points are sorted by value and offset horizontally so they spread
 * left/right of the group centre axis without colliding.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface SwarmPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly pointRadius?: number;
    readonly bandwidth?: number;
    readonly showN?: boolean;
    readonly showMean?: boolean;
}
export interface SwarmPlotData {
    readonly groups: readonly {
        readonly label: string;
        readonly values: readonly number[];
    }[];
    readonly testResult?: StatResult;
}
/**
 * Render a beeswarm plot with per-group colour and collision-avoidance layout.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - groups of numeric values + optional stat result
 * @param config - visual configuration
 */
export declare function renderSwarmPlot(container: HTMLElement, data: SwarmPlotData, config?: SwarmPlotConfig): void;
//# sourceMappingURL=swarm-plot.d.ts.map