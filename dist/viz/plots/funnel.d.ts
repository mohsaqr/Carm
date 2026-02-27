/**
 * Funnel chart: trapezoidal segments stacked vertically.
 * Width proportional to value. Labels on left, value + % drop on right.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface FunnelConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface FunnelStage {
    readonly label: string;
    readonly value: number;
}
export interface FunnelData {
    readonly stages: readonly FunnelStage[];
    readonly testResult?: StatResult;
}
export declare function renderFunnel(container: HTMLElement, data: FunnelData, config?: FunnelConfig): void;
//# sourceMappingURL=funnel.d.ts.map