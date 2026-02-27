/**
 * Raincloud plot: half-violin + jitter + box.
 * Premium visualization for showing full distribution shape.
 */
import type { CarmTheme } from '../themes/default.js';
import type { GroupData, StatResult } from '../../core/types.js';
export interface RaincloudConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showN?: boolean;
    readonly showMean?: boolean;
    readonly showJitter?: boolean;
}
export interface RaincloudData {
    readonly groups: readonly GroupData[];
    readonly testResult?: StatResult;
}
export declare function renderRaincloud(container: HTMLElement, data: RaincloudData, config?: RaincloudConfig): void;
//# sourceMappingURL=raincloud.d.ts.map