/**
 * Violin + Box + Jitter combo plot (ggbetweenstats style).
 * Shows distribution shape (violin), summary (box), individual points (jitter),
 * and significance brackets between groups.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult, PairwiseResult, GroupData } from '../../core/types.js';
export interface ViolinBoxConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showJitter?: boolean;
    readonly showBrackets?: boolean;
    readonly significantBracketsOnly?: boolean;
    readonly jitterWidth?: number;
    readonly violinBandwidth?: number;
    readonly showN?: boolean;
    readonly showMean?: boolean;
    readonly showMedian?: boolean;
    readonly numericP?: boolean;
}
export interface ViolinBoxData {
    readonly groups: readonly GroupData[];
    readonly testResult?: StatResult;
    readonly pairwise?: readonly PairwiseResult[];
}
/**
 * Render a violin + box + jitter plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - group data + optional stat results
 * @param config - visual configuration
 */
export declare function renderViolinBox(container: HTMLElement, data: ViolinBoxData, config?: ViolinBoxConfig): void;
//# sourceMappingURL=violin-box.d.ts.map