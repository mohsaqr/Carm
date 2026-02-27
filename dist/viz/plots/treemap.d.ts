/**
 * Treemap — rectangle area proportional to value, laid out using D3's
 * squarified treemap algorithm. Cells are labelled (truncated when small)
 * and coloured by group when provided, otherwise by index.
 * White borders separate cells. Tooltip shows label, value, and %.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface TreemapConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface TreemapChild {
    readonly label: string;
    readonly value: number;
    readonly group?: string;
}
export interface TreemapData {
    readonly children: readonly TreemapChild[];
    readonly testResult?: StatResult;
}
/**
 * Render a treemap using D3's squarified layout.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - leaf nodes with labels, values, optional groups + stat result
 * @param config - visual configuration
 */
export declare function renderTreemap(container: HTMLElement, data: TreemapData, config?: TreemapConfig): void;
//# sourceMappingURL=treemap.d.ts.map