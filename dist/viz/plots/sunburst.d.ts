/**
 * Sunburst chart — radial hierarchy visualization where each ring represents
 * one level of depth. Arc angles are proportional to node values. Center hole
 * shows total. Depth-1 nodes are listed in a right-side legend.
 * Inspired by ggstatsplot philosophy: every plot tells the full statistical story.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface SunburstNode {
    readonly name: string;
    readonly value?: number;
    readonly children?: readonly SunburstNode[];
}
export interface SunburstData {
    readonly root: SunburstNode;
    readonly testResult?: StatResult;
}
export interface SunburstConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly innerRadius?: number;
    readonly maxDepth?: number;
}
export declare function renderSunburst(container: HTMLElement, data: SunburstData, config?: SunburstConfig): void;
//# sourceMappingURL=sunburst.d.ts.map