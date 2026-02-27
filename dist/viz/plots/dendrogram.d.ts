/**
 * Dendrogram visualization for hierarchical agglomerative clustering.
 * Renders a U-shaped elbow dendrogram with cluster coloring, cut line,
 * leaf labels, and hover tooltips.
 *
 * Uses manual coordinate computation (not d3.hierarchy) because dendrograms
 * need continuous y-axis (merge height) and specific DFS leaf ordering,
 * which d3.cluster/d3.tree cannot represent correctly.
 */
import type { CarmTheme } from '../themes/default.js';
import type { HACMerge } from '../../stats/clustering.js';
export interface DendrogramData {
    readonly merges: readonly HACMerge[];
    readonly heights: readonly number[];
    readonly dendrogramOrder: readonly number[];
    readonly labels: readonly number[];
    readonly k: number;
    readonly linkage: string;
    readonly copheneticCorrelation: number;
    readonly observationLabels?: readonly string[];
}
export interface DendrogramConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly orientation?: 'vertical' | 'horizontal';
    readonly showCutLine?: boolean;
    readonly showLabels?: boolean;
    readonly colorSubtrees?: boolean;
}
/**
 * Render a dendrogram into the given container.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - merge data from runHierarchical + cutTree
 * @param config - visual configuration
 */
export declare function renderDendrogram(container: HTMLElement, data: DendrogramData, config?: DendrogramConfig): void;
//# sourceMappingURL=dendrogram.d.ts.map