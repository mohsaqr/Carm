/**
 * Hierarchical edge bundling (Holten 2006).
 *
 * Leaf nodes are arranged on a circle. Non-leaf nodes form the hierarchy used
 * to bundle edges — they are not drawn. Edges between leaf nodes are drawn as
 * cubic B-spline paths routed through the common ancestor path in the
 * hierarchy, which causes edges that share ancestry to visually cluster
 * ("bundle") together.
 *
 * Layout:
 *   d3.stratify builds the hierarchy from the flat node list.
 *   d3.cluster assigns angular positions (x) and radial depth (y).
 *   Only leaf nodes (no children) are shown as circles with radial labels.
 *   Each edge walks the path [source → LCA → target] through the tree and
 *   passes those anchor points to d3.lineRadial with d3.curveBundle.beta(β).
 *
 * Colour:
 *   Nodes are coloured by their `group` field if provided, otherwise by their
 *   immediate parent id. The same colour mapping is reused for edges (source
 *   colour at low opacity so overlapping edges remain readable).
 *
 * Interactivity:
 *   Hovering a node highlights all its edges (full opacity, thicker stroke)
 *   and dims all other edges. Clicking a node locks the highlight; clicking
 *   again unlocks it.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface BundleNode {
    readonly id: string;
    readonly parent: string;
    readonly label?: string;
    readonly group?: string;
}
export interface BundleEdge {
    readonly source: string;
    readonly target: string;
}
export interface EdgeBundlingData {
    readonly nodes: readonly BundleNode[];
    readonly edges: readonly BundleEdge[];
    readonly testResult?: StatResult;
}
export interface EdgeBundlingConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly bundlingStrength?: number;
    readonly nodeRadius?: number;
}
/**
 * Render a hierarchical edge bundling diagram.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data      - nodes (with parent hierarchy), edges between leaves, optional stat result
 * @param config    - visual configuration
 */
export declare function renderEdgeBundling(container: HTMLElement, data: EdgeBundlingData, config?: EdgeBundlingConfig): void;
//# sourceMappingURL=edge-bundling.d.ts.map