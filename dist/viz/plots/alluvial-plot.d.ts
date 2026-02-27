/**
 * Alluvial (Sankey-like) plot — nodes stacked by stage, flows as cubic bezier
 * ribbons. Every stage is a vertical column; node height is proportional to
 * total flow value. Ribbons connect adjacent stages and are coloured by source
 * node colour at reduced opacity so layering is visible.
 *
 * Layout:
 *   stages are evenly spaced across the plot width.
 *   within each stage nodes are sorted largest→smallest, stacked top→bottom
 *   with a configurable gap between them.
 *   each node tracks separate y-offsets for outgoing and incoming ribbons so
 *   ribbons stack neatly without crossing inside a node.
 *
 * Colour: all unique node labels are collected up-front and assigned a fixed
 * colour from the theme palette so the same category always has the same
 * colour regardless of which stage it appears in.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface AlluvialNode {
    readonly id: string;
    readonly stage: number;
    readonly label: string;
}
export interface AlluvialFlow {
    readonly source: string;
    readonly target: string;
    readonly value: number;
}
export interface AlluvialData {
    readonly nodes: readonly AlluvialNode[];
    readonly flows: readonly AlluvialFlow[];
    readonly stageLabels?: readonly string[];
    readonly testResult?: StatResult;
}
export interface AlluvialConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly nodePadding?: number;
    readonly nodeWidth?: number;
}
/**
 * Render an alluvial / Sankey-style plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data      - nodes, flows, optional stage labels and stat result
 * @param config    - visual configuration
 */
export declare function renderAlluvialPlot(container: HTMLElement, data: AlluvialData, config?: AlluvialConfig): void;
//# sourceMappingURL=alluvial-plot.d.ts.map