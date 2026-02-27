/**
 * Arc diagram: nodes on a horizontal baseline connected by arcs above them.
 * Short connections produce shallow arcs; distant connections produce tall arcs.
 * Useful for visualising network topology when node order carries meaning.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface ArcDiagramNode {
    readonly id: string;
    readonly label?: string;
    readonly group?: string;
}
export interface ArcDiagramEdge {
    readonly source: string;
    readonly target: string;
    readonly value?: number;
}
export interface ArcDiagramData {
    readonly nodes: readonly ArcDiagramNode[];
    readonly edges: readonly ArcDiagramEdge[];
    readonly testResult?: StatResult;
}
export interface ArcDiagramConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly nodeRadius?: number;
    readonly sortByGroup?: boolean;
}
export declare function renderArcDiagram(container: HTMLElement, data: ArcDiagramData, config?: ArcDiagramConfig): void;
//# sourceMappingURL=arc-diagram.d.ts.map