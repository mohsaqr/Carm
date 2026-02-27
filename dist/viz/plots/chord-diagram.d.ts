/**
 * Chord diagram: shows flows between groups as ribbons inside a circular layout.
 * Each group occupies an arc on the outer ring; ribbons inside represent flows.
 * Based on the D3 chord layout (d3.chord / d3.ribbon / d3.arc).
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface ChordDiagramData {
    readonly matrix: readonly (readonly number[])[];
    readonly labels: readonly string[];
    readonly testResult?: StatResult;
}
export interface ChordDiagramConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly innerRadius?: number;
    readonly padAngle?: number;
}
export declare function renderChordDiagram(container: HTMLElement, data: ChordDiagramData, config?: ChordDiagramConfig): void;
//# sourceMappingURL=chord-diagram.d.ts.map