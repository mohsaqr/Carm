/**
 * Pareto chart: descending bars + cumulative % line on dual axis.
 * Includes 80% threshold line (the "vital few" cutoff).
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface ParetoConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface ParetoData {
    readonly labels: readonly string[];
    readonly values: readonly number[];
    readonly testResult?: StatResult;
}
export declare function renderPareto(container: HTMLElement, data: ParetoData, config?: ParetoConfig): void;
//# sourceMappingURL=pareto.d.ts.map