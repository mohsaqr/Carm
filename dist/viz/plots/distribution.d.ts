/**
 * Interactive distribution explorer.
 * Shows PDF/CDF for normal, t, chi-square, F, binomial, Poisson distributions.
 * Uses jStat for accurate distribution functions.
 */
import type { CarmTheme } from '../themes/default.js';
export type DistributionName = 'normal' | 't' | 'chi-square' | 'F' | 'uniform' | 'exponential';
export interface DistributionParams {
    readonly distribution: DistributionName;
    readonly params: Readonly<Record<string, number>>;
    readonly showPDF?: boolean;
    readonly showCDF?: boolean;
    readonly highlightX?: number;
}
export interface DistributionConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export declare function renderDistribution(container: HTMLElement, params: DistributionParams, config?: DistributionConfig): void;
//# sourceMappingURL=distribution.d.ts.map