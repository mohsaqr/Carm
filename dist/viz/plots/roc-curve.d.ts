/**
 * ROC curve: TPR vs FPR with AUC annotation.
 * Includes diagonal reference line, filled area under curve.
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface ROCCurveConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface ROCCurveData {
    readonly fpr: readonly number[];
    readonly tpr: readonly number[];
    readonly auc: number;
    readonly testResult?: StatResult;
}
export declare function renderROCCurve(container: HTMLElement, data: ROCCurveData, config?: ROCCurveConfig): void;
//# sourceMappingURL=roc-curve.d.ts.map