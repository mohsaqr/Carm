/**
 * 4-panel regression diagnostics plot:
 * 1. Residuals vs Fitted
 * 2. QQ of residuals
 * 3. Scale-Location (sqrt|residuals| vs fitted)
 * 4. Residuals vs Leverage
 */
import type { CarmTheme } from '../themes/default.js';
import type { RegressionResult } from '../../core/types.js';
export interface ResidualPanelConfig {
    readonly title?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export declare function renderResidualPanel(container: HTMLElement, result: RegressionResult, leverage: readonly number[], config?: ResidualPanelConfig): void;
//# sourceMappingURL=residual-panel.d.ts.map