/**
 * Regression analysis module.
 * Simple and multiple OLS, logistic regression, polynomial regression,
 * diagnostics (R², AIC, BIC, VIF, residual plots).
 */
import type { RegressionResult } from '../core/types.js';
/**
 * Simple linear regression: y = β₀ + β₁·x
 */
export declare function linearRegression(x: readonly number[], y: readonly number[], ciLevel?: number): RegressionResult;
/**
 * Multiple linear regression: y = β₀ + β₁·x₁ + ... + βₖ·xₖ
 * `predictors`: named columns { name: values[] }
 */
export declare function multipleRegression(y: readonly number[], predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>, ciLevel?: number): RegressionResult;
/**
 * Polynomial regression: y = β₀ + β₁·x + β₂·x² + ... + βₖ·xᵏ
 */
export declare function polynomialRegression(x: readonly number[], y: readonly number[], degree: number, ciLevel?: number): RegressionResult;
/**
 * Binary logistic regression via IRLS (iteratively reweighted least squares).
 * Outcome y must be 0/1.
 *
 * Implements Fisher scoring matching R's glm.fit exactly:
 *   Working response: z = η + (y - μ) / w   where w = μ(1-μ)
 *   Solve WLS: (X'WX) β_new = X'Wz
 *   Equivalently: (Xw'Xw) β_new = Xw' zw   with Xw = √W·X, zw = √W·z
 *
 * Cross-validated with R:
 * > glm(y ~ x1 + x2, family = binomial, data = df)
 *
 * Cross-validated with Python:
 * > import statsmodels.api as sm
 * > sm.GLM(y, sm.add_constant(X), family=sm.families.Binomial()).fit()
 */
export declare function logisticRegression(y: readonly number[], predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>, ciLevel?: number, maxIter?: number, tol?: number): RegressionResult;
export interface RegressionDiagnostics {
    readonly leverage: readonly number[];
    readonly cooksDistance: readonly number[];
    readonly standardizedResiduals: readonly number[];
    readonly vif: readonly number[];
}
/**
 * Compute regression diagnostics.
 * Returns leverage (hat values), Cook's distance, standardized residuals, VIF.
 */
export declare function regressionDiagnostics(result: RegressionResult, predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>): RegressionDiagnostics;
//# sourceMappingURL=regression.d.ts.map