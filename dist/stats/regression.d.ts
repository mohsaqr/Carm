/**
 * Regression analysis module.
 * Simple and multiple OLS, logistic regression, polynomial regression,
 * diagnostics (R², AIC, BIC, VIF, residual plots).
 */
import { Matrix } from '../core/matrix.js';
import type { RegressionResult, OrdinalRegressionResult } from '../core/types.js';
/**
 * Compute residual sum of squares from a design matrix and response.
 * RSS = y'y - y'X(X'X)⁻¹X'y = ||y - X·β̂||²
 * Used by two-way ANOVA and ANCOVA for Type II/III SS via model comparison.
 */
export declare function computeRSS(X: Matrix, y: readonly number[]): number;
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
/**
 * Poisson regression via IRLS (iteratively reweighted least squares).
 * Link: log(μ) = Xβ,  Variance: V(μ) = μ
 * Outcome y must be non-negative (ideally integer counts).
 *
 * Implements Fisher scoring matching R's glm(family = poisson):
 *   Working response: z = η + (y - μ) / μ
 *   Weights: W = diag(μ)
 *   Solve WLS: (X'WX) β_new = X'Wz
 *
 * Cross-validated with R:
 * > glm(y ~ x1 + x2, family = poisson, data = df)
 * > coef(mod); summary(mod)$coefficients; AIC(mod); deviance(mod)
 */
export declare function poissonRegression(y: readonly number[], predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>, ciLevel?: number, maxIter?: number, tol?: number): RegressionResult;
/**
 * Quasi-Poisson regression.
 * Thin wrapper around poissonRegression that estimates the dispersion parameter
 * phi = Pearson chi^2/(n-p) and scales standard errors by sqrt(phi).
 * AIC/BIC are NaN (not defined for quasi-likelihood).
 *
 * Cross-validated with R:
 * > glm(y ~ x, family = quasipoisson, data = df)
 * > summary(mod)$dispersion
 */
export declare function quasiPoissonRegression(y: readonly number[], predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>, ciLevel?: number, maxIter?: number, tol?: number): RegressionResult & {
    readonly dispersion: number;
};
/**
 * Negative binomial regression via IRLS.
 * Link: log(mu) = X*beta, Variance: V(mu) = mu + mu^2/theta
 * Uses outer loop for theta: method-of-moments init, Newton steps on profile logLik.
 *
 * Cross-validated with R:
 * > library(MASS)
 * > glm.nb(y ~ x1 + x2, data = df)
 */
export declare function negativeBinomialRegression(y: readonly number[], predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>, ciLevel?: number, maxIter?: number, tol?: number): RegressionResult & {
    readonly theta: number;
};
/**
 * Ordinal logistic regression (proportional odds / cumulative logit model).
 * P(Y <= j | x) = logistic(alpha_j - x'beta), j = 1,...,J-1
 * Uses Fisher scoring (Newton-Raphson on the logLik).
 *
 * Cross-validated with R:
 * > library(MASS)
 * > polr(factor(y) ~ x1 + x2, method = "logistic", data = df)
 */
export declare function ordinalLogisticRegression(y: readonly number[], // integer categories 1..J
predictors: ReadonlyArray<{
    name: string;
    values: readonly number[];
}>, ciLevel?: number, maxIter?: number, tol?: number): OrdinalRegressionResult;
//# sourceMappingURL=regression.d.ts.map