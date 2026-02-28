/**
 * Linear Mixed Models (LMM) via REML or ML.
 * Model: y = Xβ + Zb + ε
 *   b ~ N(0, G), ε ~ N(0, σ²_e · I)
 * G is the random-effects covariance (intercept + optional slopes).
 *
 * Parameterization: log-Cholesky for G/σ²_e (lme4 style).
 * Optimization: multi-start Nelder-Mead on profiled log-likelihood.
 *
 * Supports:
 *   - Random intercepts (default)
 *   - Random slopes via Cholesky parameterization of G
 *   - REML (default) or ML estimation
 *   - Nakagawa R² (marginal and conditional)
 *   - Model comparison via likelihood ratio test (compareLMM)
 *
 * Cross-validated with R:
 * > lme4::lmer(y ~ x + (1|group), data = df, REML = TRUE)
 * > lme4::lmer(y ~ x + (1 + x|group), data = df, REML = TRUE)
 */
import type { LMMResult } from '../core/types.js';
export interface LMMInput {
    readonly outcome: readonly number[];
    readonly fixedPredictors: Readonly<Record<string, readonly number[]>>;
    readonly groupId: readonly (string | number)[];
    readonly randomSlopes?: readonly string[];
    readonly method?: 'REML' | 'ML';
    readonly ciLevel?: number;
}
export interface LMMComparison {
    readonly chiSq: number;
    readonly df: number;
    readonly pValue: number;
    readonly preferred: 'model1' | 'model2';
    readonly warning?: string;
}
/**
 * Fit a linear mixed model with random intercepts and optional random slopes.
 *
 * Cross-validated with R lme4:
 * > mod <- lmer(y ~ x + (1|group), data = df, REML = TRUE)
 * > mod <- lmer(y ~ x + (1 + x|group), data = df, REML = TRUE)
 * > fixef(mod); VarCorr(mod); logLik(mod); AIC(mod)
 * > MuMIn::r.squaredGLMM(mod)
 */
export declare function runLMM(input: LMMInput): LMMResult;
/**
 * Compare two LMM models via likelihood ratio test (LRT).
 * Both models should be fitted with ML for valid comparison.
 * The model with more parameters is treated as the "full" model.
 *
 * Cross-validated with R:
 * > anova(mod1, mod2)
 *
 * @returns LRT statistic, df, p-value, and which model is preferred
 */
export declare function compareLMM(model1: LMMResult, model2: LMMResult): LMMComparison;
/**
 * Compute BLUPs (Best Linear Unbiased Predictors) — the random intercepts.
 * b_hat = σ²_b Z'V^{-1}(y - Xβ)
 */
export declare function computeBLUPs(input: LMMInput, result: LMMResult): ReadonlyArray<{
    group: string | number;
    blup: number;
}>;
//# sourceMappingURL=mixed.d.ts.map