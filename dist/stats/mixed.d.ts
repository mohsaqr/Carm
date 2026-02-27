/**
 * Linear Mixed Models (LMM) via REML.
 * Model: y = Xβ + Zb + ε
 *   b ~ N(0, σ²_b · I), ε ~ N(0, σ²_e · I)
 * Random intercepts + optional random slopes.
 * Optimization: Nelder-Mead on REML profile log-likelihood.
 *
 * Cross-validated with R:
 * > lme4::lmer(y ~ x + (1|group), data = df, REML = TRUE)
 */
import type { LMMResult } from '../core/types.js';
export interface LMMInput {
    readonly outcome: readonly number[];
    readonly fixedPredictors: Readonly<Record<string, readonly number[]>>;
    readonly groupId: readonly (string | number)[];
    readonly randomSlopes?: readonly string[];
    readonly ciLevel?: number;
}
/**
 * Fit a linear mixed model with random intercepts (and optionally random slopes).
 *
 * Cross-validated with R lme4:
 * > mod <- lmer(y ~ x + (1|group), data = df, REML = TRUE)
 * > fixef(mod)
 * > VarCorr(mod)
 * > icc(mod)
 */
export declare function runLMM(input: LMMInput): LMMResult;
/**
 * Compute BLUPs (Best Linear Unbiased Predictors) — the random intercepts.
 * b_hat = σ²_b Z'V^{-1}(y - Xβ)
 */
export declare function computeBLUPs(input: LMMInput, result: LMMResult): ReadonlyArray<{
    group: string | number;
    blup: number;
}>;
//# sourceMappingURL=mixed.d.ts.map