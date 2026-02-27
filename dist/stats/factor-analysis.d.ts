/**
 * Factor Analysis module for Carm.
 * Exploratory Factor Analysis (EFA), Confirmatory Factor Analysis (CFA),
 * and psychometric diagnostics (KMO, Bartlett, MAP, Parallel Analysis).
 *
 * All functions are pure — no DOM, no D3, no side effects.
 * Uses splitmix32 PRNG for deterministic reproducibility.
 *
 * Algorithms adapted from:
 * - Iterated PAF: Gorsuch (1983), Factor Analysis (2nd ed.)
 * - ML extraction: Jöreskog (1967), Psychometrika 32(4)
 * - Varimax: Kaiser (1958), Psychometrika 23:187-200
 * - Oblimin/Promax: Jennrich (2002), Psychometrika 67(1)
 * - CFA: Bollen (1989), Structural Equations with Latent Variables
 * - Fit indices: Browne & Cudeck (1993), Hu & Bentler (1999)
 */
import type { FADiagnostics, FAResult, CFAResult } from '../core/types.js';
export interface EFAOptions {
    readonly nFactors?: number;
    readonly extraction?: 'paf' | 'ml';
    readonly rotation?: 'varimax' | 'oblimin' | 'promax' | 'quartimin' | 'geomin' | 'none';
    readonly geominDelta?: number;
    readonly seed?: number;
    readonly maxIter?: number;
    readonly tol?: number;
    readonly randomStarts?: number;
    readonly variableNames?: readonly string[];
}
export interface CFAOptions {
    readonly maxIter?: number;
    readonly tol?: number;
    readonly variableNames?: readonly string[];
    readonly factorNames?: readonly string[];
}
export interface FADiagnosticsOptions {
    readonly seed?: number;
    readonly parallelIterations?: number;
}
/**
 * Compute factor analysis diagnostics: KMO, Bartlett's test,
 * Velicer's MAP, and parallel analysis.
 *
 * Cross-validated with R:
 * > library(psych)
 * > KMO(data)
 * > cortest.bartlett(cor(data), n = nrow(data))
 * > fa.parallel(data, fm = "ml", fa = "fa")
 */
export declare function runFADiagnostics(data: readonly (readonly number[])[], options?: FADiagnosticsOptions): FADiagnostics;
/**
 * Exploratory Factor Analysis with extraction (PAF or ML) and rotation.
 *
 * Cross-validated with R:
 * > library(psych)
 * > fa(data, nfactors = 3, fm = "ml", rotate = "promax")
 * > fa(data, nfactors = 3, fm = "minres", rotate = "varimax")
 */
export declare function runEFA(data: readonly (readonly number[])[], options?: EFAOptions): FAResult;
/**
 * Confirmatory Factor Analysis via ML estimation.
 *
 * Model is specified as a record: { F1: [0, 1, 2], F2: [3, 4, 5] }
 * where keys are factor names and values are 0-indexed item indices.
 *
 * Cross-validated with R:
 * > library(lavaan)
 * > model <- 'F1 =~ x1 + x2 + x3\n F2 =~ x4 + x5 + x6'
 * > fit <- cfa(model, data)
 * > standardizedSolution(fit)
 * > fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli", "srmr"))
 */
export declare function runCFA(data: readonly (readonly number[])[], model: Readonly<Record<string, readonly number[]>>, options?: CFAOptions): CFAResult;
//# sourceMappingURL=factor-analysis.d.ts.map