# Session Handoff — 2026-03-01

## Completed This Session

### Logistic GLMM (Generalized Linear Mixed Model)

1. **Types** (`src/core/types.ts`): Added `GLMMFixedEffect` interface with z-tests, odds ratios (OR, orCI), and `GLMMResult` interface with latent-scale ICC, deviance, family/link fields.

2. **APA Formatter** (`src/core/apa.ts`): Added `formatGLMM()` producing `ICC (latent) = 0.42, AIC = 234.5, Deviance = 223.5`.

3. **Helper Exports** (`src/stats/mixed.ts`): Exported `buildCholFactor`, `buildExtendedZ`, `buildA` for reuse by GLMM.

4. **Core Implementation** (`src/stats/glmm.ts` — NEW): ~400 lines implementing:
   - `runGLMM()` — Laplace-approximated logistic GLMM
   - PIRLS inner loop (augmented normal equations with step-halving)
   - Nelder-Mead outer loop over log-Cholesky θ parameters
   - Multi-start initialization (θ ∈ {-4, -2, -1, 0, 1, 2})
   - Fixed effects: Wald z-tests, odds ratios with CIs
   - Variance components, random correlations
   - Latent-scale ICC = σ²_b / (σ²_b + π²/3)
   - AIC, BIC, deviance, log-likelihood

5. **Tests**: 34 new tests across 2 test files:
   - `tests/stats/glmm.test.ts`: 14 unit tests
   - `tests/stats/glmm-numerical-equivalence.test.ts`: 20 R cross-validation checks

6. **R cross-validation**: R script + JSON fixture:
   - `validation/r-reference/glmm-ref.R` → `tests/fixtures/glmm-ref.json`

7. **aistatia Integration**:
   - `src/data/types.ts`: Added `glmm` kind to `AppRunResult` union
   - `src/analysis/registry.ts`: Added `logistic-glmm` test entry
   - `src/analysis/runner.ts`: Added runner case for `logistic-glmm`
   - `src/views/wizard-defs.ts`: Added wizard config (reuses LMM group col, random slopes, CI)
   - `src/views/results.ts`: Added GLMM rendering with ICC/AIC tags
   - `src/views/layout.ts`: Added icon for `logistic-glmm`
   - `src/results/tables.ts`: Added `renderGLMMSummary()` with OR table
   - `src/results/apa-banner.ts`: Added GLMM APA string
   - `src/results/plot-panel.ts`: Added GLMM coefplot/forest rendering

## Current State

- **All 1434 tests pass** (1400 pre-existing + 34 new), zero regressions
- TypeScript compiles cleanly in both JStats and aistatia
- GLMM fixed effects match R lme4::glmer() within tolerances (0.5 for coefficients, 1.5 for variance, 5 for logLik)
- aistatia builds and renders GLMM results with OR table and coefficient plot

## Key Decisions

- **lme4 parameterization**: Used `log|Λ'Z'WZΛ + I|` for log-determinant (not `log|Z'WZ + G⁻¹|`) to avoid θ bias. This was the critical algorithmic fix.
- **z-tests** (not t-tests) for fixed effects — standard for GLMMs
- **Latent-scale ICC**: σ²_b / (σ²_b + π²/3) — matches lme4 convention
- **No REML**: GLMMs only support ML estimation
- **Wider tolerances**: GLMM R cross-validation uses wider tolerances than LMM due to Laplace approximation differences between implementations
- **Reused LMM infrastructure**: shared `buildCholFactor`, `buildExtendedZ`, multi-start strategy, wizard options (group col, random slopes)

## Open Issues

- None identified

## Next Steps

- Consider adding adaptive Gauss-Hermite quadrature (AGQ) as an alternative to Laplace for better accuracy with few groups
- Consider adding Poisson GLMM (count outcomes)
- Consider adding prediction/diagnostic functions (GLMM residuals, random effects BLUPs)
- Visualization: ICC forest plot, random effects caterpillar plot

## Context

- Node.js/TypeScript project with Vitest testing framework
- R required for generating cross-validation fixtures (lme4 package needed)
- All stats functions are pure (no DOM), all in `src/stats/`
- aistatia app in separate repo at `/Users/mohammedsaqr/Documents/Github/aistatia/`
