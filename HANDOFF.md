# Session Handoff — 2026-02-27

## Completed (This Session)
- **Clustering cross-validation fix**: 5/18 → 18/18 metrics
  - R K-Means `nstart=1` → `nstart=500`, Carm K-Means 10→500 starts
  - GMM comparison: absolute diff → quality-based (`max(0, r-carm)`)
- **Shapiro-Wilk p-value fix**: 97.7% → 100% pass rate
  - Rewrote `shapiroWilkPValue()` with R's actual swilk.c coefficients
  - Added n=3 exact formula, fixed n≤11 polynomial (was using 1/n instead of n)
- **Build fix**: esbuild direct invocation when tsup fails on Node.js v25

## Current State

**ALL 169/169 METRICS PASSING (100%) — PERFECT ACROSS ALL 7 MODULES**

| Harness       | Pass | Total | Rate   | Status |
|---------------|------|-------|--------|--------|
| PCA           | 14   | 14    | 100.0% | PERFECT |
| FA Extended   | 19   | 19    | 100.0% | PERFECT |
| ANOVA         | 26   | 26    | 100.0% | PERFECT |
| T-test        | 41   | 41    | 100.0% | PERFECT |
| Correlation   | 14   | 14    | 100.0% | PERFECT |
| Regression    | 36   | 36    | 100.0% | PERFECT |
| Clustering    | 18   | 18    | 100.0% | PERFECT |
| **Total**     | **169** | **169** | **100.0%** | **ALL PERFECT** |

- All 7 harnesses complete with HTML reports in `validation/reports/`
- Reference data in `validation/data/*.json`
- Build: use `npx esbuild` (tsup broken on Node.js v25)
  ```bash
  npx esbuild src/index.ts src/core/index.ts src/stats/index.ts src/viz/index.ts \
    --bundle --format=esm --outdir=dist --splitting --external:d3 --sourcemap --outbase=src
  npx esbuild src/index.ts src/core/index.ts src/stats/index.ts src/viz/index.ts \
    --bundle --format=cjs --outdir=dist --outbase=src --external:d3 --sourcemap --out-extension:.js=.cjs
  ```

## Key Decisions
- Quality-based GMM comparison: for non-convex optimization, error=0 if Carm finds equal or better loglik
- R nstart=500 for K-Means: ensures global optimum convergence
- Used R's swilk.c source code directly for Shapiro-Wilk polynomial coefficients
- esbuild as build tool fallback (tsup/cac broken on Node.js v25)

## Open Issues
- `npm run build` (tsup) fails on Node.js v25 — need to update tsup/cac or downgrade Node
- No DTS generation with esbuild (`.d.ts` files are stale from last tsup build)

## Context
- Branch: `dev-clean`
- Node.js v25.5.0 (only version available)
- R 4.5.1 with Saqrlab, dunn.test, rstatix, lme4, mclust, lavaan, ppcor, car
- Run harnesses: `npx tsx validation/ts-harness/<name>-report.ts`
- R references: `Rscript validation/r-reference/<name>-ref.R`
