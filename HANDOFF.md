# Session Handoff — 2026-02-27

## Completed (This Session)
- **Build system replaced**: tsup → esbuild+tsc (`build.mjs`)
  - `npm run build` now works on Node.js v25 (tsup/consola CJS interop was broken)
  - ESM + CJS bundles with code splitting, plus `tsc --emitDeclarationOnly` for `.d.ts` + declaration maps
  - All 12 expected output files generated: `{index,core/index,stats/index,viz/index}.{js,cjs,d.ts}`
- **Vitest upgraded**: 3.2.4 → 4.0.18
  - Fixes picomatch `scan is not a function` crash on Node.js v25 (CJS interop breakage in tinyglobby)
  - All 496 tests pass across 17 suites
- **Stats improvements committed**:
  - `core/math.ts`: `normalSurvival`, `erfc`, `pKendallExact` (Best & Gipps 1974), `pSpearmanExact` (AS89 + Edgeworth), `ptukeyApprox` (Copenhaver & Holland 1988)
  - `stats/correlation.ts`: exact p-values for Kendall's tau and Spearman's rho
  - `stats/pca.ts`: PCA computation enhancements
  - `stats/post-hoc.ts`: Tukey HSD using proper studentized range distribution
  - `stats/regression.ts`: regression diagnostics improvements
- **Aistatia integration verified**: symlink is live, `tsc --noEmit` + `vite build` both clean

## Current State

**ALL 496 UNIT TESTS PASSING (17/17 suites)**

| Suite                          | Tests | Status  |
|--------------------------------|-------|---------|
| core/math                      | 41    | PASS    |
| core/matrix                    | 18    | PASS    |
| stats/analyze                  | 28    | PASS    |
| stats/clustering               | 58    | PASS    |
| stats/clustering-engagement    | 16    | PASS    |
| stats/clustering-xval          | 44    | PASS    |
| stats/comparison               | 25    | PASS    |
| stats/correlation              | 21    | PASS    |
| stats/dbscan                   | 22    | PASS    |
| stats/descriptive              | 29    | PASS    |
| stats/hac                      | 29    | PASS    |
| stats/mixed                    | 20    | PASS    |
| stats/pca                      | 14    | PASS    |
| stats/preprocess               | 18    | PASS    |
| stats/regression               | 24    | PASS    |
| large_sample                   | 47    | PASS    |
| stress                         | 42    | PASS    |
| **Total**                      | **496** | **ALL PASS** |

**169/169 R cross-validation metrics also passing (7 harnesses)**

- Build: `npm run build` (runs `node build.mjs`)
- Tests: `npm test` (vitest 4.0.18)
- Aistatia: `file:../JStats` symlink, builds clean

## Key Decisions
- Replaced tsup with `build.mjs` (esbuild + tsc) rather than patching node_modules — permanent fix
- Kept tsup in devDependencies and `tsup.config.ts` for reference if upstream fixes land
- Upgraded vitest to v4 (major version bump) — no test API changes needed, all tests pass as-is
- Quality-based GMM comparison for cross-validation (non-convex optimization)

## Open Issues
- None blocking. All build, test, and integration pipelines are green.
- LTA cross-validation still pending (depmixS4 is non-deterministic)
- LMM Satterthwaite df approximation is crude

## Context
- Branch: `main`
- Node.js v25.5.0 (only version available)
- R 4.5.1 with Saqrlab, dunn.test, rstatix, lme4, mclust, lavaan, ppcor, car
- Build: `npm run build` → `node build.mjs` (esbuild ESM+CJS + tsc DTS)
- Tests: `npm test` → vitest 4.0.18
- Validation harnesses: `npx tsx validation/ts-harness/<name>-report.ts`
- R references: `Rscript validation/r-reference/<name>-ref.R`
- Aistatia: `../aistatia`, depends on carm via `file:../JStats` symlink
