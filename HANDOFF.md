# Session Handoff — 2026-02-26

## Completed
- **FA cross-validation: 100/100 synthetic datasets PASS** against R's `psych::fa()`
- **Real dataset verified**: 525×31 survey data (LOC/CCA/ER/FSI/TW scales) matches R for k=3,4,5,6 — loadings to 4 decimal places
- **Root cause fixed**: `psych::fa(rotate="promax")` dispatches via `psych::kaiser()` wrapper which adds outer Kaiser normalization before Promax. Without this, promax target Q is computed on different-scale loadings, causing MAE=0.055 on dataset 86.
- **Three major algorithmic improvements in `src/stats/factor-analysis.ts`**:
  1. Promax rotation: outer Kaiser normalization (normalize rows by communality → varimax+promax → denormalize)
  2. Varimax rotation: SVD-based algorithm matching R's `stats::varimax()` with polar decomposition via eigendecomposition of B'B
  3. ML extraction: hybrid Jöreskog gradient descent + Nelder-Mead polish on R's concentrated ML objective

## Current State
- `npx tsc --noEmit`: 0 errors
- `npx tsup`: Build success
- Cross-validation: 100/100 synthetic + real dataset verified
- Debug logging removed, clean production build
- Changes NOT committed

## Key Decisions
- **Outer Kaiser normalization in promax**: Matches `psych::fa`'s `kaiser()` wrapper exactly. This is the correct behavior since `psych::fa` is the standard R reference.
- **Eigendecomposition-based polar decomposition**: Avoids Carm's buggy `Matrix.svd()` for small k×k matrices. Computes B'B eigenvalues/vectors instead.
- **Hybrid ML extraction**: Jöreskog gradient is reliable for convergence, Nelder-Mead polishes to match R's exact concentrated ML objective.

## Open Issues
- CFA cross-validation against lavaan not yet done
- No vitest unit tests for FA (only R cross-validation scripts)
- aistatia FA integration plan exists but not yet implemented

## Next Steps
1. Write vitest tests encoding the R-verified expected values
2. CFA cross-validation against lavaan
3. Implement aistatia FA integration (plan in `cosmic-nibbling-boole.md`)
4. Edge case tests: single factor, Heywood cases, perfect correlation

## Context
- Carm dir: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/JStats`
- Aistatia dir: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/aistatia`
- Branch: `dev-clean`
- Build: `npx tsup`
- Type check: `npx tsc --noEmit`
- Cross-validation: `cd aistatia && npx tsx tmp/fa-crossval-report.ts` (100 datasets)
- Real dataset: `~/Downloads/rraw_dataaw_data.csv` (525×31)
- R reference data: `aistatia/tmp/fa-crossval-data.json`, `aistatia/tmp/fa-real-crossval-ref.json`
