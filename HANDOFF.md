# Session Handoff — 2026-02-28

## Completed This Session

### EFA Geomin Rotation: Model-Implied Standardization + Column Reflection
Implemented the planned lavaan-matching changes:

1. **Model-implied standardization in `runEFA()`** (~line 1804): Computes `ovVar_i = h²_i + ψ_i` per variable from ML extraction, standardizes loadings by `sqrt(ovVar)` before rotation, de-standardizes after. Matches lavaan's `std.ov=TRUE` default.

2. **Post-rotation column reflection in `gpfOblq()`** (~line 908): After computing Phi, flips factor columns where `sum(loadings) < 0`, including Phi row/column sign flips. Matches lavaan's default post-rotation processing.

3. **Removed `criterionGeominCov()`**: Unused function was causing `noUnusedLocals` build error. Was added as utility in previous session but never called.

4. **Updated `tmp/baseline-compare.ts`**: Rewrote with proper psychometric comparison metrics:
   - Tucker's congruence coefficient (φ) per factor
   - Communality comparison (rotation-invariant)
   - Factor correlation matrix (Φ) comparison
   - Reproduced correlation matrix comparison

### Results
- **784/784 vitest** pass
- **200/200 geomin cross-validation** pass
- **rraw dataset**: Tucker's φ all > 0.95 (min=0.975), mean primary Δ = 2.57%
- The model-implied standardization is nearly a no-op for correlation matrix input (ovVar ≈ 1.0), confirming the 2.5% diff is inherent to different GPA basins

## Current State
- **JStats build:** `node build.mjs` — clean
- **JStats tests:** `npx vitest run` — **784/784 pass**
- **Geomin cross-validation:** `npx tsx validation/ts-harness/fa-geomin-crossval.ts` — **200/200 pass**
- **lavaan comparison:** ~2.57% mean primary loading diff on rraw_dataaw_data (inherent to near-degenerate optima)

## Key Decisions
- **Kept both changes** (standardization + reflection) despite minimal impact on rraw dataset. They are algorithmically correct and match lavaan's procedure. On other datasets they could make a difference.
- **Removed `criterionGeominCov`** rather than exporting it — it was unused and the underscore prefix didn't suppress `noUnusedLocals`.
- **The 2.5% diff is definitively confirmed** as a basin selection issue. Tucker's φ > 0.95 for all factors demonstrates the solutions are psychometrically equivalent.

## Open Issues
- **lavaan matching on rraw_dataaw_data**: ~2.5% primary loading diff remains. All algorithmic alignment attempts have been exhausted. The diff is in factor 5 (weakest factor, most rotationally indeterminate).
- **Manual UI test still needed**: Aistatia EFA wizard preset UI not yet tested in browser.

## Next Steps
1. Test the Aistatia EFA wizard preset UI in the browser
2. Consider adding `rotationScale: 'correlation' | 'covariance'` option to `EFAOptions`
3. Any remaining items from the broader project roadmap

## Context
- JStats (Carm) repo: `/Users/mohammedsaqr/Documents/Github/JStats`
- Aistatia: `/Users/mohammedsaqr/Documents/Github/aistatia/`
- Build JStats: `node build.mjs`
- Test JStats: `NODE_OPTIONS='--max-old-space-size=4096' npx vitest run` (784 tests)
- Build Aistatia: `cd ../aistatia && npx vite build`
- Dev Aistatia: `cd ../aistatia && npx vite`
- Geomin validation: `npx tsx validation/ts-harness/fa-geomin-crossval.ts` (200/200)
- lavaan comparison: `npx tsx tmp/baseline-compare.ts` (requires rraw_dataaw_data.csv in ~/Downloads)
