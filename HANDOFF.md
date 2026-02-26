# Session Handoff — 2026-02-26

## Completed
- **Random starts for GPFoblq rotation** — matches lavaan's GPA×30 strategy
  - `randomOrthogonalMatrix(k, rng)`: Haar-distributed via Modified Gram-Schmidt QR + sign correction
  - `gpfOblq()`: now accepts optional `Tinit` and returns criterion value `f`
  - `gpfOblqWithRandomStarts()`: orchestrator — T=I first, then N-1 random starts, picks lowest f
  - `applyRotation()`: delegates oblique methods through multi-start when `randomStarts > 1`
  - `runEFA()`: reads `randomStarts` from `EFAOptions`, threads to `applyRotation`
- **Default `randomStarts=50`** — empirically validated on two real datasets:
  - rraw dataset (525×31, 5 factors): 50 starts achieves MAE=0.000001 vs lavaan (perfect match)
  - Teacher burnout (876×23, 4 factors): 10 starts sufficient; 50 starts gives extra safety margin
  - Timing: 50 starts ~3s, 100 starts ~4.6s, 1000 starts ~31s (linear scaling)
- **Full cross-validation**: Promax 100/100, Geomin 100/100, Real dataset 6/6, Diagnostics 100/100
- **Comprehensive technical reports** (validation/ folder, NEVER delete):
  - `validation/FA-TECHNICAL-REPORT.md` (1,085 lines): Complete EFA/CFA coverage — ML/PAF extraction, all 6 rotations, GPFoblq engine, random starts, CFA, fit indices, diagnostics, numerical infrastructure, cross-validation results, loading comparisons, API reference, 18 references
  - `validation/RANDOM-STARTS-REPORT.md` (692 lines): Focused random starts report — Haar matrices, multi-start GPFoblq, empirical validation, convergence analysis, timing benchmarks
  - `validation/GMM-TECHNICAL-REPORT.md`: GMM/clustering technical report — EM algorithm, covariance models, K-Means++, numerical tricks, cross-validation against R mclust
- **All committed and pushed to main**:
  - `fa3c31b` — default randomStarts=50, cross-validation passing
  - `a2141ac` — FA technical reports
  - `13ff4a9` — GMM/clustering technical report (HEAD)

## Current State
- Branch: `main` at `13ff4a9`
- `npx tsc --noEmit`: 0 errors
- `npx tsup`: Build success
- `npx vitest run`: 496/496 pass
- Cross-validation: 100/100 promax, 100/100 geomin, real dataset pass, diagnostics pass
- All technical reports committed to validation/
- Unstaged changes: CHANGES.md, HANDOFF.md, LEARNINGS.md (session documentation updates)

## Key Decisions
- Default randomStarts=50 (was 100, reduced after empirical testing showed 50 sufficient for both datasets)
- Start 0 is always T=I (backward compatible, deterministic)
- Random orthogonal matrices: Haar distribution via QR of N(0,1) matrix + sign correction
- Sign correction: Q[:,j] *= sign(R[j][j]) — standard QR convention
- When Tinit != I, compute initial L = A * inv(T)' (not just copy A)

## Open Issues
- Oblimin/quartimin cross-validation against R not yet done (uses same GPFoblq, should work)
- CFA cross-validation against lavaan not yet done
- No vitest unit tests for FA (only R cross-validation scripts in validation/)
- `validation/FA-TECHNICAL-REPORT.html` is untracked (rendered version of the MD report)

## Next Steps
1. Write vitest tests encoding R-verified expected values for FA
2. Oblimin/quartimin cross-validation
3. CFA cross-validation against lavaan
4. Edge case tests: single factor, Heywood cases, perfect correlation
5. Consider adding more modules to Carm (SEM, mediation, IRT, etc.)

## Context
- Carm dir: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/JStats`
- Branch: `main`
- Build: `npx tsup`
- Type check: `npx tsc --noEmit`
- Tests: `npx vitest run`
- Validation: `npx tsx validation/ts-harness/fa-full-report.ts`
- Reports: `validation/FA-TECHNICAL-REPORT.md`, `validation/RANDOM-STARTS-REPORT.md`, `validation/GMM-TECHNICAL-REPORT.md`
