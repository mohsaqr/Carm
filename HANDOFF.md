# Session Handoff — 2026-02-23

## Completed
- **Shapiro-Wilk** (`src/stats/descriptive.ts`): Full AS R94 algorithm (Royston 1995). `shapiroWilk([1..10])` → W=0.9728, p=0.9177 matching R exactly.
- **LMM scale≈0 GLS bug** (`src/stats/mixed.ts`): Fixed crash/wrong-slope when ICC=0. Fix: detect `scale < 1e-10` → `VinvScaled = I` (OLS).
- **LMM multi-start optimizer** (`src/stats/mixed.ts`): 5 starts [-4,-2,0,2,4], keep best.
- **LMM logLik/AIC constant** (`src/stats/mixed.ts`): Added `−½(n−p)(1+log(2π))` so logLik matches R lme4 (diff was 422.8 for n=300,p=2).
- **Large-sample tests** (`tests/large_sample.test.ts`): 47 tests, R cross-validated, n=50–300.
- **100-dataset stress tests** (`tests/stress.test.ts`): 42 tests, seeded LCG + Box-Muller, n=30–114, all mathematical invariants and statistical bias checks.
- **LMM verified correct** (`tests/stats/mixed.test.ts`): R lme4 re-run on the 2-group fixture confirmed ICC=0.977 (not 0.4 as previously wrongly stated). Exact-value tests added: ICC>0.95, logLik≈-12.399, AIC≈32.797, fixef=(2.1,1.2). Wrong ground truth corrected in all docs.

## Current State
- **All 278/278 tests passing** across 10 test files.
- Build clean (`npm run build`).
- All stats modules implemented and cross-validated vs R.
- All viz modules implemented (12 plot types).

| Test file | Tests |
|-----------|-------|
| core/matrix | 18 |
| core/math | 38 |
| stats/descriptive | 29 |
| stats/comparison | 25 |
| stats/correlation | 21 |
| stats/regression | 24 |
| stats/mixed | 20 |
| stats/pca | 14 |
| large_sample | 47 |
| stress | 42 |

## Key Decisions
- **roundTo(x,6) rounding**: Derived invariants (iqr=q3-q1, variance=sd², se=sd/√n) tested at `toBeCloseTo(x,5)` not 10, because values are independently rounded.
- **Pearson CI at r≈1**: Fisher-z CI upper bound can be 0.9999830 while stored r=1.0. Test uses `+1e-4` tolerance on upper bound; not a real bug.
- **Profiled REML (1D)**: Parameterize by log(ψ); σ²_e profiled out analytically. Scale≈0 directly sets `VinvScaled=I` (GLS=OLS).

## Open Issues
- **LMM SEs**: Satterthwaite df approximation crude (`df = n−p−nGroups+1`). lme4 uses proper Satterthwaite. Test verifies finite SEs but not exact values.
- **Random slopes**: `randomSlopes` field in interface, not yet implemented.

## Next Steps
1. **Demo app** (`demo/`): Build interactive Vite demo showcasing all analyses.
2. **tnadesktop integration**: `"carm": "github:mohsaqr/Carm"` in tnadesktop, one new tab.
3. `git commit` when ready (ask user first per CLAUDE.md).

## Context
- Working directory: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/Carm`
- Branch: `dev-clean`
- `npm run build` → `dist/`, `NODE_OPTIONS='--max-old-space-size=4096' npx vitest run` → 278/278 (45s, LMM stress is slow)
- `node tmp/lmm-verify.mjs` → 5 deterministic datasets all correct
