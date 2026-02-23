# Session Handoff — 2026-02-23

## Completed
- **analyze() dispatch layer** (`src/stats/analyze.ts`): High-level field-based API. Pass `NumericField`/`GroupField` objects → automatically selects the right test via Shapiro-Wilk normality routing.
- **New types** (`src/core/types.ts`): `FieldType`, `NumericField`, `GroupField`, `Field`, `AnalyzeOptions`, `AnalysisResult` appended at end of file.
- **Export updated** (`src/stats/index.ts`): `export * from './analyze.js'` added.
- **Tests** (`tests/stats/analyze.test.ts`): 28 tests covering all dispatch paths (t-test, MW, paired t, Wilcoxon, ANOVA+Tukey, KW+Dunn, chi-square/Fisher, descriptive-only), error throws, and `detectFieldType`.

## Current State
- **All 309/309 tests passing** across 11 test files.
- Build clean (`npm run build`).
- All stats modules implemented and cross-validated vs R.

| Test file | Tests |
|-----------|-------|
| core/matrix | 18 |
| core/math | 41 |
| stats/descriptive | 29 |
| stats/comparison | 25 |
| stats/correlation | 21 |
| stats/regression | 24 |
| stats/mixed | 20 |
| stats/pca | 14 |
| stats/analyze | 28 |
| large_sample | 47 |
| stress | 42 |

## Key Decisions
- **n < 3 treated as normal**: `checkNormality` returns W=1, p=1 for groups with fewer than 3 observations (SW undefined). n > 50 also skip SW (CLT).
- **exactOptionalPropertyTypes**: Return object uses `...(posthoc !== undefined && { posthoc })` spread to conditionally include optional fields.
- **forceTest bypasses normality routing entirely**: `selectTest` returns `forceTest` immediately if set.
- **tukeyHSD only for ANOVA**: Post-hoc is added only when dispatching to `one-way-anova` or `kruskal-wallis`; not for 2-group tests.
- **Fisher fallback**: Auto-applied only for 2×2 tables with expected count < 5 when `forceTest` is not set.

## Open Issues
- **LMM SEs**: Satterthwaite df approximation crude (`df = n−p−nGroups+1`). lme4 uses proper Satterthwaite. Test verifies finite SEs but not exact values.
- **Random slopes**: `randomSlopes` field in interface, not yet implemented.

## Next Steps
1. **Demo app** (`demo/`): Build interactive Vite demo showcasing all analyses, including `analyze()` as the primary entry point.
2. **tnadesktop integration**: `"carm": "github:mohsaqr/Carm"` in tnadesktop, one new tab.
3. `git commit` when ready (ask user first per CLAUDE.md).

## Context
- Working directory: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/JStats`
- Package name: `carm` (v0.1.0)
- Branch: `dev-clean`
- `npm run build` → `dist/`
- `NODE_OPTIONS='--max-old-space-size=4096' npx vitest run` → 309/309 (42s, LMM stress is slow)
