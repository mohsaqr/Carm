# Session Handoff — 2026-02-27

## Completed (This Session)
- **Numerical equivalence tests** for all 14 new statistical methods against R reference values
- **R reference script** (`tests/fixtures/new-methods-ref.R`) generates fixture JSON with 24 test cases
- **R fixture** (`tests/fixtures/new-methods-ref.json`) contains full-precision reference values
- **Vitest file** (`tests/stats/numerical-equivalence.test.ts`) — 24 tests, all passing
- **5 bug fixes** found during equivalence testing:
  1. Binomial test boundary crash when k=n (logGamma(0) error)
  2. Binomial two-sided p-value algorithm (exact enumeration vs simple doubling)
  3. McNemar's Yates correction (|b-c|-1, not |b-c|-0.5)
  4. Two-way ANOVA Type III SS (treatment coding, not effects coding, to match car::Anova)
  5. Factor level sorting (alphabetical to match R's default)

## Current State
- **784/784 tests passing** across 23 suites
- **Build clean** (`tsc --noEmit` passes)
- All 14 new methods are R-equivalent with appropriate tolerances
- PRNG shared from `core/prng.ts` — clustering.ts and factor-analysis.ts both import from there

## Key Decisions
- **Binomial two-sided p-value**: Implemented R's exact algorithm (enumerate all outcomes, sum those with P ≤ P(observed)) instead of simple 2*min approach. More complex but matches R exactly.
- **McNemar correction**: Changed to |b-c|-1 to match R. Some textbooks use 0.5 but R uses 1.
- **Type III SS**: Uses treatment (dummy) coding same as Type II, but with Type III model comparison strategy (keeping interaction when testing main effects). This matches R's car::Anova(type=3) default behavior.
- **Factor levels sorted alphabetically**: Matches R's `factor()` default. Previously used insertion order which gave different dummy coding.

## Open Issues
- BCa bootstrap acceleration factor uses jackknife which is O(n²) — could be slow for large n
- Ordinal logistic and negative binomial convergence tolerance may need tuning for edge cases
- Two-way ANOVA only supports balanced designs well; unbalanced designs not extensively tested

## Next Steps
1. Consider adding formal cross-validation harness (like existing validation/ infrastructure) for the 14 new methods
2. Add visualization renderers for new methods (two-way interaction plot, forest plot for ordinal)
3. Consider multinomial logistic regression (natural extension of ordinal)
4. Consider robust regression (M-estimation, Huber weights)

## Context
- Git repo: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/JStats`
- Branch: `dev-clean`
- Node.js v25, TypeScript strict mode, Vitest 4.0.18
- R 4.5.1 for cross-validation
