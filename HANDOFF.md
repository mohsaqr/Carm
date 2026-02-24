# Session Handoff — 2026-02-24

## Completed
- **Configurable plot annotations**: All 8 plot types support conditional rendering of n-labels, median, mean diamond, jitter, brackets, outliers, density overlay, normal curve, CI band, equation, rug, legend, counts, percentages, values, significance.
- **Numeric p-values on brackets**: `numericP` flag — shows `p = .025` instead of `***` stars. Default true.
- **Mean diamond marker + median value labels**: White-filled diamond at group mean, numeric label next to median line and mean diamond.
- **Redesigned plot toolbar**: Modern pill toggles, font family/size selects, subtitle toggle, SVG/PNG export.
- **Floating gear settings panel**: ⚙ icon opens popover with sliders (jitter width, bandwidth, point size, opacity, bins, plot height) and bold/halo text styling toggles. Per-plot-type control registry.
- **Editable plot title**: Click SVG title → inline text input → commit on blur.
- **Git repos set up**: Carm pushed to `mohsaqr/Carm` (dev-clean). Aistatia initialized and pushed to `mohsaqr/aistatia` (private, main).

## Current State
- **JStats (Carm)**: All 309/309 tests passing. Build clean. Pushed to `origin/dev-clean`.
- **Aistatia**: `tsc --noEmit` clean. `npm run build` clean. Pushed to `origin/main`.
- **Carm dist in aistatia**: node_modules/carm/dist has the latest bracket/annotation changes.

## Key Decisions
- **Two-tier toolbar**: Pill bar (quick booleans) + gear popover (sliders + text styling). Avoids overcrowding.
- **defaultValue=0 = "auto"**: Bandwidth, bins, plotHeight use 0 to mean "let renderer decide".
- **pointOpacity flows through CarmTheme**: Overridden in `buildTheme()`, affects all point renders.
- **Bold/halo applied post-render**: `font-weight: 600` on SVG text, `paint-order: stroke` for halo.

## Open Issues
- **p-values may be hidden by stale localStorage**: If user toggled `showBrackets` off in a previous session, it persists. Clear with: `Object.keys(localStorage).filter(k => k.startsWith('aistatia-')).forEach(k => localStorage.removeItem(k))`
- **LMM SEs**: Satterthwaite df approximation still crude.
- **Random slopes**: Not yet implemented.

## Next Steps
1. Investigate p-values disappearing issue — likely stale localStorage.
2. Visual verification of gear panel with real data.
3. Consider adding "Reset defaults" button to gear panel.

## Context
- JStats dir: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/JStats`
- Aistatia dir: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/aistatia`
- Carm package: `carm` (v0.1.0), branch `dev-clean`
- Aistatia: branch `main`, private repo `mohsaqr/aistatia`
- `npm run build` (JStats) → `dist/` via tsup
- `npm run build` (aistatia) → `dist/` via vite
- Tests: `NODE_OPTIONS='--max-old-space-size=4096' npx vitest run` → 309/309
