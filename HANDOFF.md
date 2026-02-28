# Session Handoff — 2026-02-28

## Completed This Session

### 1. EFA Rotation Method Presets & Advanced Controls (Aistatia UI)
Implemented the full preset system across 5 files in `aistatia/src/`:

| File | What changed |
|------|-------------|
| `state.ts` | Added `faRotationPreset` ('lavaan'), `faRandomStarts` (50), `faMaxIter` (1000), `faTol` (1e-6) to AppState interface + defaults |
| `views/wizard-defs.ts` | Added `FA_ROTATION_PRESET` (select: lavaan/psych/custom), `FA_RANDOM_STARTS` (number 0–500), `FA_MAX_ITER` (number 100–50000), `FA_TOL` (select: 1e-4 to 1e-7). Updated `WizardOptionSelect.id` type union. Inserted into EFA wizard options array. |
| `views/modal.ts` | Added fields to `ModalLocal` + `ModalConfig`. Init in `openModal()`. Value resolution in `renderOptionEl()` for selects + numbers. Change handlers for all 4 new fields. Added `isOptionVisible()` function + `GPF_ROTATIONS` set for conditional display. Filtering applied in `buildOptionsHtml()`. |
| `main.ts` | Added 4 new fields to `setState()` call in `handleOpenModal()` |
| `analysis/runner.ts` | Added `LAVAAN_PRESET` + `PSYCH_PRESET` constants. Resolve `rotOpts` from preset, spread into `runEFA()` call. |

**Preset values:**

| Parameter | lavaan | psych/GPArotation | Custom |
|-----------|--------|-------------------|--------|
| geominDelta | 0.001 | 0.01 | user-configurable |
| randomStarts | 30 | 0 | 0–500 |
| maxIter | 10000 | 1000 | 100–50000 |
| tol | 1e-5 | 1e-5 | 1e-4 to 1e-7 |

**Visibility rules:**
- Preset selector: shown only for GPF oblique rotations (geomin, oblimin, quartimin)
- Custom controls (delta, starts, iter, tol): shown only when preset = "Custom"
- Geomin delta: shown only when preset = "Custom" AND rotation = geomin

**No Carm changes** — Carm already accepts all these EFAOptions, only the UI wiring was needed.

### 2. Fixed Google Drive → Local Folder Migration Issues
Both repos moved from Google Drive to `/Users/mohammedsaqr/Documents/Github/`. Google Drive sync corrupted `node_modules`:
- Symlinks in `.bin/` became plain text files (21 bytes, just the relative path)
- Native binaries (`@rollup/rollup-darwin-arm64`, `@esbuild/darwin-arm64`) got `com.apple.quarantine` xattr from Chrome, causing macOS Gatekeeper `dlopen()` failures
- Execute permissions (`+x`) stripped from binaries

**Fix:** `rm -rf node_modules && npm install` in both repos.

### 3. Updated Documentation
- `HANDOFF.md`: Updated context with new repo paths and migration notes
- `LEARNINGS.md`: Added Google Drive `node_modules` corruption details and new repo paths

## Current State
- **JStats build:** `node build.mjs` — clean
- **JStats tests:** `npx vitest run` — **784/784 pass**
- **Aistatia tsc:** `npx tsc --noEmit` — **0 errors**
- **Aistatia build:** `npx vite build` — **clean** (601 modules, 1.32s)
- **Manual testing NOT YET DONE** — need to open the app and verify the EFA wizard preset UI works end-to-end

## Key Decisions
- Preset values are **constants in `runner.ts`** (not in state) — only "Custom" reads from state fields
- Visibility filtering lives in `modal.ts` via `isOptionVisible()`, not in `wizard-defs.ts`
- `onCommit` uses `...rest` spread which automatically passes new fields — no explicit field listing needed
- Default preset is "lavaan" (most common in academic use)

## Open Issues
- **Manual UI test needed:** Open Aistatia dev server, go to EFA wizard → pick geomin rotation → verify preset selector appears → switch lavaan/psych/custom → verify controls show/hide correctly
- **Geomin cross-validation not re-run** this session (200/200 from last session, no Carm changes)

## Next Steps
1. **Manual test the preset UI** — `cd aistatia && npx vite` → open browser → load data → EFA wizard → geomin → test all 3 presets
2. Consider adding preset info to the EFA results output (e.g., "Rotation: geomin (lavaan defaults)" in the summary)
3. Any remaining items from the broader project roadmap

## Context
- JStats (Carm) repo: `/Users/mohammedsaqr/Documents/Github/JStats`
- Aistatia: `/Users/mohammedsaqr/Documents/Github/aistatia/`
- Build JStats: `node build.mjs`
- Test JStats: `NODE_OPTIONS='--max-old-space-size=4096' npx vitest run` (784 tests)
- Build Aistatia: `cd ../aistatia && npx vite build`
- Dev Aistatia: `cd ../aistatia && npx vite`
- Validation: `npx tsx validation/ts-harness/fa-geomin-crossval.ts` (200/200)
- After moving from Google Drive: must `rm -rf node_modules && npm install` (Drive corrupts symlinks + quarantines native binaries)
