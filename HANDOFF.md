# Session Handoff — 2026-02-25

## Completed
- **DBSCAN clustering** (Phase 3): Full implementation in Carm `clustering.ts` — eps-neighborhood BFS, core/border/noise classification, silhouette scores (excluding noise), k-distance plot. 20 tests cross-validated against R `dbscan::dbscan()`.
- **Hierarchical Agglomerative Clustering** (Phase 4): Lance-Williams recurrence for 4 linkages (single/complete/average/ward). `cutTree()` and `cutTreeHeight()` via union-find. Cophenetic correlation. 29 tests cross-validated against R `hclust()` — heights match to 1e-10.
- **Preprocessing module** (Phase 1): `preprocessData()` with center/standardize/log/sqrt methods. `inverseTransform()` for back-transformation. Refactored PCA to use it. 15 tests cross-validated against R `scale()`.
- **Shared utilities** (Phase 2): `euclideanDistMatrix()` (N×N Float64Array), `silhouetteScores()` (excluding noise labels).
- **Aistatia integration** (Phase 5): Full UI flow for DBSCAN and Hierarchical — wizards, runner dispatch, profile tables, APA banners, copy-for-paper, narrative, plots. Also added preprocess option to existing GMM/KMeans wizards. 12 files modified, `tsc --noEmit` clean.
- **HAC test fix**: Corrected R reference heights that were wrong for 19 of 29 values. Extracted correct values from `r_hac_reference.json`.

## Current State
- **Carm**: 496/496 tests passing across 17 files. Branch `dev-clean`. Build clean.
- **Aistatia**: `tsc --noEmit` clean, 0 errors. DBSCAN and Hierarchical fully wired.
- **Not committed**: All changes are local in both repos.

## Key Decisions
- **ClusterPlotData adapter**: Unified plot rendering for 4 clustering methods (GMM, KMeans, DBSCAN, Hierarchical) via a single adapter interface instead of duplicating code.
- **Preprocessing on data, stats on original**: Preprocessing transforms the data matrix for clustering, but per-cluster means/SDs are computed on original data for interpretability.
- **Ward.D2**: Squared Euclidean internally, heights = sqrt(merge distance). Matches R convention.
- **DBSCAN noise = -1**: Following the convention of -1 for noise, 0-indexed clusters. R uses 0 for noise, 1-indexed.
- **Dendrogram**: Added as a plot type placeholder in wizard-defs. Actual D3 rendering not yet implemented.

## Open Issues
- **Dendrogram D3 visualization**: Placeholder exists in wizard-defs and plot-panel but no actual rendering code. Needs a D3 dendrogram renderer.
- **LCA/LTA UI**: Carm module ready, but aistatia UI needs binary variable detection and time-series input.
- **Cluster label stability**: GMM/KMeans label ordering depends on initialization.
- **Not yet committed**: Changes in both repos are local only.

## Next Steps
1. Implement D3 dendrogram renderer for hierarchical clustering visualization.
2. Commit both repos (pending user approval).
3. Test in aistatia: load data → DBSCAN → verify noise + clusters in browser.
4. Test in aistatia: load data → Hierarchical (Ward) → verify cluster assignment in browser.
5. Consider silhouette visualization (bar chart of per-point silhouette scores).
6. LCA/LTA UI integration.

## Context
- Carm dir: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/JStats`
- Aistatia dir: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/aistatia`
- Branch: `dev-clean`
- Build: `npm run build` (Carm) → copy dist to aistatia `node_modules/carm/dist/`
- Tests: `npx vitest run` → 496/496 (Carm)
- Plan file: `/Users/mohammedsaqr/.claude/plans/cosmic-nibbling-boole.md`
