## 2026-02-26 (Factor Analysis — Geomin rotation + GPFoblq)

### Random starts for GPFoblq (lavaan GPA×30 strategy)
- lavaan uses 30 random orthogonal starting matrices (Haar-distributed via QR of N(0,1) matrix with sign correction), picks the solution with lowest criterion value. GPArotation uses single start from T=I.
- Carm now supports `randomStarts` option (default: 50). With 50 random starts, it matches lavaan's GPA×30 behavior on all tested datasets. Start 0 is always T=I (deterministic).
- Random orthogonal matrix generation: fill k×k with N(0,1), Modified Gram-Schmidt QR, then sign-correct Q columns by sign(R diagonal). Matches both `GPArotation::Random.Start()` and `lavaan::lav_matrix_rotate_gen()`.
- When Tinit ≠ I, initial L = A × inv(T)' (not just A). The gradient formula is the same general form.
- Start 0 is always T=I (deterministic, reproducible). Random starts use the PRNG seeded from options.seed.

### Teacher burnout cross-validation: randomStarts=30 achieves exact lavaan match
- eps=0.001 (lavaan default): single start MAE=0.089, maxErr=0.446, Tucker's φ=0.903 → randomStarts=30 MAE=0.000001, maxErr=0.000005, φ=1.000000. 100% MAE reduction, all 115 cells within 0.05.
- eps=0.01: single start already matches lavaan perfectly (MAE=0.000003). Criterion surface has a single basin — no random starts needed.
- Performance: single start ~500ms, 30 starts ~1100ms (2× slower, still fast for 438×23 dataset).
- The local optimum issue with eps=0.001 is concentrated on EE and DE items (EE MAE=0.199, DE MAE=0.097 with single start).
- `lavInspect(fit, "std")$lambda` gives standardized loadings matching the summary output. `lavInspect(fit, "est")$lambda` gives unstandardized — different values.

### lavaan vs GPArotation geomin: different optimization strategies
- lavaan uses GPA(30) — 30 random rotation starts, picks global optimum. GPArotation (and Carm with randomStarts=1) use single start from T=I. With randomStarts=30, Carm now matches lavaan exactly (MAE < 1e-5 on teacher burnout 438×23). lavaan's epsilon=0.001 defaults also differ from GPArotation's delta=0.01.

### geomin delta=0.001 convergence
- With delta=0.001 (lavaan default), GPFoblq may not converge in 1000 iterations on complex datasets (23 vars, 5 factors). delta=0.01 converges reliably. Both produce valid solutions.

### GPFoblq gradient uses L (rotated), not A (unrotated)
- R's `GPArotation::GPFoblq` computes gradient as `G = -t(t(L) %*% VgQ$Gq %*% solve(Tmat))` = `-(L' × Gq × T^{-1})'`.
- Note: uses **L** (current rotated loadings) and **T^{-1}**, NOT A (original unrotated) and inv(T').
- The two formulas are only equivalent when T = I (first iteration). After step 1, using A gives the wrong gradient direction.
- This caused 70/100 geomin cross-validation failures (all k≥3 cases) despite our unrotated loadings matching R 100/100.
- The initial gradient formula (using A and solve(t(Tmat))) in R is a shortcut valid only at T=I. The loop formula (using L and solve(Tmat)) is the correct general form.

### GPFoblq always takes a step after line search
- R's GPFoblq updates `Tmat <- Tmatt; f <- VgQt$f; G <- ...` AFTER the inner line search loop, regardless of whether the Armijo condition was satisfied.
- Our initial implementation only updated inside the `if (improvement > ...)` block, causing the algorithm to stall when no step size satisfied Armijo.
- Fix: declare Tmatt/Lnew/vgQnew outside the inner loop, always update after it.

### Canonical sign convention for eigenvector determinism
- After ML or PAF extraction, for each factor column, ensure the element with the largest absolute value is positive.
- Matches LAPACK's dsyev convention used by R. Essential for non-rotation-invariant methods (geomin, oblimin) where the starting orientation determines which local optimum is found.
- Not needed for varimax/promax (rotation-invariant), but harmless to apply universally.

## 2026-02-26 (Factor Analysis — Promax cross-validation fixes)

### psych::fa promax rotation uses outer Kaiser normalization
- `psych::fa(rotate="promax")` does NOT call `stats::promax()` or `psych::Promax()` directly. It dispatches via `psych::kaiser(loadings, rotate="Promax")`, which adds outer Kaiser normalization:
  1. Compute communalities h² = rowSums(L²)
  2. Normalize: L_norm = L / sqrt(h²) (rows → unit norm)
  3. Call Promax(L_norm) (which internally calls varimax on the already-normalized loadings)
  4. Denormalize: result = rotated × sqrt(h²)
- This changes the promax target Q nonlinearly because Q = V⊙|V|^(m-1) is computed on unit-norm rows instead of original-scale rows.
- `psych::Promax(L)` ≡ `stats::promax(L)` (MAE=0), but `psych::fa(promax)` ≠ both due to the kaiser wrapper.
- Source: `psych:::fac()` line 555-557: `promax = { pro <- kaiser(loadings, rotate="Promax", ...) }`

### SVD-based varimax using eigendecomposition of B'B
- R's `stats::varimax()` uses SVD of the B matrix for polar decomposition. The Carm `Matrix.svd()` is unreliable for small k×k matrices. Fix: compute polar decomposition via eigendecomposition of B'B → get V and σ², then T = B·V·diag(1/σ)·V'.

### ML extraction hybrid approach
- Pure R-style concentrated gradient descent (FAfn/FAgr) doesn't converge as well as R's L-BFGS-B. Solution: Phase 1 uses Jöreskog gradient descent (reliable convergence), Phase 2 polishes with Nelder-Mead on R's exact concentrated ML objective. This achieves matching uniquenesses (max Δ ~5e-5) across 100 synthetic datasets.

### R's factanal initialization formula
- R initializes uniquenesses as `(1 - 0.5*k/p) / diag(R^{-1})`, not `1 - 1/diag(R^{-1})` (which computes communality). Both converge to the same optimum in practice.

### aistatia npm install gotcha
- When aistatia installs carm via `npm install ../JStats`, it copies the dist/ files (not symlinks). After rebuilding carm, must run `npm install` again in aistatia to pick up changes. Stale node_modules/carm was the cause of mysterious "debug output not appearing" during development.

## 2026-02-25 (Factor Analysis — EFA + CFA)

### Matrix immutability in Carm
- Carm's `Matrix` class has NO `set()` method — it's immutable. The FA.ts reference code uses `L.set(i, j, v)` everywhere, which won't compile. Pattern: use mutable `number[][]` arrays for iterative algorithms (PAF, ML, CFA gradient descent), then `Matrix.fromArray()` at boundaries when you need matrix operations (multiply, inverse, eigen, logDet).

### CFA standard errors via numerical Hessian
- The Fisher information approach requires the Hessian of the Wishart ML discrepancy w.r.t. all free parameters. Analytic Hessian is complex (involves ∂Σ⁻¹/∂θ terms). Numerical central finite differences (h=1e-4) work well and give SE estimates that are reasonable for moderate n. The Hessian is symmetric so only upper triangle needs computation.

### RMSEA 90% CI
- Exact RMSEA CI requires non-central chi-square quantile inversion, which is computationally expensive. The Steiger (1990) approximation using ±1.645×√(2df) on the NCP works well for moderate df. More exact methods (e.g., bisection on non-central chi-square CDF) can be added later if needed.

### Parallel analysis eigenvalue comparison
- Compare observed eigenvalues to the 95th percentile of simulated eigenvalues (not the mean). Using the mean can over-retain factors. Store all simulated eigenvalues, sort per position, and take the 95th percentile index.

### Float64Array and noUncheckedIndexedAccess
- When using `Float64Array` under strict `noUncheckedIndexedAccess`, the `+=` operator fails because `arr[i]` returns `number | undefined`. Use `arr[i] = arr[i]! + value` instead of `arr[i] += value`.

## 2026-02-25 (Dendrogram visualization)

### Dendrogram coordinate computation
- d3.hierarchy/d3.cluster is wrong for dendrograms — it assumes integer depth levels and uniform leaf spacing. Real dendrograms need continuous y-axis (merge height) and DFS leaf ordering from the HAC. Manual computation: leafX from dendrogramOrder index, internal nodeY from merge heights, internal nodeX = midpoint of children.
- U-shaped elbow path: `M ax,ay V my H bx V by` — from child A up to merge height, across to child B, down.

### Subtree coloring propagation
- Leaves get their cluster label directly. Internal nodes: if both children share the same cluster, inherit it; otherwise mark as -1 (mixed). Mixed subtrees get the neutral `theme.axisLine` color. This naturally colors below-cut subtrees by cluster and above-cut merges in gray.

### Cut line placement
- For K clusters from n observations: the cut goes between merge index `n-k-1` (last merge kept) and `n-k` (first merge cut). The midpoint of their heights gives a clean visual separation.

## 2026-02-25 (DBSCAN + HAC + preprocessing)

### HAC R reference heights
- When hard-coding R reference values from JSON, **always extract programmatically** — never type-transcribe. In this session, the first 10 ward heights were correct but the remaining 19 were wrong (likely copy-paste from a different dataset or prior run). Debug was painful because the values were plausible but slightly off.
- HAC is fully deterministic — heights must match R to 1e-10. Any deviation means the reference data is wrong, not the algorithm.
- Ward.D2 last 3 heights for the 30-observation test data: 2.9423, 21.3337, 38.8347. The wrong values (5.86, 6.14, 8.20) were dramatically different — always sanity-check scale.

### DBSCAN label conventions
- R dbscan: 0 = noise, 1-indexed clusters. Carm: -1 = noise, 0-indexed clusters. Conversion: `rLabel === 0 ? -1 : rLabel - 1`.
- Silhouette scores exclude noise points (label === -1). The mean silhouette is computed only over non-noise points.

### ClusterPlotData adapter pattern
- When 4+ clustering result kinds (GMM, KMeans, DBSCAN, Hierarchical) need the same plot rendering (profile bar, scatter), use an adapter interface (`ClusterPlotData`) with a `toClusterPlotData(runResult)` converter. This avoids duplicating plot code 4x and handles the discriminated union narrowing in one place.

### Preprocessing integration
- Added `preprocess?` option to all clustering wizards (GMM, KMeans, DBSCAN, Hierarchical). Preprocessing runs on the data matrix before the clustering algorithm, but per-cluster means/SDs are computed on the **original** (unpreprocessed) data for interpretability.
- `preprocessData` returns `{ data, colMeans, colSDs, method, centered, scaled }` — the `data` field is the transformed matrix.

### Lance-Williams recurrence for HAC
- Ward coefficients use dynamic cluster sizes: `α_i = (n_i + n_k)/N_t`, `α_j = (n_j + n_k)/N_t`, `β = -n_k/N_t`, `γ = 0` where `N_t = n_i + n_j + n_k`.
- Ward.D2: squared Euclidean internally, final heights = `sqrt(merge_distance)`. This matches R `hclust(method="ward.D2")`.

### TypeScript narrowing with discriminated unions
- When a function returns early for multiple `kind` values, TypeScript doesn't always narrow the remaining type. Adding explicit `if (runResult.kind === 'clustering') return` guards before accessing `.result` on a standard result kind is needed to satisfy the compiler.

## 2026-02-25 (cluster finder wizard)

### Normalized entropy (mclust convention)
- The `diagnostics.entropy` field must be **normalized**: `1 - E / (N * log(K))`, range [0,1], higher = better separation. This matches mclust's `1 + sum(probs * log(probs)) / (n * log(K))`.
- The **raw** entropy `E = -sum(z * log(z))` is only used internally for ICL = BIC + 2E. Never expose it directly.
- With well-separated data, posteriors are near one-hot → rawE ≈ 0 → normalizedE ≈ 1.0.
- Vite pre-bundling cache (`node_modules/.vite/`) must be cleared after rebuilding Carm, otherwise the dev server serves stale code.

### Elbow heuristic for KMeans
- The elbow method uses the maximum second derivative of inertia: `d2 = inertia[i-1] - 2*inertia[i] + inertia[i+1]`. Requires at least 3 entries in the range to compute. Falls back to first K if fewer than 3.

### clusterMinK/clusterMaxK state routing
- Same pattern as `clusterK` and `clusterModel`: top-level state fields routed through ModalLocal → ModalConfig → setState, not through plotConfig. The modal number input handler needs a special-case lookup for `curVal` (reading from `local.clusterMinK` / `local.clusterMaxK` instead of `local.plotConfig[opt]`).

### fitGMMRange/fitKMeansRange error handling
- These functions silently skip failed fits (e.g. singular covariance for large K relative to data). Only throws if ALL fits fail. This is important because high-K fits can be unstable with small datasets.

## 2026-02-25 (clustering UI in aistatia)

### Wizard option routing for non-plotConfig fields
- `clusterK` and `clusterModel` are top-level state fields, not PlotConfig fields. The modal's `data-opt` change handler must special-case them (like `ciLevel` and `pAdjMethod`) rather than routing through `local.plotConfig[opt]`.
- WizardOptionSelect `id` union must include the new option IDs (added `'clusterModel'`).
- WizardOptionNumber with `id: 'clusterK'` flows through the standard number-input codepath but with a special-case read from `local.clusterK` instead of `local.plotConfig.clusterK`.

### ClusterDiagnostics type
- `ClusterDiagnostics` does not have a `model` field — the covariance model name is not stored in the diagnostics struct. It's embedded in the `formatted` string instead.

### Scatter plot cluster coloring
- Carm's `renderScatterStats` renders points as `<circle>` elements in data order. Post-render coloring by cluster label works by iterating circles and matching index to `labels[i]`.
- Must use requestAnimationFrame polling because D3 renders asynchronously via `import('d3').then(...)`.

## 2026-02-25 (clustering module)

### GMM cross-validation with R mclust
- K-Means++ initialization (splitmix32 PRNG) finds slightly different local optima than mclust's hierarchical init. Multi-seed search (9+ seeds) closes the gap to ΔLL ≈ 0.004 for VVV. EII and VVI reach near-exact match (ΔLL ≈ 0.0001 and 0.0000).
- mclust returns *negative* BIC (higher = better). Convert to standard BIC: `bic = -fit$bic`.
- mclust's `attr(fit$bic, "df")` stores the DF — useful for cross-validation but not always well-documented.

### LCA smoothing: MLE vs Beta(1,1)
- poLCA uses raw MLE for rho (item-response probabilities): `sumX / sumW`.
- Beta(1,1) smoothing `(sumX + 1) / (sumW + 2)` biases rho toward 0.5. For a 2-class/5-item dataset (n=100), this causes a ΔLL gap of ~0.478 — unacceptable for R equivalence.
- Fix: use raw MLE with minimal floor to prevent log(0): `Math.max(Math.min(sumX / sumW, 1 - 1e-10), 1e-10)`. After fix: ΔLL = 0.0000000076 (exact to 10 decimals).

### poLCA data format
- poLCA requires data coded as 1/2 (not 0/1). The R reference script does `lca_data_raw + 1`. Our TS implementation works with 0/1 directly.
- `probs[[d]][k, 2]` gives P(value=2) = P(original=1), which is what we store as rho.

### K-Means exact match with R
- When both R and TS converge to the same optimum (which K-Means++ usually does for well-separated clusters), inertia matches to 10+ decimal places.
- R's `kmeans()` with `algorithm = "Lloyd"` matches our Lloyd's implementation exactly.

### LTA cross-validation
- depmixS4 is non-deterministic (different runs give different results even with set.seed). Cannot use for cross-validation.
- Alternatives: seqHMM package, hmmlearn (Python), or manual hand-computation for small cases.

### TypeScript noUncheckedIndexedAccess with +=
- With `noUncheckedIndexedAccess: true`, `arr[i] += x` fails because `arr[i]` is `T | undefined`. Must use `arr[i]! += x`.
- This applies to all compound assignment operators: `+=`, `-=`, `*=`, `++`.

## 2026-02-24 (floating settings panel + bold/halo)

### Gear popover architecture
- Two-tier toolbar: pill bar (quick boolean toggles + selects) stays in `.plot-bar`, gear popover (sliders + text styling) floats absolutely positioned from the bar.
- `position: relative` on `.plot-bar` is set dynamically in JS when the popover is attached, so it only applies to bars that have gear controls.
- Slider controls with `defaultValue: 0` use "auto" display for bandwidth/bins/plotHeight — 0 means "let the renderer decide". Non-zero values pass through to the renderer config.
- `pointOpacity` flows through `CarmTheme.pointOpacity` (overridden in `buildTheme()`). All other slider values pass directly to renderer config objects via conditional spread.
- Bold labels: `font-weight: 600` on all `svg text` except `.plot-title` (which is already bold). Applied post-render.
- Text halo: SVG `paint-order: stroke` + white stroke + `stroke-linejoin: round`. Creates a clean outline behind text for readability over data.
- Close-on-outside-click uses `document.addEventListener('click')` — must check `!popover.contains(target)` and `target !== gearBtn` to avoid immediate close on open.

## 2026-02-24 (configurable annotations + toolbar redesign)

### Annotation toggle pattern
- All new config booleans use `!== false` pattern (e.g. `if (config.showN !== false)`), which means `undefined` defaults to `true`. Only `showMean` uses `if (config.showMean)` to default to `false`.
- Mean diamond marker uses SVG `<polygon>` with 4 points (top, right, bottom, left) rather than a rotated `<rect>` — simpler SVG and no transform needed.
- When copying JStats dist/ to aistatia's node_modules/carm/dist/, TypeScript picks up the new types immediately without needing `npm install`. This is because `carm` is referenced via GitHub branch, not a versioned registry.
- For categorical re-renders (bar/pie/lollipop), the context map stores `catVarName` and `catFreqTable` separately because these plots go through `renderCategoricalPlot` instead of the standard `renderPlot` dispatcher.

### Toolbar design: pills > checkboxes
- Pill-shaped buttons with `.on` CSS class toggle look much better than native checkbox inputs for boolean toggles. Use `border-radius: 10px`, subtle border, blue background when on.
- For font control: construct a `CarmTheme` override from the user's selection and pass it as the `theme` field in each renderer's config. The renderer uses `theme.fontFamily` in D3 `.attr()` calls, so the theme must be set BEFORE rendering (CSS custom properties from `applyTheme` won't override inline SVG attributes).
- To show/hide the statistics subtitle after render, add a CSS class (`.plot-subtitle`) to the subtitle text element in Carm's `addSubtitle()`, then toggle `style.display = 'none'` post-render. This avoids modifying every renderer's internal logic.
- Merging toggle bar + export buttons into one unified toolbar reduces DOM complexity and looks cleaner than separate bars.

### Numeric p-values vs stars
- `formatBracketP(p, numeric)` — when `numeric=true`, shows `p = .025` or `p < .001`. When false, shows `***`/`**`/`*`/`ns`. Default to numeric since researchers need exact values.
- For `exactOptionalPropertyTypes: true`, passing `config.numericP` (type `boolean | undefined`) directly to a `boolean` field errors. Use conditional spread: `...(config.numericP !== undefined && { numericP: config.numericP })`.

### Editable title/subtitle
- Custom title/subtitle stored as `customTitle`/`customSubtitle` in PlotConfig. Applied post-render by finding `.plot-title`/`.plot-subtitle` SVG text elements and setting `textContent`.
- Live text update on `input` event (fast, just replaces text content in DOM). Full re-render on `blur` (for title that affects layout).
- `escapeAttr()` needed when injecting user text into HTML attribute values (title/subtitle default values in input fields).

## 2026-02-23 (analyze dispatch)

### analyze() auto-routing via Shapiro-Wilk
- Groups with n < 3 are treated as normal (SW undefined for n<3). This means: if you want deterministic routing to t-test in tests, use n=2 per group.
- Groups with n > 50 also skip SW (CLT applies — SW is overly sensitive at large n) → treated as normal.
- `exactOptionalPropertyTypes: true` in tsconfig: must use `...(x !== undefined && { x })` spread pattern to conditionally include optional object properties; cannot assign `undefined` to an optional field.
- `tukeyHSD(groups, msWithin, dfWithin, ciLevel)` requires ANOVA output — run `oneWayANOVA` first, then extract `.msWithin` and `.dfWithin` from the `ANOVAResult`.
- `fisherExactTest(a, b, c, d)` takes 4 cell values (not a matrix). Extract from `contingencyTable().table[i][j]`.
- Auto-fallback chi-square → Fisher only applies for 2×2 tables with any expected count < 5.

## 2026-02-23

### Shapiro-Wilk AS R94 Algorithm
- Our simplified halfNormalQuantiles + phi formula for `an` was wrong: phi formula gives `an ≈ 0.868` for n=10, but AS R94 polynomial correction gives `an ≈ 0.573` (matches R). The difference is ~50% and cascades into a completely wrong p-value (0.09 vs 0.92).
- The full AS R94 algorithm (Royston 1995) stores only the upper-half coefficients `a[0..nn2-1]` (all positive), then computes W as `(Σ a[i] * (sorted[n-1-i] - sorted[i]))² / SST`. No full antisymmetric array needed.
- Key step: save `a0orig` and `a1orig` BEFORE overwriting with polynomial corrections — needed to compute `fac` for scaling remaining inner coefficients.
- Coefficient arrays SW_C1 and SW_C2 are ascending-degree (constant term first). Use `swPoly(c, x)` with Horner's method: iterate from last index downward.
- n=3,4,5: use exact tabulated coefficients (not polynomial). n=6: correct only outermost (no a[2] correction). n>5: correct two outer coefficients.
- p-value function (`shapiroWilkPValue`) was already correct; only the W computation needed fixing.

### LMM AIC/logLik constant (found in large-sample testing)
- Our profiled REML minimized `−½{(n−p)log(σ²_e) + log|V_ψ| + log|X'V_ψ⁻¹X|}`, omitting `−½(n−p)(1+log(2π))`.
- lme4 includes this constant in reported logLik/AIC/BIC.
- Fix: `logLik = profiled_reml − 0.5*(n-p)*(1+log(2π))`.
- For n=300, p=2: constant = 0.5*298*(1+log(2π)) = 422.8. Confirmed: -245.43 - 422.84 = -668.27 = R's logLik. AIC then matches R to <1 unit.
- ANOVA effect size: `oneWayANOVA` returns ω² (omega-squared), NOT η² (eta-squared). ω² = (SS_B − df_B·MS_W)/(SS_T + MS_W). R's `eta²` is larger (SS_B/SS_T). Both correct — different estimands.

### LMM / REML
- **R lme4 match confirmed (2026-02-23)**: Our implementation matches R lme4 exactly on the 2-group fixture y=[1..5,6..9,12]. Verified: fixef=(2.1,1.2), σ²_b=14.511, σ²_e=0.343, ICC=0.977, logLik=-12.399, AIC=32.797. The earlier "R ground truth" of ICC=0.4092 was wrong/fabricated — ICC=0.977 is structurally correct (SS_between=72.9, SS_within=2.4, naive ICC≈0.968).
- **scale≈0 bug (FIXED 2026-02-23)**: When `sigmab2 ≈ 0` (ICC ≈ 0), `scale = sigmab2/sigmae2 ≈ 0`, so `1/scale = Infinity`. Fix: check `scale < 1e-10` and directly set `VinvScaled = Matrix.identity(n)` → GLS correctly equals OLS.
- **Multi-start optimizer (2026-02-23)**: Replaced single-start with 5 starts [-4,-2,0,2,4]. Prevents getting stuck in local optima.
- Verified on 5 deterministic datasets: ICC=0 correctly falls back to OLS; ICC>0.9/0.3/0.8 cases correctly estimated.
- Large-sample benchmark (n=300, 10 groups): logLik matches R lme4 to <1 unit after adding the `−½(n-p)(1+log(2π))` constant.

### Nelder-Mead sync bug (fixed in previous session)
- After sorting `fvals`, the `simplex` array must be re-indexed to match. Failure to do this breaks vertex-function mapping after iteration 1 and prevents convergence.
- Fix: after sort, immediately reassign `simplex[i] = sorted[i]` for all i.

### SVD truncation (fixed in previous session)
- For m×n matrix, SVD must return k=min(m,n) singular values. Returning n values instead of min(m,n) breaks rank-deficient or non-square cases.

### normalCDF / normalQuantile precision
- A&S 7.1.26 polynomial has ~1e-9 error at z=0. `normalCDF(0)` returns 0.5 + 5e-10, not exactly 0.5. Test tolerance should be precision=6, not 10.
- BSM normalQuantile (A&S §26.2.16) has max error ~4.5e-4. Round-trip `normalCDF(normalQuantile(p))` accurate to ~3 decimal places, not 4.

### Kendall tau p-value
- Our implementation uses normal approximation; R uses exact distribution for small n. For n=5, normal approx gives p≈0.07, R exact gives p≈0.10. Both indicate non-significance. Accept normal approx.

### Pearson CI at r≈±1 (found in 100-dataset stress test)
- When x is exactly 1,2,…,n (integer sequence) and the noise is very small, r can round to exactly 1.0 (stored as `roundTo(r, 6) = 1.0`).
- The Fisher-z CI upper bound is `tanh(atanh(r) + z/√(n−3))`. At r=1.0, `atanh(1)=Infinity`, so the CI is computed at the pre-rounded r≈0.9999830 instead, giving ci[1]≈0.9999830 < 1.0.
- Net effect: stored r=1.0 but ci[1]=0.9999830 → naive `r ≤ ci[1]` fails by ~1.7e-5.
- Fix in tests: use tolerance `1e-4` in the upper-bound check. The discrepancy is a rounding artefact, not a real logical error. Production code does not need changing.
- This affects ~1 of 100 stress datasets (dataset 19), only when the regression has near-perfect fit.

### stress.test.ts — tolerance for derived quantities
- `roundTo(x, 6)` rounds independently for each quantity. So `iqr`, `variance`, `se` are rounded separately from `q3-q1`, `sd^2`, `sd/√n`. Differences can be up to ~1e-6.
- `toBeCloseTo(x, 10)` fails because tolerance is `5e-11`. Use `toBeCloseTo(x, 5)` (tolerance `5e-6`) for derived-quantity consistency checks.

### Seeded LCG + Box-Muller for deterministic stress tests
- Numerical Recipes LCG: `s = (1664525·s + 1013904223) mod 2^32`. Box-Muller for N(μ,σ).
- All 100 stress datasets are deterministic and reproducible across platforms.
- LCG seeded per dataset (i×prime+offset) to ensure independence.

### ground_truth.json corrections (from previous session)
- sd was 2.0 (wrong, used population SD); correct Bessel-corrected value is 2.13809.
- BH p-adjustment formula: sort descending, multiply by n/rank, cumulative min, restore order.
