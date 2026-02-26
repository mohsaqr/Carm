# Carm — Complete Architecture & Engineering Reference

> **Purpose of this document**: If you are an LLM (or a human engineer) picking up Carm for the first time — after a crash, a context reset, or a handoff — this document tells you everything you need to know to continue working. It is the single source of truth for the architecture, every module, every trick, every convention, and every cross-validation target. Read this first, then read `CARM-PROMPT.md` for code-level rules.

---

## Table of Contents

1. [What is Carm](#1-what-is-carm)
2. [Philosophy & Standards](#2-philosophy--standards)
3. [Repository Layout](#3-repository-layout)
4. [Build, Test, Validate](#4-build-test-validate)
5. [TypeScript Strictness Contract](#5-typescript-strictness-contract)
6. [Core Layer — The Numerical Engine](#6-core-layer--the-numerical-engine)
7. [Stats Layer — 14 Pure Computation Modules](#7-stats-layer--14-pure-computation-modules)
8. [Viz Layer — 42 Publication-Quality Plots](#8-viz-layer--42-publication-quality-plots)
9. [Deterministic PRNG (splitmix32)](#9-deterministic-prng-splitmix32)
10. [Numerical Tricks & Stability Patterns](#10-numerical-tricks--stability-patterns)
11. [Cross-Validation Protocol & R Equivalence](#11-cross-validation-protocol--r-equivalence)
12. [Validation Infrastructure](#12-validation-infrastructure)
13. [Module Deep Dives](#13-module-deep-dives)
14. [Known Limitations & Open Issues](#14-known-limitations--open-issues)
15. [Commit & Workflow Rules](#15-commit--workflow-rules)
16. [Quick Reference Tables](#16-quick-reference-tables)

---

## 1. What is Carm

**Carm** is an enterprise-grade, zero-dependency TypeScript statistical analysis and visualization library. It runs entirely client-side — no R, no Python, no server. Every distribution, every hypothesis test, every matrix decomposition is implemented from first principles. The sole runtime dependency is D3.js (peer, lazy-loaded), used only by the visualization layer.

**Vital statistics:**

| Metric | Value |
|--------|-------|
| Language | TypeScript (ES2022, strict mode) |
| Total source | ~18,700 LOC |
| Core (math/types/matrix) | 1,490 LOC across 5 files |
| Stats (pure computation) | 6,957 LOC across 14 files |
| Viz (D3 renderers) | 9,728 LOC across 42 plot files + 5 component files |
| Tests | 496/496 passing (Vitest), 17 test files |
| Cross-validation | 200+ synthetic datasets vs R, 6+ real datasets |
| Build | tsup → ESM + CJS dual output with .d.ts |
| Package | `carm` v0.1.0, subpath exports for `/stats`, `/viz`, `/core` |
| Repository | Git, folder `JStats/`, branch `main` |

**Git path**: `/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/My Drive/Git/JStats`

---

## 2. Philosophy & Standards

### 2.1 Uncompromising R Equivalence

Carm does not approximate R. It matches R. Every function is cross-validated against the canonical R implementation with numerical tolerances documented per category:

| Category | Tolerance | Rationale |
|----------|-----------|-----------|
| Deterministic algorithms (eigenvalues, KMO, chi-square) | 1e-8 | IEEE-754 rounding only |
| ML extraction (uniquenesses, communalities) | MAE ≤ 0.05 | Different optimizer (Jöreskog vs L-BFGS-B) |
| Rotation loadings (promax, varimax) | MAE ≤ 0.001 | Same algorithm, same result |
| Geomin/oblimin with random starts=50 | MAE ≤ 0.001 | Global optimization matches lavaan |
| EM log-likelihood (GMM vs mclust) | ±2.0 | Different initialization (K-Means++ vs hierarchical) |
| LCA vs poLCA | 1e-10 | Exact MLE, same algorithm |
| K-Means vs R kmeans (Lloyd) | 1e-10 | Identical algorithm, same init |
| HAC heights vs R hclust | 1e-10 | Deterministic Lance-Williams |

### 2.2 Zero External Math Dependencies

No jStat. No simple-statistics. No numeric.js. Every mathematical function is implemented from scratch in `core/math.ts` and `core/matrix.ts`:

- **Gamma/Beta**: Lanczos approximation (Numerical Recipes §6.1)
- **Normal CDF/quantile**: Acklam rational approximation (~1e-9 error)
- **t/F/chi-square**: Via incomplete beta/gamma continued fractions
- **Matrix inverse**: LU with partial pivoting
- **Eigendecomposition**: Jacobi rotations for symmetric matrices
- **SVD**: Golub-Kahan bidiagonalization
- **Cholesky**: Standard lower-triangular factorization
- **Optimization**: Nelder-Mead simplex

**Why?** Full control over precision, edge cases, and testing. No version conflicts, no supply chain risk, no mysterious numerical discrepancies.

### 2.3 Every Result Has `formatted`

Every statistical function returns a structured `readonly` object with a `formatted: string` field containing the APA 7th edition result string. This string is ready for plot subtitles, clipboard copy, or publication:

```
t(48) = 2.31, p = .025, d = 0.65, 95% CI [0.08, 1.22]
F(2, 87) = 4.21, p = .018, ω² = 0.09
χ²(3) = 8.46, p = .037, V = 0.29
```

The `formatP()` function drops leading zeros: `p = .025`, not `p = 0.025`. This is APA convention.

### 2.4 Deterministic Reproducibility

Every stochastic algorithm (GMM, K-Means, LCA, LTA, random starts for rotation, parallel analysis, bootstrap) uses a seeded **splitmix32** PRNG. Default seed: 42. Same input + same seed = identical output, always. Every public function with randomness accepts `seed?: number` in its options.

---

## 3. Repository Layout

```
JStats/                          ← Git root (the project name is Carm, folder is still JStats)
├── src/
│   ├── core/                    ← 1,490 LOC — Pure math, types, formatting
│   │   ├── types.ts       (303) ← THE contract — ALL interfaces defined here
│   │   ├── math.ts        (547) ← Distributions, Gamma/Beta, Nelder-Mead, adjustPValues
│   │   ├── matrix.ts      (414) ← Matrix class: LU, Cholesky, SVD, eigen, inverse
│   │   ├── apa.ts         (221) ← APA 7th ed. formatters (formatTTest, formatANOVA, ...)
│   │   └── index.ts         (5) ← Re-exports
│   │
│   ├── stats/                   ← 6,957 LOC — Pure statistical computation (NO DOM)
│   │   ├── descriptive.ts (287) ← Mean, SD, skewness, kurtosis, Shapiro-Wilk (AS R94)
│   │   ├── comparison.ts  (476) ← t-tests, ANOVA, Mann-Whitney, Kruskal-Wallis, Friedman
│   │   ├── correlation.ts (282) ← Pearson, Spearman, Kendall, partial, correlationMatrix
│   │   ├── regression.ts  (379) ← OLS, multiple, polynomial, logistic (IRLS)
│   │   ├── clustering.ts (1850) ← GMM (6 models), LCA, LTA (Baum-Welch), K-Means, DBSCAN, HAC
│   │   ├── factor-analysis.ts (1923) ← EFA (ML/PAF), CFA, 6 rotations, diagnostics
│   │   ├── pca.ts         (156) ← PCA via SVD on correlation matrix
│   │   ├── mixed.ts       (297) ← LMM (profiled REML, random intercepts)
│   │   ├── post-hoc.ts    (221) ← Tukey HSD, Dunn, Games-Howell
│   │   ├── effect-size.ts (188) ← Cohen's d, Hedges' g, η², ω², rank-biserial, Cramér's V
│   │   ├── frequency.ts   (287) ← Contingency tables, chi-square, Fisher exact
│   │   ├── preprocess.ts  (181) ← Standardize, center, log, sqrt transforms
│   │   ├── analyze.ts     (416) ← Auto-dispatch: detects field types → routes to correct test
│   │   └── index.ts        (14) ← Re-exports all stats modules
│   │
│   ├── viz/                     ← ~10,000 LOC — D3 visualization renderers
│   │   ├── plots/          (42 files, 9,728 LOC total)
│   │   │   ├── violin-box.ts    (306) ← Violin + box + jitter + brackets
│   │   │   ├── scatter-stats.ts (159) ← Scatter + regression + CI band + marginals
│   │   │   ├── histogram.ts     (177) ← Histogram + density + normal curve
│   │   │   ├── fa-plot.ts       (600) ← 6 FA visualizations (scree, loadings, path, ...)
│   │   │   ├── dendrogram.ts    (230) ← HAC dendrogram (custom layout, not d3.hierarchy)
│   │   │   └── ... (37 more)
│   │   ├── components/     (4 files)
│   │   │   ├── annotations.ts   ← Title, subtitle, caption, regression equation
│   │   │   ├── brackets.ts      ← Significance brackets with p-values/stars
│   │   │   ├── axis.ts          ← Styled axes with labels
│   │   │   └── tooltip.ts       ← Hover tooltips (DOM-based, not SVG)
│   │   ├── themes/
│   │   │   └── default.ts       ← CarmTheme interface, CARM_PALETTE (8 colorblind-safe colors)
│   │   └── export.ts            ← SVG/PNG export (300 DPI)
│   │
│   └── index.ts                 ← Main entry point
│
├── tests/                       ← 17 test files, 496 tests
│   ├── stats/                   ← Per-module tests (clustering, comparison, correlation, ...)
│   ├── core/                    ← Math and matrix tests
│   ├── large_sample.test.ts     ← Cross-validation vs R (n=200-300)
│   ├── stress.test.ts           ← 100-dataset stress suite
│   └── r_*_reference.json       ← R-generated reference data
│
├── validation/                  ← PERMANENT — NEVER delete this folder
│   ├── VALIDATION-STRATEGY.md   ← Master plan for all validation rounds
│   ├── FA-TECHNICAL-REPORT.md   ← 1,085-line EFA/CFA engineering report
│   ├── RANDOM-STARTS-REPORT.md  ← 692-line random starts validation
│   ├── GMM-TECHNICAL-REPORT.md  ← Clustering engineering report
│   ├── data/                    ← JSON reference datasets (100 synthetic + real)
│   ├── r-reference/             ← R scripts that generate reference values
│   ├── ts-harness/              ← TypeScript validation harnesses
│   └── reports/                 ← Generated HTML cross-validation reports
│
├── tmp/                         ← Scratch scripts, comparison outputs (not committed)
├── CARM-PROMPT.md               ← System prompt for LLMs writing Carm code — READ THIS
├── CLAUDE.md                    ← Project-level coding rules
├── HANDOFF.md                   ← Session state for continuity
├── LEARNINGS.md                 ← Accumulated engineering insights
├── CHANGES.md                   ← Human-readable changelog
├── package.json                 ← npm config, build scripts
├── tsconfig.json                ← TypeScript strict config
└── tsup.config.ts               ← Build config (ESM + CJS dual output)
```

### Layer Discipline (CRITICAL)

```
core/  →  stats/  →  viz/
```

- **core/** is pure math. No stats logic, no DOM.
- **stats/** imports core/. Never imports viz/. No DOM, no D3, no side effects. Pure functions.
- **viz/** imports both core/ and stats/. Never imported by stats/.
- All shared interfaces live in `core/types.ts` — the single contract.

---

## 4. Build, Test, Validate

```bash
# From repo root:
npx tsup                                    # Build → dist/ (ESM + CJS + .d.ts)
npx tsc --noEmit                            # Type check (0 errors required)
npx vitest run                              # Run all 496 tests
npx tsx validation/ts-harness/fa-full-report.ts  # FA cross-validation → HTML report
```

### Package Exports

```json
{
  ".":       { "import": "./dist/index.js",       "require": "./dist/index.cjs" },
  "./stats": { "import": "./dist/stats/index.js", "require": "./dist/stats/index.cjs" },
  "./viz":   { "import": "./dist/viz/index.js",   "require": "./dist/viz/index.cjs" },
  "./core":  { "import": "./dist/core/index.js",  "require": "./dist/core/index.cjs" }
}
```

Consumers can import just the stats engine (`import { fitGMM } from 'carm/stats'`) without pulling in D3.

---

## 5. TypeScript Strictness Contract

**tsconfig.json** enforces maximum strictness:

```json
{
  "strict": true,
  "noUncheckedIndexedAccess": true,
  "exactOptionalPropertyTypes": true,
  "noUnusedLocals": true,
  "noUnusedParameters": true
}
```

**What this means in practice:**

| Rule | Effect | Pattern |
|------|--------|---------|
| `noUncheckedIndexedAccess` | `arr[i]` is `T \| undefined` | Use `arr[i]!` after bounds check |
| `noUncheckedIndexedAccess` | `arr[i] += x` won't compile | Use `arr[i]! += x` or `arr[i] = arr[i]! + x` |
| `exactOptionalPropertyTypes` | Cannot assign `undefined` to optional field | Use conditional spread: `...(x !== undefined && { x })` |
| `strict` | All returns and params have explicit types | No `any` — use `unknown` with guards |
| All results | Use `readonly` arrays and tuples | `readonly number[]`, `readonly [number, number]` |

---

## 6. Core Layer — The Numerical Engine

### 6.1 `core/math.ts` (547 LOC)

The mathematical foundation. Everything statistical builds on this.

**Special Functions:**

| Function | Algorithm | Citation |
|----------|-----------|----------|
| `logGamma(z)` | Lanczos approximation + reflection formula | Press et al., NR §6.1 |
| `gamma(z)` | `exp(logGamma(z))` | |
| `logBeta(a, b)` | `logΓ(a) + logΓ(b) - logΓ(a+b)` | |
| `incompleteBeta(x, a, b)` | Lentz continued fraction | NR §6.4 |
| `incompleteGamma(a, x)` | Series + continued fraction | NR §6.2 |

**Distribution Functions:**

| Function | Uses | Returns |
|----------|------|---------|
| `normalCDF(z)` | Error function polynomial | P(Z ≤ z) |
| `normalQuantile(p)` | Acklam rational approximation (~1e-9) | z such that P(Z ≤ z) = p |
| `tDistCDF(t, df)` | incompleteBeta | P(T ≤ t) |
| `tDistQuantile(p, df)` | Bisection (100 iters, tol 1e-10) | t such that P(T ≤ t) = p |
| `tDistPValue(t, df)` | Two-tailed via tDistCDF | p-value |
| `fDistCDF(f, df1, df2)` | incompleteBeta | P(F ≤ f) |
| `fDistPValue(f, df1, df2)` | Upper-tail | p-value |
| `chiSqCDF(x, df)` | incompleteGamma | P(χ² ≤ x) |
| `chiSqPValue(x, df)` | Upper-tail | p-value |
| `chiSqQuantile(p, df)` | Bisection | Quantile |

**Utilities:**

| Function | Notes |
|----------|-------|
| `mean(x)`, `median(x)`, `variance(x)`, `sd(x)`, `se(x)` | Bessel-corrected (n-1) |
| `quantile(x, p)` | Linear interpolation between sorted values |
| `rank(x)` | Average rank with tie-breaking |
| `sortAsc(x)` | Stable merge sort |
| `cov(x, y)` | Bessel-corrected covariance |
| `roundTo(x, decimals)` | Consistent rounding |
| `adjustPValues(pvals, method)` | Bonferroni, Holm, BH/FDR, BY |

**Optimization:**

| Function | Algorithm | Used By |
|----------|-----------|---------|
| `nelderMead(fn, x0, opts)` | Simplex direct search | ML extraction, LMM REML |

### 6.2 `core/matrix.ts` (414 LOC)

Immutable matrix class. All data stored as flat row-major `readonly number[]`.

**Construction:**

| Method | Purpose |
|--------|---------|
| `new Matrix(rows, cols, data?)` | From flat array |
| `Matrix.fromArray(2D)` | From 2D array |
| `Matrix.identity(n)` | n×n identity |
| `Matrix.zeros(rows, cols)` | Zero matrix |

**Operations:**

| Method | Returns | Algorithm |
|--------|---------|-----------|
| `get(i, j)` | `number` | Bounds-checked |
| `toArray()` | `number[][]` | Convert to 2D |
| `transpose()` | `Matrix` | Row-major transpose |
| `multiply(B)` | `Matrix` | O(n³) DGEMM-style |
| `scale(s)` | `Matrix` | Scalar × each element |
| `add(B)`, `subtract(B)` | `Matrix` | Element-wise |
| `inverse()` | `Matrix` | LU with partial pivoting |
| `pseudoInverse(tol?)` | `Matrix` | Moore-Penrose via SVD |
| `cholesky()` | `Matrix` | Lower-triangular L: A = LL' |
| `logDet()` | `number` | Via Cholesky (2·Σlog(L_ii)) |
| `svd()` | `{U, S, V}` | Golub-Kahan + Jacobi |
| `eigen()` | `{values, vectors}` | Jacobi rotations (symmetric) |
| `trace()` | `number` | Sum of diagonal |
| `norm()` | `number` | Frobenius norm |

**Key design**: Matrix is **immutable**. No `set(i,j,v)` method. Iterative algorithms (EM, gradient descent) use mutable `number[][]` arrays internally, then call `Matrix.fromArray()` at boundaries when matrix operations are needed.

### 6.3 `core/types.ts` (303 LOC)

**The contract.** Every interface used across layers is defined here.

**Universal result types:**

```typescript
interface StatResult {
  readonly testName: string
  readonly statistic: number
  readonly df: number | readonly [number, number]
  readonly pValue: number
  readonly effectSize: EffectSize
  readonly ci: readonly [number, number]
  readonly ciLevel: number
  readonly n: number
  readonly formatted: string                    // ← APA string, ALWAYS present
}

interface EffectSize {
  readonly value: number
  readonly name: string                         // "Cohen's d", "η²", "Cramér's V"
  readonly interpretation: EffectInterpretation  // 'negligible'|'small'|'medium'|'large'
}
```

**Domain-specific result types** (all extend or follow the same `readonly` + `formatted` pattern):

- `DescriptiveResult` — 20+ fields: mean, median, sd, variance, skewness, kurtosis, shapiroWilk, ci, ...
- `ANOVAResult extends StatResult` — adds msWithin, dfWithin, msGroup, dfGroup, groups
- `RegressionResult` — coefficients[], r2, adjR2, fStatistic, aic, bic, residuals, fitted
- `FAResult` — loadings[][], uniqueness[], communalities[], factorCorrelations, fit: FactorFit, eigenvalues
- `CFAResult extends FAResult` — adds parameterEstimates with SEs and z-values
- `FADiagnostics` — kmo, kmoPerItem, bartlett, mapSuggested, parallelEigenvalues, parallelSuggested
- `FactorFit` — chiSq, df, pValue, rmsea, rmseaCI, cfi, tli, srmr, aic, bic
- `GMMResult` — weights, means, covariances, posteriors, labels, diagnostics: ClusterDiagnostics
- `ClusterDiagnostics` — converged, iterations, logLikelihood, df, aic, bic, icl, entropy, avepp, formatted
- `LMMResult` — fixedEffects, varianceComponents, icc, logLik, aic, bic
- `PairwiseResult` — group1, group2, meanDiff, se, statistic, pValue, pValueAdj, ci, significant
- `AnalysisResult` — test, outcome, predictor, result, descriptives, posthoc, normality

### 6.4 `core/apa.ts` (221 LOC)

APA 7th edition formatting. Key functions:

| Function | Output Example |
|----------|---------------|
| `formatP(p)` | `"p = .025"` or `"p < .001"` (drops leading zero) |
| `formatStat(v, dec)` | `"2.31"` |
| `formatCI([lo, hi])` | `"[0.08, 1.22]"` |
| `formatTTest(t, df, p, d, ci)` | `"t(48) = 2.31, p = .025, d = 0.65, 95% CI [0.08, 1.22]"` |
| `formatANOVA(F, df1, df2, p, ω²)` | `"F(2, 87) = 4.21, p = .018, ω² = 0.09"` |
| `formatChiSq(χ², df, p, V)` | `"χ²(3) = 8.46, p = .037, V = 0.29"` |
| `formatCorrelation(r, df, p, ci)` | `"r(48) = 0.54, p = .003, 95% CI [0.28, 0.73]"` |
| `formatCFAFit(fit)` | `"χ²(24) = 28.42, p = .241, RMSEA = 0.042 [0.000, 0.085], CFI = 0.987, ..."` |

---

## 7. Stats Layer — 14 Pure Computation Modules

All functions are pure: data in → result out. No DOM, no D3, no side effects.

### 7.1 `descriptive.ts` (287 LOC)

| Function | Returns | Algorithm |
|----------|---------|-----------|
| `describe(x, ciLevel?)` | `DescriptiveResult` | Full univariate summary (20 fields) |
| `trimmedMean(x, alpha?)` | `number` | α-trimmed (default 5%) |
| `skewness(x)` | `number` | Adjusted Fisher-Pearson: `n/((n-1)(n-2)) × Σ((x-m)/s)³` |
| `kurtosis(x)` | `number` | Excess kurtosis (Fisher-Pearson adjusted) |
| `ciMean(x, ciLevel?)` | `[lo, hi]` | t-based: `m ± t_{α/2,n-1} × se` |
| `shapiroWilk(x)` | `{statistic, pValue}` | Full AS R94 (Royston 1995) — matches R exactly |

**Shapiro-Wilk implementation details**: Uses the complete Royston (1995) polynomial correction algorithm. For n=3..5, uses exact tabulated coefficients. For n≥6, expected normal quantiles → Modified Gram-Schmidt QR → polynomial corrections via C1/C2 arrays → W statistic → p-value via Royston's transformation. Valid for n=3..5000.

### 7.2 `comparison.ts` (476 LOC)

| Function | Test | Effect Size | Notes |
|----------|------|-------------|-------|
| `tTestIndependent(x1, x2, equalVar?, ci?, alt?)` | Welch's / Student's t | Cohen's d | Default: Welch's (unequal variance) |
| `tTestPaired(x1, x2, ci?)` | Paired t | Cohen's d (paired) | Requires equal-length arrays |
| `oneWayANOVA(groups)` | F-test | **ω²** (omega-squared, NOT η²) | Returns msWithin/dfWithin for post-hoc |
| `mannWhitneyU(x1, x2, alt?)` | Mann-Whitney U | Rank-biserial r | Normal approx with tie correction |
| `wilcoxonSignedRank(x1, x2)` | Wilcoxon signed-rank | Rank-biserial r | Exact p for n≤20, normal approx n>20 |
| `kruskalWallis(groups)` | Kruskal-Wallis H | η²_H | Tie-corrected H statistic |
| `friedmanTest(data)` | Friedman χ² | Kendall's W | Repeated measures |

**Important**: `oneWayANOVA` returns **ω² (omega-squared)**, which is a bias-corrected effect size. R's default `eta²` is larger. Both are correct — they measure different things. ω² = (SS_b - df_b × MS_w) / (SS_total + MS_w).

### 7.3 `correlation.ts` (282 LOC)

| Function | Method | CI | Notes |
|----------|--------|-----|-------|
| `pearsonCorrelation(x, y, ci?, alt?)` | Pearson r | Fisher z-transform | Supports one-tailed |
| `spearmanCorrelation(x, y, ci?)` | Spearman ρ | Pearson on ranks | t-approximation for p |
| `kendallTau(x, y)` | Kendall τ-b | Normal approx | Tie-corrected |
| `partialCorrelation(x, y, controls)` | Partial r | Fisher z | Residualize both on controls |
| `correlationMatrix(data, method?)` | Pairwise r/ρ/τ | — | Returns r[][] and p[][] |

### 7.4 `regression.ts` (379 LOC)

| Function | Model | Notes |
|----------|-------|-------|
| `linearRegression(x, y, ci?)` | Simple OLS | β = (X'X)⁻¹X'y |
| `multipleRegression(y, predictors[], ci?)` | Multiple OLS | Design matrix [1, x₁, x₂, ...] |
| `polynomialRegression(x, y, degree, ci?)` | Polynomial OLS | Design matrix [1, x, x², ..., xᵈ] |
| `logisticRegression(y, predictors[], ci?)` | Binary logistic | IRLS (Newton-Raphson on Fisher information) |

All regression returns include: coefficients (estimate, SE, t, p, CI), R², adjR², F-statistic, AIC, BIC, residuals, fitted values.

### 7.5 `clustering.ts` (1,850 LOC) — THE CLUSTERING ENGINE

**Six algorithms**, all with deterministic PRNG:

#### fitGMM — Gaussian Mixture Models

```typescript
fitGMM(data: number[][], options?: {
  k?: number, model?: 'VVV'|'EEE'|'VVI'|'EEI'|'VII'|'EII',
  seed?: number, tol?: number, maxIter?: number, regCovar?: number
}): GMMResult
```

**EM algorithm:**
1. **Init**: K-Means++ (seeded PRNG) → initial means, equal weights, empirical covariance
2. **E-step**: Compute responsibilities z_ik = w_k × φ(x_i; μ_k, Σ_k) / Σ w_l φ(...) — **all in log-space** to prevent underflow
3. **M-step**: Update weights, means, covariances per constraint model
4. **Convergence**: |ΔLL| < tol

**Six covariance models** (mclust naming):

| Model | Constraints | Free Params | Description |
|-------|-----------|-------------|-------------|
| VVV | None | k × d(d+1)/2 | Full per-cluster |
| EEE | Σ_1 = ... = Σ_k | d(d+1)/2 | Single shared full |
| VVI | Σ_k = diag(σ²_{k1}, ..., σ²_{kd}) | k × d | Diagonal per-cluster |
| EEI | Σ_k = diag(σ²_1, ..., σ²_d) | d | Single shared diagonal |
| VII | Σ_k = σ²_k × I | k | Spherical per-cluster |
| EII | Σ_k = σ² × I | 1 | Single shared spherical |

**Numerical trick**: MVN log-PDF via eigendecomposition. Instead of computing Σ⁻¹ directly: `log φ(x; μ, Σ) = -d/2·log(2π) - ½·Σlog(λ_j) - ½·Σ((V'(x-μ))_j / √λ_j)²`. Avoids matrix inversion entirely.

#### fitLCA — Latent Class Analysis

```typescript
fitLCA(data: number[][], options?: { k?, seed?, tol?, maxIter? }): LCAResult
```

EM on binary (0/1) data. Parameters: rho[k][d] = P(item_d = 1 | class k), weights[k]. Uses raw MLE (no Beta(1,1) smoothing — this is critical for matching poLCA exactly). Floor clamping at 1e-10.

#### fitLTA — Latent Transition Analysis

```typescript
fitLTA(data: number[][][], options?: { k?, seed?, tol?, maxIter? }): LTAResult
```

Hidden Markov LCA. Baum-Welch algorithm **entirely in log-space**:
- Forward: `log α_t[s]` accumulates via logSumExp
- Backward: `log β_t[s]`
- Posteriors (γ), transition posteriors (ξ): all log-space
- Viterbi decoding for most-likely state sequence

#### runKMeans — K-Means

```typescript
runKMeans(data: number[][], options?: { k?, seed?, maxIter?, tol? }): KMeansResult
```

Lloyd's algorithm with K-Means++ initialization. Empty-cluster re-seeding (assign farthest point).

#### runDBSCAN

```typescript
runDBSCAN(data: number[][], options?: { eps?, minPts? }): DBSCANResult
```

Density-based clustering. Labels: -1 = noise, 0+ = cluster index.

#### runHierarchical — HAC

```typescript
runHierarchical(data: number[][], options?: { k?, linkage? }): HACResult
```

Lance-Williams recurrence formula. Four linkage methods: single, complete, average, ward.D2.

**Ward.D2 detail**: Internally uses squared Euclidean distances; final heights = √(merge distance). Coefficients: α_i = (n_i + n_k)/N_t, α_j = (n_j + n_k)/N_t, β = -n_k/N_t, γ = 0.

#### Helper functions:

| Function | Purpose |
|----------|---------|
| `predictGMM(model, newData)` | Posterior probabilities for new observations |
| `predictKMeans(model, newData)` | Assign new points to nearest centroid |
| `findBestGMM(data, kRange?, models?)` | Grid search K × model, return lowest BIC |
| `fitGMMRange(data, kRange, model?)` | Fit across K values, return sorted array |
| `fitKMeansRange(data, kRange)` | Fit K-Means across K values |
| `euclideanDistMatrix(data)` | N×N pairwise distance matrix |
| `silhouetteScores(data, labels)` | Silhouette per point (excludes noise) |
| `kDistancePlot(data, k)` | k-distance values for DBSCAN eps selection |
| `cutTree(merges, heights, n, k)` | Cut HAC at k clusters |

### 7.6 `factor-analysis.ts` (1,923 LOC) — THE CROWN JEWEL

This is the largest and most sophisticated module. It implements EFA with two extraction methods, six rotation methods (including a full gradient projection engine for oblique rotations with random starts), CFA with ML estimation and numerical Hessian standard errors, and four psychometric diagnostics.

#### Public API:

```typescript
runEFA(data: number[][], options?: EFAOptions): FAResult
runCFA(model: Record<string, number[]>, data: number[][], options?: CFAOptions): CFAResult
runFADiagnostics(data: number[][], options?: FADiagnosticsOptions): FADiagnostics
```

#### EFAOptions (all optional):

| Field | Default | Description |
|-------|---------|-------------|
| `nFactors` | Auto (parallel analysis) | Number of factors to extract |
| `extraction` | `'ml'` | `'ml'` (Maximum Likelihood) or `'paf'` (Principal Axis Factoring) |
| `rotation` | `'promax'` | `'varimax'`, `'promax'`, `'geomin'`, `'oblimin'`, `'quartimin'`, `'none'` |
| `randomStarts` | **50** | Number of random orthogonal starting matrices for oblique rotations |
| `geominDelta` | 0.01 | Geomin δ (smoothing parameter) |
| `seed` | 42 | PRNG seed for random starts and parallel analysis |
| `maxIter` | 1000 | Maximum iterations for rotation |
| `tol` | 1e-6 | Convergence tolerance |

#### ML Extraction — Two-Phase Optimization

**Phase 1: Jöreskog gradient descent** on uniquenesses Ψ

```
Init: Ψ = (1 - 0.5k/d) / diag(R⁻¹)    ← R's factanal formula
Loop:
  Eigen(Ψ⁻¹/² R Ψ⁻¹/²) → λ, V
  L = √Ψ × V × diag(√max(λ-1, 0))
  Σ = LL' + Ψ
  Gradient = diag(Σ⁻¹(Σ-R)Σ⁻¹)
  Ψ ← clamp(Ψ - lr × grad, [0.005, 0.995])
  lr uses cosine annealing
```

**Phase 2: Nelder-Mead polish** on R's concentrated ML objective

```
F(Ψ) = -Σ_{j>k} (log λ_j - λ_j) + k - d
```

This hybrid achieves matching uniquenesses (max Δ ~5e-5) vs R across 100 synthetic datasets.

#### PAF Extraction — Iterated Principal Axis

```
Init: h² = SMC (Squared Multiple Correlation) = 1 - 1/R⁻¹_ii
Loop:
  R_reduced = R with diagonal replaced by h²
  Eigen(R_reduced) → top k eigenvalues/vectors
  L = V × diag(√λ)
  h²_new = rowSums(L²), clamped [0.001, 0.9999]
  Converge when max|Δh²| < tol
```

#### Varimax — SVD-Based Orthogonal Rotation

Kaiser normalization → polar decomposition via eigendecomposition of B'B (avoids buggy SVD for small k×k). Matches R's `stats::varimax()` exactly.

#### Promax — Orthogonal-to-Oblique with Outer Kaiser

**Critical detail**: `psych::fa(rotate="promax")` wraps the promax call in `psych::kaiser()`, which adds **outer Kaiser normalization**:

```
1. h² = rowSums(L²)
2. L_norm = L / √h²              ← normalize rows by communality
3. V = varimax(L_norm)
4. Q = V ⊙ sign(V)|V|^(power-1)  ← promax target (on NORMALIZED loadings)
5. U = (V'V)⁻¹V'Q                ← rotation matrix
6. rotated_norm = V × U
7. rotated = rotated_norm × √h²   ← denormalize
```

Without this outer Kaiser, the promax target Q is computed on different-scale loadings, causing MAE=0.055 on difficult datasets.

#### GPFoblq — General Oblique Gradient Projection

```typescript
function gpfOblq(A, criterionFn, maxIter, tol, Tinit?):
  { rotated, T, Phi, f }
```

The **workhorse** for geomin, oblimin, and quartimin. Algorithm (Bernaards & Jennrich 2005):

```
T = Tinit ?? I_k
L = A × (T⁻¹)'                           ← rotated loadings
Loop:
  (f, Gq) = criterionFn(L)               ← criterion value + gradient w.r.t. L
  G = -(L' × Gq × T⁻¹)'                 ← gradient w.r.t. T
  Gp = G - T × diag(colSums(T ⊙ G))     ← oblique projection
  Line search (Armijo):
    al = 2·al
    for s = 0..10:
      X = T - al·Gp
      X columns normalized to unit length
      L_new = A × (X⁻¹)'
      f_new = criterionFn(L_new).f
      if improvement > 0.5·‖Gp‖²·al: break
      al /= 2
    T = X (ALWAYS update, even if Armijo not satisfied — matches R exactly)
  Converge when ‖Gp‖ < tol
Return: rotated = L, Phi = T'T, f = criterion value
```

**Two critical bugs fixed during development**:
1. Gradient used A (unrotated) instead of L (rotated) — only equivalent at T=I
2. T was only updated inside Armijo success block — R always updates after line search

#### Multi-Start Random Rotation

```typescript
function gpfOblqWithRandomStarts(loadings, criterionFn, maxIter, tol, k, randomStarts, seed):
  { rotated, Phi }
```

- **Start 0**: T = I (always, deterministic, backward-compatible)
- **Starts 1..randomStarts-1**: Haar-distributed random orthogonal matrices
- Keep solution with **lowest criterion value**
- Default: 50 starts. Matches lavaan's GPA×30 strategy.

**Random orthogonal matrix generation** (`randomOrthogonalMatrix(k, rng)`):

```
1. Fill k×k matrix M with N(0,1) from splitmix32 PRNG
2. Modified Gram-Schmidt QR: M = QR
3. Sign correction: Q[:,j] *= sign(R[j][j])
→ Produces Haar-distributed random orthogonal matrix
```

This matches both `GPArotation::Random.Start()` and `lavaan::lav_matrix_rotate_gen()`.

#### Why Random Starts Matter

With geomin δ=0.01 (GPArotation default), the criterion surface is smooth — single start finds the global minimum. With δ=0.001 (lavaan default), multiple local minima exist:

| Dataset | δ | 1 start MAE | 50 starts MAE |
|---------|---|-------------|---------------|
| rraw 525×31, k=5 | 0.001 | 0.043 | 0.000001 |
| Burnout 876×23, k=4 | 0.001 | 0.033 | 0.000002 |
| rraw 525×31, k=5 | 0.01 | 0.000003 | 0.000001 |

#### Criterion Functions

**Geomin** (Browne 2001):

```
f = Σ_i [Π_j (λ²_ij + δ)]^(1/k)
∂f/∂λ_ij = (2/k) × λ_ij/(λ²_ij + δ) × [Π_m (λ²_im + δ)]^(1/k)
```

**Oblimin** (Jennrich 2002):

```
f = (1/4) Σ_i Σ_{j≠m} λ²_ij λ²_im - (γ/4d) Σ_j Σ_m (Σ_i λ²_ij)(Σ_i λ²_im)
γ = 0 → quartimin
```

#### Canonical Sign Convention

After extraction (ML or PAF), for each factor column: if the largest absolute loading is negative, flip the entire column. This matches LAPACK's dsyev convention used by R. Essential for non-rotation-invariant methods (geomin, oblimin).

#### FA Diagnostics

| Diagnostic | Algorithm |
|------------|-----------|
| **KMO** | `Σr²_ij / (Σr²_ij + Σp²_ij)` where p = partial correlation |
| **Bartlett** | `-(n-1-(2p+5)/6) × log|R|` with df = p(p-1)/2 |
| **MAP** (Velicer 1976) | Extract 0..d-1 components, compute avg squared off-diagonal of residual correlation, return k at minimum |
| **Parallel Analysis** (Horn 1965) | 100 simulations of random normal data, compare observed eigenvalues to 95th percentile |

#### Fit Indices

| Index | Formula | Notes |
|-------|---------|-------|
| **χ²** | (n-1-(2p+5)/6-2k/3) × F_ml | Bartlett correction for EFA |
| **RMSEA** | √(max(χ²-df,0) / (df(n-1))) | Browne & Cudeck 1993 |
| **RMSEA CI** | Steiger approximation via NCP | 90% confidence interval |
| **CFI** | 1 - NCP_model / NCP_null | Bentler 1990, clamped [0,1] |
| **TLI** | (χ²_null/df_null - χ²/df) / (χ²_null/df_null - 1) | Can exceed 1.0 |
| **SRMR** | √(mean((r_obs - r_implied)²)) | On correlation scale |

### 7.7 `mixed.ts` (297 LOC) — Linear Mixed Models

```typescript
runLMM(input: LMMInput): LMMResult
```

Random-intercept LMM via profiled REML:

```
y = Xβ + Zb + ε
b ~ N(0, σ²_b I), ε ~ N(0, σ²_e I)
```

**Algorithm**:
1. Profile out σ²_e analytically (single parameter: logψ = log(σ²_b/σ²_e))
2. Multi-start Nelder-Mead: 5 starting points [logψ = -4, -2, 0, 2, 4]
3. Inverse via Woodbury: V⁻¹ = I - Z·D⁻¹·Z' where D = Z'Z + (1/ψ)I
4. GLS: β = (X'V⁻¹X)⁻¹X'V⁻¹y
5. ICC = σ²_b / (σ²_b + σ²_e)
6. logLik includes normalizing constant: -½(n-p)(1+log(2π))

**Known fix**: When ψ ≈ 0 (ICC ≈ 0), 1/ψ = ∞. Detect and set VinvScaled = I (GLS = OLS).

Cross-validated against R lme4. Matches logLik to <1 unit for n=300.

### 7.8 Other Modules (Brief)

**`post-hoc.ts` (221 LOC)**: Tukey HSD (studentized range), Games-Howell (Welch-type SE + Satterthwaite df), Dunn (rank-based + p-adjustment). All return `PairwiseResult[]`.

**`effect-size.ts` (188 LOC)**: Cohen's d (pooled SD), Hedges' g (bias-corrected), η² (SS_b/SS_total), ω² (bias-corrected), rank-biserial (1-2U/(n₁n₂)), Cramér's V (√(χ²/(n·(min(R,C)-1)))). All with interpretation functions.

**`frequency.ts` (287 LOC)**: frequencyTable (counts/relative/cumulative), contingencyTable (cross-tab), chiSquareTest (with optional Yates correction), fisherExactTest (hypergeometric PMF for 2×2), goodnessOfFit.

**`pca.ts` (156 LOC)**: PCA via SVD on correlation matrix. Returns loadings, scores, eigenvalues, variance explained.

**`preprocess.ts` (181 LOC)**: preprocessData (none/center/standardize/log/sqrt). Returns transformed data + column means/SDs for inverse transform.

**`analyze.ts` (416 LOC)**: High-level auto-dispatch. Detects field types → checks normality (Shapiro-Wilk per group, n∈[3,50]) → routes to parametric or non-parametric test → runs post-hoc for 3+ groups. Returns structured AnalysisResult.

---

## 8. Viz Layer — 42 Publication-Quality Plots

### 8.1 Architecture

Every plot follows the same pattern:

```typescript
export function renderViolinBox(
  container: HTMLElement,
  data: ViolinBoxData,
  config: ViolinBoxConfig = {}
): void {
  import('d3').then(d3 => renderViolinBoxD3(d3, container, data, config))
}
```

D3 is lazy-loaded as a peer dependency. The stats layer has zero visualization cost.

### 8.2 Theme System

```typescript
interface CarmTheme {
  readonly background: string      // '#ffffff'
  readonly text: string            // Foreground
  readonly textMuted: string       // Secondary text
  readonly gridLine: string        // Grid color
  readonly colors: readonly string[] // 8-color palette
  readonly fontFamily: string      // System fonts
  readonly fontSize: number        // 12px
  readonly pointOpacity: number    // 0.55
  readonly violinOpacity: number   // 0.72
  // ... margins, font sizes, etc.
}
```

**CARM_PALETTE**: 8 colorblind-safe colors (Tableau 10 inspired): blue, orange, red, teal, green, mauve, pink, umber.

### 8.3 Complete Plot Inventory

| # | Plot | File | Purpose |
|---|------|------|---------|
| 1 | Violin + Box + Jitter | `violin-box.ts` | Group comparison (ggbetweenstats style) |
| 2 | Scatter + Regression | `scatter-stats.ts` | Correlation with CI band + marginals |
| 3 | Histogram + Density | `histogram.ts` | Distribution with normal curve overlay |
| 4 | Bar + Significance | `bar-stats.ts` | Categorical with brackets |
| 5 | Correlogram | `correlogram.ts` | Correlation matrix heatmap |
| 6 | Coefficient Plot | `coef-plot.ts` | Dot-and-whisker with CI |
| 7 | Q-Q Plot | `qq-plot.ts` | Normality check with confidence band |
| 8 | Residual Panel | `residual-panel.ts` | 4-panel regression diagnostics |
| 9 | Raincloud | `raincloud.ts` | Half-violin + jitter + box |
| 10 | Boxplot | `boxplot.ts` | Traditional box-and-whisker |
| 11 | Density | `density.ts` | Density plot with area fill |
| 12 | Strip Plot | `strip-plot.ts` | Jittered points only |
| 13 | Swarm Plot | `swarm-plot.ts` | Beeswarm force-direct layout |
| 14 | Dot Plot | `dot-plot.ts` | Cleveland dot plot |
| 15 | Lollipop | `lollipop.ts` | Dots on stems |
| 16 | Line Chart | `line-chart.ts` | Time series |
| 17 | Area Chart | `area-chart.ts` | Stacked area |
| 18 | Grouped Bar | `grouped-bar.ts` | Multiple series grouped/stacked |
| 19 | Bubble Chart | `bubble-chart.ts` | 3D scatter (x, y, size) |
| 20 | Pareto | `pareto.ts` | Bars + cumulative line |
| 21 | Funnel | `funnel.ts` | Stage-wise reduction |
| 22 | Pie/Donut | `pie-chart.ts` | Proportional |
| 23 | Waffle | `waffle-chart.ts` | 100-square grid |
| 24 | Sparkline | `sparkline.ts` | Small inline chart |
| 25 | Forest Plot | `forest-plot.ts` | Odds ratio / CI visualization |
| 26 | ROC Curve | `roc-curve.ts` | Receiver operating characteristic |
| 27 | Radar | `radar-chart.ts` | Spider chart |
| 28 | Pair Plot | `pair-plot.ts` | Scatterplot matrix |
| 29 | Parallel Coords | `parallel-coords.ts` | Multivariate comparison |
| 30 | Mosaic | `mosaic-plot.ts` | Contingency table |
| 31 | Treemap | `treemap.ts` | Hierarchical rectangles |
| 32 | Sunburst | `sunburst.ts` | Radial treemap |
| 33 | Marimekko | `marimekko.ts` | Variable-width stacked bar |
| 34 | Chord Diagram | `chord-diagram.ts` | Flow/connection |
| 35 | Arc Diagram | `arc-diagram.ts` | Hierarchical network |
| 36 | Alluvial/Sankey | `alluvial-plot.ts` | Flow diagram |
| 37 | Edge Bundling | `edge-bundling.ts` | Hierarchical edge bundling |
| 38 | Dendrogram | `dendrogram.ts` | HAC tree (custom layout) |
| 39 | PCA Plot | `pca-plot.ts` | Score + loading plots |
| 40 | Distribution Explorer | `distribution.ts` | Interactive PDF/CDF |
| 41 | Mixed Plot | `mixed-plot.ts` | Composite multi-axis |
| 42 | FA Plot | `fa-plot.ts` | 6 types: scree, loadings heatmap, path diagram, communality bar, factor correlation, fit dashboard |

### 8.4 Components

| Component | File | Purpose |
|-----------|------|---------|
| Annotations | `annotations.ts` | Title + subtitle + caption + regression equation |
| Brackets | `brackets.ts` | Significance brackets with p-values/stars |
| Axis | `axis.ts` | Styled axes with labels |
| Tooltip | `tooltip.ts` | Hover tooltips (DOM element, not SVG) |
| Export | `export.ts` | SVG download + PNG at 300 DPI |

---

## 9. Deterministic PRNG (splitmix32)

```typescript
class PRNG {
  private state: number
  constructor(seed: number) { this.state = seed >>> 0 }
  next(): number {
    this.state = (this.state + 0x9E3779B9) | 0
    let t = this.state ^ (this.state >>> 16)
    t = Math.imul(t, 0x21F0AAAD)
    t = t ^ (t >>> 15)
    t = Math.imul(t, 0x735A2D97)
    t = t ^ (t >>> 15)
    return (t >>> 0) / 4294967296   // uniform [0, 1)
  }
}
```

**Properties**: Period 2³², uniform distribution, deterministic. Default seed: 42.

**Normal variates**: Box-Muller transform: `√(-2·log(u)) × cos(2πv)` where u, v ~ Uniform(0,1).

**Used by**: GMM (K-Means++ init), LCA, LTA, K-Means, factor analysis (random starts, parallel analysis), bootstrap.

---

## 10. Numerical Tricks & Stability Patterns

### 10.1 Log-Space Arithmetic

All probabilistic models (GMM, LCA, LTA) compute in log-space to prevent underflow:

```typescript
function logSumExp(arr: ArrayLike<number>): number {
  let max = -Infinity
  for (let i = 0; i < arr.length; i++) if (arr[i]! > max) max = arr[i]!
  if (max === -Infinity) return -Infinity
  let sum = 0
  for (let i = 0; i < arr.length; i++) sum += Math.exp(arr[i]! - max)
  return max + Math.log(sum)
}
```

Pattern: Compute log-probabilities → logSumExp for normalization → exp only at the end.

### 10.2 Eigendecomposition for MVN Log-PDF

Instead of computing Σ⁻¹ directly:

```
Σ = V diag(λ) V'
log φ(x; μ, Σ) = -d/2·log(2π) - ½·Σlog(λ_j) - ½·Σ((V'(x-μ))_j)² / λ_j
```

Avoids matrix inversion, numerically stable even for near-singular covariances.

### 10.3 Float64Array for Hot Accumulation Loops

Covariance computation, EM responsibilities, and other inner loops use `Float64Array` for cache-friendly sequential access.

### 10.4 Regularization for Singularity Prevention

GMM covariances: `Σ_k += regCovar × I` (default regCovar = 1e-6). Prevents singular covariance during EM when a cluster shrinks.

### 10.5 Death Penalty in Optimization

Nelder-Mead boundary constraints: when uniquenesses Ψ_i leave [0.005, 0.995], return `Infinity` as the objective value. Forces the simplex back into the feasible region.

### 10.6 Woodbury Matrix Identity for LMM

```
V⁻¹ = I - Z·(Z'Z + (1/ψ)I)⁻¹·Z'
```

O(n·g²) instead of O(n³) for n observations, g groups.

### 10.7 Cosine-Annealed Learning Rate

ML extraction Phase 1 uses cosine annealing:

```
lr = lr_min + (lr_0 - lr_min) × 0.5 × (1 + cos(π·t/T))
```

Prevents oscillation in later iterations while allowing large initial steps.

---

## 11. Cross-Validation Protocol & R Equivalence

### 11.1 The Protocol

For every statistical function:

1. Write the TypeScript implementation
2. Create identical test data (deterministic or from reference JSON)
3. Run the same analysis in R (document the R code in a comment block)
4. Compare results numerically with documented tolerance
5. If results differ: investigate, document why in LEARNINGS.md
6. The test must encode the R-verified expected values

### 11.2 R Equivalence Map

| Carm Function | R Package::function | Verified? | Notes |
|---------------|---------------------|-----------|-------|
| `shapiroWilk` | `stats::shapiro.test` | Yes | AS R94, exact match |
| `tTestIndependent` | `stats::t.test` | Yes | Welch and Student |
| `oneWayANOVA` | `stats::aov` | Yes | Returns ω² not η² |
| `mannWhitneyU` | `stats::wilcox.test` | Yes | Tie correction |
| `kruskalWallis` | `stats::kruskal.test` | Yes | |
| `pearsonCorrelation` | `stats::cor.test` | Yes | Fisher z CI |
| `spearmanCorrelation` | `stats::cor.test(method="spearman")` | Yes | |
| `kendallTau` | `stats::cor.test(method="kendall")` | Yes | Normal approx |
| `chiSquareTest` | `stats::chisq.test` | Yes | |
| `fisherExactTest` | `stats::fisher.test` | Yes | Hypergeometric |
| `linearRegression` | `stats::lm` | Yes | R², AIC, BIC |
| `logisticRegression` | `stats::glm(family=binomial)` | Yes | IRLS |
| `runEFA(extraction='paf')` | `psych::fa(fm="pa")` | Yes | 100/100 synthetic |
| `runEFA(extraction='ml')` | `psych::fa(fm="ml")` / `stats::factanal` | Yes | 100/100 synthetic |
| `runEFA(rotation='varimax')` | `stats::varimax` | Yes | SVD-based |
| `runEFA(rotation='promax')` | `psych::fa(rotate="promax")` | Yes | With outer Kaiser |
| `runEFA(rotation='geomin')` | `GPArotation::geominQ` / lavaan | Yes | 100/100 + random starts |
| `runEFA(rotation='oblimin')` | `GPArotation::oblimin` | Yes | Same GPFoblq engine |
| `runFADiagnostics` | `psych::KMO`, `stats::bartlett.test` | Yes | 100/100 |
| `fitGMM` | `mclust::Mclust` | Yes | 6 covariance models |
| `fitLCA` | `poLCA::poLCA(nrep=1)` | Yes | 10-decimal match |
| `runKMeans` | `stats::kmeans(algorithm="Lloyd")` | Yes | 10-decimal match |
| `runHierarchical` | `stats::hclust` | Yes | 4 linkages, 1e-10 |
| `runDBSCAN` | `dbscan::dbscan` | Yes | Labels: -1=noise |
| `runLMM` | `lme4::lmer(REML=TRUE)` | Yes | logLik <1 unit |
| `runPCA` | `stats::prcomp` | Yes | SVD-based |

### 11.3 What's NOT Yet Cross-Validated

- CFA vs lavaan (estimation + standard errors)
- Oblimin/quartimin rotations (use same GPFoblq, should work)
- LTA (no deterministic R equivalent — depmixS4 is non-deterministic)

---

## 12. Validation Infrastructure

**Location**: `validation/` (PERMANENT — never delete)

```
validation/
├── VALIDATION-STRATEGY.md        ← Master plan for all validation rounds
├── FA-TECHNICAL-REPORT.md        ← 1,085 lines: complete EFA/CFA engineering report
├── RANDOM-STARTS-REPORT.md       ← 692 lines: random starts validation
├── GMM-TECHNICAL-REPORT.md       ← Clustering engineering report
├── data/
│   ├── fa-crossval-data.json     ← 100 synthetic datasets + R promax results
│   ├── fa-geomin-ref.json        ← R geomin results for 100 datasets
│   ├── fa-real-crossval-ref.json ← R results for real 525×31 dataset
│   └── ...
├── r-reference/
│   ├── fa-promax-ref.R           ← R code generating promax reference
│   ├── fa-geomin-ref.R           ← R code generating geomin reference
│   └── ...
├── ts-harness/
│   ├── fa-full-report.ts         ← Main FA cross-validation → HTML report
│   └── fa-geomin-crossval.ts     ← Geomin-specific validation
└── reports/
    └── fa-crossval-report.html   ← Generated report
```

**Current pass rates**: Promax 100/100, Geomin 100/100, Real dataset 6/6, Diagnostics 100/100.

**To run validation**: `npx tsx validation/ts-harness/fa-full-report.ts` (from repo root)

---

## 13. Module Deep Dives

### 13.1 How GMM EM Works in Carm

```
Input: data[n×d], k clusters, covariance model
1. K-Means++ init:
   - Pick random point as first centroid (seeded PRNG)
   - For each subsequent centroid: pick proportional to distance² from nearest existing
2. Initial params:
   - weights = 1/k
   - means = K-Means++ centroids
   - covariances = empirical covariance per cluster + regularization
3. EM loop:
   E-step:
     For each i,j: log_resp[i][j] = log(w_j) + mvnLogPdf(x_i, μ_j, Σ_j)
     Normalize: resp[i][j] = exp(log_resp[i][j] - logSumExp(log_resp[i]))
   M-step:
     N_j = Σ_i resp[i][j]
     w_j = N_j / n
     μ_j = Σ_i resp[i][j] × x_i / N_j
     Σ_j = (Σ_i resp[i][j] × (x_i-μ_j)(x_i-μ_j)') / N_j + regCovar×I
     Apply covariance constraint (VVV/EEE/VVI/EEI/VII/EII)
   LL = Σ_i logSumExp(log_resp[i])
   If |LL - LL_prev| < tol: converge
4. Compute diagnostics: BIC, AIC, ICL, entropy, AvePP
5. Labels = argmax(resp[i])
```

### 13.2 How GPFoblq With Random Starts Works

```
Input: A (p×k unrotated loadings), criterionFn (geomin/oblimin/quartimin),
       randomStarts=50, seed=42

1. Run gpfOblq(A, criterionFn, maxIter, tol)  // Start 0: T = I (deterministic)
   → best = {rotated, Phi, f}

2. rng = new PRNG(seed)
   For s = 1 to randomStarts-1:
     T_rand = randomOrthogonalMatrix(k, rng)
     result = gpfOblq(A, criterionFn, maxIter, tol, T_rand)
     If result.f < best.f:
       best = result

3. Return best.rotated, best.Phi
```

### 13.3 How LTA Baum-Welch Works (All Log-Space)

```
Data: sequences[n][T][M] (n subjects, T timepoints, M binary items)

Forward (per subject):
  log_α[0][s] = log(π[s]) + Σ_m log_emission(o[0][m]; ρ[s][m])
  log_α[t][s] = Σ_m log_emission(o[t][m]; ρ[s][m]) + logSumExp_p(log_α[t-1][p] + log(τ[p][s]))

Backward (per subject):
  log_β[T-1][s] = 0
  log_β[t][s] = logSumExp_l(log(τ[s][l]) + Σ_m log_emission(o[t+1][m]; ρ[l][m]) + log_β[t+1][l])

Posteriors:
  γ_t[s] = exp(log_α[t][s] + log_β[t][s] - logSumExp_s(log_α[T-1][s]))
  ξ_t[s][l] = exp(log_α[t][s] + log(τ[s][l]) + emission[t+1][l] + log_β[t+1][l] - LL)

M-step:
  π[s] = (Σ_i γ[0][s] + 1) / (n + k)         // Dirichlet(1) smoothing
  τ[s][l] = Σ_i Σ_t ξ_t[s][l] / Σ_i Σ_t γ_t[s]
  ρ[s][m] = Σ_i Σ_t γ_t[s] × o[t][m] / Σ_i Σ_t γ_t[s]

Viterbi: Dynamic programming for most-likely state path
```

---

## 14. Known Limitations & Open Issues

1. **CFA not cross-validated against lavaan** — estimation works, but SEs from numerical Hessian not yet compared
2. **Oblimin/quartimin cross-validation** — uses same GPFoblq engine as geomin, should match, but not formally tested with 100 datasets
3. **No vitest unit tests for FA** — only R cross-validation scripts in validation/
4. **LTA cross-validation** — no deterministic R equivalent (depmixS4 is non-deterministic); uses hand-computed examples
5. **Matrix SVD** — Jacobi-based, occasionally unreliable for small k×k matrices. Eigendecomposition-based polar decomposition used as workaround in varimax.
6. **No IRT, SEM (beyond CFA), mediation analysis, or Bayesian methods** yet
7. **Visualization is async** (D3 dynamic import) — `import('d3').then(...)` pattern means plots render on next microtask

---

## 15. Commit & Workflow Rules

1. **Never push to remote without explicit user permission**
2. **Never delete `validation/` folder or its contents** — it's permanent infrastructure
3. **Build before commit**: `npx tsup && npx tsc --noEmit && npx vitest run`
4. **After any change to stats/**: run validation if applicable
5. **After multi-file tasks**: update CHANGES.md, HANDOFF.md, LEARNINGS.md
6. **Read CARM-PROMPT.md** before writing any code for Carm — it contains the exact coding rules
7. **Read HANDOFF.md** at the start of every session — it has current state

---

## 16. Quick Reference Tables

### 16.1 All Public Functions

```typescript
// Descriptive
describe(x, ciLevel?) → DescriptiveResult
trimmedMean(x, alpha?) → number
skewness(x) → number
kurtosis(x) → number
shapiroWilk(x) → {statistic, pValue}

// Comparison
tTestIndependent(x1, x2, equalVar?, ci?, alt?) → StatResult
tTestPaired(x1, x2, ci?) → StatResult
oneWayANOVA(groups) → ANOVAResult
mannWhitneyU(x1, x2, alt?) → StatResult
wilcoxonSignedRank(x1, x2) → StatResult
kruskalWallis(groups) → StatResult
friedmanTest(data) → StatResult

// Correlation
pearsonCorrelation(x, y, ci?, alt?) → StatResult
spearmanCorrelation(x, y, ci?) → StatResult
kendallTau(x, y) → StatResult
partialCorrelation(x, y, controls) → StatResult
correlationMatrix(data, method?) → CorrelationMatrix

// Regression
linearRegression(x, y, ci?) → RegressionResult
multipleRegression(y, predictors[], ci?) → RegressionResult
polynomialRegression(x, y, degree, ci?) → RegressionResult
logisticRegression(y, predictors[], ci?) → LogisticResult

// Clustering
fitGMM(data, options?) → GMMResult
predictGMM(model, newData) → posteriors
findBestGMM(data, kRange?, models?) → GMMResult
fitGMMRange(data, kRange, model?) → GMMRangeEntry[]
fitLCA(data, options?) → LCAResult
fitLTA(data, options?) → LTAResult
runKMeans(data, options?) → KMeansResult
predictKMeans(model, newData) → labels
fitKMeansRange(data, kRange) → KMeansRangeEntry[]
runDBSCAN(data, options?) → DBSCANResult
runHierarchical(data, options?) → HACResult
euclideanDistMatrix(data) → number[][]
silhouetteScores(data, labels) → number[]
kDistancePlot(data, k) → number[]
cutTree(merges, heights, n, k) → number[]

// Factor Analysis
runEFA(data, options?) → FAResult
runCFA(model, data, options?) → CFAResult
runFADiagnostics(data, options?) → FADiagnostics

// PCA
runPCA(data, nComponents?, options?) → PCAResult

// Mixed Models
runLMM(input) → LMMResult

// Post-hoc
tukeyHSD(groups, msWithin, dfWithin, ci?, method?) → PairwiseResult[]
dunnTest(groups, method?, ci?) → PairwiseResult[]
gamesHowellTest(groups, ci?) → PairwiseResult[]

// Effect Size
cohensD(x1, x2) → EffectSize
cohensDPaired(diffs) → EffectSize
hedgesG(x1, x2) → EffectSize
etaSquared(ssBetween, ssTotal) → EffectSize
omegaSquared(ssBetween, ssTotal, dfBetween, msWithin) → EffectSize
rankBiserial(U, n1, n2) → EffectSize
cramerV(chiSq, n, minDim) → EffectSize

// Frequency
frequencyTable(data) → FrequencyRow[]
contingencyTable(x, y) → ContingencyTable
chiSquareTest(observed, expected?) → FrequencyTestResult
fisherExactTest(a, b, c, d) → StatResult
goodnessOfFit(observed, probabilities) → FrequencyTestResult

// Preprocessing
preprocessData(data, options?) → PreprocessResult
inverseTransform(data, params) → number[][]

// Auto-dispatch
analyze(outcome, predictor?, options?) → AnalysisResult

// Core
normalCDF(z), normalQuantile(p), tDistCDF(t, df), tDistQuantile(p, df),
tDistPValue(t, df), fDistCDF(f, df1, df2), fDistPValue(f, df1, df2),
chiSqCDF(x, df), chiSqPValue(x, df), chiSqQuantile(p, df),
logGamma(z), gamma(z), incompleteBeta(x, a, b), incompleteGamma(a, x),
mean(x), median(x), variance(x), sd(x), se(x), quantile(x, p),
rank(x), sortAsc(x), cov(x, y), roundTo(x, dec),
adjustPValues(pvals, method), nelderMead(fn, x0, opts),
Matrix class (see §6.2)
```

### 16.2 Key Files to Read When Starting

| Priority | File | Why |
|----------|------|-----|
| 1 | `CARM-PROMPT.md` | Exact coding rules for writing Carm code |
| 2 | `HANDOFF.md` | Current project state |
| 3 | `LEARNINGS.md` | Accumulated engineering insights and gotchas |
| 4 | `src/core/types.ts` | The interface contract |
| 5 | `validation/VALIDATION-STRATEGY.md` | Testing standards |

### 16.3 Build Commands

| Command | Purpose |
|---------|---------|
| `npx tsup` | Build ESM + CJS + .d.ts to dist/ |
| `npx tsc --noEmit` | Type check (must be 0 errors) |
| `npx vitest run` | Run all 496 tests |
| `npx vitest run tests/stats/clustering.test.ts` | Run specific test file |
| `npx tsx validation/ts-harness/fa-full-report.ts` | FA cross-validation report |
| `npx tsx tmp/some-script.ts` | Run scratch scripts |

---

*This document was generated on 2026-02-26 from a complete exploration of the Carm codebase at commit `dc56b9d` on branch `main`. Line counts, test counts, and cross-validation results reflect the state of the code at that moment.*
