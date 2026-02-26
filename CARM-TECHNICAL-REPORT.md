# Carm Technical Report — Full Library Architecture

> **Library**: Carm — Statistical analysis and visualization for TypeScript/JavaScript
> **Source**: `src/` — 67 TypeScript files, 18,897 lines
> **Tests**: `tests/` — 17 test files, 4,456 lines, 496 passing tests
> **Build**: tsup (ESM + CJS dual output) + TypeScript 5.7 (maximum strictness)
> **Runtime dependencies**: Zero (D3 is a peer dependency, loaded dynamically)
> **Commit**: `d8cb3bb` on `main`

---

## Table of Contents

1. [Design Thesis](#1-design-thesis)
2. [Three-Layer Architecture](#2-three-layer-architecture)
3. [Module Map and Dependency Graph](#3-module-map-and-dependency-graph)
4. [Package Design and Distribution](#4-package-design-and-distribution)
5. [TypeScript Strictness Regime](#5-typescript-strictness-regime)
6. [Core Layer: Mathematics from Scratch](#6-core-layer-mathematics-from-scratch)
7. [Core Layer: Matrix Algebra](#7-core-layer-matrix-algebra)
8. [Core Layer: Type Contract](#8-core-layer-type-contract)
9. [Core Layer: APA Formatting](#9-core-layer-apa-formatting)
10. [Stats Layer: Descriptive Statistics](#10-stats-layer-descriptive-statistics)
11. [Stats Layer: Group Comparisons](#11-stats-layer-group-comparisons)
12. [Stats Layer: Correlations](#12-stats-layer-correlations)
13. [Stats Layer: Regression](#13-stats-layer-regression)
14. [Stats Layer: Frequency Analysis](#14-stats-layer-frequency-analysis)
15. [Stats Layer: Effect Sizes and Post-Hoc Tests](#15-stats-layer-effect-sizes-and-post-hoc-tests)
16. [Stats Layer: PCA](#16-stats-layer-pca)
17. [Stats Layer: Clustering (GMM, LCA, LTA, K-Means, DBSCAN, HAC)](#17-stats-layer-clustering)
18. [Stats Layer: Factor Analysis (EFA, CFA)](#18-stats-layer-factor-analysis)
19. [Stats Layer: Linear Mixed Models](#19-stats-layer-linear-mixed-models)
20. [Stats Layer: Auto-Dispatch (analyze)](#20-stats-layer-auto-dispatch)
21. [Stats Layer: Preprocessing](#21-stats-layer-preprocessing)
22. [Viz Layer: Architecture and Theme System](#22-viz-layer-architecture-and-theme-system)
23. [Viz Layer: Shared Components](#23-viz-layer-shared-components)
24. [Viz Layer: 42 Plot Types](#24-viz-layer-42-plot-types)
25. [Viz Layer: FA Path Diagram (Flagship Visualization)](#25-viz-layer-fa-path-diagram)
26. [Viz Layer: Export Pipeline](#26-viz-layer-export-pipeline)
27. [Deterministic Reproducibility (splitmix32 PRNG)](#27-deterministic-reproducibility)
28. [Cross-Validation Infrastructure](#28-cross-validation-infrastructure)
29. [Public API Reference](#29-public-api-reference)
30. [Engineering Decisions: Problems, Solutions, and Optimizations](#30-engineering-decisions)
31. [Mathematical Tricks That Made It Possible](#31-mathematical-tricks)

---

## 1. Design Thesis

Carm is an enterprise-grade statistical library that runs entirely in the browser with **zero runtime dependencies**. Its core thesis:

> **Every statistical method implemented from first principles — no jStat, no simple-statistics, no BLAS/LAPACK. Every distribution function, matrix decomposition, and special function built from scratch in TypeScript.**

This is not a wrapper around existing libraries. Every line of mathematical code was written, validated against R, and optimized for browser execution. The motivation:

1. **Zero dependency risk** — no supply chain, no version conflicts, no abandoned packages
2. **Auditable numerics** — every algorithm readable in a single codebase
3. **Browser-first** — no Node.js native modules, no WASM compilation step
4. **Deterministic** — seeded PRNG ensures identical results across runs

---

## 2. Three-Layer Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     Carm (18,897 LoC)                           │
│                                                                 │
│  ┌──────────────── viz/ (10,433 LoC) ─────────────────────┐    │
│  │  42 D3 plot renderers + 4 components + themes + export  │    │
│  │  import('d3').then(d3 => ...)   [peer dependency]       │    │
│  └─────────────────────┬───────────────────────────────────┘    │
│                        │ consumes                               │
│  ┌─────────────── stats/ (7,481 LoC) ────────────────────┐     │
│  │  13 modules: descriptive, comparison, correlation,     │     │
│  │  regression, frequency, effect-size, post-hoc, pca,    │     │
│  │  mixed, clustering, factor-analysis, analyze, preprocess│     │
│  └─────────────────────┬───────────────────────────────────┘    │
│                        │ consumes                               │
│  ┌──────────────── core/ (1,489 LoC) ────────────────────┐     │
│  │  types.ts (303L) — all interfaces                      │     │
│  │  math.ts  (547L) — distributions, special functions    │     │
│  │  matrix.ts(414L) — Matrix class, decompositions        │     │
│  │  apa.ts   (221L) — APA 7 formatting                    │     │
│  └────────────────────────────────────────────────────────┘     │
└─────────────────────────────────────────────────────────────────┘
```

**Layer discipline is enforced by convention:**
- `core/` imports nothing from `stats/` or `viz/`
- `stats/` imports only from `core/`
- `viz/` imports from `core/` and `stats/` (for type definitions and `linearRegression` for scatter+fit)
- No circular dependencies exist

---

## 3. Module Map and Dependency Graph

### Core (4 files, 1,489 lines)

| File | Lines | Purpose | Dependencies |
|------|-------|---------|-------------|
| `types.ts` | 303 | All shared interfaces | None |
| `math.ts` | 547 | Distributions, special functions, optimization | None |
| `matrix.ts` | 414 | Matrix class with decompositions | None |
| `apa.ts` | 221 | APA 7 string formatting | `math.ts` (roundTo) |

### Stats (13 files, 7,481 lines)

| File | Lines | Purpose | Core deps | Stats deps |
|------|-------|---------|-----------|-----------|
| `descriptive.ts` | 287 | Descriptive stats, Shapiro-Wilk | math, types | — |
| `comparison.ts` | 476 | t-tests, ANOVA, Mann-Whitney, Wilcoxon, Kruskal-Wallis, Friedman | math, types, apa | — |
| `effect-size.ts` | 188 | Cohen's d, Hedges' g, η², ω², rank-biserial | math, types, apa | — |
| `frequency.ts` | 287 | Chi-square, Fisher's exact, contingency tables | math, types, apa | — |
| `post-hoc.ts` | 221 | Tukey HSD, Games-Howell, Dunn's test | math, types | — |
| `correlation.ts` | 282 | Pearson, Spearman, Kendall, partial, matrix | math, types, apa | regression |
| `regression.ts` | 379 | OLS, logistic, polynomial, diagnostics | math, matrix, types, apa | — |
| `pca.ts` | 156 | PCA via SVD, varimax rotation | matrix, types | — |
| `mixed.ts` | 297 | Linear mixed models (REML) | math, matrix, types, apa | — |
| `clustering.ts` | 1,850 | GMM, LCA, LTA, K-Means, DBSCAN, HAC | math, matrix, types | — |
| `factor-analysis.ts` | 1,923 | EFA (PAF, ML), CFA, rotations, diagnostics | math, matrix, types, apa | — |
| `analyze.ts` | 416 | Auto-routing dispatcher | types | descriptive, comparison, frequency, effect-size, post-hoc |
| `preprocess.ts` | 181 | Center, standardize, log, sqrt | types | — |

### Viz (50 files, 10,433 lines)

| Category | Files | Lines | Contents |
|----------|-------|-------|----------|
| Themes | 1 | 116 | Light/dark theme, CSS properties, color scales |
| Components | 4 | 462 | Annotations, brackets, tooltip, axis helpers |
| Plots | 42 | 9,785 | All visualization types |
| Export | 1 | 70 | SVG + PNG export |

---

## 4. Package Design and Distribution

### Dual ESM/CJS Output

```json
{
  "exports": {
    ".":      { "import": "./dist/index.js",      "require": "./dist/index.cjs" },
    "./stats": { "import": "./dist/stats/index.js", "require": "./dist/stats/index.cjs" },
    "./viz":   { "import": "./dist/viz/index.js",   "require": "./dist/viz/index.cjs" },
    "./core":  { "import": "./dist/core/index.js",  "require": "./dist/core/index.cjs" }
  }
}
```

Three import paths:
- `import { pearsonCorrelation, renderScatterStats } from 'carm'` — everything
- `import { pearsonCorrelation } from 'carm/stats'` — stats only (no D3 dependency)
- `import { renderScatterStats } from 'carm/viz'` — viz only

### D3 as Peer Dependency

D3 is declared as `peerDependencies`, not `dependencies`. Every viz function loads D3 via dynamic `import('d3').then(...)`. This means:

1. Importing `carm/stats` adds zero bytes of D3 to your bundle
2. D3 is only loaded when a plot is actually rendered
3. The consumer controls the D3 version (>=7.0.0)

---

## 5. TypeScript Strictness Regime

Carm uses the maximum possible TypeScript strictness:

```json
{
  "strict": true,
  "noUncheckedIndexedAccess": true,
  "noUnusedLocals": true,
  "noUnusedParameters": true,
  "exactOptionalPropertyTypes": true
}
```

**`noUncheckedIndexedAccess`**: Every array access `arr[i]` returns `T | undefined`. Code must use `arr[i]!` after validation or handle the `undefined` case. This prevents out-of-bounds access bugs.

**`exactOptionalPropertyTypes`**: Optional properties like `{ x?: number }` cannot be assigned `undefined` directly — must use conditional spread `...(x !== undefined && { x })`. This catches subtle bugs where `undefined` is accidentally passed as a value.

All function return types are `readonly` — results are immutable:
```typescript
readonly loadings: readonly (readonly number[])[]
readonly eigenvalues: readonly number[]
```

---

## 6. Core Layer: Mathematics from Scratch

**File**: `src/core/math.ts` (547 lines)

Every distribution function and special function in Carm is implemented from scratch. No external math library is used.

### Special Functions

| Function | Algorithm | Reference |
|----------|-----------|-----------|
| `logGamma(z)` | Lanczos approximation (g=7, 9 coefficients) | Numerical Recipes §6.1 |
| `gamma(z)` | `exp(logGamma(z))` | |
| `betaFn(a,b)` | `exp(logGamma(a) + logGamma(b) - logGamma(a+b))` | |
| `incompleteBeta(x,a,b)` | Lentz continued fraction with symmetry optimization | NR §6.4 |
| `incompleteGamma(a,x)` | Series for x<a+1, continued fraction complement for x≥a+1 | NR §6.2 |
| `erf(x)` | A&S 7.1.26 polynomial (5 terms, ~1.5e-7 accuracy) | Abramowitz & Stegun |

### Distribution Functions

| Distribution | CDF | Quantile | P-value |
|-------------|-----|----------|---------|
| Normal | `erf`-based | Peter Acklam's rational approximation (~1.15e-9) | 2-tailed via CDF |
| t | `incompleteBeta`-based | Bisection (100 iterations) | 2-tailed |
| F | `incompleteBeta`-based | — | 1-tailed |
| χ² | `incompleteGamma`-based | Bisection (100 iterations) | 1-tailed |

### Optimization

| Function | Algorithm | Parameters |
|----------|-----------|-----------|
| `nelderMead(fn, x0, opts)` | Standard simplex: reflect/expand/contract/shrink | α=1, γ=2, ρ=0.5, σ=0.5; 5% perturbation initialization |

### P-Value Adjustment

| Method | Algorithm |
|--------|-----------|
| Bonferroni | `min(p × m, 1)` |
| Holm | Step-down: `min(p × (m-rank+1), 1)` with running max |
| BH | Step-up: `min(p × m/rank, 1)` with running min (descending) |
| BY | Like BH with harmonic correction `c(m) = Σ 1/i` |

### Descriptive Utilities

`mean()`, `variance()`, `sd()`, `se()`, `median()`, `quantile()` (R type=7), `rank()` (average-tie), `cov()`, `clamp()`, `roundTo()`

---

## 7. Core Layer: Matrix Algebra

**File**: `src/core/matrix.ts` (414 lines)

### Matrix Class

Row-major flat `Float64Array` storage. Immutable — all operations return new matrices.

| Operation | Algorithm | Complexity | Notes |
|-----------|-----------|-----------|-------|
| `multiply(B)` | i-k-j loop order (cache-friendly) | O(n²m) | Sequential access for B's columns |
| `cholesky()` | Cholesky-Banachiewicz | O(n³/3) | Throws if not SPD |
| `inverse()` | Gauss-Jordan with partial pivoting | O(n³) | Singularity check at 1e-12 |
| `eigen()` | Classical Jacobi (largest off-diagonal pivot) | O(100·n²) max | Convergence at 1e-12 |
| `svd()` | One-sided Jacobi | O(200·n²) max | Convergence at 1e-15·√(αβ) |
| `pseudoInverse()` | SVD-based with threshold truncation | O(n³) | tol = 1e-10 × max(σ) |
| `logDet()` | Cholesky: `2·Σ log(L_ii)` | O(n²) | Avoids overflow of `det()` |

### Static Constructors

`Matrix.fromArray(data)`, `Matrix.identity(n)`, `Matrix.zeros(r,c)`, `Matrix.colVec(v)`, `Matrix.rowVec(v)`

### Arithmetic

`add(B)`, `subtract(B)`, `scale(c)`, `transpose()`, `trace()`, `diagonal()`

---

## 8. Core Layer: Type Contract

**File**: `src/core/types.ts` (303 lines)

All shared interfaces live here. No file in `types.ts` imports from any other source file.

### Universal Result Types

| Interface | Fields | Used by |
|-----------|--------|---------|
| `StatResult` | testName, statistic, df, pValue, effectSize, ci, ciLevel, n, formatted | All tests |
| `DescriptiveResult` | n, mean, median, mode, trimmedMean, sd, se, variance, min, max, range, iqr, q1, q3, skewness, kurtosis, ci, shapiroWilk, formatted | `describe()` |
| `RegressionResult` | coefficients[], r2, adjR2, fStatistic, fDf, fPValue, aic, bic, residuals, fitted, formatted | All regression |
| `PairwiseResult` | group1, group2, meanDiff, se, statistic, pValue, pValueAdj, ci, significant | Post-hoc |
| `EffectSize` | value, name, interpretation | All effect sizes |
| `AnalysisResult` | test, outcome, predictor, result, descriptives, posthoc, normality | `analyze()` |

### Factor Analysis Types

| Interface | Fields |
|-----------|--------|
| `FAResult` | loadings[][], uniqueness[], communalities[], factorCorrelations[][], fit, eigenvalues[], nFactors, rotation, extraction, variableNames[], formatted |
| `CFAResult` | extends FAResult + parameterEstimates (loadings, uniquenesses, factorCovariances with SE/z/p) + model |
| `FADiagnostics` | kmo, kmoPerItem[], bartlett, mapSuggested, parallelEigenvalues[], parallelSimulated[], parallelSuggested |
| `FactorFit` | chiSq, df, pValue, rmsea, rmseaCI, cfi, tli, srmr, aic, bic |

### Clustering Types

| Interface | Fields |
|-----------|--------|
| `GMMResult` | labels[], posteriors[][], means[][], covariances[][][], weights[], diagnostics |
| `KMeansResult` | labels[], centroids[][], inertia, diagnostics |
| `LCAResult` | labels[], posteriors[][], classProbabilities[], itemProbabilities[][], diagnostics |
| `LTAResult` | labels[], posteriors[][], initialProbs[], transitionMatrix[][], itemProbabilities[][], diagnostics |
| `ClusterDiagnostics` | converged, iterations, logLikelihood, df, aic, bic, icl, entropy, avepp, formatted |

### Data Types

| Type | Definition |
|------|-----------|
| `FieldType` | `'numeric' \| 'binary' \| 'categorical' \| 'ordinal'` |
| `Field` | `NumericField \| GroupField` |
| `PAdjMethod` | `'bonferroni' \| 'holm' \| 'BH' \| 'BY' \| 'none'` |
| `CovarianceModel` | `'VVV' \| 'EEE' \| 'VVI' \| 'EEI' \| 'VII' \| 'EII'` |
| `DataMatrix` | `readonly (readonly number[])[]` |

---

## 9. Core Layer: APA Formatting

**File**: `src/core/apa.ts` (221 lines)

Every statistical result includes a `formatted` string in APA 7 style. The formatting functions are centralized here.

### Formatters

| Function | Output example |
|----------|---------------|
| `formatP(p)` | `p < .001` or `p = .025` (no leading zero) |
| `formatStat(v, decimals)` | `2.31` (2dp default) |
| `formatCI(ci, decimals)` | `[0.08, 1.22]` |
| `formatDF(df)` | `48` (integer) or `23.41` (Welch) or `2, 87` (F-test) |
| `formatTTest(...)` | `t(48) = 2.31, p = .025, d = 0.65, 95% CI [0.08, 1.22]` |
| `formatANOVA(...)` | `F(2, 87) = 4.21, p = .018, η² = 0.09` |
| `formatChiSq(...)` | `χ²(3) = 8.46, p = .037, V = 0.29` |
| `formatCorrelation(...)` | `r(48) = 0.54, p = .003, 95% CI [0.28, 0.73]` |
| `formatRegression(...)` | `R² = 0.42, adj. R² = 0.40, F(3, 96) = 23.2, p < .001` |
| `formatMannWhitney(...)` | `W = 423, p = .031, r = 0.27` |
| `formatKruskalWallis(...)` | `H(2) = 12.3, p = .002, η²_H = 0.18` |
| `formatLMM(...)` | `ICC = 0.42, AIC = 1234.5, BIC = 1256.3, logLik = -612.3` |
| `formatCFAFit(fit)` | `χ²(24) = 28.42, p = .241, RMSEA = 0.042 [0.000, 0.085], CFI = 0.987, TLI = 0.982, SRMR = 0.038` |

### Effect Size Interpreters

| Function | Thresholds |
|----------|-----------|
| `interpretCohensD(d)` | negligible < 0.2 < small < 0.5 < medium < 0.8 < large |
| `interpretEtaSq(η²)` | negligible < 0.01 < small < 0.06 < medium < 0.14 < large |
| `interpretR(r)` | negligible < 0.1 < small < 0.3 < medium < 0.5 < large < 0.7 < very large |
| `interpretCramerV(V, df)` | df-adjusted thresholds (Cohen 1988) |
| `interpretEffect(value, [s,m,l])` | Generic with custom thresholds |

---

## 10. Stats Layer: Descriptive Statistics

**File**: `src/stats/descriptive.ts` (287 lines)

### `describe(x, ciLevel?)` → `DescriptiveResult`

Computes 17 descriptive statistics in one pass:
- Central tendency: mean, median, mode, 5% trimmed mean
- Dispersion: sd, se, variance, min, max, range, IQR, Q1, Q3
- Shape: skewness (adjusted Fisher-Pearson), excess kurtosis
- CI: t-based confidence interval for the mean
- Normality: embedded Shapiro-Wilk test

### Shapiro-Wilk Test (AS R94, Royston 1995)

The most complex function in descriptive.ts. Implementation follows the AS R94 algorithm exactly:

1. **Coefficient computation** (3 paths): n=3 (exact), n=4-5 (polynomial), n≥6 (normal quantile + polynomial correction)
2. **Polynomial corrections** for outermost coefficients (SW_C1, SW_C2 arrays)
3. **W statistic**: W = (Σ aᵢ x_(i))² / Σ (xᵢ - x̄)²
4. **P-value**: Royston approximation branching on n≤11 (exact polynomial) vs n>11 (transformed normal)

Valid for n = 3 to 5000.

---

## 11. Stats Layer: Group Comparisons

**File**: `src/stats/comparison.ts` (476 lines)

### Parametric Tests

| Function | Test | df | Effect Size |
|----------|------|----|----|
| `tTestIndependent(x1, x2, equalVar?, ci?)` | Welch's (default) or Student's | Welch-Satterthwaite | Cohen's d |
| `tTestPaired(x1, x2, ci?)` | Paired t-test | n-1 | Cohen's d (paired) |
| `oneWayANOVA(groups)` | F-test | [k-1, N-k] | η² (eta-squared) |

### Non-Parametric Tests

| Function | Test | Distribution | Effect Size |
|----------|------|-------------|-------------|
| `mannWhitneyU(x1, x2, alt?)` | Mann-Whitney U | Normal approx with tie correction | Rank-biserial r |
| `wilcoxonSignedRank(x1, x2)` | Wilcoxon signed-rank | Exact DP (n≤20) or normal approx | Rank-biserial r |
| `kruskalWallis(groups)` | Kruskal-Wallis H | χ² with tie correction | η²_H |
| `friedmanTest(data)` | Friedman χ² | χ² distribution | Kendall's W |

### Key Implementation Details

**Wilcoxon exact distribution** (n ≤ 20): Dynamic programming builds the full null distribution. `dist[w]` = number of sign assignments of ranks 1..n yielding W⁺ = w. Space O(n·maxW), time O(n·maxW).

**Mann-Whitney tie correction**: Variance adjusted by `1 - Σ(t³-t)/(n³-n)` where t = tie group sizes.

**Kruskal-Wallis**: `H = 12/(N(N+1)) · Σ(R²ⱼ/nⱼ) - 3(N+1)`, divided by tie correction `C = 1 - Σ(t³-t)/(N³-N)`.

---

## 12. Stats Layer: Correlations

**File**: `src/stats/correlation.ts` (282 lines)

| Function | Method | CI method | p-value |
|----------|--------|-----------|---------|
| `pearsonCorrelation(x, y, ci?)` | Standard formula | Fisher z-transform | t-test, df=n-2 |
| `spearmanCorrelation(x, y, ci?)` | Pearson on ranks | Fisher z-transform | t-test on rₛ |
| `kendallTau(x, y, ci?)` | O(n²) concordance counting | Normal approx | Normal approx |
| `partialCorrelation(x, y, controls)` | OLS residualization | — | t-test |
| `correlationMatrix(data, labels?, method?)` | Pairwise computation | Per-pair CI | Per-pair p |

### Fisher z-Transform CI

```
z = 0.5 · ln((1+r)/(1-r))
SE_z = 1/√(n-3)
CI_z = z ± z_α/2 · SE_z
CI_r = tanh(CI_z)     ← guarantees [-1, 1]
```

### Kendall τ-b

Handles ties via: `τ = (concordant - discordant) / √((n₀ - tiesX)(n₀ - tiesY))` where `n₀ = n(n-1)/2`.

---

## 13. Stats Layer: Regression

**File**: `src/stats/regression.ts` (379 lines)

### OLS Engine

All regression routes through a single OLS engine: `β = (X'X)⁻¹X'y`

| Function | Design matrix |
|----------|--------------|
| `linearRegression(x, y)` | `[1, x]` |
| `multipleRegression(y, predictors)` | `[1, x₁, x₂, ...]` |
| `polynomialRegression(x, y, degree)` | `[1, x, x², ..., xᵈ]` |

### Logistic Regression (IRLS)

`logisticRegression(y, predictors, ci?, maxIter?, tol?)` — Binary outcomes via Fisher scoring:

1. Initialize β = 0
2. Compute μ = logistic(Xβ)
3. Weights W = diag(μ(1-μ)), clamped to min 1e-10
4. Working response z = Xβ + (y-μ)/(μ(1-μ))
5. Update β via weighted LS: `(X'WX)⁻¹X'Wz`
6. Converge when max|Δβ| < tol

McFadden pseudo-R² = 1 - logLik/logLik_null. P-values from Wald z-test.

### Regression Diagnostics

`regressionDiagnostics(result, predictors)`:
- **Leverage**: `h_ii = [X(X'X)⁻¹X']_ii`
- **Standardized residuals**: `e_i / √(MSE(1-h_i))`
- **Cook's D**: `r²_i · h_i / (p · MSE · (1-h_i)²)`
- **VIF**: `1/(1-R²_j)` from regressing each predictor on all others

---

## 14. Stats Layer: Frequency Analysis

**File**: `src/stats/frequency.ts` (287 lines)

| Function | Test | Effect Size |
|----------|------|-------------|
| `chiSquareTest(observed, yates?)` | χ² independence | Cramér's V |
| `fisherExactTest(a, b, c, d)` | Fisher's exact (2×2) | Odds ratio + CI |
| `goodnessOfFit(observed, expected?)` | χ² goodness-of-fit | Cohen's w |
| `frequencyTable(data)` | Frequency table | — |
| `phiCoefficient(a, b, c, d)` | φ coefficient | — |

### Fisher's Exact Test

Log-space hypergeometric PMF: `logP = logC(K,k) + logC(N-K,n-k) - logC(N,n)` via `logFactorial`. Two-tailed p = sum of all tables with P ≤ P_observed.

Haldane-Anscombe correction: adds 0.5 to all cells when any cell is 0 (for OR computation).

### Cramér's V

`V = √(χ²/(n·(min(R,C)-1)))` with df-adjusted interpretation thresholds.

---

## 15. Stats Layer: Effect Sizes and Post-Hoc Tests

### Effect Sizes (`stats/effect-size.ts`, 188 lines)

| Function | Formula | CI |
|----------|---------|-----|
| `cohensD(x1, x2)` | d = (M₁-M₂) / SD_pooled | Normal approx (Hedges & Olkin) |
| `hedgesG(x1, x2)` | g = d × J, J = 1 - 3/(4df-1) | Same with bias correction |
| `etaSquared(SSb, SSt)` | η² = SSb/SSt | — |
| `omegaSquared(SSb, SSt, dfb, MSw)` | ω² = (SSb-dfb·MSw)/(SSt+MSw) | max(0, ...) guard |
| `rankBiserial(U, n1, n2)` | r = 1 - 2U/(n₁n₂) | — |
| `etaSquaredKW(H, k, n)` | η²_H = (H-k+1)/(n-k) | — |

### Post-Hoc Tests (`stats/post-hoc.ts`, 221 lines)

| Function | Assumptions | p-adjustment |
|----------|-------------|--------------|
| `tukeyHSD(groups, MSw, dfw)` | Equal variances | Bonferroni-approximated studentized range |
| `gamesHowell(groups)` | Unequal variances | Per-pair Welch df + Bonferroni |
| `dunnTest(groups, method?)` | Non-parametric (ranks) | User-specified (Bonferroni/Holm/BH/BY) |

**Why Bonferroni approximation for Tukey**: The true studentized range distribution requires `qtukey()` — a special function not available without tables or complex computation. The Bonferroni-adjusted t approximation is conservative but adequate and avoids implementing another special function.

---

## 16. Stats Layer: PCA

**File**: `src/stats/pca.ts` (156 lines)

### `runPCA(data, nComponents?, scale?)`

1. Standardize (default) or center data
2. Scale by `1/√(n-1)` — ensures singular values² = eigenvalues of covariance matrix
3. SVD of scaled data → singular values, left/right singular vectors
4. Eigenvalues = σ²; loadings = right singular vectors; scores = data × loadings

### `varimaxRotation(loadings, maxIter?, tol?)`

Kaiser (1958) pairwise algorithm:
1. Kaiser normalize (divide rows by communality)
2. For each pair (p,q): compute optimal rotation angle via `atan2(Y,X)/4`
3. Apply Givens rotation
4. Denormalize
5. Converge when total angle change < tol

---

## 17. Stats Layer: Clustering

**File**: `src/stats/clustering.ts` (1,850 lines)

The largest stats module. Implements 6 clustering algorithms.

### Gaussian Mixture Model (GMM)

`fitGMM(data, options)` — Full EM with mclust-style covariance constraints.

**6 covariance models** (matching R's mclust):

| Model | Shape | Volume | Orientation | Parameters per component |
|-------|-------|--------|-------------|------------------------|
| VVV | Variable | Variable | Variable | d(d+1)/2 (full) |
| EEE | Equal | Equal | Equal | d(d+1)/2 (shared) |
| VVI | Variable | Variable | Axis-aligned | d (diagonal) |
| EEI | Equal | Equal | Axis-aligned | d (shared diagonal) |
| VII | Variable | Spherical | — | 1 (σ²) |
| EII | Equal | Spherical | — | 1 (shared σ²) |

**EM algorithm**:
- E-step: log-space responsibilities via logSumExp trick
- M-step: weighted updates of means, covariances, weights
- Convergence: |ΔlogLik| < 1e-6
- Regularization: regCovar = 1e-6 on diagonal

**K-Means++ initialization**: First centroid random, subsequent proportional to D² from nearest existing centroid.

### K-Means

`fitKMeans(data, options)` — Lloyd's algorithm with K-Means++ init and empty-cluster re-seeding.

### DBSCAN

`runDBSCAN(data, eps, minPts)` — Density-based clustering. Noise label = -1.

### Hierarchical Agglomerative Clustering (HAC)

`runHierarchical(data, options)` — 5 linkage methods: single, complete, average, Ward's, centroid.

`cutTree(merges, heights, n, k)` — Cut dendrogram at k clusters.

`silhouetteScores(data, labels)` — Per-point silhouette width for cluster validation.

### Latent Class Analysis (LCA)

`fitLCA(data, options)` — Bernoulli mixture model for binary data. Beta(1,1) prior smoothing.

### Latent Transition Analysis (LTA)

`fitLTA(data, options)` — Hidden Markov model for sequences of binary items. Full Baum-Welch in log-space (forward-backward algorithm).

### Model Comparison

`fitGMMRange(data, kRange?, models?)` — Grid search over K × covariance model, returns best by BIC.
`fitKMeansRange(data, kRange?)` — Elbow method via inertia.

### Cluster Diagnostics

Every clustering result includes `ClusterDiagnostics`:
- `entropy`: Normalized classification entropy (0 = perfect, 1 = random)
- `avepp`: Average posterior probability per assigned class
- `aic`, `bic`, `icl`: Information criteria (ICL = BIC + 2·entropy penalty)
- `df`: Model degrees of freedom (varies by covariance model)

---

## 18. Stats Layer: Factor Analysis

**File**: `src/stats/factor-analysis.ts` (1,923 lines)

The largest and most algorithmically complex file in Carm.

### EFA: `runEFA(data, options)`

**Two extraction methods**:

| Method | Algorithm | Lines | Matches R |
|--------|-----------|-------|-----------|
| PAF | Iterative eigendecomposition of reduced R | ~100 | `psych::fa(fm='pa')` |
| ML | Two-phase: Jöreskog gradient + Nelder-Mead polish | ~150 | `stats::factanal()` |

**Six rotation methods**:

| Rotation | Type | Algorithm | Random starts |
|----------|------|-----------|--------------|
| varimax | Orthogonal | Kaiser pairwise angle | No (deterministic) |
| promax | Oblique | Varimax → target → Procrustes | No (deterministic) |
| oblimin | Oblique | GPFoblq gradient descent (γ=0) | 50 (Haar random) |
| quartimin | Oblique | GPFoblq gradient descent (γ=0, quartimin criterion) | 50 (Haar random) |
| geomin | Oblique | GPFoblq with geometric mean criterion | 50 (Haar random) |
| none | — | No rotation | — |

**Random starts**: Start 0 = identity (deterministic reference). Starts 1-49 = Haar-distributed random orthogonal matrices via QR of random Gaussian matrix with sign correction.

### CFA: `runCFA(model, options)`

User specifies factor structure as `Record<string, number[]>` (factor name → item indices).

1. Initialize: loadings = 0.7 on specified paths, uniquenesses = 0.5
2. Gradient descent on ML discrepancy: F = log|Σ| + tr(Σ⁻¹S) - log|S| - p
3. Convergence: max|gradient| < tol
4. Standard errors via numerical Hessian (central finite differences, h=1e-4)
5. Triple fallback: inverse → pseudo-inverse → default SE = 0.05

### Fit Indices

| Index | Formula | Good fit |
|-------|---------|---------|
| χ² | (n-1-correction) × F_ML | p > 0.05 |
| RMSEA | √(max(0, (χ²/df - 1)/(n-1))) | < 0.06 |
| CFI | 1 - max(0, ncp_model) / max(ncp_null, ncp_model, 1) | > 0.95 |
| TLI | ((χ²_null/df_null) - (χ²_model/df_model)) / ((χ²_null/df_null) - 1) | > 0.95 |
| SRMR | √(mean((r_obs - r_implied)²)) | < 0.08 |
| AIC | χ² + 2·df_params | Lower is better |
| BIC | χ² + log(n)·df_params | Lower is better |

### Diagnostics: `runFADiagnostics(data, options)`

| Test | Purpose | Threshold |
|------|---------|-----------|
| KMO | Sampling adequacy | > 0.6 acceptable |
| Bartlett's sphericity | Tests R = I | p < 0.05 |
| MAP | Velicer's minimum average partial | nFactors suggestion |
| Parallel analysis | Monte Carlo eigenvalue comparison | nFactors suggestion |

---

## 19. Stats Layer: Linear Mixed Models

**File**: `src/stats/mixed.ts` (297 lines)

### `runLMM(input)` → `LMMResult`

Random intercept model: `y = Xβ + Zb + ε` where `b ~ N(0, σ²_b I)`, `ε ~ N(0, σ²_e I)`.

**Profiled REML**: Optimize over single parameter `log(ψ)` where `ψ = σ²_b/σ²_e`:

1. `V_ψ⁻¹` via Woodbury identity: `I - Z(Z'Z + ψ⁻¹I)⁻¹Z'` — reduces O(n³) to O(q³)
2. `log|V_ψ|` via matrix determinant lemma: `q·log(ψ) + log|D|`
3. Profile out σ²_e analytically
4. 5-point multi-start Nelder-Mead: logψ ∈ {-4, -2, 0, 2, 4}

### `computeBLUPs(input, result)` → `{group, blup}[]`

Empirical Bayes: `b̂_j = ψ/(1 + ψ·n_j) · Σ_{i∈j} e_i`

### Variance Components

- ICC = σ²_b / (σ²_b + σ²_e)
- AIC = -2·logLik + 2·(p+2) where p = number of fixed effects
- Satterthwaite df = max(1, n - p - nGroups + 1) (simplified)

---

## 20. Stats Layer: Auto-Dispatch

**File**: `src/stats/analyze.ts` (416 lines)

### `analyze(outcome, predictor?, options?)` → `AnalysisResult`

Automatic test selection based on field types and normality:

```
outcome type    predictor type    normality    selected test
─────────────   ──────────────    ─────────    ─────────────
numeric         —                 —            describe-only
numeric         binary            all normal   t-test (paired/independent)
numeric         binary            any non-norm mann-whitney / wilcoxon
numeric         categorical (3+)  all normal   one-way-anova + post-hoc
numeric         categorical (3+)  any non-norm kruskal-wallis + dunn
binary/cat      binary/cat        —            chi-square (or Fisher if expected<5)
```

### `detectFieldType(values)` → `FieldType`

1. All finite numbers with only {0,1} → `'binary'`
2. All finite numbers → `'numeric'`
3. Exactly 2 unique values → `'binary'`
4. Otherwise → `'categorical'`

---

## 21. Stats Layer: Preprocessing

**File**: `src/stats/preprocess.ts` (181 lines)

### `preprocessData(data, options?)` → `PreprocessResult`

| Method | Transform | Inverse |
|--------|-----------|---------|
| `'none'` | Pass-through | — |
| `'center'` | `x - mean` | `x + mean` |
| `'standardize'` | `(x - mean) / sd` | `x·sd + mean` |
| `'log'` | `ln(x)` (requires x > 0) | `exp(x)` |
| `'sqrt'` | `√x` (requires x ≥ 0) | `x²` |

Zero-variance columns get SD = 1 to prevent division by zero.

---

## 22. Viz Layer: Architecture and Theme System

**Files**: `viz/themes/default.ts` (116L), `viz/index.ts` (57L)

### Dynamic D3 Import Pattern

Every single one of the 42 plot functions uses:

```typescript
export function renderPlot(container: HTMLElement, data: Data, config: Config = {}): void {
  import('d3').then(d3 => renderPlotD3(d3, container, data, config))
}
```

D3 is loaded once (cached by the module system), but the `import()` call is non-blocking — the plot renders asynchronously.

### Theme System

```typescript
interface CarmTheme {
  // Colors
  colors: readonly string[]     // 8-color palette
  background: string            // '#ffffff'
  surface: string               // '#f8f9fa'
  text: string                  // '#1a1a2e'
  textMuted: string             // '#6c757d'
  textAnnotation: string        // '#495057'
  gridLine: string              // '#eaeef3'
  axisLine: string              // '#c4cdd6'
  // Typography
  fontFamily: string            // System sans-serif stack
  fontFamilyMono: string        // Monospace stack
  fontSize: number              // 12
  fontSizeSmall: number         // 11
  fontSizeTitle: number         // 16
  // Spacing
  marginTop: number             // 58
  marginRight: number           // 32
  marginBottom: number          // 84
  marginLeft: number            // 64
  // Opacity
  pointOpacity: number          // 0.55
  violinOpacity: number         // 0.72
  ciOpacity: number             // 0.15
}
```

`applyTheme(container, theme)` injects CSS custom properties: `--js-bg`, `--js-text`, `--js-color-0` through `--js-color-7`, etc.

### Default Palette (8 colors)

Tableau 10 variant, colorblind-safe:
1. Cornflower blue `#5b8bd6`
2. Amber orange `#e8993e`
3. Soft crimson `#d65b5b`
4. Dusty teal `#6bb5a0`
5. Forest green `#5bab5e`
6. Dusty mauve `#a87dba`
7. Rose pink `#d6839b`
8. Warm umber `#b09e82`

---

## 23. Viz Layer: Shared Components

### Annotations (`components/annotations.ts`, 182 lines)

| Function | Purpose |
|----------|---------|
| `addSubtitle(svg, title, subtitle, W, theme)` | Title + APA subtitle with accent bar |
| `addCaption(svg, caption, W, H, theme)` | Bottom-left italic caption |
| `addRegressionEquation(svg, reg, x, y, theme)` | Monospace pill with ŷ = b₀ + b₁x, R² |
| `addNLabel(g, text, x, y, theme)` | Sample size label below group |
| `addStatBadge(g, lines, x, y, theme)` | Multi-line stat box |

### Brackets (`components/brackets.ts`, 129 lines)

Greedy level-assignment algorithm for non-overlapping significance brackets:
1. Sort brackets shortest first
2. Assign each to the lowest level where it doesn't overlap existing brackets
3. Render: horizontal lines at each level with vertical connectors
4. P-value label centered above each bracket

### Tooltip (`components/tooltip.ts`, 64 lines)

Floating `<div>` with smart viewport clamping. Appears on mouseover, fades on mouseout.

### Export (`viz/export.ts`, 70 lines)

- `exportSVG(container, filename)` — XML serialization with xmlns
- `exportPNG(container, filename)` — Canvas at 300 DPI (scale 300/96 ≈ 3.125×)

---

## 24. Viz Layer: 42 Plot Types

### Group 1: Core Statistical Plots (12)

| Plot | File | Lines | Purpose |
|------|------|-------|---------|
| Violin + Box | `violin-box.ts` | 306 | KDE violin + box + jitter + significance brackets |
| Scatter Stats | `scatter-stats.ts` | 159 | Scatter + regression line + CI band + equation |
| Histogram | `histogram.ts` | 128 | Bins + optional density + normal curve |
| Bar Stats | `bar-stats.ts` | 120 | Categorical bars with count/% labels |
| Correlogram | `correlogram.ts` | 155 | Heatmap with RdBu diverging color + significance stars |
| Coef Plot | `coef-plot.ts` | 118 | Coefficient dot-and-whisker with CI |
| Q-Q Plot | `qq-plot.ts` | 122 | Quantile-quantile with KS confidence band |
| Residual Panel | `residual-panel.ts` | 144 | 2×2 diagnostic grid |
| Raincloud | `raincloud.ts` | 179 | Half-violin + box + jitter (premium) |
| Mixed Plot | `mixed-plot.ts` | 102 | BLUP caterpillar for LMM |
| PCA Plot | `pca-plot.ts` | 203 | Biplot + scree + loadings heatmap |
| Distribution | `distribution.ts` | 266 | Interactive PDF/CDF explorer |

### Group 2: Distribution & Comparison (7)

| Plot | File | Lines | Purpose |
|------|------|-------|---------|
| Density | `density.ts` | 171 | Multi-series KDE with rug |
| Boxplot | `boxplot.ts` | 205 | Standalone box + whisker + outliers |
| Lollipop | `lollipop.ts` | 133 | Dot + line from baseline |
| Dot Plot | `dot-plot.ts` | 153 | Dot position per category |
| Grouped Bar | `grouped-bar.ts` | 179 | Multi-series categorical bars |
| Line Chart | `line-chart.ts` | 151 | Connected line series |
| Bubble Chart | `bubble-chart.ts` | 181 | Scatter with size dimension |

### Group 3: Ranking & Composition (7)

| Plot | File | Lines | Purpose |
|------|------|-------|---------|
| Pareto | `pareto.ts` | 202 | Bars + cumulative % line |
| Funnel | `funnel.ts` | 188 | Stage-by-stage narrowing |
| Pie Chart | `pie-chart.ts` | 160 | Slices with donut option |
| Area Chart | `area-chart.ts` | 218 | Stacked/baseline area |
| Forest Plot | `forest-plot.ts` | 232 | CI bars for effect sizes |
| ROC Curve | `roc-curve.ts` | 172 | Sensitivity vs 1-specificity |
| Strip Plot | `strip-plot.ts` | 158 | Jittered points (no violin) |

### Group 4: Multivariate & Categorical (8)

| Plot | File | Lines | Purpose |
|------|------|-------|---------|
| Swarm Plot | `swarm-plot.ts` | 245 | Force-directed non-overlapping jitter |
| Mosaic Plot | `mosaic-plot.ts` | 306 | Area-proportional tiles |
| Pair Plot | `pair-plot.ts` | 216 | Matrix of scatter/density |
| Radar Chart | `radar-chart.ts` | 205 | Radial axes (spider web) |
| Parallel Coords | `parallel-coords.ts` | 188 | Multiple parallel axes |
| Treemap | `treemap.ts` | 196 | Hierarchical area-proportional rectangles |
| Waffle Chart | `waffle-chart.ts` | 159 | Grid squares by proportion |
| Sparkline | `sparkline.ts` | 126 | Inline mini time-series |

### Group 5: Hierarchical & Proportional (2)

| Plot | File | Lines | Purpose |
|------|------|-------|---------|
| Sunburst | `sunburst.ts` | 256 | Radial partition tree |
| Marimekko | `marimekko.ts` | 346 | Width-varied stacked bar |

### Group 6: Networks & Flows (4)

| Plot | File | Lines | Purpose |
|------|------|-------|---------|
| Chord Diagram | `chord-diagram.ts` | 244 | Circular network with ribbons |
| Arc Diagram | `arc-diagram.ts` | 264 | Linear nodes with curved arcs |
| Alluvial Plot | `alluvial-plot.ts` | 403 | Sankey-style flow diagram |
| Edge Bundling | `edge-bundling.ts` | 409 | Force-directed bundled curves |

### Group 7: Specialized (2)

| Plot | File | Lines | Purpose |
|------|------|-------|---------|
| Dendrogram | `dendrogram.ts` | 311 | Hierarchical clustering tree |
| FA Plot | `fa-plot.ts` | 1,349 | Factor analysis (6 sub-types) |

---

## 25. Viz Layer: FA Path Diagram (Flagship Visualization)

**File**: `src/viz/plots/fa-plot.ts` (1,349 lines)

The largest and most sophisticated visualization in Carm. Renders 6 sub-types via `FAPlotType`:

| Sub-type | Purpose | Complexity |
|----------|---------|-----------|
| `'path'` | CFA/EFA path diagram (AMOS-style) | ~600 lines |
| `'scree'` | Scree + parallel analysis + Kaiser line | ~150 lines |
| `'loadings'` | Factor loadings heatmap | ~200 lines |
| `'communality'` | Communality bar chart | ~100 lines |
| `'factor-correlation'` | Factor inter-correlation matrix | ~100 lines |
| `'fit-indices'` | Fit indices dashboard (CFI, RMSEA, SRMR, TLI) | ~100 lines |

### Path Diagram Architecture

The path diagram is a fully parametrized SEM-style visualization with **30+ configurable style parameters** (`FAPathStyle`):

**Visual elements**:
- Factor ellipses with radial gradients and optional glow filter
- Item rectangles with linear gradients and colored accent bars
- Error circles with self-loop arcs
- Loading arrows as Bézier curves with width/opacity proportional to magnitude
- Loading labels as halo-text pills
- Factor covariance arcs (nested on left side)
- Group brackets (EFA only)
- Fit indices card

**SVG features used**:
- `<filter>` — Gaussian blur glow, drop shadows (medium + small)
- `<radialGradient>` — Factor ellipses (55% white center → color edge)
- `<linearGradient>` — Item rectangles
- `<marker>` — Custom arrowheads per factor + covariance markers
- `paint-order: stroke` — White text halo for readability

**Layout algorithm**:
1. Factors positioned vertically, evenly spaced
2. Items positioned right of their primary factor, stacked with `itemGap`
3. Errors positioned right of items with `errorSpan`
4. Covariance arcs on left with `covArcReserve` spacing
5. Cross-loadings rendered as dashed arrows (EFA only, above threshold)

---

## 26. Viz Layer: Export Pipeline

### SVG Export

```
SVG element → XMLSerializer → add xmlns → Blob('image/svg+xml') → download
```

### PNG Export (Publication Quality)

```
SVG element → XMLSerializer → encodeURIComponent → data:image/svg+xml
  → Image.src → Canvas (300/96 × upscale) → Canvas.toBlob('image/png') → download
```

Scale factor 300/96 ≈ 3.125× ensures 300 DPI output for print publication.

---

## 27. Deterministic Reproducibility

### splitmix32 PRNG

Every stochastic function in Carm (clustering, factor analysis random starts, parallel analysis) uses a seeded splitmix32 generator:

```typescript
function splitmix32(seed: number): () => number {
  return () => {
    seed += 0x9E3779B9
    let t = seed
    t = Math.imul(t ^ (t >>> 15), t | 1)
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61)
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}
```

- Full period: 2³²
- Default seed: 42
- All functions accept `seed?` in their options

**Box-Muller** for normal variates (used in Haar matrix generation, parallel analysis):

```typescript
u1 = rng(), u2 = rng()
z = sqrt(-2·ln(u1)) · cos(2π·u2)
```

**Haar-distributed random orthogonal matrices** (used in FA random starts):
1. Generate k×k matrix of standard normals
2. QR decomposition via modified Gram-Schmidt
3. Sign correction: `Q[:,j] *= sign(R[j,j])` to ensure Haar measure

---

## 28. Cross-Validation Infrastructure

### Test Suite (4,456 lines, 496 tests)

| Test file | Lines | Tests | Validated against |
|-----------|-------|-------|------------------|
| `math.test.ts` | 214 | ~30 | R `pt()`, `pf()`, `pchisq()`, `qnorm()` |
| `matrix.test.ts` | 190 | ~25 | R `eigen()`, `svd()`, `chol()` |
| `descriptive.test.ts` | 123 | ~15 | R `shapiro.test()`, `e1071::skewness()` |
| `comparison.test.ts` | 157 | ~20 | R `t.test()`, `wilcox.test()`, `kruskal.test()` |
| `correlation.test.ts` | 139 | ~15 | R `cor.test()`, `ppcor::pcor()` |
| `regression.test.ts` | 143 | ~15 | R `lm()`, `glm()` |
| `analyze.test.ts` | 281 | ~30 | R full pipeline |
| `mixed.test.ts` | 128 | ~15 | R `lme4::lmer()` |
| `pca.test.ts` | 103 | ~12 | R `prcomp()`, `psych::principal()` |
| `clustering.test.ts` | 702 | ~40 | R `mclust`, `poLCA`, `seqHMM` |
| `clustering-xval.test.ts` | 440 | ~30 | R `mclust` (deterministic) |
| `clustering-engagement.test.ts` | 180 | ~12 | Real dataset validation |
| `dbscan.test.ts` | 230 | ~15 | R `dbscan::dbscan()` |
| `hac.test.ts` | 312 | ~20 | R `hclust()` |
| `preprocess.test.ts` | 193 | ~15 | R `scale()` |
| `stress.test.ts` | 576 | ~30 | Edge cases, large inputs |
| `large_sample.test.ts` | 345 | ~20 | N=10,000+ performance |

### Validation Reports

6 comprehensive technical reports in `validation/`:
- `FA-TECHNICAL-REPORT.md` (1,629 lines)
- `GMM-TECHNICAL-REPORT.md` (1,536 lines)
- `CORE-MATH-TECHNICAL-REPORT.md` (1,566 lines)
- `DESCRIPTIVE-COMPARISON-TECHNICAL-REPORT.md` (1,338 lines)
- `CORRELATION-REGRESSION-TECHNICAL-REPORT.md` (1,194 lines)
- `MIXED-MODEL-TECHNICAL-REPORT.md` (925 lines)

Total: 8,188 lines of validation documentation.

---

## 29. Public API Reference

### Core Exports

```typescript
// Math
logGamma, gamma, betaFn, logBeta, incompleteBeta, incompleteGamma
normalCDF, normalQuantile, tDistCDF, tDistPValue, tDistQuantile
fDistCDF, fDistPValue, chiSqCDF, chiSqPValue, chiSqQuantile
nelderMead, adjustPValues
mean, variance, sd, se, median, quantile, rank, cov, clamp, roundTo

// Matrix
class Matrix { fromArray, identity, zeros, multiply, inverse, pseudoInverse,
  cholesky, eigen, svd, logDet, transpose, trace, diagonal, add, subtract, scale }
solveLinear

// Types
StatResult, DescriptiveResult, RegressionResult, PairwiseResult, EffectSize,
AnalysisResult, FrequencyTestResult, FAResult, CFAResult, FADiagnostics, FactorFit,
GMMResult, KMeansResult, LCAResult, LTAResult, ClusterDiagnostics,
LMMResult, PCAResult, Field, FieldType, PAdjMethod, CovarianceModel

// APA
formatP, formatStat, formatCI, formatDF, formatTTest, formatANOVA,
formatChiSq, formatCorrelation, formatRegression, formatMannWhitney,
formatKruskalWallis, formatLMM, formatCFAFit,
interpretCohensD, interpretEtaSq, interpretR, interpretCramerV, interpretEffect
```

### Stats Exports

```typescript
// Descriptive
describe, trimmedMean, skewness, kurtosis, ciMean, shapiroWilk

// Comparison
tTestIndependent, tTestPaired, oneWayANOVA, mannWhitneyU,
wilcoxonSignedRank, kruskalWallis, friedmanTest

// Effect Size
cohensD, cohensDPaired, hedgesG, etaSquared, omegaSquared,
rankBiserial, rankBiserialWilcoxon, etaSquaredKW, cohensDCI

// Post-Hoc
tukeyHSD, gamesHowell, dunnTest

// Correlation
pearsonCorrelation, spearmanCorrelation, kendallTau,
partialCorrelation, correlationMatrix

// Regression
linearRegression, multipleRegression, polynomialRegression,
logisticRegression, regressionDiagnostics

// Frequency
frequencyTable, contingencyTable, chiSquareTest, fisherExactTest,
phiCoefficient, goodnessOfFit

// PCA
runPCA, varimaxRotation, screeData

// Mixed
runLMM, computeBLUPs

// Clustering
fitGMM, predictGMM, fitGMMRange, fitKMeans, fitKMeansRange,
fitLCA, fitLTA, runDBSCAN, runHierarchical, cutTree, silhouetteScores

// Factor Analysis
runEFA, runCFA, runFADiagnostics (aliased as computeFADiagnostics)

// Analyze
analyze, detectFieldType

// Preprocess
preprocessData, inverseTransform
```

### Viz Exports

```typescript
// Theme
DEFAULT_THEME, DARK_THEME, applyTheme, getColor, themeColorScale

// Components
addSubtitle, addCaption, addRegressionEquation, addNLabel, addStatBadge
renderBrackets, formatBracketP, totalBracketHeight
showTooltip, hideTooltip, formatTooltipRow

// 42 Plot Renderers
renderViolinBox, renderScatterStats, renderHistogram, renderBarStats,
renderCorrelogram, renderCoefPlot, renderQQPlot, renderResidualPanel,
renderRaincloud, renderMixedPlot, renderPCAPlot, renderDistribution,
renderDensity, renderBoxplot, renderLollipop, renderDotPlot,
renderGroupedBar, renderLineChart, renderBubbleChart,
renderPareto, renderFunnel, renderPieChart, renderAreaChart,
renderForestPlot, renderROCCurve, renderStripPlot,
renderSwarmPlot, renderMosaicPlot, renderPairPlot,
renderRadarChart, renderParallelCoords, renderTreemap,
renderWaffleChart, renderSparkline,
renderSunburst, renderMarimekko,
renderChordDiagram, renderArcDiagram, renderAlluvialPlot, renderEdgeBundling,
renderDendrogram, renderFAPlot

// Export
exportSVG, exportPNG
```

---

## 30. Engineering Decisions: Problems, Solutions, and Optimizations

### 30.1 Zero External Math Dependencies

**Problem**: JavaScript's math ecosystem is fragmented. jStat has API inconsistencies. simple-statistics lacks advanced distributions. Both have abandoned maintenance periods.

**Root cause**: Statistical libraries need precise numerical computation (15-digit accuracy for p-values, stable matrix decompositions). Third-party libraries may change algorithms between minor versions, breaking cross-validation against R.

**Solution**: Implement every distribution function, special function, and matrix decomposition from scratch in `core/math.ts` and `core/matrix.ts`.

**Why this over alternatives**: Using jStat would save ~500 lines but introduce a dependency that could change behavior in patch releases. Using WASM (e.g., compiled Fortran LAPACK) would add a build step and break in environments without WASM support.

**Result**: 961 lines of pure math. Complete control over numerical behavior. Zero supply chain risk. Every algorithm is readable, auditable, and tested against R.

### 30.2 D3 as Peer Dependency with Dynamic Import

**Problem**: D3 is ~200KB. Consumers who only need statistics shouldn't pay for visualization code.

**Root cause**: Traditional bundling includes all transitive dependencies. Even tree-shaking can't eliminate D3 if any viz function is imported.

**Solution**: Declare D3 as `peerDependencies`. Every viz function uses `import('d3').then(...)` for lazy loading. Export stats and viz through separate entry points (`carm/stats`, `carm/viz`).

**Why this over alternatives**: Bundling D3 would force all consumers to include it. Making D3 optional via `try/require` would lose TypeScript types. Dynamic import is the standard ESM pattern for optional dependencies.

**Result**: `import { pearsonCorrelation } from 'carm/stats'` adds zero bytes of D3. `import { renderViolinBox } from 'carm/viz'` loads D3 only when the first plot renders.

### 30.3 Immutable Results with `readonly`

**Problem**: Statistical results must not be accidentally mutated after computation. A consumer modifying `result.loadings[0][1] = 999` would silently corrupt data.

**Root cause**: JavaScript arrays and objects are mutable by default. TypeScript's `readonly` is compile-time only — but it catches the vast majority of accidental mutations.

**Solution**: All result interfaces use `readonly` arrays and properties. Every function returns frozen-at-the-type-level results.

**Why this over alternatives**: `Object.freeze()` at runtime has performance cost and doesn't work recursively without a deep freeze utility. TypeScript's `readonly` is zero-cost at runtime and catches mutations at compile time.

**Result**: `result.loadings.push(...)` is a compile error. No runtime overhead.

### 30.4 `StatResult` as Universal Contract

**Problem**: Every test returns different fields. A t-test has df; a chi-square has a contingency table; ANOVA has SS decomposition. How do you build a generic results viewer?

**Root cause**: Statistical tests are inherently diverse. But downstream consumers (APA formatting, plotting, narrative generation) need a common interface.

**Solution**: `StatResult` with: testName, statistic, df, pValue, effectSize, ci, ciLevel, n, formatted. Every test fills all fields. Test-specific results extend `StatResult` (e.g., `FrequencyTestResult extends StatResult` adds `table`).

**Why this over alternatives**: A union type would require consumers to pattern-match on every test type. A minimal interface (just pValue) would lose information. The rich `StatResult` gives consumers everything they need for display without knowing the test type.

**Result**: Any `StatResult` can be displayed as an APA string, plotted with significance brackets, or interpreted by the narrative engine — without knowing whether it came from a t-test or chi-square.

### 30.5 Row-Major Flat Array for Matrix Storage

**Problem**: JavaScript's `number[][]` (array of arrays) has poor cache locality for matrix multiplication. Each row is a separate heap allocation.

**Root cause**: Modern CPUs use cache lines (64 bytes = 8 doubles). Accessing `A[i][k]` then `A[i][k+1]` may cause cache misses if rows are scattered in memory.

**Solution**: `Matrix` uses a single flat `Float64Array` with row-major layout: element (i,j) at index `i*cols + j`. The multiply method uses i-k-j loop order for sequential access of the second matrix's elements.

**Why this over alternatives**: Column-major (Fortran convention) would optimize `A'x` but we more often multiply `AB`. Nested arrays are simpler but 2-4× slower for matrix multiply due to cache misses.

**Result**: Matrix multiply is ~2× faster than nested arrays for typical sizes (p < 100). `Float64Array` also prevents silent coercion of non-numeric values.

### 30.6 Lanczos over Stirling for logGamma

**Problem**: The gamma function is the foundation of t, F, and chi-square distributions. Stirling's approximation gives ~8 digits of precision — not enough for p-values near machine epsilon.

**Root cause**: Stirling's `log Γ(z) ≈ (z-½)log(z) - z + ½log(2π) + 1/(12z) - ...` is asymptotic — it works well for large z but poorly for z < 10.

**Solution**: Lanczos approximation with g=7 and 9 coefficients from Numerical Recipes. Gives ~15 digits of precision for all z > 0.5. Reflection formula `Γ(z)Γ(1-z) = π/sin(πz)` handles z < 0.5.

**Why this over alternatives**: Spouge's formula is simpler but less accurate. The exact coefficients from NR have been validated against Mathematica to 15 digits.

**Result**: `logGamma(0.001)` matches R's `lgamma(0.001)` to 14 significant figures.

### 30.7 Two-Phase ML Extraction for Factor Analysis

**Problem**: R's `factanal()` uses L-BFGS-B (bounded quasi-Newton) for ML factor extraction. Implementing L-BFGS-B from scratch is ~400 lines of complex code.

**Root cause**: ML factor extraction optimizes over uniquenesses Ψ (constrained to [0.005, 0.995]). Nelder-Mead is unconstrained and slow in high dimensions. L-BFGS-B is the right tool but prohibitively complex to implement.

**Solution**: Phase 1 uses Jöreskog gradient descent with cosine-annealed learning rate and momentum for fast convergence to the vicinity. Phase 2 uses Nelder-Mead with death penalty (1000× for out-of-bounds) for final polish.

**Why this over alternatives**: Nelder-Mead alone is too slow (hundreds of iterations). Gradient descent alone oscillates near the optimum. The two-phase approach combines the speed of gradient methods with the robustness of simplex methods.

**Result**: Matches R's `factanal()` to 4 decimal places on loadings. Convergence in ~50+20 iterations vs ~200 for Nelder-Mead alone.

### 30.8 50 Random Starts for Oblique Rotations

**Problem**: Oblique rotations (oblimin, quartimin, geomin) have multiple local optima. A single starting point may converge to a suboptimal rotation.

**Root cause**: The GPFoblq criterion is non-convex. Different starting rotation matrices lead to different solutions. The global optimum is not guaranteed by any single run.

**Solution**: 50 random starts. Start 0 = identity (deterministic reference). Starts 1-49 = Haar-distributed random orthogonal matrices. Keep the solution with the best (lowest) criterion value.

**Why this over alternatives**: Tested empirically on rraw (525×31, 5 factors) and teacher burnout (876×23, 4 factors). 10 starts are sufficient for teacher burnout but rraw needs 30+. 50 gives a safety margin. Timing: 50 starts ≈ 3s, acceptable for interactive use.

**Result**: Consistent global optimum on all tested datasets. Matches lavaan's GPA strategy (which uses 30 starts by default).

### 30.9 Woodbury Identity for LMM

**Problem**: The marginal covariance `V = σ²_eI + σ²_bZZ'` is n×n. Inverting it directly is O(n³) — prohibitive for n > 10,000.

**Root cause**: In a random-intercept model with q groups, Z is n×q and ZZ' is n×n. But q << n (typically 10-100 groups vs 1,000-100,000 observations).

**Solution**: Woodbury identity: `V⁻¹ = I/σ²_e - Z(σ²_eI_q + σ²_bZ'Z)⁻¹Z'/σ²_e`. This requires inverting a q×q matrix instead of n×n. Similarly, `log|V|` via matrix determinant lemma avoids the n×n determinant.

**Why this over alternatives**: Cholesky of V (n×n) is O(n³/3). Sparse Cholesky (exploiting block structure) is complex to implement. Woodbury is exact, simple, and reduces to O(q³ + nq²).

**Result**: LMM with n=10,000 and q=50 runs in ~200ms (vs ~30s for naive V inversion).

### 30.10 Log-Space EM for GMM and LTA

**Problem**: GMM responsibilities `γ_ik = π_k · N(x_i|μ_k,Σ_k) / Σ_j π_j · N(x_i|μ_j,Σ_j)` involve products of small probabilities that underflow to 0.

**Root cause**: For d=20 dimensions, `N(x|μ,Σ)` can be as small as 10⁻¹⁰⁰. Multiplying by weights and summing in linear space produces 0/0 = NaN.

**Solution**: Compute everything in log-space. `log γ_ik = log π_k + log N(x_i|μ_k,Σ_k) - logSumExp_j(log π_j + log N(x_i|μ_j,Σ_j))`. The logSumExp trick: `logΣ exp(a_j) = max(a) + log(Σ exp(a_j - max(a)))`.

**Why this over alternatives**: Scaling (dividing by max) works for simple cases but still underflows when the gap between log-probabilities exceeds ~700 (the log of the double max). Log-space with logSumExp handles arbitrary gaps.

**Result**: Stable EM convergence for d=30 dimensions. No NaN propagation even with widely separated clusters.

### 30.11 splitmix32 Over Math.random()

**Problem**: `Math.random()` is non-deterministic — different results on every page load. This makes clustering and FA results non-reproducible.

**Root cause**: `Math.random()` uses implementation-specific PRNG (V8 uses xorshift128+) with an automatically generated seed. There is no standard API to seed it.

**Solution**: splitmix32 — a fast 32-bit PRNG with known properties (full period 2³², passes BigCrush). Default seed = 42. Every stochastic function accepts `seed?` in options.

**Why this over alternatives**: xoshiro128** is higher quality but more complex. Mersenne Twister has 2.5KB of state. splitmix32 is 6 lines of code, has full period, and is fast enough for statistical use (not cryptographic).

**Result**: `fitGMM(data, { seed: 42 })` produces identical results on every call, every browser, every platform. Cross-validation against R is possible because results are deterministic.

### 30.12 Bartlett Correction: EFA Yes, CFA No

**Problem**: Chi-square test for factor model fit uses different scaling in R's `psych::fa()` (EFA) vs lavaan (CFA). Applying the wrong correction produces 5-20% discrepancies.

**Root cause**: EFA uses Bartlett's (1950) correction: `χ² = (n-1-(2p+5)/6-2k/3) × F`. CFA (lavaan) uses uncorrected: `χ² = (n-1) × F`. This is not documented in most textbooks — discovered by reading R source code.

**Solution**: Apply Bartlett correction in `runEFA()`, omit it in `runCFA()`. Both produce correct chi-square values matching their respective R implementations.

**Result**: EFA chi-square matches `psych::fa()` to 4 decimal places. CFA chi-square matches lavaan to 4 decimal places.

### 30.13 Heywood Case Clamping

**Problem**: PAF communalities can exceed 1.0 (Heywood case). This makes the reduced correlation matrix non-positive-definite, crashing eigendecomposition.

**Root cause**: When items are highly correlated, the iterative PAF procedure can produce communality estimates > 1.0 — a mathematical impossibility but a numerical artifact.

**Solution**: Clamp communalities to [0.001, 0.9999]. The 0.001 floor prevents items from being completely excluded. The 0.9999 ceiling prevents singularity.

**Why this over alternatives**: Setting a hard boundary at 1.0 allows the matrix to be exactly singular (determinant = 0). The 0.0001 buffer prevents this while keeping the communality effectively at the boundary.

**Result**: PAF converges on all tested datasets including those with r > 0.95 between items.

### 30.14 Fisher's Exact Test in Log-Space

**Problem**: Fisher's exact test computes hypergeometric probabilities: `P = C(K,k)·C(N-K,n-k)/C(N,n)`. For large contingency tables, these combinatorials overflow.

**Root cause**: `C(100,50)` ≈ 10²⁹ — far beyond double precision max (~1.8×10³⁰⁸ but with only 15 digits of precision). Products of large combinatorials produce Infinity.

**Solution**: Compute in log-space: `logP = logC(K,k) + logC(N-K,n-k) - logC(N,n)` where `logC(n,k) = logFactorial(n) - logFactorial(k) - logFactorial(n-k)`.

**Result**: Fisher's exact test works for contingency tables up to N ≈ 10,000 without overflow.

### 30.15 Tukey HSD via Bonferroni Approximation

**Problem**: True Tukey HSD requires the studentized range distribution quantile function `qtukey()`. Implementing this requires either precomputed tables or a complex integral approximation.

**Root cause**: The studentized range distribution has no closed-form CDF. R uses Algorithm AS 190 (Lund & Lund, 1983) with iterative numerical integration — ~200 lines of Fortran.

**Solution**: Conservative Bonferroni-adjusted t-test approximation: `p_adj = min(1, k(k-1)/2 × p_unadjusted)`. This is conservative (wider CIs, larger p-values) but correct.

**Why this over alternatives**: Implementing AS 190 would be ~200 lines of code for a function used only in post-hoc tests. The Bonferroni approximation is a standard alternative used in many software packages.

**Result**: Tukey HSD p-values are conservative (larger than exact). This is acceptable because it never produces false positives — only occasional false negatives.

---

## 31. Mathematical Tricks That Made It Possible

### 31.1 Lanczos Approximation for logGamma

**Why needed**: `Γ(z)` is the foundation of all distribution CDFs (t, F, χ²). Without it, no hypothesis testing.

**The trick**: Transform `Γ(z)` into a sum of rational functions plus an exponential:
```
Γ(z+1) = √(2π) · (z + g + ½)^(z+½) · e^(-(z+g+½)) · (c₀ + c₁/(z+1) + c₂/(z+2) + ...)
```
With g=7 and 9 specific coefficients, this gives 15 digits of precision for all z > 0.5.

**Implementation**: `core/math.ts` lines 1-40

**Impact**: One function replaces tables, infinite products, and series expansions. All subsequent distribution functions build on it.

### 31.2 Lentz's Continued Fraction for incompleteBeta

**Why needed**: `I_x(a,b)` (regularized incomplete beta) is the CDF of the t, F, and beta distributions. The standard integral has no closed form.

**The trick**: Express `I_x(a,b)` as a continued fraction `a₁/(b₁ + a₂/(b₂ + ...))` and evaluate via Lentz's modified algorithm. The key insight: maintain two running values C and D, updating them per term, and the ratio converges to the fraction's value.

**Implementation**: `core/math.ts` lines 50-100. FPMIN = 1e-30 prevents 0/0. Convergence at ε = 3e-7.

**Impact**: One function handles t-distribution CDF, F-distribution CDF, and beta distribution CDF. Without it, Carm would need separate numerical methods for each.

### 31.3 logSumExp for Stable Mixture Model EM

**Why needed**: GMM E-step computes `γ_ik ∝ exp(log π_k + log N(...))`. Direct exponentiation overflows/underflows.

**The trick**: `log(Σ exp(aᵢ)) = max(a) + log(Σ exp(aᵢ - max(a)))`. By subtracting the max first, the largest term becomes exp(0) = 1, and all other terms are ≤ 1.

**Implementation**: `stats/clustering.ts` E-step

**Impact**: GMM converges stably in 20-50+ dimensions. Without this, EM fails for d > 10.

### 31.4 Woodbury Identity for LMM

**Why needed**: `V = σ²_eI + σ²_bZZ'` is n×n. Inverting it is O(n³). For n=10,000, this is ~20 seconds.

**The trick**: `V⁻¹ = (1/σ²_e)(I - Z(Z'Z + (σ²_e/σ²_b)I)⁻¹Z')`. The inner matrix is q×q (number of groups). For q=50, inversion is ~0.1ms.

**Implementation**: `stats/mixed.ts` lines 100-150

**Impact**: LMM scales to n=100,000 observations with q=100 groups in seconds.

### 31.5 Concentrated ML Objective for Factor Analysis

**Why needed**: ML extraction optimizes over both loadings Λ (p×k matrix) and uniquenesses Ψ (p-vector). That's pk+p parameters.

**The trick**: For fixed Ψ, the optimal Λ can be computed analytically from the eigendecomposition of Ψ⁻¹/²RΨ⁻¹/². So the objective becomes a function of Ψ only: `F(Ψ) = -Σ(log λⱼ - λⱼ) + k - d`. This halves the optimization dimension from pk+p to just p.

**Implementation**: `stats/factor-analysis.ts` lines 494-513

**Impact**: ML extraction for 31 items × 5 factors optimizes over 31 parameters instead of 186. Converges in 50 iterations instead of 500.

### 31.6 Cosine Annealing for ML Gradient Descent

**Why needed**: Fixed learning rate either oscillates near the optimum (too large) or converges too slowly (too small).

**The trick**: `lr = lr_min + (lr₀ - lr_min) × 0.5 × (1 + cos(π × iter/maxIter))`. Starts at lr₀ = 0.02, smoothly decays to lr_min = 0.001 following a cosine curve.

**Implementation**: `stats/factor-analysis.ts` Phase 1 gradient loop

**Impact**: Fast early convergence (high LR) + stable final convergence (low LR). Removes the need to hand-tune the learning rate.

### 31.7 Eigendecomposition for Varimax Polar Factor

**Why needed**: Varimax rotation needs the polar decomposition of a k×k matrix. The standard approach is SVD, but Carm's SVD is one-sided Jacobi — optimized for tall matrices, not small square ones.

**The trick**: For a k×k matrix B, compute B'B → eigendecompose → singular values = √(eigenvalues) → U = BV·diag(1/σ) → polar factor T = UV'. This uses the existing `Matrix.eigen()` method (well-tested) instead of SVD.

**Implementation**: `stats/factor-analysis.ts` varimax section

**Impact**: Numerically stable varimax for any k ≤ 30. Matches R's `varimax()` to machine precision.

### 31.8 Modified Gram-Schmidt for Haar QR

**Why needed**: Random starts for oblique rotations require Haar-distributed random orthogonal matrices. Naive random matrix + classical Gram-Schmidt loses orthogonality for k > 5.

**The trick**: Modified Gram-Schmidt reorthogonalizes against already-processed columns at each step, maintaining `|Q'Q - I| < 1e-14` even for k = 10. Plus sign correction `Q[:,j] *= sign(R[j,j])` ensures Haar measure.

**Implementation**: `stats/factor-analysis.ts` random start generation

**Impact**: 50 random starts produce truly uniform-on-SO(k) rotations. Without this, rotation starts would be biased toward certain orientations.

### 31.9 Dynamic Programming for Wilcoxon Exact Distribution

**Why needed**: For n ≤ 20, the exact Wilcoxon signed-rank distribution gives precise p-values. The brute-force approach (enumerate all 2ⁿ sign assignments) is O(2ⁿ) — too slow for n > 15.

**The trick**: Build `dist[w]` iteratively. Start with rank 1: dist[0] = 1, dist[1] = 1. For each subsequent rank r, update: `dist[w] += dist[w-r]` (choosing + for rank r). Total states: O(n × maxW) where maxW = n(n+1)/2.

**Implementation**: `stats/comparison.ts` wilcoxonSignedRank exact path

**Impact**: Exact p-values for n ≤ 20 in O(n²) time instead of O(2ⁿ). For n=20: ~200 operations vs ~1,000,000.

### 31.10 Cholesky Log-Determinant

**Why needed**: `log|V|` appears in every likelihood function (LMM, GMM, FA). Computing `det(V)` directly overflows for d > 30.

**The trick**: If V = LL' (Cholesky), then `log|V| = 2·Σ log(L_ii)`. This is O(n²) for the Cholesky + O(n) for the sum, and never overflows because it stays in log-space.

**Implementation**: `core/matrix.ts` `logDet()` method

**Impact**: Log-determinant works for matrices up to d = 1000+. Without this, GMM with d > 30 dimensions would fail.

---

## References

- Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences* (2nd ed.). Erlbaum.
- Lanczos, C. (1964). A precision approximation of the gamma function. *SIAM J. Numer. Anal.*, 1, 86–96.
- Press, W., Teukolsky, S., Vetterling, W., & Flannery, B. (2007). *Numerical Recipes* (3rd ed.). Cambridge.
- Abramowitz, M. & Stegun, I. (1972). *Handbook of Mathematical Functions*. Dover.
- Acklam, P. (2004). An algorithm for computing the inverse normal CDF.
- Royston, P. (1995). Remark AS R94: A remark on Algorithm AS 181. *Appl. Statist.*, 44, 547–551.
- Kaiser, H.F. (1958). The varimax criterion for analytic rotation in factor analysis. *Psychometrika*, 23, 187–200.
- Jennrich, R.I. (2002). A simple general method for oblique rotation. *Psychometrika*, 67, 7–19.
- Browne, M.W. & Cudeck, R. (1993). Alternative ways of assessing model fit. In *Testing Structural Equation Models*, 136–162.
- Hu, L. & Bentler, P.M. (1999). Cutoff criteria for fit indexes. *Struct. Eq. Modeling*, 6, 1–55.
- Henderson, C.R. (1950). Estimation of genetic parameters. *Annals Math. Statist.*, 21, 309–310.
- Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models using lme4. *J. Stat. Software*, 67, 1–48.
- Scrucca, L., Fop, M., Murphy, T.B., & Raftery, A.E. (2016). mclust 5: Clustering, classification, and density estimation. *The R Journal*, 8, 289–317.
- Woodbury, M.A. (1950). Inverting modified matrices. *Statistical Research Group*, Memo Rep. 42.
