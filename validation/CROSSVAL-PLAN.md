# Comprehensive Cross-Validation Plan: Carm vs R (Saqrlab-Powered)

## Context

Carm has 69 exported stats functions across 13 modules, but only Factor Analysis has rigorous cross-validation (100/100 promax, 100/100 geomin). Every other module lacks systematic R equivalence testing. This plan uses the `Saqrlab` R package's `simulate_data()` function to generate thousands of random datasets across all 6 data types and cross-validate every Carm stats function against R reference implementations. The goal: prove laser-sharp numerical equivalence across the entire library.

## Architecture

Each validation round follows the established pattern:
```
R script (Saqrlab → compute refs → JSON)  →  TS harness (Carm → compare → HTML)
```

### File Layout
```
validation/
├── r-reference/
│   ├── ttest-ref.R              # Round 1 (1000 datasets)
│   ├── correlation-ref.R        # Round 2 (1000 datasets)
│   ├── anova-ref.R              # Round 3 (500 datasets)
│   ├── regression-ref.R         # Round 4 (500 datasets)
│   ├── pca-ref.R                # Round 5 (500 datasets)
│   ├── clustering-ref.R         # Round 6 (200 datasets)
│   └── fa-extended-ref.R        # Round 7 (200 datasets)
├── ts-harness/
│   ├── ttest-report.ts
│   ├── correlation-report.ts
│   ├── anova-report.ts
│   ├── regression-report.ts
│   ├── pca-report.ts
│   ├── clustering-report.ts
│   └── fa-extended-report.ts
├── data/
│   ├── ttest-ref.json           # ~30-50 MB for 1000 datasets
│   ├── correlation-ref.json
│   ├── anova-ref.json
│   ├── regression-ref.json
│   ├── pca-ref.json
│   ├── clustering-ref.json
│   └── fa-extended-ref.json
└── reports/
    ├── ttest-report.html
    ├── correlation-report.html
    ├── anova-report.html
    ├── regression-report.html
    ├── pca-report.html
    ├── clustering-report.html
    └── fa-extended-report.html
```

### Saqrlab API
```r
library(Saqrlab)
d <- simulate_data("ttest", seed = 42)          # columns: group (factor A/B), score (numeric)
d <- simulate_data("anova", seed = 42)          # columns: group (factor G1..Gk), score (numeric)
d <- simulate_data("correlation", seed = 42)    # columns: x1..xp (numeric, p=4-7)
d <- simulate_data("clusters", seed = 42)       # columns: x1..xd (numeric), true_cluster (int)
d <- simulate_data("factor_analysis", seed = 42)# columns: x1..xp (numeric); attrs: n_factors, loadings
d <- simulate_data("prediction", seed = 42)     # columns: y, x1-x4 (numeric), cat1, cat2 (factor)
```
- `seed` controls both reproducibility AND structural parameters (n, effect_size, n_groups, etc.)
- Parameters can be overridden via `...`: `simulate_data("ttest", seed=42, n=200, effect_size=0.8)`

### Tolerance Tiers
| Category | Tolerance | Rationale |
|---|---|---|
| Deterministic (eigenvalues, corr matrix, sums) | 1e-8 | Pure floating-point arithmetic |
| Descriptive stats (mean, sd, skew, kurtosis) | 1e-8 | Closed-form formulas |
| Hypothesis test statistics (t, F, chi2, U, W) | 1e-4 | Same formula, minor float differences |
| p-values | 1e-3 | Distribution function approximations differ |
| Confidence intervals | 1e-3 | Quantile function + stat tolerance |
| Effect sizes (d, g, eta2, omega2) | 1e-4 | Closed-form |
| Regression coefficients | 1e-3 | Matrix inversion differences |
| Logistic regression | 1e-2 | Iterative IRLS convergence |
| EM-based (GMM log-lik) | 2.0 | Different local optima |
| EM-based (BIC) | 5.0 | LL tolerance amplified by penalty |
| Rotated loadings (FA) | MAE ≤ 0.05 | Sign/permutation indeterminacy |

### HTML Report Template
All reports use the established dark GitHub-inspired theme (from `fa-full-report.ts`):
- Summary cards at top: total datasets, pass/fail counts
- Aggregate table: metric, threshold, mean/median/p95/max error, pass rate
- Per-dataset detail table: every dataset with per-metric color coding
- Color coding: green (≤ threshold), yellow (≤ 3x threshold), red (> 3x)
- Timestamp and metadata in footer

---

## Round 1: T-Test & Descriptive Suite (1000 datasets)

**Saqrlab type**: `"ttest"` — 2-group data with `group` (A/B) and `score`

### R Reference Script (`validation/r-reference/ttest-ref.R`)

For each seed 1..1000:
```r
d <- simulate_data("ttest", seed = i)
g <- split(d$score, d$group)
x1 <- g[[1]]; x2 <- g[[2]]
```

Compute and save:

| Category | R Function | Carm Function | Fields |
|---|---|---|---|
| **Descriptive (per group)** | `mean()`, `median()`, `sd()`, `var()` | `mean`, `median`, `sd`, `variance` | value |
| | `mean(x, trim=0.1)` | `trimmedMean(x, 0.1)` | value |
| | Custom `skewness/kurtosis` (see note) | `skewness`, `kurtosis` | value |
| | `se = sd/sqrt(n)` | `se` | value |
| | `t.test(x)$conf.int` | `ciMean` | lower, upper |
| | `shapiro.test(x)` | `shapiroWilk` | statistic, pValue |
| **T-tests** | `t.test(x1, x2, var.equal=FALSE)` | `tTestIndependent(x1, x2, false)` | statistic, pValue, df, ci |
| | `t.test(x1, x2, var.equal=TRUE)` | `tTestIndependent(x1, x2, true)` | statistic, pValue, df, ci |
| | `t.test(x1, x2, paired=TRUE)` | `tTestPaired(x1[1:n], x2[1:n])` | statistic, pValue, df, ci |
| | `t.test(x1, x2, alt="less")` | `tTestIndependent(x1, x2, false, 0.95, 'less')` | statistic, pValue |
| | `t.test(x1, x2, alt="greater")` | `tTestIndependent(x1, x2, false, 0.95, 'greater')` | statistic, pValue |
| **Nonparametric** | `wilcox.test(x1, x2)` | `mannWhitneyU(x1, x2)` | statistic, pValue |
| | `wilcox.test(x1, x2, alt="less")` | `mannWhitneyU(x1, x2, 'less')` | statistic, pValue |
| | `wilcox.test(x1[1:n], x2[1:n], paired=TRUE)` | `wilcoxonSignedRank(x1, x2)` | statistic, pValue |
| **Effect sizes** | `effsize::cohen.d(x1, x2)` | `cohensD(x1, x2)` | value, ci |
| | `effsize::cohen.d(x1, x2, hedges=TRUE)` | `hedgesG(x1, x2)` | value, ci |
| | Rank-biserial from U | `rankBiserial(U, n1, n2)` | value |
| | NCC CI for d | `cohensDCI(d, n1, n2)` | lower, upper |
| **Frequency** | `chisq.test(table(d$group))` | `goodnessOfFit(observed)` | statistic, pValue |
| | `table(d$group)` | `frequencyTable(d$group)` | counts, percentages |
| | `chisq.test(matrix, correct=FALSE)` | `chiSquareTest(matrix, false)` | statistic, pValue, df |
| | `chisq.test(matrix, correct=TRUE)` | `chiSquareTest(matrix, true)` | statistic, pValue, df |
| | `fisher.test(matrix)` | `fisherExactTest(a, b, c, d)` | pValue |

**Note on skewness/kurtosis**: R `e1071::skewness(x, type=2)` matches Carm's formula (sample excess). Verify formula alignment before generating refs.

**Note on paired tests**: `simulate_data("ttest")` may have unequal group sizes. For paired tests, truncate to `min(n1, n2)` observations.

**Note on frequency tests**: Create a 2×2 contingency table by median-splitting `score` within each group.

### TS Harness (`validation/ts-harness/ttest-report.ts`)

Import from `carm`: `mean`, `median`, `sd`, `variance`, `se`, `trimmedMean`, `skewness`, `kurtosis`, `ciMean`, `shapiroWilk`, `describe`, `tTestIndependent`, `tTestPaired`, `mannWhitneyU`, `wilcoxonSignedRank`, `cohensD`, `hedgesG`, `rankBiserial`, `cohensDCI`, `goodnessOfFit`, `frequencyTable`, `chiSquareTest`, `fisherExactTest`

Per dataset: compute each Carm function, compare against R ref. Report sections:
1. Descriptive Statistics (8 metrics)
2. T-Tests (5 variants: Welch, Student, paired, one-sided less, one-sided greater)
3. Nonparametric Tests (3 variants)
4. Effect Sizes (4 metrics)
5. Frequency Tests (4 tests)

---

## Round 2: Correlation Suite (1000 datasets)

**Saqrlab type**: `"correlation"` — multivariate data with `x1..xp` (p=4-7)

### R Reference Script (`validation/r-reference/correlation-ref.R`)

For each seed 1..1000:
```r
d <- simulate_data("correlation", seed = i)
```

Compute and save:

| Category | R Function | Carm Function | Fields |
|---|---|---|---|
| **Pearson** (all pairs) | `cor.test(x, y, method="pearson")` | `pearsonCorrelation(x, y)` | statistic (t), pValue, ci, effectSize (r) |
| **Spearman** (all pairs) | `cor.test(x, y, method="spearman")` | `spearmanCorrelation(x, y)` | statistic (rho), pValue |
| **Kendall** (all pairs) | `cor.test(x, y, method="kendall")` | `kendallTau(x, y)` | statistic (tau), pValue |
| **Partial** | `ppcor::pcor.test(x, y, controls)` | `partialCorrelation(x, y, controls)` | statistic, pValue, effectSize |
| **Correlation matrix** | `cor(d)` | `correlationMatrix(data, labels, 'pearson')` | full matrix (p×p) |
| | `cor(d, method="spearman")` | `correlationMatrix(data, labels, 'spearman')` | full matrix (p×p) |
| | `cor(d, method="kendall")` | `correlationMatrix(data, labels, 'kendall')` | full matrix (p×p) |

**Partial correlation**: For each dataset, pick `(x1, x2)` as the pair and `x3..xp` as controls. Also test a 3-variable partial: `(x1, x3)` controlling for `x2`.

### TS Harness (`validation/ts-harness/correlation-report.ts`)

Report sections:
1. Pearson Correlation (per-pair: r, t-stat, p-value, CI)
2. Spearman Correlation (per-pair: rho, p-value)
3. Kendall Tau (per-pair: tau, p-value)
4. Partial Correlation (2 configurations)
5. Correlation Matrix (Pearson/Spearman/Kendall: element-wise MAE)

---

## Round 3: ANOVA & Post-Hoc Suite (500 datasets)

**Saqrlab type**: `"anova"` — k-group data with `group` (G1..Gk) and `score`

### R Reference Script (`validation/r-reference/anova-ref.R`)

For each seed 1..500:
```r
d <- simulate_data("anova", seed = i)
groups <- split(d$score, d$group)
```

Compute and save:

| Category | R Function | Carm Function | Fields |
|---|---|---|---|
| **ANOVA** | `summary(aov(score ~ group, data=d))` | `oneWayANOVA(groups)` | statistic (F), pValue, df1, df2 |
| **Effect sizes** | `eta2 = SS_B / SS_T` | `etaSquared(ssB, ssT)` | value |
| | `omega2 formula` | `omegaSquared(ssB, ssT, dfB, msW)` | value |
| **Kruskal-Wallis** | `kruskal.test(score ~ group, data=d)` | `kruskalWallis(groups)` | statistic (H), pValue, df |
| | `eta2_H = (H - k + 1) / (n - k)` | `etaSquaredKW(H, k, n)` | value |
| **Tukey HSD** | `TukeyHSD(aov(...))` | `tukeyHSD(groups, msW, dfW)` | diff, pValue, ci (per pair) |
| **Games-Howell** | `rstatix::games_howell_test(d, score ~ group)` | `gamesHowell(groups)` | diff, pValue, ci (per pair) |
| **Dunn test** | `dunn.test::dunn.test(d$score, d$group, method="bonferroni")` | `dunnTest(groups, 'bonferroni')` | statistic (z), pValue (per pair) |
| | Same with `method="holm"` | `dunnTest(groups, 'holm')` | statistic, pValue (per pair) |
| | Same with `method="bh"` | `dunnTest(groups, 'bh')` | statistic, pValue (per pair) |
| **Friedman** | `friedman.test(matrix)` | `friedmanTest(data)` | statistic, pValue |
| **LMM** | `lme4::lmer(score ~ 1 + (1\|group), data=d)` | `runLMM(input)` | fixedEffects, randomVariance, residualVariance, icc |

**Friedman note**: Requires repeated measures. Reshape anova data: take first `min_group_size` from each group to build a balanced matrix.

**LMM note**: Treat `group` as random effect. The LMM intercept-only model gives ICC and variance components.

**Post-hoc pair matching**: Sort pairwise results by group-pair name to ensure consistent ordering for comparison.

### TS Harness (`validation/ts-harness/anova-report.ts`)

Report sections:
1. One-Way ANOVA (F, p, df)
2. Effect Sizes (eta², omega², eta²_KW)
3. Kruskal-Wallis (H, p)
4. Tukey HSD (per-pair: diff, p, CI)
5. Games-Howell (per-pair: diff, p, CI)
6. Dunn Test (3 p-adjustment methods)
7. Friedman Test (chi², p)
8. Linear Mixed Model (ICC, variance components)

---

## Round 4: Regression Suite (500 datasets)

**Saqrlab type**: `"prediction"` — `y`, `x1`-`x4` (numeric), `cat1`, `cat2` (factor)

### R Reference Script (`validation/r-reference/regression-ref.R`)

For each seed 1..500:
```r
d <- simulate_data("prediction", seed = i)
```

Compute and save:

| Category | R Function | Carm Function | Fields |
|---|---|---|---|
| **Simple linear** | `lm(y ~ x1, data=d)` | `linearRegression(x1, y)` | coefficients (intercept, slope), se, t, p, r², adjR², fStatistic, ci |
| **Multiple** | `lm(y ~ x1 + x2 + x3 + x4, data=d)` | `multipleRegression(y, predictors)` | coefficients, se, t, p, r², adjR², fStatistic, aic, bic |
| **Polynomial** | `lm(y ~ poly(x1, 2, raw=TRUE), data=d)` | `polynomialRegression(x1, y, 2)` | coefficients, r², adjR² |
| **Logistic** | `glm(ybin ~ x1+x2+x3+x4, binomial)` | `logisticRegression(ybin, predictors)` | coefficients, se, z, p, aic |
| **Diagnostics** | `hatvalues()`, `cooks.distance()`, `vif()` | `regressionDiagnostics(result, preds)` | leverage, cooksD, vif |

**Logistic note**: Create binary outcome by median-splitting `y`: `ybin = as.integer(y > median(y))`.

**VIF**: Use `car::vif(model)` for R reference.

**Diagnostics**: Compare first 10 leverage values and Cook's distances (not all n, to keep JSON small), plus all VIF values.

### TS Harness (`validation/ts-harness/regression-report.ts`)

Report sections:
1. Simple Linear Regression (coefficients, R², F-test)
2. Multiple Regression (coefficients, R², adj-R², AIC, BIC)
3. Polynomial Regression (coefficients, R²)
4. Logistic Regression (coefficients, AIC)
5. Regression Diagnostics (leverage, Cook's D, VIF)

---

## Round 5: PCA Suite (500 datasets)

**Saqrlab type**: `"correlation"` — multivariate continuous data (ideal for PCA)

### R Reference Script (`validation/r-reference/pca-ref.R`)

For each seed 1..500:
```r
d <- simulate_data("correlation", seed = i)
pca <- prcomp(d, center = TRUE, scale. = TRUE)
```

Compute and save:

| Category | R Function | Carm Function | Fields |
|---|---|---|---|
| **Eigenvalues** | `pca$sdev^2` | `runPCA(data).eigenvalues` | array |
| **Variance explained** | `summary(pca)$importance[2,]` | `runPCA(data).varianceExplained` | array |
| **Cumulative var** | `summary(pca)$importance[3,]` | `runPCA(data).cumulativeVariance` | array |
| **Loadings** | `pca$rotation` | `runPCA(data).loadings` | matrix (with sign alignment) |
| **Scores** (first 5 rows) | `pca$x[1:5, ]` | `runPCA(data).scores[0:4]` | matrix subset |
| **Varimax rotation** | `varimax(pca$rotation[,1:k])` | `varimaxRotation(loadings)` | rotatedLoadings matrix |

**Sign indeterminacy**: PCA components can flip sign. Align by dot-product: if dot(R_col, Carm_col) < 0, flip Carm_col.

**Varimax**: Extract first k=min(3, p) components and apply varimax. Compare rotated loadings after sign alignment.

### TS Harness (`validation/ts-harness/pca-report.ts`)

Report sections:
1. Eigenvalues (element-wise MAE)
2. Variance Explained (element-wise MAE)
3. Loadings (per-component MAE after sign alignment)
4. Scores (first 5 rows, after sign alignment)
5. Varimax Rotation (rotated loadings MAE)

---

## Round 6: Clustering Suite (200 datasets)

**Saqrlab type**: `"clusters"` — `x1..xd` features + `true_cluster` label

### R Reference Script (`validation/r-reference/clustering-ref.R`)

For each seed 1..200:
```r
d <- simulate_data("clusters", seed = i)
features <- d[, grep("^x", names(d))]
k <- max(d$true_cluster)
```

Compute and save:

| Category | R Function | Carm Function | Fields |
|---|---|---|---|
| **K-Means** | `kmeans(features, k, nstart=25, algorithm="Lloyd")` | `runKMeans(data, {k, maxIter, seed})` | centroids, labels, withinSS, totalSS |
| **GMM (VVV)** | `mclust::Mclust(features, G=k, modelNames="VVV")` | `fitGMM(data, {k, model:'VVV'})` | means, logLik, bic, labels |
| **GMM (VVI)** | `mclust::Mclust(features, G=k, modelNames="VVI")` | `fitGMM(data, {k, model:'VVI'})` | means, logLik, bic, labels |
| **GMM (EEI)** | `mclust::Mclust(features, G=k, modelNames="EEI")` | `fitGMM(data, {k, model:'EEI'})` | means, logLik, bic, labels |
| **HAC (ward.D2)** | `hclust(dist(features), "ward.D2")` | `runHierarchical(data, {linkage:'ward'})` | merge order, heights, cutree labels |
| **HAC (complete)** | `hclust(dist(features), "complete")` | `runHierarchical(data, {linkage:'complete'})` | merge order, heights |
| **HAC (single)** | `hclust(dist(features), "single")` | `runHierarchical(data, {linkage:'single'})` | merge order, heights |
| **HAC (average)** | `hclust(dist(features), "average")` | `runHierarchical(data, {linkage:'average'})` | merge order, heights |
| **Silhouette** | `cluster::silhouette(labels, dist(features))` | `silhouetteScores(data, labels)` | mean silhouette, per-point scores (sample) |

**K-Means note**: Use `algorithm="Lloyd"` in R to match Carm's standard Lloyd's algorithm. Use the same `nstart=1` with fixed centers from K-Means++ (or compare final centroids sorted by first dimension, since label assignment differs).

**GMM label matching**: Use permutation-based label alignment (same pattern as FA factor matching but for cluster labels). Compare centroids sorted by first dimension.

**HAC merge comparison**: Compare merge heights (sorted) rather than exact merge order (which can differ for tied distances).

**Silhouette**: Compare mean silhouette score (scalar) and per-point scores for first 20 points.

### TS Harness (`validation/ts-harness/clustering-report.ts`)

Report sections:
1. K-Means (centroid MAE, within-SS, total-SS)
2. GMM — VVV (means MAE, logLik diff, BIC diff)
3. GMM — VVI (means MAE, logLik diff, BIC diff)
4. GMM — EEI (means MAE, logLik diff, BIC diff)
5. HAC — 4 Linkages (merge height MAE)
6. Silhouette (mean score diff, per-point MAE)

---

## Round 7: Factor Analysis Extended Suite (200 datasets)

**Saqrlab type**: `"factor_analysis"` — `x1..xp` items with known factor structure

### R Reference Script (`validation/r-reference/fa-extended-ref.R`)

For each seed 1..200:
```r
d <- simulate_data("factor_analysis", seed = i)
nf <- attr(d, "n_factors")
```

Compute and save:

| Category | R Function | Carm Function | Fields |
|---|---|---|---|
| **PAF + oblimin** | `psych::fa(d, nf, fm="pa", rotate="oblimin")` | `runEFA(data, {extraction:'paf', rotation:'oblimin'})` | loadings, communalities, Phi |
| **PAF + quartimin** | `psych::fa(d, nf, fm="pa", rotate="quartimin")` | `runEFA(data, {extraction:'paf', rotation:'quartimin'})` | loadings, communalities |
| **ML + oblimin** | `psych::fa(d, nf, fm="ml", rotate="oblimin")` | `runEFA(data, {extraction:'ml', rotation:'oblimin'})` | loadings, fit indices |
| **ML + varimax** | `psych::fa(d, nf, fm="ml", rotate="varimax")` | `runEFA(data, {extraction:'ml', rotation:'varimax'})` | loadings |
| **CFA** | `lavaan::cfa(model, data=d)` | `runCFA(data, model)` | fit: chiSq, df, pValue, rmsea, cfi, tli, srmr; loadings |
| **Diagnostics** | `psych::KMO(cor(d))`, `cortest.bartlett(cor(d), n)` | `runFADiagnostics(data)` | kmo, bartlett |

**CFA model**: Build from Saqrlab's `attr(d, "loadings")` — assign items to factors based on which column has the highest loading (simple structure). Feed as `{ F1: [0,1,2], F2: [3,4,5], ... }`.

**Factor matching**: Same exhaustive permutation × sign search from existing `fa-full-report.ts`.

### TS Harness (`validation/ts-harness/fa-extended-report.ts`)

Report sections:
1. PAF + Oblimin (loading MAE, communality MAE)
2. PAF + Quartimin (loading MAE)
3. ML + Oblimin (loading MAE, fit indices)
4. ML + Varimax (loading MAE)
5. CFA (fit indices: χ², RMSEA, CFI, TLI, SRMR; loading MAE)
6. Diagnostics (KMO, Bartlett)

---

## Implementation Order

1. **Round 1** (ttest-suite, 1000) — covers the most Carm functions (20+), immediate breadth
2. **Round 2** (correlation-suite, 1000) — second-most functions, builds confidence
3. **Round 3** (anova-suite, 500) — post-hoc tests, LMM
4. **Round 4** (regression-suite, 500) — regression family
5. **Round 5** (pca-suite, 500) — PCA + varimax
6. **Round 6** (clustering-suite, 200) — GMM/K-Means/HAC (CPU-intensive)
7. **Round 7** (fa-extended, 200) — extends existing FA coverage

Each round: write R script → run → write TS harness → run → verify report → fix any failures → commit.

## R Dependencies

```r
# All must be installed:
library(Saqrlab)      # Data generation
library(jsonlite)     # JSON output (digits=12)
library(effsize)      # Cohen's d, Hedges' g
library(e1071)        # Skewness, kurtosis (type=2)
library(ppcor)        # Partial correlation
library(rstatix)      # Games-Howell
library(dunn.test)    # Dunn's test
library(car)          # VIF
library(lme4)         # LMM
library(mclust)       # GMM
library(cluster)      # Silhouette
library(psych)        # FA, KMO, Bartlett
library(lavaan)       # CFA
library(GPArotation)  # Oblimin, quartimin
```

## Verification

After all 7 rounds:
1. `Rscript validation/r-reference/ttest-ref.R` → JSON written
2. `npx tsx validation/ts-harness/ttest-report.ts` → HTML report
3. Repeat for all 7 rounds
4. `npx tsc --noEmit` — still clean
5. `npx vitest run` — existing 496 tests still pass
6. All HTML reports show green (pass rates ≥ 95% for each metric)
7. Any failures: investigate, fix Carm or adjust tolerance, re-run

**Master command** (runs all 7):
```bash
for suite in ttest correlation anova regression pca clustering fa-extended; do
  Rscript validation/r-reference/${suite}-ref.R && npx tsx validation/ts-harness/${suite}-report.ts
done
```

## Existing Validation (DO NOT MODIFY)
The following files are already validated and committed — do not touch them:
- `validation/r-reference/fa-promax-ref.R` (100 datasets)
- `validation/r-reference/fa-geomin-ref.R` (100 datasets)
- `validation/r-reference/fa-real-ref.R` (real dataset)
- `validation/ts-harness/fa-full-report.ts` (promax + geomin + diagnostics + real)
- `validation/data/fa-crossval-data.json`, `fa-geomin-ref.json`, `fa-geomin-real-ref.json`, `fa-real-crossval-ref.json`
- `validation/reports/fa-crossval-report.html`
- `validation/FA-TECHNICAL-REPORT.md`, `validation/RANDOM-STARTS-REPORT.md`, `validation/GMM-TECHNICAL-REPORT.md`
