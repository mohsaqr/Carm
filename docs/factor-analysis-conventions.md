# Factor Analysis Conventions: Carm's Approach and Differences with lavaan

## Overview

This document describes how Carm implements Exploratory Factor Analysis (EFA), the algorithmic conventions it follows, where it matches lavaan exactly, and where it diverges. It also covers how other packages (psych, GPArotation, factor_analyzer) differ from each other and from lavaan.

The short version: Carm matches R's GPArotation to 1e-6 on 200/200 synthetic datasets. Against lavaan, it matches on most datasets but shows ~2.5% loading differences on some, due entirely to different random number generators exploring different basins in a multi-modal criterion surface.

---

## 1. ML Extraction

### What it does

Maximum Likelihood factor extraction minimizes the Wishart discrepancy function:

```
F_ml = log|Sigma| + tr(Sigma^-1 S) - log|S| - p
```

where S is the observed correlation matrix and Sigma = L L' + Psi is the model-implied correlation matrix (L = loadings, Psi = diagonal uniquenesses).

### Cross-package comparison

| | Carm | lavaan | psych | factanal | factor_analyzer |
|---|---|---|---|---|---|
| Algorithm | Newton-Raphson on uniquenesses | Custom Newton | Calls `factanal()` | Newton-Raphson | EM or minres |
| Input | Correlation matrix | Covariance or correlation | Correlation (converts) | Correlation | Correlation |
| Heywood bounds | Psi clamped to [0.005, 0.995] | Bounded + warnings | Bounded at 0.005 | Bounded | Bounded |
| Chi-square | Bartlett correction | Uncorrected (CFA) / Bartlett (EFA) | Bartlett | Uncorrected | Uncorrected |

### Do these differences matter?

**No, not for loadings.** All implementations maximize the same likelihood. Unrotated loadings agree to < 1e-6 across packages. The only differences are:

- **Bartlett correction**: affects chi-square statistic, not loadings.
- **Heywood handling**: only matters for degenerate solutions (communality > 1).
- **Covariance vs correlation input**: produces different unstandardized estimates but identical standardized loadings and fit indices.

### Sign ambiguity in unrotated loadings

ML extraction produces loadings that are unique up to column sign flips. Factor column j can be multiplied by -1 without changing the likelihood. Different packages may produce different sign conventions for the unrotated solution.

For sign-invariant rotations (varimax, quartimax, oblimin, geomin), this is irrelevant — the rotation compensates. For target rotation, the sign convention must be aligned between the unrotated loadings and the target matrix.

---

## 2. Rotation

Rotation is where packages differ most. The differences fall into three categories: the rotation algorithm, default parameters, and pre/post-processing.

### 2.1 The GPA Algorithm

All modern packages use Gradient Projection Algorithm (GPA) from Bernaards & Jennrich (2005). Carm implements both variants:

**GPFoblq (oblique)** — used for oblimin, quartimin, geomin, target:
```
Parameterization:  L = A * inv(T)'
Gradient:          G = -(L' * Gq * T^-1)'
Projection:        Gp = G - T * diag(colSums(T . G))
Step update:       Column normalization
Factor correlations: Phi = T' * T
```

**GPForth (orthogonal)** — used for quartimax:
```
Parameterization:  L = A * T
Gradient:          G = A' * Gq
Projection:        Gp = G - T * sym(T'G),  sym(M) = (M + M')/2
Step update:       SVD polar decomposition, T = U * V'
Factor correlations: Phi = I (always identity)
```

Both use Armijo backtracking line search with the same step-size logic as R's GPArotation package.

**Carm matches GPArotation exactly** (200/200 synthetic datasets, MAE < 1e-6) when using the same starting matrix (T = I) and parameters.

### 2.2 Rotation Criteria

| Criterion | Formula | Gradient | Type | Reference |
|---|---|---|---|---|
| Varimax | Maximize variance of squared loadings per factor | Kaiser (1958) varimax algorithm | Orthogonal | Kaiser (1958) |
| Quartimax | f = -(1/4) sum(L^4) | Gq = -L^3 | Orthogonal | Neuhaus & Wrigley (1954) |
| Oblimin/Quartimin | f = (1/4) sum_i sum_{j!=m} L^2_ij L^2_im | See Jennrich (2002) | Oblique | Jennrich (2002) |
| Geomin | f = sum_i prod_j (L^2_ij + delta)^(1/k) | See Browne (2001) | Oblique | Browne (2001) |
| Target | f = sum_ij w^2_ij (L_ij - H_ij)^2 | Gq = 2 w^2 (L - H) | Oblique | Browne (1972) |
| Promax | Varimax -> power target -> Procrustes | n/a (closed-form) | Oblique | Hendrickson & White (1964) |

### 2.3 Default Parameters — The Biggest Source of Cross-Package Disagreement

| Parameter | Carm | lavaan | psych | GPArotation | factor_analyzer |
|---|---|---|---|---|---|
| Geomin delta | **0.001** | **0.001** | 0.01 | **0.01** | 0.01 |
| Random starts | 50 | **100** | 1 | **1** | 1 |
| Default rotation | promax | geomin (SEM context) | oblimin | n/a | promax |
| Max iterations | 1000 | 10000 | 1000 | 500 | 1000 |

**Geomin delta** is the most impactful parameter. It controls the penalty for near-zero loadings:
- delta = 0.01 (GPArotation/psych default): smoother criterion surface, fewer local optima
- delta = 0.001 (lavaan/Carm default): sharper criterion surface, more local optima, stronger simple structure

Changing delta from 0.01 to 0.001 can change which basin GPA converges to, producing visibly different loading patterns on the same data.

### 2.4 Pre/Post-Processing — Where lavaan Diverges from Everyone Else

This is the most subtle source of differences. The rotation pipeline differs:

```
GPArotation / psych / factor_analyzer:
  unrotated loadings --> GPA rotation --> done

lavaan:
  unrotated loadings
    --> model-implied standardization (divide by sqrt(h^2 + psi))
    --> GPA rotation (on standardized loadings)
    --> de-standardize (multiply back by sqrt(h^2 + psi))
    --> column reflection (flip columns with negative loading sum)
    --> done

Carm (matching lavaan):
  unrotated loadings
    --> model-implied standardization
    --> GPA rotation
    --> de-standardize
    --> column reflection
    --> done
```

#### Model-implied standardization

lavaan standardizes loadings before rotation:

```
ov_var_i = h^2_i + psi_i    (model-implied variance per variable)
L_std_i  = L_i / sqrt(ov_var_i)
```

For correlation matrix input, ov_var is approximately 1.0 for all variables (typical deviation < 0.001). This makes the standardization nearly a no-op numerically. However, it changes the gradient path through the criterion surface by applying slightly different weights to each variable's loadings, which can shift which basin GPA enters.

Carm implements this standardization to match lavaan's pipeline.

#### Column reflection

After rotation, lavaan flips factor columns where the sum of loadings is negative:

```
for each factor j:
  if sum(L[,j]) < 0:
    L[,j] = -L[,j]
    Phi[j,] = -Phi[j,]  (off-diagonals)
    Phi[,j] = -Phi[,j]  (off-diagonals)
```

This is a sign convention choice that improves interpretability (primary loadings tend to be positive). It does not affect optimality — both sign choices are valid.

GPArotation and psych do not perform column reflection. Carm does, matching lavaan.

**Exception**: Carm disables column reflection for target rotation. Reflection can move loadings away from the target and invalidates the criterion value used for multi-start selection, since the target criterion is not sign-invariant.

---

## 3. Multi-Start Strategy

### Why multiple starts matter

The geomin (and oblimin, target) criterion surfaces have multiple local optima. GPA is a descent algorithm — it always converges to a local minimum, but which one depends on the starting rotation matrix T.

### Cross-package comparison

| | Carm | lavaan | GPArotation | psych |
|---|---|---|---|---|
| Start 0 | T = I (identity) | T = I | T = I | T = I |
| Start 1 | T from QR(A') | T from random | - | - |
| Start 2 | T from varimax | T from random | - | - |
| Starts 3+ | Haar-random orthogonal | Haar-random orthogonal | - | - |
| Total starts | 50 (default) | 100 (default) | 1 | 1 |
| RNG | splitmix32 (deterministic, seed=42) | R Mersenne Twister | n/a | n/a |

### Carm's start strategy

1. **T = I** (always): Matches GPArotation's default. Deterministic, reproducible.
2. **T from QR(A')**: Produces lower-triangular initial loadings, matching lavaan's Cholesky-parameterized starting point.
3. **T from varimax**: Varimax finds an orthogonal simple-structure solution that is often close to the oblique optimum. Uses it as a warm start for GPA.
4. **Random orthogonal matrices**: Haar-distributed (uniform over the orthogonal group) via Modified Gram-Schmidt QR with sign correction.

The best solution (lowest criterion value) across all starts is returned.

### Why Carm and lavaan can find different basins

Even with the same algorithm, delta, standardization, and reflection, different random number generators produce different starting matrices:

```
Carm:   splitmix32(42)       -> start matrices S1, S2, S3, ..., S50
lavaan: Mersenne Twister(?)  -> start matrices T1, T2, T3, ..., T100
```

On a criterion surface with near-degenerate optima (multiple basins with similar criterion values), the specific set of starting points determines which basin is found. This is the source of the ~2.5% loading difference observed on some datasets.

Both solutions are valid — Tucker's congruence coefficient exceeds 0.95 for all factors, confirming psychometric equivalence. Carm sometimes finds a mathematically better (lower) criterion value than lavaan.

---

## 4. Rotation Methods Available

| Method | Carm | lavaan | psych | GPArotation | factor_analyzer |
|---|---|---|---|---|---|
| None | Y | Y | Y | - | Y |
| Varimax | Y | Y | Y | Y | Y |
| Quartimax | Y | Y | Y | Y | Y |
| Oblimin/Quartimin | Y | Y | Y | Y | Y |
| Promax | Y | Y | Y | Y | Y |
| Geomin | Y | Y | Y | Y | N |
| Target (oblique) | Y | Y | Y | Y | N |
| Target (orthogonal) | N | Y | Y | Y | N |
| Equamax | N | N | Y | Y | N |
| Simplimax | N | N | Y | Y | N |
| Bentler's invariant | N | N | N | Y | N |

---

## 5. CFA Parameterization

Carm's CFA uses a specific parameterization that differs from lavaan's default:

| | Carm | lavaan (default) | lavaan (std.lv=TRUE) |
|---|---|---|---|
| Input | Correlation matrix | Covariance matrix | Covariance matrix |
| Factor identification | Factor variance = 1 | First loading = 1 (marker variable) | Factor variance = 1 |
| Estimated loadings | All freely estimated | First per factor fixed to 1 | All freely estimated |

**Effect on results**: These parameterizations identify the same model. The standardized solution (`std.all` loadings), fit indices (chi-square, RMSEA, CFI, TLI, SRMR), and model-implied correlations are identical regardless of parameterization.

The only differences are in the unstandardized estimates and their standard errors, which reflect the arbitrary scale of the latent variables. For most practical purposes (reporting standardized loadings and fit indices), the choice is irrelevant.

To get exact unstandardized agreement with lavaan's default, one would need to:
1. Fit on the covariance matrix (not correlation)
2. Use the marker-variable constraint (first loading fixed to 1)

This is not currently implemented because Carm reports standardized results, making the parameterization choice transparent to the user.

---

## 6. Numerical Equivalence Summary

### Where Carm matches exactly (< 1e-6)

- **GPArotation**: 200/200 synthetic datasets, geomin with delta=0.01 and delta=0.001
- **GPArotation quartimax**: 100/100 synthetic datasets via GPForth
- **GPArotation target**: 93/100 synthetic datasets (7 failures due to GPA basin differences)
- **psych promax**: Matches psych::fa(rotate="promax") via kaiser normalization
- **All fit indices**: chi-square, RMSEA, CFI, TLI, SRMR match R to 4+ decimal places

### Where Carm diverges from lavaan (~2-3% on some datasets)

- **Cause**: Different RNG producing different starting matrices -> different basins found
- **Magnitude**: Mean absolute loading difference ~2.5% on affected factors
- **Validation**: Tucker's congruence > 0.95 for all factors (psychometric equivalence)
- **Not fixable without**: Reimplementing R's Mersenne Twister with identical seed state

### What would NOT help close the gap

- Changing CFA parameterization (std.lv=TRUE vs marker variable) — affects unstandardized CFA estimates only, not EFA rotation
- Using covariance instead of correlation matrix — same standardized solution
- Changing Bartlett correction — affects chi-square only, not loadings
- Column reflection convention — sign choice, not optimality

---

## 7. References

- Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and software for arbitrary rotation criteria in factor analysis. *Educational and Psychological Measurement*, 65(5), 676-696.
- Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis. *Multivariate Behavioral Research*, 36(1), 111-150.
- Browne, M. W. (1972). Oblique rotation to a partially specified target. *British Journal of Mathematical and Statistical Psychology*, 25(2), 207-212.
- Jennrich, R. I. (2002). A simple general method for oblique rotation. *Psychometrika*, 67(1), 7-27.
- Kaiser, H. F. (1958). The varimax criterion for analytic rotation in factor analysis. *Psychometrika*, 23(3), 187-200.
- Neuhaus, J. O., & Wrigley, C. (1954). The quartimax method: An analytic approach to orthogonal simple structure. *British Journal of Statistical Psychology*, 7(2), 81-91.
