# Factor Analysis in Carm: Technical Report

## A Complete Implementation of EFA, CFA, and Psychometric Diagnostics in TypeScript

**Date:** 2026-02-26
**Version:** Carm 1.0
**Module:** `src/stats/factor-analysis.ts` (1,923 lines)
**Dependencies:** Zero external math libraries — all linear algebra, distributions, and optimization from scratch
**Validation:** 100/100 synthetic promax, 100/100 synthetic geomin, 6/6 real dataset, 100/100 diagnostics
**Cross-validated against:** R psych 2.4.6, R GPArotation 2024.3-1, R lavaan 0.6-19, R factanal (stats)

---

## Table of Contents

1. [Architecture and Design Principles](#1-architecture-and-design-principles)
2. [Extraction: Maximum Likelihood](#2-extraction-maximum-likelihood)
3. [Extraction: Principal Axis Factoring](#3-extraction-principal-axis-factoring)
4. [Rotation: Varimax (Orthogonal)](#4-rotation-varimax-orthogonal)
5. [Rotation: Promax (Oblique)](#5-rotation-promax-oblique)
6. [Rotation: The GPFoblq Engine](#6-rotation-the-gpfoblq-engine)
7. [Rotation Criteria: Geomin, Oblimin, Quartimin](#7-rotation-criteria-geomin-oblimin-quartimin)
8. [Random Starting Matrices](#8-random-starting-matrices)
9. [Confirmatory Factor Analysis (CFA)](#9-confirmatory-factor-analysis-cfa)
10. [Fit Indices](#10-fit-indices)
11. [Psychometric Diagnostics](#11-psychometric-diagnostics)
12. [Numerical Infrastructure](#12-numerical-infrastructure)
13. [Cross-Validation Methodology](#13-cross-validation-methodology)
14. [Cross-Validation Results](#14-cross-validation-results)
15. [Detailed Loading Comparisons](#15-detailed-loading-comparisons)
16. [Public API Reference](#16-public-api-reference)
17. [References](#17-references)

---

## 1. Architecture and Design Principles

### 1.1 Zero-Dependency Mathematics

Carm implements all mathematical operations from scratch in TypeScript:

- **Matrix algebra** (`core/matrix.ts`): inverse, pseudo-inverse, eigen-decomposition (QR algorithm with Wilkinson shifts), Cholesky, determinant, trace, SVD
- **Distributions** (`core/math.ts`): chi-square CDF/inverse, normal CDF/inverse, F distribution, t distribution
- **Optimization**: Nelder-Mead simplex, gradient descent with Armijo line search, modified Gram-Schmidt QR

No jStat, no simple-statistics, no numeric.js. This eliminates version conflicts, bundle size concerns, and behavioral surprises from third-party code.

### 1.2 Deterministic Reproducibility

Every stochastic operation uses a seeded splitmix32 PRNG (default seed = 42). Given identical inputs and options, `runEFA` produces bit-identical output across runs, platforms, and JavaScript engines.

### 1.3 Layer Discipline

```
core/         Pure math: Matrix, distributions, PRNG
  ↓
stats/        Pure statistics: factor-analysis.ts, clustering.ts, ...
  ↓
viz/          D3-based visualization (consumes stats, never reverse)
```

`factor-analysis.ts` lives in `stats/` — it has no DOM access, no D3, no side effects. Every function is a pure transformation from data to results.

### 1.4 Type Safety

TypeScript strict mode with `noUncheckedIndexedAccess`. All array accesses use `arr[i]!` after validation. All return types are `readonly`. All options are `Readonly<>`. The result interfaces enforce immutability:

```typescript
interface FAResult {
  readonly loadings: readonly (readonly number[])[]
  readonly standardizedLoadings: readonly (readonly number[])[]
  readonly uniqueness: readonly number[]
  readonly communalities: readonly number[]
  readonly factorCorrelations: readonly (readonly number[])[]
  readonly fit: FactorFit
  readonly eigenvalues: readonly number[]
  readonly nFactors: number
  readonly rotation: string
  readonly extraction: string
  readonly variableNames: readonly string[]
  readonly factorNames: readonly string[]
  readonly formatted: string    // APA-formatted result string
}
```

---

## 2. Extraction: Maximum Likelihood

### 2.1 The ML Discrepancy Function

Given observed correlation matrix **R** (p × p) and model-implied covariance **Σ** = **LL**ᵀ + **Ψ**, the ML discrepancy is:

```
F_ML(Ψ) = log|Σ| + tr(Σ⁻¹R) − log|R| − p
```

where **L** is the p × k loading matrix and **Ψ** = diag(ψ₁, ..., ψ_p) is the diagonal uniqueness matrix. The parameters are the p uniquenesses; loadings are derived from the spectral decomposition of the reduced correlation matrix.

### 2.2 Two-Phase Optimization

**Phase 1: Jöreskog Gradient Descent** (Jöreskog, 1967)

This provides a fast warm start. For given uniquenesses **Ψ**:

1. Form the scaled matrix:
   ```
   R* = Ψ^(−1/2) R Ψ^(−1/2)
   ```

2. Eigendecompose R* = VΛVᵀ, keep top k components:
   ```
   L = Ψ^(1/2) × V_k × diag(√max(λⱼ − 1, 0))
   ```

3. Compute gradient:
   ```
   ∂F/∂ψᵢ = (Σ⁻¹(Σ − R)Σ⁻¹)ᵢᵢ
   ```

4. Update with momentum (β = 0.8) and cosine-annealed learning rate (0.02 → 0.001):
   ```
   vᵢ ← β vᵢ + (1 − β) ∂F/∂ψᵢ
   ψᵢ ← ψᵢ − α vᵢ
   ```

5. Clamp: ψᵢ ∈ [0.005, 0.995] (prevents Heywood cases and singular Ψ).

Initialization uses R's `factanal()` convention:
```
ψᵢ = (1 − 0.5k/p) / R⁻¹ᵢᵢ
```

**Phase 2: Nelder-Mead Polish**

The Jöreskog solution provides an excellent starting simplex. Nelder-Mead then minimizes R's concentrated negative log-likelihood:

```
f(Ψ) = −Σⱼ₌ₖ₊₁ᵖ (log λⱼ − λⱼ) + k − p
```

where λⱼ are eigenvalues of Ψ^(−1/2) R Ψ^(−1/2), and only the p − k smallest eigenvalues contribute (the k largest are used for loading extraction).

Parameters: up to 5000p iterations, tolerance 1e-10, bounds [0.005, 0.995] with penalty 1000 for violations.

### 2.3 Canonical Sign Convention

After extraction, each factor column is sign-corrected to match LAPACK's `dsyev` convention:

```
for f = 1, ..., k:
    i* = argmax_i |L[i,f]|
    if L[i*,f] < 0:
        L[:,f] ← −L[:,f]
```

This ensures the unrotated loadings have a deterministic orientation, so that rotation algorithms (starting from T = I) begin from the same point as R's `factanal()` and `eigen()`.

---

## 3. Extraction: Principal Axis Factoring

### 3.1 Algorithm (Gorsuch, 1983)

PAF is an iterative eigendecomposition of the reduced correlation matrix:

1. **Initialize communalities** from squared multiple correlations:
   ```
   h²ᵢ = max(0.01, 1 − 1/R⁻¹ᵢᵢ)
   ```
   Falls back to h²ᵢ = 0.5 if R is singular.

2. **Iterate** until max |Δh²| < tolerance:
   - Replace R's diagonal with current h²: R* = R, R*ᵢᵢ = h²ᵢ
   - Eigendecompose R*: top k eigenvalues λ and eigenvectors V
   - Loadings: Lᵢⱼ = Vᵢⱼ √λⱼ
   - Update: h²ᵢ = Σⱼ L²ᵢⱼ

3. Apply canonical sign convention (same as ML).

### 3.2 When to Use PAF vs ML

| Criterion | ML | PAF |
|-----------|-----|-----|
| Distributional assumption | Multivariate normality | None |
| Fit indices | Full suite (χ², RMSEA, CFI, TLI) | Limited |
| Heywood cases | Can occur; clamped | Less common |
| Small samples | May not converge | More robust |
| Cross-validation target | R `factanal()` | R `psych::fa(fm="pa")` |

---

## 4. Rotation: Varimax (Orthogonal)

### 4.1 The Varimax Criterion (Kaiser, 1958)

Varimax maximizes the variance of squared loadings within each factor:

```
V(Λ) = Σⱼ₌₁ᵏ [ (1/p) Σᵢ λ⁴ᵢⱼ − ((1/p) Σᵢ λ²ᵢⱼ)² ]
```

Equivalently, it seeks the orthogonal rotation T that makes each column of ΛT have loadings that are either large or near-zero (not intermediate).

### 4.2 Implementation (R's `stats::varimax` Algorithm)

1. **Kaiser normalization:** Divide each row by its communality norm:
   ```
   scᵢ = √(Σⱼ L²ᵢⱼ)
   xᵢⱼ = Lᵢⱼ / scᵢ
   ```

2. **Iterative SVD-based rotation:**
   For each iteration:
   - Apply current rotation: z = x × T
   - Compute varimax criterion matrix:
     ```
     Bᵣ꜀ = Σᵢ xᵢᵣ (z³ᵢ꜀ − zᵢ꜀ × (Σᵢ z²ᵢ꜀)/p)
     ```
   - Polar decomposition of B via eigendecomposition of BᵀB:
     ```
     BᵀB = VΣ²Vᵀ
     T_new = B V Σ⁻¹ Vᵀ
     ```
   - Convergence: sum of singular values stabilizes.

3. **Denormalize:** Multiply rows by scᵢ.

This matches R's `stats::varimax()` exactly — not the textbook algorithm, but R's specific SVD-via-eigendecomposition implementation.

---

## 5. Rotation: Promax (Oblique)

### 5.1 Algorithm (Hendrickson & White, 1964; psych::kaiser)

Promax converts an orthogonal varimax solution into an oblique one by fitting to a simplified target:

1. **Outer Kaiser normalization** (critical for matching `psych::fa(rotate="promax")`):
   ```
   h²ᵢ = Σⱼ L²ᵢⱼ
   L_normᵢⱼ = Lᵢⱼ / √h²ᵢ
   ```
   This is applied **before** the varimax step, not after.

2. **Varimax** on normalized loadings: V = rotateVarimax(L_norm)

3. **Power target** (κ = 4):
   ```
   Hᵢⱼ = sign(Vᵢⱼ) × |Vᵢⱼ|^κ
   ```

4. **Least-squares Procrustes fit:**
   ```
   U = (VᵀV)⁻¹ VᵀH
   ```

5. **Rescaling** for unit factor variances:
   ```
   dⱼ = diag((UᵀU)⁻¹)ⱼⱼ
   U ← U × diag(√d)
   ```

6. **Oblique loadings:** L_rot = V × U, then denormalize by √h²ᵢ.

7. **Factor correlations:** Φ = T⁻¹(T⁻¹)ᵀ where T = T_varimax × U.

### 5.2 Critical Detail: Outer vs Inner Kaiser

R's `stats::promax()` and `psych::fa(rotate="promax")` differ in Kaiser normalization. Carm implements `psych::fa()`'s approach (outer Kaiser — normalize before promax), which produces different results from `stats::promax()` (which does its own internal normalization). Cross-validation confirms exact match with `psych::fa()`.

---

## 6. Rotation: The GPFoblq Engine

### 6.1 The Oblique Rotation Problem

Given unrotated loadings **A** (p × k), find non-singular **T** (k × k) minimizing:

```
min_T  f(A (T⁻¹)ᵀ)
```

subject to each column of **T** having unit Euclidean norm. The rotated loadings are:

```
L = A (T⁻¹)ᵀ
```

and the factor correlation matrix is:

```
Φ = TᵀT
```

Unlike orthogonal rotation (where T is restricted to the orthogonal group), oblique rotation allows correlated factors — the columns of T need not be orthogonal.

### 6.2 The GPFoblq Algorithm (Bernaards & Jennrich, 2005)

Carm implements GPFoblq exactly as in R's GPArotation package:

**Gradient computation:**

The gradient of f with respect to T is:

```
G = −(Lᵀ Gq T⁻¹)ᵀ
```

where Gq = ∂f/∂L is the criterion gradient with respect to the rotated loadings.

**Critical:** This formula uses L (the rotated loadings) and T⁻¹, not A and (T⁻¹)ᵀ. Getting this wrong is the most common implementation bug.

**Oblique projection:**

Project G onto the tangent space of the constraint manifold (unit-norm columns of T):

```
Gp = G − T × diag(diag(TᵀG))
```

This removes the component that would change column norms.

**Armijo line search:**

```
α ← 2α
for i = 0, ..., 10:
    X ← T − α × Gp
    Normalize columns of X
    L_new ← A (X⁻¹)ᵀ
    if f(L) − f(L_new) > 0.5 ‖Gp‖²_F α:  break
    α ← α/2
```

**Unconditional update** (matching R exactly): After the line search, T, L, f, and G are always updated — even if the Armijo condition was never satisfied. This is a deliberate design choice in GPArotation that prevents the algorithm from stalling.

**Convergence:** ‖Gp‖_F < ε (default ε = 10⁻⁶).

**Output:** Rotated loadings L, transformation T, factor correlations Φ = TᵀT, final criterion value f.

---

## 7. Rotation Criteria: Geomin, Oblimin, Quartimin

### 7.1 Geomin (Yates, 1987; Browne, 2001)

The geomin criterion encourages simple structure by penalizing the geometric mean of squared loadings per variable:

```
f(L) = Σᵢ₌₁ᵖ  [ Πⱼ₌₁ᵏ (λ²ᵢⱼ + δ) ]^(1/k)
```

The gradient:

```
∂f/∂λᵢⱼ = (2/k) × λᵢⱼ / (λ²ᵢⱼ + δ) × [ Πₘ (λ²ᵢₘ + δ) ]^(1/k)
```

**Log-space computation:** To prevent underflow with small loadings and small δ, the product is computed as:

```
proᵢ = exp( (1/k) Σⱼ log(λ²ᵢⱼ + δ) )
```

**The δ parameter:** The regularization constant δ controls the criterion surface topology:

| δ | Source | Surface | Starts needed |
|---|--------|---------|---------------|
| 0.01 | GPArotation default | Smooth, typically unimodal | 1 |
| 0.001 | lavaan default | Rougher, multiple local minima | 10–50 |

This is confirmed empirically: with δ = 0.01, single-start Carm matches lavaan perfectly. With δ = 0.001, 50 random starts are needed.

### 7.2 Oblimin (Jennrich, 2002)

The oblimin family of criteria is parameterized by γ:

```
f(L) = (1/4) Σᵢ Σⱼ≠ₘ λ²ᵢⱼ λ²ᵢₘ − (γ/4p) Σⱼ Σₘ≠ⱼ (Σᵢ λ²ᵢⱼ)(Σᵢ λ²ᵢₘ)
```

The gradient:

```
∂f/∂λᵢⱼ = λᵢⱼ ( Σₘ≠ⱼ λ²ᵢₘ − (γ/p) Σₘ λ²ᵢₘ )
```

### 7.3 Quartimin

Setting γ = 0 in oblimin gives quartimin, which minimizes:

```
f(L) = (1/4) Σᵢ Σⱼ≠ₘ λ²ᵢⱼ λ²ᵢₘ
```

Quartimin seeks to minimize the sum of products of squared loadings across factor pairs — variables should load on as few factors as possible.

### 7.4 Criterion Dispatch

All three criteria share the same GPFoblq engine. The criterion function is passed as a parameter:

```typescript
// Oblimin/Quartimin
const criterionFn = (L: number[][]) => criterionOblimin(L, gamma)

// Geomin
const criterionFn = (L: number[][]) => criterionGeomin(L, delta)

// Same optimizer for all
gpfOblqWithRandomStarts(loadings, criterionFn, maxIter, tol, k, randomStarts, seed)
```

---

## 8. Random Starting Matrices

### 8.1 The Problem: Local Optima in Oblique Rotation

The geomin and oblimin criterion surfaces are non-convex. Starting from T = I (the identity), gradient descent converges to whichever local minimum the identity falls into. For datasets where the criterion surface has competing basins of attraction, this single-start approach can find suboptimal solutions.

Empirical evidence from two real-world datasets:

| Dataset | n × p | k | Single-start MAE vs lavaan | 50-start MAE vs lavaan |
|---------|-------|---|---------------------------|------------------------|
| Teacher Burnout | 876 × 23 | 4 | 0.033 | 0.000002 |
| rraw (5-subscale) | 525 × 31 | 5 | 0.043 | 0.000001 |

The errors are not uniformly distributed — they concentrate on specific factors where the identity start converges to the wrong basin.

### 8.2 Haar-Distributed Random Orthogonal Matrices

Random starting matrices are drawn from the Haar measure on O(k) — the unique distribution on orthogonal matrices that is invariant under multiplication by any orthogonal matrix. This is the "uniform distribution" on rotations.

**Algorithm** (matching GPArotation::Random.Start() and lavaan::lav_matrix_rotate_gen()):

1. Fill k × k matrix M with independent N(0,1) entries.

2. Modified Gram-Schmidt QR decomposition:
   ```
   for j = 1, ..., k:
       for prev = 1, ..., j−1:
           Q[:,j] ← Q[:,j] − ⟨Q[:,j], Q[:,prev]⟩ Q[:,prev]
       R[j,j] ← ‖Q[:,j]‖
       Q[:,j] ← Q[:,j] / R[j,j]
   ```

3. Sign correction (critical for Haar uniformity):
   ```
   for j = 1, ..., k:
       if R[j,j] < 0:  Q[:,j] ← −Q[:,j]
   ```

**Mathematical justification** (Stewart, 1980; Mezzadri, 2007): If M has i.i.d. N(0,1) entries and M = QR is the QR decomposition with R having positive diagonal, then Q follows the Haar measure on O(k). The sign correction maps the standard QR (which may have negative R diagonal entries) to the unique canonical form.

### 8.3 The Multi-Start Strategy

```
Algorithm: GPFoblq with Multiple Random Starts

Input:  A (p×k unrotated loadings), criterion f, S starts, seed
Output: (L*, Φ*) with lowest criterion value

1.  best ← GPFoblq(A, f, T₀ = Iₖ)              // deterministic start
2.  rng ← splitmix32(seed)
3.  for s = 1, ..., S−1:
4.      T_rand ← randomOrthogonalMatrix(k, rng)  // Haar-distributed
5.      result ← GPFoblq(A, f, T₀ = T_rand)
6.      if result.f < best.f:  best ← result
7.  return (best.L, best.Φ)
```

**Design decisions:**

- **Start 0 = identity always:** With `randomStarts: 1`, output matches GPArotation exactly.
- **Selection by criterion value f**, not loading similarity — this is the mathematically correct measure of rotation quality.
- **Default: 50 starts.** Sufficient for both test datasets; lavaan uses 30.

### 8.4 Empirical Calibration

**Teacher Burnout (876 × 23, k = 4, geomin δ = 0.001):**

| Starts | Time | All 92 cells within 0.001 of lavaan |
|--------|------|--------------------------------------|
| 1 | 0.4s | No — DE items off by 0.12–0.21 |
| 10 | 0.5s | Yes (92/92) |
| 50 | 1.2s | Yes |

**rraw (525 × 31, k = 5, geomin δ = 0.001):**

| Starts | Time | All 155 cells within 0.001 of lavaan |
|--------|------|---------------------------------------|
| 1 | 1.4s | No — FSI items off by 0.11–0.15 |
| 30 | 2.4s | No — ER items still off by 0.07 |
| 50 | 3.0s | Yes (155/155) |
| 100 | 4.9s | Yes |

Runtime scales linearly at ~0.03s per start per dataset.

---

## 9. Confirmatory Factor Analysis (CFA)

### 9.1 Model Specification

CFA models are specified as a mapping from factor names to item indices:

```typescript
const model = {
  TE:  [0, 1, 2, 3, 4],     // Teacher Efficacy → items 0-4
  EE:  [5, 6, 7, 8, 9],     // Emotional Exhaustion → items 5-9
  DE:  [10, 11, 12],         // Depersonalization → items 10-12
  RPA: [13, 14, 15, 16, 17]  // Reduced Personal Accomplishment → items 13-17
}
```

### 9.2 ML Estimation

The implied covariance matrix is:

```
Σ = ΛΦΛᵀ + Θ
```

where Λ is loadings (p × k, zero for non-specified paths), Φ is factor correlations (k × k), and Θ = diag(ψ₁, ..., ψ_p).

Optimization minimizes:

```
F_ML = log|Σ| + tr(Σ⁻¹S) − log|S| − p
```

via gradient descent with Armijo backtracking:

**Gradients:**

Let Δ = Σ⁻¹(Σ − S)Σ⁻¹. Then:

```
∂F/∂Λ = 2ΔΛΦ           (only for specified loadings)
∂F/∂Φᵢⱼ = ΛᵀΔΛᵢⱼ      (off-diagonal only; diagonal = 1)
∂F/∂ψᵢ = Δᵢᵢ
```

**Bounds enforcement:** ψᵢ ∈ [0.001, 0.995], Φᵢⱼ ∈ [−0.99, 0.99].

**Convergence:** max|gradient| < tolerance.

### 9.3 Standard Errors

Standard errors are computed via numerical differentiation of the Hessian:

1. Collect all free parameters into vector θ (specified loadings + uniquenesses + off-diagonal Φ).

2. Compute Hessian via central finite differences (h = 10⁻⁴):
   ```
   H_ii = (f(θ + heᵢ) − 2f(θ) + f(θ − heᵢ)) / h²
   H_ij = (f(θ + heᵢ + heⱼ) − f(θ + heᵢ − heⱼ) − f(θ − heᵢ + heⱼ) + f(θ − heᵢ − heⱼ)) / (4h²)
   ```

3. Information matrix: I = ((n−1)/2) × H

4. Standard errors: SE = √diag(I⁻¹)

5. z-statistics: z = estimate / SE, p-value from normal distribution.

### 9.4 CFA Output

```typescript
interface CFAResult extends FAResult {
  readonly parameterEstimates: {
    readonly loadings: readonly (readonly ParameterEstimate[])[]
    readonly uniquenesses: readonly ParameterEstimate[]
    readonly factorCovariances: readonly (readonly ParameterEstimate[])[]
  }
  readonly model: Readonly<Record<string, readonly number[]>>
}

interface ParameterEstimate {
  readonly estimate: number
  readonly se: number
  readonly z: number
  readonly pValue: number
  readonly stdAll: number    // STDYX standardized estimate
}
```

---

## 10. Fit Indices

### 10.1 Chi-Square Test

For EFA (with Bartlett, 1950, correction matching `psych::fa()`):

```
χ² = [n − 1 − (2p + 5)/6 − 2k/3] × F_ML
```

For CFA (no Bartlett correction, matching lavaan):

```
χ² = (n − 1) × F_ML
```

Degrees of freedom:

```
df = p(p+1)/2 − nFreeParams
```

where nFreeParams = pk + p (EFA with oblique rotation) or model-implied (CFA).

### 10.2 RMSEA (Browne & Cudeck, 1993)

```
RMSEA = √( max(χ² − df, 0) / (df × (n − 1)) )
```

90% confidence interval via Steiger (1990):

```
NCP_lower = max(χ² − df − 1.645√(2df), 0)
NCP_upper = max(χ² − df + 1.645√(2df), 0)
RMSEA_CI = [√(NCP_lower / (df(n−1))), √(NCP_upper / (df(n−1)))]
```

### 10.3 CFI (Bentler, 1990)

```
CFI = 1 − max(0, χ² − df) / max(χ² − df, χ²₀ − df₀)
```

where χ²₀ and df₀ are from the null (independence) model. Clamped to [0, 1].

### 10.4 TLI / NNFI (Tucker & Lewis, 1973)

```
TLI = (χ²₀/df₀ − χ²/df) / (χ²₀/df₀ − 1)
```

Not clamped — can exceed 1 for well-fitting models.

### 10.5 SRMR

```
SRMR = √( (2 / p(p+1)) × Σᵢ≥ⱼ (rᵢⱼ − σ̂ᵢⱼ)² )
```

where rᵢⱼ and σ̂ᵢⱼ are standardized (correlation-scale) observed and implied values.

### 10.6 Information Criteria

```
AIC = χ² + 2 × nFreeParams
BIC = χ² + nFreeParams × log(n)
```

---

## 11. Psychometric Diagnostics

### 11.1 Kaiser-Meyer-Olkin (KMO) Sampling Adequacy

KMO measures whether partial correlations are small relative to raw correlations — a prerequisite for factor analysis:

```
KMO = Σᵢ<ⱼ r²ᵢⱼ / (Σᵢ<ⱼ r²ᵢⱼ + Σᵢ<ⱼ a²ᵢⱼ)
```

where rᵢⱼ are correlations and aᵢⱼ are anti-image correlations:

```
aᵢⱼ = −R⁻¹ᵢⱼ / √(R⁻¹ᵢᵢ R⁻¹ⱼⱼ)
```

Per-item MSA is computed analogously for each variable's row.

KMO interpretation: ≥ 0.90 marvelous, ≥ 0.80 meritorious, ≥ 0.70 middling, ≥ 0.60 mediocre, < 0.50 unacceptable.

### 11.2 Bartlett's Test of Sphericity

Tests H₀: R = I (correlation matrix is identity — variables are uncorrelated):

```
χ² = −(n − 1 − (2p + 5)/6) × log|R|
df = p(p − 1)/2
```

A significant result (p < 0.05) indicates correlations exist and factor analysis is appropriate.

### 11.3 Parallel Analysis (Horn, 1965)

Monte Carlo method for determining the number of factors:

1. Generate `iterations` (default 100) random normal datasets of size n × p.
2. Compute eigenvalues of each random correlation matrix.
3. For each eigenvalue position, compute the 95th percentile across iterations.
4. Number of factors = count of observed eigenvalues exceeding their simulated 95th percentile.

### 11.4 Velicer's MAP (Minimum Average Partial, 1976)

1. For k = 0, 1, ..., p−1 components:
   - Extract k components from R
   - Compute partial correlation matrix (R with k components removed)
   - Average the squared off-diagonal elements
2. Return k with the smallest average squared partial correlation.

MAP tends to recommend fewer factors than parallel analysis — it identifies the point where extracting additional factors no longer reduces the average partial correlation.

---

## 12. Numerical Infrastructure

### 12.1 splitmix32 PRNG

```
state ← state + 0x9E3779B9         (golden ratio constant)
z ← state
z ← (z ⊕ (z >> 16)) × 0x45D9F3B
z ← (z ⊕ (z >> 16)) × 0x45D9F3B
z ← z ⊕ (z >> 16)
return z / 2³²
```

Period: 2³². Quality: passes BigCrush. Speed: ~1 ns per call in V8.

Normal deviates via Box-Muller:
```
z = √(−2 ln u₁) × cos(2π u₂),  u₁, u₂ ~ Uniform(0,1)
```

### 12.2 Matrix Operations (core/matrix.ts)

| Operation | Algorithm | Complexity |
|-----------|-----------|------------|
| Multiply | Standard triple loop | O(n³) |
| Inverse | Gauss-Jordan elimination | O(n³) |
| Pseudo-inverse | SVD-based | O(n³) |
| Eigendecomposition | QR algorithm + Wilkinson shifts | O(n³) per iteration |
| Cholesky | Standard lower-triangular | O(n³/3) |
| Determinant | Via Cholesky or eigen | O(n³) |
| Log-determinant | log of eigenvalues or Cholesky | O(n³) |

All operations use `Float64Array` for accumulator loops to minimize floating-point drift.

### 12.3 Robustness

- **Near-singular matrices:** `inverse()` throws → caught → `pseudoInverse()` used as fallback.
- **Heywood cases:** Uniquenesses clamped to [0.005, 0.995] throughout optimization.
- **Numerical log-determinant:** If Cholesky fails (non-positive-definite), falls back to eigenvalue-based computation with negative eigenvalues treated as ε.
- **Input validation:** Every public function checks minimum observations (≥3), minimum variables (≥2), nFactors < p.

---

## 13. Cross-Validation Methodology

### 13.1 Protocol

```
┌──────────────────────────────────────────┐
│ R generates ground truth                 │
│   psych::fa() + GPArotation::geominQ()   │
│   → JSON with 10-digit precision         │
└──────────────┬───────────────────────────┘
               │
┌──────────────▼───────────────────────────┐
│ Carm runs identical analysis             │
│   runEFA(data, { same options })         │
│   → loading matrix, fit indices          │
└──────────────┬───────────────────────────┘
               │
┌──────────────▼───────────────────────────┐
│ Factor matching                          │
│   Exhaustive search: k! × 2^k combos    │
│   Minimize total absolute error          │
└──────────────┬───────────────────────────┘
               │
┌──────────────▼───────────────────────────┐
│ Comparison metrics                       │
│   MAE, max error, pass/fail threshold    │
│   → HTML report with color-coded tables  │
└──────────────────────────────────────────┘
```

### 13.2 Factor Matching Algorithm

Factor analysis has inherent indeterminacy: factors can be permuted and sign-flipped without changing the model. To compare loadings across software, we search over all k! permutations and 2^k sign patterns:

```
for each permutation π of {1,...,k}:
    for each factor j, set sign_j = sign(Σᵢ R_ref[i,j] × L_carm[i,π(j)])
    MAE = (1/pk) Σᵢ Σⱼ |R_ref[i,j] − sign_j × L_carm[i,π(j)]|
    keep if MAE < best
```

For k ≤ 6 (120 permutations × 64 signs = 7,680 combinations), this is exhaustive and fast. For k > 6, Hungarian algorithm with sign search would be needed.

### 13.3 Thresholds

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| Eigenvalues | 1e-8 | Deterministic — should be machine precision |
| KMO | 1e-8 | Deterministic |
| Bartlett χ² | 0.01 | Deterministic with minor correction differences |
| Loading MAE | 0.05 | Accounts for different optimizers + numerical paths |
| Communality MAE | 0.05 | Derived from loadings |
| χ² | 10.0 | Different Bartlett corrections |
| RMSEA | 0.02 | Derived from χ² |
| CFI | 0.02 | Sensitive to null model differences |
| TLI | 0.02 | Sensitive to null model differences |
| SRMR | 0.01 | Direct residual computation |

### 13.4 Synthetic Data Generation

100 datasets generated in R with controlled factor structure:

```r
# For each dataset:
n ~ Uniform({100, 150, 200, 300, 500})
k ~ Uniform({2, 3, 4, 5})
p ~ k × Uniform({3, 4, 5})
loading_strength ~ Uniform(0.5, 0.9)

# True loading matrix: simple structure with noise
L_true[i,j] = loading_strength  if item i belongs to factor j
L_true[i,j] ~ Uniform(0, 0.15)  otherwise

# True uniquenesses:
Ψ_true[i] ~ Uniform(0.2, 0.6)

# Generate data:
Σ = L_true Φ L_true' + Ψ_true
X ~ MVN(0, Σ), n × p
```

---

## 14. Cross-Validation Results

### 14.1 Summary

| Module | Rotation | Datasets | Pass Rate |
|--------|----------|----------|-----------|
| EFA ML | Promax | 100 synthetic | **100/100** |
| EFA ML | Geomin (δ=0.01) | 100 synthetic | **100/100** |
| EFA ML | Promax | Real 525×31 (k=3,4,5,6) | **4/4** |
| EFA ML | Geomin | Real 525×31 (k=3,5) | **2/2** |
| Diagnostics | — | 100 synthetic | **100/100** |

### 14.2 Promax Aggregate Statistics (100 datasets)

| Metric | Mean | Median | 95th pct | Max | Threshold | Pass |
|--------|------|--------|----------|-----|-----------|------|
| Eigenvalue MAE | 2.1e-14 | 1.5e-14 | 5.3e-14 | 1.7e-13 | 1e-8 | 100% |
| KMO |Δ| | 8.2e-15 | 5.6e-15 | 2.1e-14 | 5.7e-14 | 1e-8 | 100% |
| Bartlett χ² |Δ| | 3.0e-12 | 1.4e-12 | 1.2e-11 | 4.4e-11 | 0.01 | 100% |
| Loading MAE | 2.9e-5 | 1.7e-5 | 8.1e-5 | 1.8e-4 | 0.05 | 100% |
| Communality MAE | 1.7e-5 | 1.0e-5 | 5.0e-5 | 1.1e-4 | 0.05 | 100% |
| RMSEA |Δ| | 5.1e-4 | 3.2e-4 | 1.5e-3 | 2.8e-3 | 0.02 | 100% |
| CFI |Δ| | 2.6e-10 | 1.1e-10 | 7.8e-10 | 1.8e-9 | 0.02 | 100% |
| TLI |Δ| | 7.2e-4 | 4.4e-4 | 2.1e-3 | 3.7e-3 | 0.02 | 100% |
| SRMR |Δ| | 1.3e-3 | 9.8e-4 | 3.5e-3 | 6.7e-3 | 0.01 | 100% |

Eigenvalues, KMO, and Bartlett match at machine precision (10⁻¹⁴). Loading MAE peaks at 1.8×10⁻⁴, three orders of magnitude below the 0.05 threshold.

### 14.3 Lavaan Equivalence (δ = 0.001, 50 starts)

**rraw dataset (525 × 31, 5 factors):**

| Metric | Value |
|--------|-------|
| Overall MAE | 0.000001 |
| Max single-cell error | 0.000005 |
| Cells within 0.001 | 155/155 (100%) |
| Factor correlation max |Δ| | 0.000009 |

Per-factor:

| Factor | MAE | Max Error |
|--------|-----|-----------|
| LOC | 0.000001 | 0.000003 |
| CCA | 0.000002 | 0.000005 |
| ER | 0.000001 | 0.000003 |
| FSI | 0.000001 | 0.000002 |
| TW | 0.000001 | 0.000002 |

**Teacher Burnout (876 × 23, 4 factors):**

| Metric | Value |
|--------|-------|
| Overall MAE | 0.000002 |
| Cells within 0.001 | 92/92 (100%) |

---

## 15. Detailed Loading Comparisons

### 15.1 Teacher Burnout: DE Items (Single-Start Failure)

**Factor 2 (Emotional Exhaustion) — DE cross-loadings:**

| Item | lavaan | Carm (1 start) | Carm (50 starts) |
|------|--------|----------------|-------------------|
| DE1 | −0.230 | −0.016 | −0.230 |
| DE2 | −0.002 | +0.208 | −0.002 |
| DE3 | +0.028 | +0.229 | +0.028 |

**Factor 3 (Depersonalization) — DE primary loadings:**

| Item | lavaan | Carm (1 start) | Carm (50 starts) |
|------|--------|----------------|-------------------|
| DE1 | 0.753 | 0.619 | 0.753 |
| DE2 | 0.711 | 0.586 | 0.711 |
| DE3 | 0.680 | 0.558 | 0.680 |

**Factor correlations:**

| Pair | lavaan | Carm (1 start) | Carm (50 starts) |
|------|--------|----------------|-------------------|
| TE–DE | 0.612 | 0.388 | 0.612 |
| EE–DE | 0.540 | 0.211 | 0.540 |
| DE–RPA | 0.449 | 0.195 | 0.449 |

The single-start solution finds a local optimum where DE items are distributed between EE and DE factors, attenuating the DE primary loadings by 0.12–0.13 and inflating cross-loadings by 0.20–0.23. Factor correlations involving DE are underestimated by 0.22–0.33.

### 15.2 rraw Dataset: FSI and ER Factors

**FSI factor — most affected by local optima:**

| Item | lavaan | Carm (1 start) | Carm (30 starts) | Carm (50 starts) |
|------|--------|----------------|-------------------|-------------------|
| LOC2 on FSI | 0.361 | 0.213 | 0.360 | 0.361 |
| FSI4 on FSI | −0.029 | −0.166 | −0.029 | −0.029 |
| FSI6 on FSI | −0.205 | −0.331 | −0.205 | −0.205 |
| TW5 on FSI | 0.441 | 0.307 | 0.433 | 0.441 |

**ER factor — required 50 starts:**

| Item | lavaan | Carm (30 starts) | Carm (50 starts) |
|------|--------|-------------------|-------------------|
| ER2 on ER | 0.822 | 0.752 | 0.822 |
| ER3 on ER | 0.800 | 0.734 | 0.800 |
| ER5 on ER | 0.747 | 0.680 | 0.747 |
| ER6 on ER | 0.726 | 0.655 | 0.726 |

The ER factor illustrates non-monotonic convergence: 30 starts fixed FSI but made ER worse (0.071 error vs 0.059 with 1 start) because the better overall solution placed different weight on the ER–TW boundary. At 50 starts, the globally optimal basin was found.

---

## 16. Public API Reference

### 16.1 `runEFA(data, options?): FAResult`

Full exploratory factor analysis.

```typescript
const result = runEFA(data, {
  nFactors: 5,           // omit for auto-detection via parallel analysis
  extraction: 'ml',      // 'ml' | 'paf' (default: 'ml')
  rotation: 'geomin',    // 'varimax' | 'promax' | 'geomin' | 'oblimin' | 'quartimin' | 'none'
  geominDelta: 0.01,     // geomin δ parameter (default: 0.01)
  randomStarts: 50,      // for oblique rotations (default: 50)
  seed: 42,              // PRNG seed (default: 42)
  maxIter: 1000,         // max iterations (default: 1000)
  tol: 1e-6,             // convergence tolerance (default: 1e-6)
  variableNames: names,  // optional variable labels
})
```

### 16.2 `runCFA(data, model, options?): CFAResult`

Confirmatory factor analysis with specified measurement model.

```typescript
const result = runCFA(data, {
  F1: [0, 1, 2, 3, 4],
  F2: [5, 6, 7, 8, 9],
}, {
  maxIter: 1000,
  tol: 1e-6,
  variableNames: names,
  factorNames: ['Efficacy', 'Exhaustion'],
})
```

### 16.3 `runFADiagnostics(data, options?): FADiagnostics`

Pre-analysis diagnostics without extraction.

```typescript
const diag = runFADiagnostics(data, {
  seed: 42,
  parallelIterations: 100,
})
// diag.kmo, diag.bartlett, diag.parallelSuggested, diag.mapSuggested
```

---

## 17. References

1. Bartlett, M. S. (1950). Tests of significance in factor analysis. *British Journal of Statistical Psychology*, 3(2), 77–85.

2. Bentler, P. M. (1990). Comparative fit indexes in structural models. *Psychological Bulletin*, 107(2), 238–246.

3. Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and software for arbitrary rotation criteria in factor analysis. *Educational and Psychological Measurement*, 65(5), 676–696.

4. Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis. *Multivariate Behavioral Research*, 36(1), 111–150.

5. Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing model fit. In K. A. Bollen & J. S. Long (Eds.), *Testing structural equation models* (pp. 136–162). Sage.

6. Gorsuch, R. L. (1983). *Factor analysis* (2nd ed.). Erlbaum.

7. Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for rotation to oblique simple structure. *British Journal of Statistical Psychology*, 17(1), 65–70.

8. Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. *Psychometrika*, 30(2), 179–185.

9. Jennrich, R. I. (2002). A simple general method for oblique rotation. *Psychometrika*, 67(1), 7–27.

10. Jöreskog, K. G. (1967). Some contributions to maximum likelihood factor analysis. *Psychometrika*, 32(4), 443–482.

11. Kaiser, H. F. (1958). The varimax criterion for analytic rotation in factor analysis. *Psychometrika*, 23(3), 187–200.

12. Kaiser, H. F. (1974). An index of factorial simplicity. *Psychometrika*, 39(1), 31–36.

13. Mezzadri, F. (2007). How to generate random matrices from the classical compact groups. *Notices of the AMS*, 54(5), 592–604.

14. Steiger, J. H. (1990). Structural model evaluation and modification: An interval estimation approach. *Multivariate Behavioral Research*, 25(2), 173–180.

15. Stewart, G. W. (1980). The efficient generation of random orthogonal matrices with an application to condition estimators. *SIAM Journal on Numerical Analysis*, 17(3), 403–409.

16. Tucker, L. R., & Lewis, C. (1973). A reliability coefficient for maximum likelihood factor analysis. *Psychometrika*, 38(1), 1–10.

17. Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations. *Psychometrika*, 41(3), 321–327.

18. Yates, A. (1987). *Multivariate exploratory data analysis: A perspective on exploratory factor analysis*. SUNY Press.

---

## 18. Comparison with Other Software

| Feature | Carm | R psych | R lavaan | Mplus |
|---------|------|---------|----------|-------|
| ML extraction | Jöreskog + Nelder-Mead | L-BFGS-B | ML (custom) | ML (custom) |
| PAF extraction | Iterated eigen | Iterated eigen | — | — |
| Varimax | Kaiser-normalized SVD | Same | Same | Custom |
| Promax | Outer Kaiser + Procrustes | Same (psych::kaiser) | — | Custom |
| Geomin | GPFoblq | GPFoblq (via GPArotation) | GPA×30 | Custom |
| Oblimin | GPFoblq | GPFoblq (via GPArotation) | GPA×30 | — |
| Quartimin | GPFoblq | GPFoblq (via GPArotation) | GPA×30 | — |
| Random starts | 50 (default) | 1 | 30 | 0 |
| CFA | Gradient + numerical Hessian | — | Analytic derivatives | Analytic |
| CFA standard errors | Numerical Hessian | — | Expected information | Expected |
| KMO | Anti-image correlation | Same | — | — |
| Parallel analysis | Monte Carlo (100 iter) | Same | — | — |
| MAP | Velicer (1976) | Same | — | — |
| Fit indices | χ², RMSEA+CI, CFI, TLI, SRMR, AIC, BIC | Same | Same + more | Same + more |
| Deterministic | Yes (seeded PRNG) | No (R's MT) | No (R's MT) | Yes |
| Dependencies | 0 | MASS, GPArotation | — | — |
| Runtime environment | Browser / Node.js | R | R | Standalone |
| Licensing | — | GPL-3 | GPL-2 | Commercial |

---

## 19. Summary

Carm's factor analysis module provides a complete, zero-dependency TypeScript implementation of:

- **Two extraction methods** (ML and PAF) with two-phase optimization for ML
- **Six rotation methods** (varimax, promax, geomin, oblimin, quartimin, none)
- **General oblique rotation engine** (GPFoblq) with pluggable criterion functions
- **50 random starting matrices** (Haar-distributed) for robust global optimization
- **Confirmatory factor analysis** with ML estimation, standard errors, and parameter tables
- **Full fit index suite** (χ², RMSEA with 90% CI, CFI, TLI, SRMR, AIC, BIC)
- **Psychometric diagnostics** (KMO, Bartlett, parallel analysis, MAP)

The implementation achieves numerical equivalence with R's psych, GPArotation, and lavaan packages, demonstrated through:

- **100/100 synthetic datasets** for promax and geomin rotations
- **6/6 real dataset configurations** across different factor counts
- **100/100 diagnostic comparisons** at machine precision
- **Side-by-side loading equivalence** with lavaan (MAE < 0.000005) on two real-world datasets

All validation artifacts are permanently stored in `validation/` and can be regenerated on demand.
