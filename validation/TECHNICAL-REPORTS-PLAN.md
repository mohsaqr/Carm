# Plan: Comprehensive Technical Reports for All Carm Modules

## Context

Carm has one exemplary technical report — `validation/GMM-TECHNICAL-REPORT.md` (1,537 lines) — covering the clustering module with 21 sections including engineering decisions, mathematical tricks, cross-validation results, and API reference. The FA module has `validation/FA-TECHNICAL-REPORT.md` (1,085 lines) but it LACKS the critical "Engineering Decisions" and "Mathematical Tricks" sections that make the GMM report so valuable. The user wants:

1. **Extend the existing FA report** with a Section 20 (Engineering Decisions) documenting every obstacle overcome
2. **Write new reports** for all other major modules following the GMM template
3. Each report documents HOW problems were solved, not just what was implemented

## Template Structure (from GMM report)

Every report follows this structure:
1. Header: module name, file path, LOC, dependencies, validation status
2. Architecture & design principles
3. Algorithm sections (one per major algorithm)
4. Numerical stability & tricks
5. Cross-validation methodology & results (if available)
6. Public API reference
7. References (academic citations)
8. **Engineering Decisions: Problems, Solutions, and Optimizations** — the crown jewel
9. **Mathematical Tricks That Made It Possible** — pure-TS tricks replacing LAPACK/BLAS

## Files to Create/Modify

### 1. EXTEND: `validation/FA-TECHNICAL-REPORT.md` (add ~400 lines)

Add two new sections at the end of the existing report:

**Section 20: Engineering Decisions: Problems, Solutions, and Optimizations**

Document these specific decisions (all discovered from source code analysis):

| # | Problem | Solution | Lines |
|---|---------|----------|-------|
| 20.1 | Two-phase ML extraction: why not just Nelder-Mead? | Phase 1 (Jöreskog gradient with cosine-annealed LR + momentum) converges fast to vicinity; Phase 2 (Nelder-Mead) polishes to match R's `factanal()` precision. Nelder-Mead alone on d-dimensional uniqueness is too slow. | 455-553 |
| 20.2 | Cosine-annealed learning rate | `lr = lr_min + (lr0 - lr_min) * 0.5 * (1 + cos(pi * iter/maxIter))` — starts at 0.02, decays to 0.001. Prevents oscillation near convergence while maintaining fast early progress. With fixed LR, loadings oscillate and never match R. | 463 |
| 20.3 | Death penalty in Nelder-Mead | Nelder-Mead is unconstrained. Adding penalty of 1000 when uniqueness exits [0.005, 0.995] keeps it bounded without modifying the optimizer. R uses L-BFGS-B (bounded quasi-Newton) but implementing that from scratch is prohibitively complex. | 509-512 |
| 20.4 | Outer Kaiser normalization in Promax | `psych::fa()` dispatches promax via `psych::kaiser()` which normalizes by communality BEFORE rotation. Without this, results differ from R by MAE ≈ 0.03. This is NOT documented in any textbook — discovered by reading R source code. | 927-959 |
| 20.5 | Varimax SVD via eigendecomposition of B'B | Direct one-sided Jacobi SVD has accuracy issues for small k×k matrices. Computing the polar decomposition via eigendecomposition of B'B is more stable and gives identical results. | 636-673 |
| 20.6 | Log-space geomin criterion | The product `∏(L²+δ)^(1/k)` overflows/underflows for large k and extreme loadings. Computing as `exp((1/k) * Σ log(L²+δ))` prevents this entirely. | 739-742 |
| 20.7 | GPFoblq always accepts step | R's GPArotation always updates T even if Armijo condition fails after 11 backtracking iterations. Replicating this "always accept" behavior was critical for matching R exactly. Initially rejecting failed line searches caused divergence. | 871-877 |
| 20.8 | Haar random orthogonal matrices | Naive random matrix + Gram-Schmidt does NOT give Haar-distributed matrices on SO(k). Must apply sign correction `Q[:,j] *= sign(R[j,j])` after QR. Without this, rotation starts are biased and 50 starts are needed instead of 30. | 95-100 |
| 20.9 | CFA standard errors via numerical Hessian | Analytic derivatives of ML discrepancy w.r.t. all CFA parameters are complex and error-prone. Central finite differences with h=1e-4 are simpler, adequate for SE computation, and allow triple fallback (inverse → pseudo-inverse → default 0.05). | 1545-1616 |
| 20.10 | Bartlett correction: EFA yes, CFA no | R's `psych::fa()` applies Bartlett (1950) correction `(n-1-(2p+5)/6-2k/3)` to chi-square but lavaan's CFA uses uncorrected `(n-1)`. Applying the wrong correction produces chi-square discrepancies of 5-20%. | 215-220 |
| 20.11 | TLI intentionally NOT clamped to [0,1] | Many implementations clamp TLI like CFI. But TLI can legitimately exceed 1 for models that fit better than expected. Clamping would hide this information and mismatch R's `psych::fa()`. | 276 |
| 20.12 | Communality computation for oblique rotations | Oblique communalities are NOT just `Σ L²`. Must use `Σ_{f1,f2} L_{i,f1} * Φ_{f1,f2} * L_{i,f2}` to account for factor correlations. Using the orthogonal formula with oblique loadings gives systematically low communalities (error ≈ 0.05). | 1723-1735 |
| 20.13 | Float64Array accumulators for correlation matrix | Long summation loops (means, SDs) drift with standard `number[]`. Using `Float64Array` prevents accumulation errors that compound into correlation matrix discrepancies of 1e-10 → 1e-8. | 116-117 |
| 20.14 | Heywood case clamping: 0.001 to 0.9999 | PAF communalities can exceed 1.0 (Heywood case). Clamping to 0.9999 (not 1.0) prevents the reduced correlation matrix from becoming singular, which would crash eigendecomposition. The 0.001 floor prevents items from being completely excluded. | 369 |
| 20.15 | Sign convention for eigenvectors | Different eigendecomposition implementations (Intel MKL, OpenBLAS, Jacobi) produce eigenvectors with different signs. Enforcing "largest absolute element is positive" matches LAPACK's `dsyev` convention used by R's `eigen()`. Without this, loadings flip signs randomly across platforms. | 375-385 |
| 20.16 | ML uniqueness initialization matches R exactly | `start = (1 - 0.5*k/p) / R^{-1}_{ii}` — this is R's `stats::factanal` formula. Using a different initialization (e.g., 1-SMC or uniform 0.5) converges to different local optima, producing loading discrepancies of MAE ≈ 0.01. | 418-428 |
| 20.17 | Velicer MAP negative diagonal guard | After extracting k components, residual diagonal can go negative. Using `abs()` in the denominator `sqrt(|R*_{ii}| * |R*_{jj}|)` prevents NaN propagation in the partial correlation computation. | 1142-1143 |
| 20.18 | RMSEA CI via Steiger approximation | Exact CI requires non-central chi-square quantile function, which is computationally expensive. The Steiger (1990) approximation `NCP ± 1.645*sqrt(2*df)` is standard in practice and avoids implementing `qchisq(ncp=...)`. | 255-269 |
| 20.19 | CFA factor covariance constraints | Factor variances (diagonal of Φ) must be fixed at 1.0 for identification. Off-diagonal covariances are clamped to [-0.99, 0.99] to prevent degenerate solutions during gradient descent. | 1434, 1440 |
| 20.20 | 50 random starts default: empirical calibration | Tested on rraw (525×31, 5 factors) and teacher burnout (876×23, 4 factors). 10 starts sufficient for teacher burnout, but rraw needs 30+. 50 gives safety margin. Timing: 50 starts ≈ 3s, acceptable for interactive use. | validated empirically |

**Section 21: Mathematical Tricks That Made It Possible**

| # | Trick | Why needed | How |
|---|-------|-----------|-----|
| 21.1 | Eigendecomposition for varimax polar factor | Can't call LAPACK `dgesvd`. Compute B'B (k×k), eigendecompose, get singular values as `sqrt(eigenvalues)`, reconstruct U = BV·diag(1/σ), polar factor T = UV'. | 636-673 |
| 21.2 | Concentrated ML objective | Instead of optimizing over both loadings and uniquenesses, concentrate out loadings analytically. Objective becomes a function of uniquenesses only: `F(Ψ) = -Σ(log λ_j - λ_j) + k - d`. This halves the optimization dimension from p*k+p to just p. | 494-513 |
| 21.3 | Promax via OLS regression | Target = sign(V)·|V|⁴. Transformation U = (V'V)⁻¹V'Q is just least-squares. Then normalize U so Φ has unit diagonal. No iterative optimization needed — promax is purely algebraic after varimax. | 964-998 |
| 21.4 | Anti-image correlation for KMO | KMO uses the anti-image matrix `Q = -R⁻¹_{ij}/sqrt(R⁻¹_{ii}·R⁻¹_{jj})`. This requires only one matrix inversion (not per-item regressions), making it O(p³) instead of O(p⁴). | 1233-1262 |
| 21.5 | Cholesky log-determinant with eigenvalue fallback | `log|R| = 2·Σ log(L_{ii})` via Cholesky is O(p²). When R is near-singular, Cholesky fails; fall back to `Σ log(max(λ_i, 1e-15))` from eigendecomposition. | 1266-1271 |
| 21.6 | Modified Gram-Schmidt for Haar QR | Classical GS loses orthogonality for k > 5. MGS reorthogonalizes against already-processed columns, maintaining `|Q'Q - I| < 1e-14` even for k = 10. | 74-93 |

---

### 2. NEW: `validation/DESCRIPTIVE-COMPARISON-TECHNICAL-REPORT.md`

**Covers**: `descriptive.ts` (260 lines), `comparison.ts` (450 lines), `effect-size.ts` (190 lines), `frequency.ts` (290 lines), `post-hoc.ts` (220 lines)

These are grouped because they share the same data type (group comparisons) and have deep cross-references.

**Outline** (~800 lines):

1. Architecture: pure functions, zero deps, shared `StatResult` contract
2. Shapiro-Wilk: AS R94 Algorithm
   - Three-path coefficient computation (n=3, n=4-5, n≥6)
   - Polynomial correction `SW_C1/SW_C2` for outermost coefficients
   - W clamping to [0,1] and SST zero guard
   - P-value: Royston (1995) approximation branching on n≤11 vs n>11
3. Skewness & Kurtosis: adjusted Fisher-Pearson formulas matching `e1071::skewness(type=2)`
4. Independent T-Test: Welch-Satterthwaite df, one-sided conversion, zero-SE guard
5. Paired T-Test: difference-score approach
6. Mann-Whitney U: normal approximation with tie correction, exact DP for n≤20
7. Wilcoxon Signed-Rank: exact DP distribution for n≤20, continuity-corrected normal for n>20
8. Kruskal-Wallis: tie correction divisor `1 - C/(N³-N)`
9. Friedman Test: rank-sum approach for repeated measures
10. One-Way ANOVA: SS decomposition, F-test
11. Effect Sizes: Cohen's d (pooled SD), Hedges' g (bias correction J), eta²/omega² (with max(0,...) guard), rank-biserial
12. Chi-Square: Yates correction, Cramér's V
13. Fisher's Exact: log-space hypergeometric PMF, Haldane-Anscombe 0.5 correction for OR
14. Tukey HSD: Bonferroni-adjusted t approximation (NOT exact studentized range)
15. Games-Howell: Welch df + Bonferroni-adjusted approximation
16. Dunn's Test: rank means, tie correction, 3 p-adjustment methods
17. Public API Reference
18. References
19. Engineering Decisions (~15 items):
    - Why Welch is default (robust when variances differ, reduces to Student when equal)
    - Wilcoxon exact DP: building full distribution via iterative rank addition
    - Mann-Whitney: why normal approximation rather than exact tables
    - Tukey HSD: why Bonferroni approximation instead of studentized range (no `qtukey()` implementation — would require a new special function)
    - Four separate copies of normal CDF (comparison, regression, post-hoc, effect-size) to avoid circular imports
    - Cohen's d CI: Hedges & Olkin normal approximation vs exact non-central t (simpler, adequate precision, no `qt(ncp=...)` needed)
    - Omega-squared max(0,...): prevents negative estimates for negligible effects
    - Fisher log-space hypergeometric: prevents integer overflow for large contingency tables
    - Haldane-Anscombe correction: adds 0.5 to all cells when any cell is 0
    - Friedman ranking: per-row ranks, not global ranks
    - APA formatting: p < .001 threshold, no leading zero, semicolon-separated
20. Mathematical Tricks (~5 items):
    - Shapiro-Wilk polynomial corrections from AS R94
    - Wilcoxon exact via DP (O(n·maxW) instead of O(2^n) brute force)
    - Fisher's test in log-space: `logC(K,k) + logC(N-K,n-k) - logC(N,n)` via `logFactorial`
    - Erf-based normal CDF: A&S 7.1.26 polynomial (5 terms, ~1.5e-7 accuracy)

---

### 3. NEW: `validation/CORRELATION-REGRESSION-TECHNICAL-REPORT.md`

**Covers**: `correlation.ts` (260 lines), `regression.ts` (380 lines), `pca.ts` (160 lines)

Grouped because they share the linear algebra pipeline: correlation → regression → PCA.

**Outline** (~700 lines):

1. Architecture: OLS engine shared by linear/multiple/polynomial/logistic, Matrix class dependency
2. Pearson Correlation: r clamping to [-1,1], Fisher z-transform CI, perfect-correlation guard (t → ∞)
3. Spearman as Pearson-on-Ranks: average-tie ranking, why not specialized formula
4. Kendall Tau-b: O(n²) pairwise concordance, tie corrections, normal approximation p-value
5. Partial Correlation via Residualization: OLS residuals approach (generalizes to multiple controls), NOT recursive formula
6. Correlation Matrix: pairwise construction with per-pair error catching
7. OLS Engine: normal equations (X'X)⁻¹X'y, R² clamping, adjusted R², AIC/BIC log-likelihood
8. Simple/Multiple/Polynomial Regression: all route through OLS
9. Logistic Regression: IRLS (Fisher scoring)
   - Weight clamping `max(1e-10, mu*(1-mu))` to prevent singular W
   - Convergence on max coefficient change
   - McFadden pseudo-R²
   - Wald z-test vs t-test discrepancy for CI
10. Regression Diagnostics: hat matrix, standardized residuals, Cook's D, VIF via recursive regression
11. PCA via SVD: one-sided Jacobi SVD, `X/sqrt(n-1)` scaling trick, sign indeterminacy
12. Varimax Rotation: Kaiser (1958) pairwise algorithm, angle-based convergence
13. Public API Reference
14. References
15. Engineering Decisions (~12 items):
    - OLS via normal equations, not QR (simpler, existing Matrix.inverse(), adequate for p<100)
    - R² clamping: prevents negative R² from floating-point error
    - AIC/BIC: `max(ss_res, 1e-15)` prevents log(0) on perfect fit
    - Coefficient SE: `sqrt(max(0, sigma² * C_{jj}))` prevents sqrt of negative
    - Logistic IRLS weight clamping: why 1e-10, what happens without it
    - McFadden pseudo-R²: NaN guard when null log-likelihood ≈ 0
    - VIF via recursive regression: simple but O(p²·n), acceptable for p<50
    - Spearman as Pearson-on-ranks: naturally handles ties, no separate formula needed
    - Fisher z CI: guarantees bounds in [-1,1] via tanh back-transform
    - Partial correlation via residualization: handles multiple controls without recursive formula
    - PCA `X/sqrt(n-1)` scaling: ensures singular values² = eigenvalues of covariance matrix
    - SVD vs eigendecomposition for PCA: SVD is more numerically stable for n >> p
16. Mathematical Tricks (~4 items):
    - IRLS as iterated normal equations: `(X'WX)⁻¹X'Wz` where z = working response
    - Fisher z-transform: `z = 0.5·ln((1+r)/(1-r))`, SE = `1/sqrt(n-3)`, back via `tanh()`
    - Hat matrix diagonal for Cook's D: `h_ii = X_i(X'X)⁻¹X_i'`
    - Concentrated log-likelihood for AIC: `-n/2·(log(2π) + log(RSS/n) + 1)`

---

### 4. NEW: `validation/MIXED-MODEL-TECHNICAL-REPORT.md`

**Covers**: `mixed.ts` (300 lines)

Standalone because LMM is architecturally unique — the only module using profiled REML, Woodbury identity, and multi-start Nelder-Mead.

**Outline** (~500 lines):

1. Architecture: profiled REML, single-parameter optimization over log(ψ)
2. The LMM Model: `y = Xβ + Zb + ε`, random intercepts
3. Profiled REML: why optimize over `log(ψ = σ²_b/σ²_e)` instead of two variance components
4. Woodbury Identity: `V⁻¹ = I - Z·D⁻¹·Z'` where `D = Z'Z + (1/ψ)·I` — reduces O(n³) to O(q³)
5. Matrix Determinant Lemma: `log|V| = q·log(ψ) + log|D|` — avoids n×n determinant
6. Multi-Start Nelder-Mead: 5 starting values [-4,-2,0,2,4] for logψ
7. Fixed Effects & Standard Errors: GLS estimates `β = (X'V⁻¹X)⁻¹X'V⁻¹y`
8. BLUPs: closed-form `b_j = ψ/(1 + ψ·n_j) · Σ e_i`
9. ICC: `σ²_b / (σ²_b + σ²_e)` — one-way ICC(1)
10. Satterthwaite df: simplified `n - p - nGroups + 1`
11. REML normalizing constant: omitted during optimization, added back for final LL
12. Public API Reference
13. References (Henderson 1950, Bates et al. 2015, Woodbury 1950)
14. Engineering Decisions (~10 items):
    - Why profiled REML over full REML (single-parameter is more stable, avoids 2D surface)
    - Why log-parameterization for ψ (ensures positivity, improves conditioning)
    - Woodbury identity: when q << n, this is O(q³) instead of O(n³) — critical for large n
    - Multi-start: 5 values spans 8 orders of magnitude in ψ; enough for most random effect scales
    - σ²_e clamping to 1e-8: prevents zero residual variance (perfect fit edge case)
    - Singular matrix fallback: returns Infinity for negLogLik, steering optimizer away
    - Simplified Satterthwaite: less precise than lme4's full KR approximation but much simpler
    - REML vs ML: REML gives unbiased variance estimates; ML underestimates (divides by n not n-p)
    - AIC counts p+2 parameters: p fixed effects + σ²_b + σ²_e
    - Near-zero ψ fallback: when `scale < 1e-10`, random intercept is negligible → V⁻¹ ≈ I
15. Mathematical Tricks (~4 items):
    - Woodbury identity derivation and why it's exact (not an approximation)
    - Matrix determinant lemma derivation
    - BLUP as posterior mean under normal prior: `b|y ~ N(ψ·Z'V⁻¹(y-Xβ), ψ·I - ψ²·Z'V⁻¹Z)`
    - Cholesky log-determinant: `log|D| = 2·Σ log(L_{ii})`, O(q²) instead of O(q³)

---

### 5. NEW: `validation/CORE-MATH-TECHNICAL-REPORT.md`

**Covers**: `core/math.ts` (550 lines), `core/matrix.ts` (420 lines)

The foundation that all stats modules depend on. Documents every distribution, special function, and matrix algorithm.

**Outline** (~900 lines):

1. Architecture: zero external deps, all distributions from scratch, splitmix32 PRNG
2. Special Functions:
   - logGamma: Lanczos (g=7, 9 coefficients), reflection formula for z<0.5
   - incompleteBeta: Lentz continued fraction, symmetry optimization, convergence parameters
   - incompleteGamma: series for x<a+1, Lentz CF complement for x≥a+1
   - erf: A&S 7.1.26 polynomial (5 terms, ~1.5e-7)
3. Distribution Functions:
   - Normal CDF/quantile: erf-based CDF, Peter Acklam's rational approximation for quantile (~1.15e-9)
   - t-distribution: incompleteBeta-based CDF, bisection quantile
   - F-distribution: incompleteBeta-based CDF
   - Chi-square: incompleteGamma-based CDF, bisection quantile
4. Nelder-Mead Optimizer: simplex initialization, reflection/expansion/contraction/shrink
5. P-Value Adjustment: Bonferroni, Holm, BH, BY — all with enforced monotonicity
6. Descriptive Utilities: mean, median, sd, variance, quantile (R type=7), rank (average-tie)
7. Matrix Class:
   - Row-major flat storage, i-k-j multiply loop (cache-friendly)
   - Cholesky decomposition: Cholesky-Banachiewicz, positive-definite validation
   - Gauss-Jordan inverse: augmented matrix with partial pivoting (NOT LU despite docstring)
   - Jacobi eigendecomposition: classical pivot selection, O(100·n²) max iterations, 1e-12 convergence
   - One-sided Jacobi SVD: convergence `|γ| < 1e-15·sqrt(αβ)`, O(200·n²) max iterations
   - Pseudo-inverse via SVD: threshold `tol·max(σ)` for zero singular values
8. Public API Reference
9. References (Lanczos, Abramowitz & Stegun, Acklam, Press et al. Numerical Recipes, Golub & Van Loan)
10. Engineering Decisions (~15 items):
    - Lanczos vs Stirling for logGamma (Lanczos: 15 digits, Stirling: ~8 digits)
    - Reflection formula: avoids Lanczos instability near z=0 and negative integers
    - incompleteBeta symmetry trick: converges faster when `x > (a+1)/(a+b+2)`, swap a↔b and use 1-x
    - FPMIN = 1e-30 in continued fractions: prevents division by zero without affecting precision
    - Bisection for distribution quantiles: simple, guaranteed convergence, 100 iterations ≈ 2^100 precision
    - Peter Acklam vs A&S for normal quantile: Acklam has ~1.15e-9 error vs A&S ~4.5e-4 — 5 orders better
    - i-k-j loop order in matrix multiply: sequential memory access for `other` matrix, ~2x faster than i-j-k
    - Gauss-Jordan vs LU: simpler implementation, same O(n³) complexity, adequate for p<100
    - Jacobi eigendecomposition vs QR iteration: Jacobi produces eigenvectors as byproduct, QR needs separate phase
    - One-sided Jacobi SVD: simpler than Golub-Reinsch bidiagonalization, adequate for small-to-medium matrices
    - Pseudo-inverse threshold: `tol * max(σ)` zeroes out near-zero singular values, prevents noise amplification
    - Nelder-Mead simplex initialization: 5% perturbation (or 0.00025 if near zero)
    - BH monotonicity enforcement via running minimum: ensures `adj_p[i] ≤ adj_p[i+1]` for sorted p-values
    - R type=7 quantile: matches R's default, linear interpolation between order statistics
    - `Number(v)` coercion in mean(): guards against string-type bypass
11. Mathematical Tricks (~8 items):
    - Lanczos approximation: transforms `Γ(z)` into a rational function plus exponential — the key trick that makes `logGamma` possible without tables
    - Lentz's continued fraction modification: avoids 0/0 via `FPMIN`, maintains convergence by tracking `C` and `D` separately
    - incompleteBeta via continued fraction: replaces the standard integral with rapidly-converging fractions — the building block for t, F, and chi-square CDFs
    - incompleteGamma series: `γ(a,x) = e^{-x} x^a Σ (x^n / Γ(a+n+1))` — converges when x < a+1
    - Peter Acklam's normal quantile: three-region rational approximation with matched asymptotic behavior
    - Jacobi rotation: one 2×2 problem `[cos θ, -sin θ; sin θ, cos θ]` zeroes A[p,r] exactly — repeated application converges quadratically
    - Cholesky for log-determinant: `log|A| = 2Σ log(L_{ii})` avoids the numerical instability of `det(A)` for large matrices
    - Matrix determinant vs log-determinant: `det(A)` overflows for d>30; `logDet(A)` handles d>1000

## Implementation Order

1. **FA report extension** — add sections 20-21 to existing `validation/FA-TECHNICAL-REPORT.md`
2. **Core Math report** — foundation for understanding all other modules
3. **Descriptive/Comparison report** — covers most functions, references core math
4. **Correlation/Regression report** — references core math and matrix
5. **Mixed Model report** — standalone, references matrix (Woodbury, Cholesky)

Each report: read source code → write report sections → verify line numbers → save to `validation/`

## Conventions

- Line numbers reference the source file at time of writing (may drift with future changes — note the commit hash)
- Every engineering decision follows the pattern: **Problem** → **Root cause / Why it matters** → **Solution** → **Why this solution over alternatives** → **Result/Impact**
- Mathematical tricks explain WHY the trick works, not just WHAT it does
- All academic citations include year, journal, and the specific result used
- Cross-reference other Carm reports where relevant (e.g., FA report references GMM report for eigendecomposition discussion)

## Estimated Output

| Report | Status | Est. Lines | Module LOC |
|--------|--------|------------|------------|
| FA-TECHNICAL-REPORT.md (extend) | Existing + 400 new | ~1,500 | 1,923 |
| CORE-MATH-TECHNICAL-REPORT.md | New | ~900 | 970 |
| DESCRIPTIVE-COMPARISON-TECHNICAL-REPORT.md | New | ~800 | 1,410 |
| CORRELATION-REGRESSION-TECHNICAL-REPORT.md | New | ~700 | 800 |
| MIXED-MODEL-TECHNICAL-REPORT.md | New | ~500 | 300 |
| GMM-TECHNICAL-REPORT.md (existing) | Complete | 1,537 | 1,850 |
| **Total** | | **~5,937** | **7,253** |
