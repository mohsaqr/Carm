# Carm Module Development Prompt

> Use this system prompt when asking an LLM to write code for Carm or a related module.

---

You are building a module for **Carm**, a zero-dependency TypeScript statistical analysis and visualization library. Carm implements all distributions, hypothesis tests, matrix algebra, and optimization from first principles — no jStat, no simple-statistics, no external math libraries. D3 is the sole dependency, loaded dynamically as a peer dependency so the statistical core has zero runtime cost.

When writing code for Carm, you must follow its exact engineering conventions:

## 1. Architecture & Module Boundaries

Carm has three layers with strict dependency direction:

```
core/          →  Pure math, types, formatting (no stats, no DOM)
  types.ts        Central type contract — ALL interfaces live here
  math.ts         Gamma, Beta, distributions, Nelder-Mead, p-adjustment
  matrix.ts       Matrix class (LU, Cholesky, SVD, eigendecomposition)
  apa.ts          APA 7th edition formatters

stats/         →  Pure computation (imports core/, never viz/)
  descriptive.ts, comparison.ts, correlation.ts, regression.ts,
  clustering.ts, pca.ts, mixed.ts, effect-size.ts, post-hoc.ts,
  preprocess.ts, frequency.ts, analyze.ts

viz/           →  D3 rendering (imports core/ and stats/, never imported by stats/)
  plots/          40+ renderers (violin-box, scatter-stats, histogram, etc.)
  components/     Shared: axis, brackets, annotations, tooltip
  themes/         Color palettes, typography, spacing
  export.ts       SVG/PNG export
```

**Rules**: `stats/` is pure — no DOM, no D3, no side effects. `viz/` consumes `stats/`, never the reverse. All shared interfaces go in `core/types.ts`.

## 2. Result Interfaces

Every statistical function returns a structured, readonly result object. The universal pattern:

```typescript
interface StatResult {
  readonly testName: string                        // "Welch's t-test"
  readonly statistic: number                       // t, F, χ², U, etc.
  readonly df: number | readonly [number, number]  // single or [df1, df2]
  readonly pValue: number
  readonly effectSize: EffectSize                  // { value, name, interpretation }
  readonly ci: readonly [number, number]
  readonly ciLevel: number                         // 0.95
  readonly n: number
  readonly formatted: string                       // APA string, always present
}

interface EffectSize {
  readonly value: number
  readonly name: string                   // "Cohen's d", "η²", "Cramér's V"
  readonly interpretation: EffectInterpretation  // 'negligible'|'small'|'medium'|'large'
}
```

**Every result has `formatted`** — an APA 7th edition string ready for plot subtitles, reports, or clipboard. Generated via helpers in `core/apa.ts` (`formatTTest`, `formatANOVA`, `formatChiSq`, `formatRegression`, `formatP`). The `formatP` function drops leading zeros: `p = .025`, `p < .001`.

Higher-level results (clustering, regression, ANOVA) have domain-specific fields but always include a `formatted` string and follow the readonly convention.

## 3. TypeScript Strictness

Carm uses maximum TypeScript strictness:

```json
{
  "strict": true,
  "noUncheckedIndexedAccess": true,
  "exactOptionalPropertyTypes": true,
  "noUnusedLocals": true,
  "noUnusedParameters": true
}
```

**Implications you must follow:**

- `arr[i]` is `T | undefined` — use `arr[i]!` only after bounds validation, or guard with `if (arr[i] !== undefined)`.
- `arr[i]! += x` is required for compound assignment on indexed access.
- Optional properties cannot be assigned `undefined` — use conditional spread: `...(x !== undefined && { x })`.
- All return types use `readonly` arrays and tuples: `readonly number[]`, `readonly [number, number]`.
- Every function has explicit parameter and return types. No `any` — use `unknown` with type guards.

## 4. Zero-Dependency Math

Carm implements its own numerical engine in `core/math.ts` and `core/matrix.ts`:

- **Distributions**: Normal CDF/quantile (A&S approximations), t, F, chi-square, Beta, Gamma (Lanczos approximation)
- **Optimization**: Nelder-Mead simplex (used for REML in LMM, MLE in various models)
- **Linear algebra**: Matrix class with LU decomposition, Cholesky factorization, Jacobi SVD, eigendecomposition
- **P-value adjustment**: Bonferroni, Holm, Benjamini-Hochberg (BH/FDR)
- **Basics**: `mean()`, `sd()`, `variance()`, `median()`, `quantile()`, `rank()` — all with edge-case guards

**When you need a mathematical function, check `core/math.ts` first.** If it exists, use it. If not, implement it natively following the same patterns (cite the algorithm source in a comment, handle edge cases, return typed results). Never add an external dependency.

## 5. Deterministic PRNG

All stochastic algorithms use a shared splitmix32 PRNG:

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
    return (t >>> 0) / 4294967296
  }
}
```

- Default seed: `42`
- Every function with randomness accepts `seed?` in its options interface
- K-Means++, GMM initialization, LCA, LTA, bootstrap — all use this PRNG
- This guarantees reproducible results across runs and platforms

## 6. Edge Case Handling

Every public function validates inputs at entry:

```typescript
if (data.length === 0) throw new Error('fitGMM: data cannot be empty')
if (k < 1) throw new Error('fitGMM: k must be >= 1')
if (k > n) throw new Error('fitGMM: k cannot exceed n')
if (x.length !== y.length) throw new Error('pearsonCorrelation: arrays must have equal length')
if (sdX === 0 || sdY === 0) throw new Error('pearsonCorrelation: zero variance in input')
```

**Required guards:**
- Empty arrays → `throw` with function name prefix
- Zero variance → handle gracefully (throw or return degenerate result, depending on context)
- `NaN` / `Infinity` → check `isFinite()` on critical intermediate values
- Mismatched lengths → `throw`
- Log-space operations → prevent `log(0)` with floor (e.g., `Math.max(x, 1e-300)`)
- Iterative algorithms → `maxIter` + `tol` convergence with fallback if not converged

Private utility functions may use `!` non-null assertions after the outer function has validated bounds.

## 7. Dynamic D3 Imports

All visualization functions load D3 lazily:

```typescript
export function renderViolinBox(
  container: HTMLElement,
  data: ViolinBoxData,
  config: ViolinBoxConfig = {}
): void {
  import('d3').then(d3 => renderViolinBoxD3(d3, container, data, config))
}
```

D3 is a peer dependency — consumers using only the `stats/` layer pay zero bundle cost. The actual render logic lives in an inner function that receives the `d3` namespace as a parameter.

## 8. Performance Conventions

- **Typed arrays**: Use `Float64Array` for accumulation loops in hot paths (clustering, matrix operations). Convert back to `number[]` in return values.
- **Avoid O(n^2) where possible**: For distance matrices and similar, O(n^2) is unavoidable — use `Float64Array` with row-major flat storage.
- **Log-space arithmetic**: EM algorithms, LCA, LTA all work in log-space to prevent underflow. Use `logSumExp` for stable log-probability addition.
- **Early termination**: Iterative algorithms check convergence per iteration (`|LL_new - LL_old| < tol`), not just at `maxIter`.

## 9. Intelligent Dispatch

The `analyze()` function in `stats/analyze.ts` demonstrates Carm's dispatch pattern:

1. Infer field types from data (`detectFieldType`)
2. Assess preconditions (Shapiro-Wilk normality, group counts, sample sizes)
3. Route to the correct test (parametric vs non-parametric, 2-group vs k-group)
4. Apply corrections automatically (Welch for unequal variances, tie-correction for ranks)
5. Return a unified `AnalysisResult` with descriptives, test result, post-hoc, and normality info

When designing a new top-level API, follow this pattern: smart defaults, automatic precondition assessment, transparent routing with the chosen test name in the result.

## 10. Testing & Cross-Validation Protocol

Every function gets a test file in `tests/stats/`. Tests follow this structure:

### Required test categories:
1. **Known textbook examples** — cite the source in a comment
2. **R cross-validation** — run the same analysis in R, embed the R code and expected values
3. **Edge cases** — empty input, single value, zero variance, ties, perfect correlation
4. **Determinism** — same seed produces same result across runs

### R cross-validation format:
```typescript
// Cross-validated with R:
// > t.test(c(1,2,3,4,5), c(6,7,8,9,10), var.equal = FALSE)
// t = -5.0, df = 8, p-value = 0.001053
// 95% CI: [-6.306, -2.694]
// Cohen's d (effsize::cohen.d): -3.162
it('matches R Welch t-test', () => {
  const r = tTestIndependent([1,2,3,4,5], [6,7,8,9,10])
  expect(r.statistic).toBeCloseTo(-5.0, 2)
  expect(r.pValue).toBeCloseTo(0.001053, 3)
})
```

### For EM-based algorithms (GMM, LCA, LTA):
- Label order is initialization-dependent — compare means/parameters in sorted order
- Use order-invariant comparison: sort by first coordinate, then compare element-wise
- Use permutation-based classification agreement for label matching
- Allow log-likelihood tolerance of ~2.0 (different initializations converge to nearby optima)
- Try multiple seeds and take the best result to close the gap with R

### Quality-based comparison for non-convex optimization:
For algorithms with non-convex objectives (K-Means, GMM, FA rotations), Carm may find a **better** solution than R due to different initialization strategies. In cross-validation, use one-sided error:
```typescript
// WRONG: penalizes Carm for finding better solutions
error = Math.abs(carm_loglik - r_loglik)

// RIGHT: error = 0 if Carm is equal or better
error = Math.max(0, r_loglik - carm_loglik)       // for loglik (higher = better)
error = Math.max(0, carm_inertia - r_inertia)     // for inertia (lower = better)
error = Math.max(0, carm_bic - (-r_mclust_bic))   // for BIC (lower = better)
```
This applies to: GMM log-likelihood/BIC, K-Means inertia/within-SS, FA rotation criterion values.

### Multi-start strategy for reproducible cross-validation:
Non-convex algorithms require sufficient random restarts to reliably find the global optimum on both sides:
| Algorithm | R Reference | Carm Harness | Why |
|-----------|-------------|--------------|-----|
| K-Means | `nstart=500` | 500 K-Means++ starts | Both converge to same global minimum |
| GMM | mclust (hierarchical init) | 200 K-Means++ starts | Different init strategy; more starts compensate |
| FA rotations | `GPArotation` (1 start) | `randomStarts=50` (Haar matrices) | Non-convex; 50 starts matches lavaan's strategy |

### R's swilk.c polynomial conventions (Shapiro-Wilk p-value):
R's `swilk.c` uses ascending-degree polynomials evaluated with Horner's method:
```
poly(c, x) = c[0] + c[1]*x + c[2]*x² + c[3]*x³ + ...
```
This is the OPPOSITE of many textbook implementations that use descending order. Critical details:
- **n=3**: Exact formula `p = (6/π) * (asin(√W) - Φ⁻¹(0.75))`, NOT the polynomial approximation
- **n≤11**: Variable is `n` (NOT `1/n`). Coefficients: `SW_G=[-2.273, 0.459]`, `SW_C3=[0.544, -0.39978, 0.025054, -6.714e-4]`, `SW_C4=[1.3822, -0.77857, 0.062767, -0.0020322]`
- **n≥12**: Variable is `log(n)`. Coefficients: `SW_C5=[-1.5861, -0.31082, -0.083751, 0.0038915]`, `SW_C6=[-0.4803, -0.082676, 0.0030302]`
- Always verify against R's source: `src/library/stats/src/swilk.c`

### Standard tolerance levels:
| Quantity | `toBeCloseTo` precision | Absolute tolerance |
|----------|------------------------|--------------------|
| Test statistic | 2 | ~0.01 |
| p-value | 3-4 | ~0.001 |
| Effect size | 2-3 | ~0.01 |
| df (Welch) | 1 | ~0.1 |
| Log-likelihood (EM) | 0 | ~2.0 |
| Deterministic algorithms (HAC, KMeans w/ same init) | 10 | ~1e-10 |

### R reference data:
For complex modules, generate reference data with an R script (`tests/r_*_reference.R`) that outputs JSON (`tests/r_*_reference.json`). The test file reads the JSON and compares against it. **Always extract reference values programmatically from the JSON — never hand-transcribe.**

## 11. Naming Conventions

- Functions: `camelCase` with descriptive verbs — `fitGMM()`, `runDBSCAN()`, `cutTree()`, `pearsonCorrelation()`
- Interfaces: `PascalCase` — `StatResult`, `GMMOptions`, `ClusterDiagnostics`
- Constants: `UPPER_SNAKE_CASE` — `DEFAULT_ALPHA`, `MAX_ITER`
- Files: `kebab-case.ts` — `effect-size.ts`, `post-hoc.ts`
- Test files: `*.test.ts` matching source name — `clustering.test.ts`
- Error messages: prefix with function name — `'fitGMM: data cannot be empty'`

## 12. Code Style

- **Comments**: Cite algorithm sources (`// Lanczos approximation, Numerical Recipes §6.1`) and R package equivalents (`// Matches R: hclust(method="ward.D2")`). Don't comment obvious code.
- **Immutability**: Return `readonly` arrays. Don't mutate input data — copy first if mutation is needed internally.
- **No classes for data**: Use plain interfaces and functions. The only class is `Matrix` (in `core/matrix.ts`) and `PRNG` (internal to clustering).
- **Export discipline**: Export types and functions from module `index.ts`. Internal helpers stay unexported.
- **Options pattern**: Every configurable function takes an options object with defaults:
  ```typescript
  interface GMMOptions {
    readonly k: number
    readonly model?: CovarianceModel  // default: 'VVV'
    readonly seed?: number            // default: 42
    readonly maxIter?: number         // default: 200
    readonly tol?: number             // default: 1e-6
    readonly preprocess?: PreprocessMethod  // default: 'none'
  }
  ```

---

## 13. Validation: Run With Numbers and Produce HTML

**Never assume code works. Always run it with concrete numbers and inspect the output.**

After writing any module, you must:

1. **Create a test script** (`tmp/test-<module>.mjs` or `.ts`) that imports your functions, runs them with synthetic data, and prints every field of every result object.
2. **Generate a test HTML file** (`tmp/test-<module>.html`) that renders the numerical results in a readable format — tables, formatted values, pass/fail indicators. Open it in the browser for visual inspection.
3. **Compare against R**: Run the same analysis in R, capture R's output, and display both side-by-side in the HTML. Highlight any discrepancy beyond tolerance.

### HTML validation template pattern:

```typescript
// tmp/test-mymodule.mjs
import { myFunction } from '../dist/index.js'
import { writeFileSync } from 'fs'

const data = [[1.2, 3.4], [5.6, 7.8], ...]  // synthetic test data
const result = myFunction(data, { k: 3 })

// R reference (from running R script)
const rRef = { logLik: -123.45, bic: 260.12, labels: [0, 0, 1, 1, 2, 2] }

const rows = [
  ['Log-likelihood', result.diagnostics.logLikelihood.toFixed(4), rRef.logLik.toFixed(4),
   Math.abs(result.diagnostics.logLikelihood - rRef.logLik) < 2.0 ? 'PASS' : 'FAIL'],
  ['BIC', result.diagnostics.bic.toFixed(2), rRef.bic.toFixed(2),
   Math.abs(result.diagnostics.bic - rRef.bic) < 5.0 ? 'PASS' : 'FAIL'],
]

const html = `<!DOCTYPE html>
<html><head><style>
  body { font-family: system-ui; max-width: 900px; margin: 2rem auto; }
  table { border-collapse: collapse; width: 100%; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: right; }
  .pass { background: #d4edda; } .fail { background: #f8d7da; }
</style></head><body>
<h1>Test: myFunction</h1>
<table>
  <tr><th>Metric</th><th>Carm</th><th>R</th><th>Status</th></tr>
  ${rows.map(r => `<tr><td>${r[0]}</td><td>${r[1]}</td><td>${r[2]}</td>
    <td class="${r[3].toLowerCase()}">${r[3]}</td></tr>`).join('')}
</table>
</body></html>`

writeFileSync('tmp/test-mymodule.html', html)
```

Then open the HTML: `open tmp/test-mymodule.html`

**This is mandatory for every new function.** The HTML serves as a living verification artifact. Do not delete test HTML files from `tmp/` without user approval.

For visualization changes (plots, themes, annotations), the HTML must render the actual plot so the user can visually inspect it — you cannot see images, so always present the HTML for user review.

---

## Quick Reference Checklist

Before submitting code for Carm, verify:

- [ ] All interfaces defined in or imported from `core/types.ts`
- [ ] Result objects include `formatted` APA string
- [ ] All arrays in return types are `readonly`
- [ ] Input validation with function-name-prefixed error messages
- [ ] No external math dependencies — use `core/math.ts` utilities
- [ ] Stochastic functions use splitmix32 PRNG with `seed?` option
- [ ] D3 loaded via dynamic `import('d3').then(...)` in viz layer only
- [ ] `Float64Array` for hot accumulation loops
- [ ] Log-space arithmetic where underflow is possible
- [ ] Tests: textbook examples + R cross-validation + edge cases + determinism
- [ ] R code embedded in test comments
- [ ] Reference values extracted programmatically from JSON, never hand-transcribed
- [ ] `noUncheckedIndexedAccess` satisfied (no bare `arr[i]` in expressions)
- [ ] `exactOptionalPropertyTypes` satisfied (conditional spread for optional fields)
- [ ] Validation HTML generated in `tmp/` with numerical Carm vs R comparison
- [ ] HTML opened in browser and presented to user for review
- [ ] No test HTML files deleted without user approval
