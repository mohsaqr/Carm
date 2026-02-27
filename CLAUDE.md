# StatViz — Standalone Statistical Software in TypeScript

## Project Overview
A completely independent, browser-based statistical analysis and visualization application built in TypeScript. No R, no Python, no server — everything runs client-side. The software provides publication-quality, richly annotated statistical plots inspired by the ggstatsplot philosophy: every visualization includes embedded statistical details (test statistics, p-values, effect sizes, confidence intervals) directly in the plot, following APA gold standard reporting.

## Design Philosophy
- Every plot tells the full statistical story. No naked charts — every visualization is annotated with relevant test results, effect sizes, and confidence intervals directly on the plot.
- Plots must be more beautiful than ggstatsplot, not a copy. Use modern design: clean typography, generous whitespace, subtle gradients, smooth animations, and a cohesive color system.
- The user should never need a separate statistics table. The plot IS the report.
- Accessibility first: colorblind-safe palettes, screen reader labels, keyboard navigation.

## Technology Stack
- **Language**: TypeScript (strict mode, no `any` types)
- **Build**: Vite
- **Visualization**: D3.js for custom statistical plots (primary), with Plotly.js as fallback for 3D
- **Statistics engine**: jStat for distributions and hypothesis tests, simple-statistics for descriptive stats and regression, custom implementations where neither library covers the need
- **UI framework**: Vanilla TypeScript with Web Components, or React if complexity demands it
- **Testing**: Vitest for unit tests, Playwright for visual regression tests
- **Styling**: CSS custom properties for theming, no CSS framework dependency

## Statistical Modules to Implement

### 1. Descriptive Statistics
- Mean, median, mode, trimmed mean
- Variance, standard deviation, standard error
- Skewness, kurtosis (excess and regular)
- Percentiles, quartiles, IQR
- Confidence intervals for the mean (t-based and bootstrap)
- Shapiro-Wilk normality test, Anderson-Darling test
- **Visualizations**: Histogram with density overlay and normality curve, box-and-whisker with jittered data points, QQ plot with confidence band, ridgeline plot for multiple groups

### 2. Frequency Analysis
- Frequency tables (absolute, relative, cumulative)
- Cross-tabulation / contingency tables
- Chi-square test of independence, Fisher's exact test
- Cramér's V, phi coefficient, odds ratio with CI
- Goodness-of-fit test
- **Visualizations**: Annotated bar chart with counts and percentages, mosaic plot, pie/donut chart with test results in subtitle, grouped and stacked bar charts with significance brackets

### 3. Group Comparisons
- Independent samples t-test (Student's and Welch's)
- Paired samples t-test
- One-way ANOVA with post-hoc (Tukey HSD, Games-Howell)
- Repeated measures ANOVA
- Kruskal-Wallis test with Dunn's post-hoc
- Mann-Whitney U test
- Wilcoxon signed-rank test
- Friedman test
- Effect sizes: Cohen's d, Hedges' g, eta-squared, omega-squared, rank-biserial correlation
- **Visualizations**: Violin + box + jitter combo plots (like ggbetweenstats) with significance brackets, pairwise comparison annotations, effect size in subtitle, distribution shape visible. Raincloud plots as premium option.

### 4. Correlation Analysis
- Pearson correlation with CI
- Spearman rank correlation
- Kendall's tau
- Point-biserial correlation
- Partial correlation
- Correlation matrix with significance stars
- **Visualizations**: Scatter plot with regression line, CI band, marginal distributions, and correlation stats in subtitle (like ggscatterstats). Correlogram/heatmap with significance indicators and hierarchical clustering.

### 5. Regression Analysis
- Simple linear regression
- Multiple linear regression
- Logistic regression (binary)
- Polynomial regression
- Diagnostics: residual plots, leverage, Cook's distance, VIF
- R², adjusted R², AIC, BIC
- Coefficient table with CI and p-values
- **Visualizations**: Scatter with fitted line and prediction/confidence bands, residual diagnostic panel (residuals vs fitted, QQ of residuals, scale-location, leverage), coefficient dot-and-whisker plot (like ggcoefstats) with CI bars and labels.

### 6. Distribution Fitting
- Normal, t, chi-square, F, binomial, Poisson, exponential, uniform
- PDF, CDF, quantile function visualization
- Parameter estimation (MLE)
- Goodness-of-fit overlay on histogram
- **Visualizations**: Interactive distribution explorer with adjustable parameters, overlay of theoretical vs empirical distributions.

## Visualization Design Standards

### Every Statistical Plot Must Include
1. **Title**: Clear, descriptive title of the analysis
2. **Subtitle**: Full statistical result string, e.g. `t(48) = 2.31, p = .025, d = 0.65, 95% CI [0.08, 1.22]`
3. **Axis labels**: Clear, with units where applicable
4. **Annotations on the plot**: Significance brackets for pairwise comparisons, regression equations on scatter plots, distribution parameters on density curves
5. **Caption**: Data source, sample size, method notes
6. **Legend**: Only when needed, never redundant with axis

### Visual Style Guide
- **Color palette**: A curated set of 8-10 colors that are colorblind-safe (Okabe-Ito or similar). Define as CSS custom properties and D3 scale.
- **Typography**: System font stack with monospace for statistical annotations. Statistical results use a slightly smaller, muted font. Titles are bold and prominent.
- **Data points**: Always show individual data points where feasible (jitter for overlapping). Points are semi-transparent (opacity 0.4-0.6) so density is visible.
- **Distributions**: Violin/density shapes use soft fills with visible outlines. Box plots inside violins use contrasting but harmonious color.
- **Significance brackets**: Clean bracket lines connecting compared groups, with p-value or significance stars. Brackets stack neatly without overlap.
- **Confidence intervals**: Shaded bands (not just lines) with low opacity.
- **Grid**: Minimal. Light horizontal gridlines only. No vertical gridlines unless scatter plot.
- **Margins**: Generous. Plots should breathe. Never cramped.
- **Animation**: Subtle enter transitions for data points and lines. No gratuitous animation.
- **Export**: Every plot exportable as SVG and PNG at publication resolution (300 DPI).

### Reference Style (Improve Upon, Don't Copy)
ggstatsplot's approach: violin + box + jitter with stats in subtitle. Our improvements:
- Better typography and spacing
- Smoother violin shapes
- Raincloud option (half-violin + jitter + box)
- Interactive tooltips showing individual data point values
- Responsive design that works at any viewport size
- Dark mode support
- Animated transitions when switching between groups or analyses

## Architecture

### File Structure
```
src/
├── stats/                  # Pure statistical computation (no DOM)
│   ├── descriptive.ts      # Mean, median, variance, etc.
│   ├── frequency.ts        # Frequency tables, chi-square
│   ├── comparison.ts       # t-tests, ANOVA, nonparametric
│   ├── correlation.ts      # Pearson, Spearman, partial
│   ├── regression.ts       # Linear, logistic, diagnostics
│   ├── distributions.ts    # PDF, CDF, quantile functions
│   ├── effect-size.ts      # Cohen's d, eta-squared, etc.
│   ├── post-hoc.ts         # Tukey, Dunn, Games-Howell
│   └── utils.ts            # Shared math utilities
├── viz/                    # Visualization layer (D3-based)
│   ├── plots/
│   │   ├── violin-box.ts   # Violin + box + jitter combo
│   │   ├── scatter-stats.ts # Scatter with regression + marginals
│   │   ├── histogram.ts    # Histogram with density + stats
│   │   ├── bar-stats.ts    # Bar chart with significance
│   │   ├── correlogram.ts  # Correlation matrix heatmap
│   │   ├── coef-plot.ts    # Coefficient dot-and-whisker
│   │   ├── qq-plot.ts      # QQ plot with confidence band
│   │   ├── residual-panel.ts # Regression diagnostics
│   │   ├── raincloud.ts    # Raincloud plots
│   │   └── distribution.ts # Interactive distribution explorer
│   ├── components/
│   │   ├── axis.ts         # Custom axis rendering
│   │   ├── brackets.ts     # Significance brackets
│   │   ├── annotations.ts  # Stats text, equations, labels
│   │   ├── legend.ts       # Legend component
│   │   └── tooltip.ts      # Interactive tooltips
│   ├── themes/
│   │   ├── default.ts      # Light theme
│   │   ├── dark.ts         # Dark theme
│   │   └── print.ts        # Print/export optimized
│   └── export.ts           # SVG/PNG export
├── data/
│   ├── parser.ts           # CSV, JSON, Excel import
│   ├── transform.ts        # Filtering, grouping, reshaping
│   └── validation.ts       # Type checking, missing values
├── ui/
│   ├── app.ts              # Main application shell
│   ├── data-panel.ts       # Data import and preview
│   ├── analysis-panel.ts   # Analysis selection and configuration
│   ├── output-panel.ts     # Results display
│   └── settings.ts         # User preferences
├── types/
│   └── index.ts            # All shared TypeScript interfaces
└── tests/
    ├── stats/              # Mirror of src/stats/ with test files
    ├── viz/                # Visual regression tests
    └── fixtures/           # Test datasets
```

### Key Architecture Rules
- **stats/ is pure**: No DOM, no D3, no side effects. Pure functions that take data in and return results. Every function independently testable.
- **viz/ consumes stats/**: Visualization functions call stats functions, never the other way around.
- **types/ is the contract**: All interfaces defined centrally. A StatResult type carries the test name, statistic, p-value, effect size, CI, and formatted APA string.
- **No circular dependencies**: Use dependency injection where needed.

## Coding Standards

### TypeScript Strictness
- `strict: true` in tsconfig, no exceptions
- No `any` types. Use `unknown` with type guards if type is truly unknown.
- All functions have explicit return types
- All function parameters have explicit types
- Use discriminated unions for result types (e.g. `{ type: 'parametric', ... } | { type: 'nonparametric', ... }`)
- Prefer `readonly` arrays and properties where mutation is not needed

### Naming
- Functions: `camelCase`, descriptive verbs — `computeMean()`, `runTTest()`, `renderViolin()`
- Types/Interfaces: `PascalCase` — `StatResult`, `ViolinPlotConfig`, `AnalysisOptions`
- Constants: `UPPER_SNAKE_CASE` — `DEFAULT_ALPHA`, `COLOR_PALETTE`
- Files: `kebab-case.ts`
- Test files: `*.test.ts` matching source file name

### Statistical Computation Rules
- Always validate inputs at function boundaries: check for empty arrays, NaN, Infinity, mismatched lengths.
- Return structured results, never just a number. Every stat function returns an object with at least: `{ statistic, pValue, effectSize, ci, formatted }` where `formatted` is the APA-style string.
- Use existing jStat/simple-statistics functions where they exist and are correct. Verify their output against known values before trusting.
- When implementing from scratch, cite the formula source in a comment.
- All p-values should support multiple comparison correction (Bonferroni, Holm, BH/FDR).
- Handle edge cases: n < 3, zero variance, perfect correlation, singular matrices.

### Testing Requirements
- Every stats function gets a test file with:
  - Known textbook examples with expected values (cite the textbook/source)
  - Edge cases: empty input, single value, identical values, very large/small numbers
  - Cross-validation against R or Python (document the R/Python code used to generate expected values in a comment)
  - Numerical tolerance: use `toBeCloseTo()` with appropriate decimal places
- Every visualization gets:
  - A unit test verifying the correct SVG elements are generated
  - A visual regression test capturing a screenshot
  - A test with synthetic data verifying annotations show correct values
- Minimum 90% code coverage for stats/. No merging below this.
- Run all tests after every change. Fix failures before moving on.

## Core Principle (Inherited from Global Rules)
- NEVER assume a jStat or simple-statistics function returns what you expect. Run it with synthetic data first, inspect the output, then use it.
- After writing any module, create synthetic test data, run it, inspect the structure, verify correctness against R or Python output.
- Always create tests. No function ships without tests. No exceptions.

## Cross-Validation Protocol
For every statistical function implemented:
1. Write the TypeScript implementation.
2. Create identical test data.
3. Run the same analysis in R (document the R code in a comment block in the test file).
4. Compare results numerically with tolerance (typically 1e-6 for most stats, 1e-4 for iterative methods).
5. If results differ, investigate and document why (different default parameters, different algorithms, etc.) in LEARNINGS.md.
6. The test must encode the R-verified expected values.

Example test comment:
```typescript
// Cross-validated with R:
// > t.test(c(1,2,3,4,5), c(6,7,8,9,10), var.equal = FALSE)
// t = -5.0, df = 8, p-value = 0.001053
// 95% CI: [-6.306, -2.694]
// Cohen's d (effsize::cohen.d): -3.162
```

## Numerical Equivalence Testing (MANDATORY)
Every new statistical implementation MUST have a dedicated numerical equivalence test against R. This is non-negotiable — no implementation is considered complete without it.

### Process
1. Write an R script that runs the equivalent analysis on synthetic test data and outputs all key quantities (statistic, df, p-value, effect size, CI, coefficients, etc.) to a JSON file in `tests/fixtures/`.
2. Write a Vitest test that loads the JSON fixture and compares each value from the Carm implementation against R's output using `toBeCloseTo()` with appropriate tolerance.
3. The R script and its invocation must be documented in a comment block at the top of the test file.
4. Cover at minimum: the main test statistic, degrees of freedom, p-value, effect size, and confidence intervals. For regression models: all coefficients, SEs, AIC/BIC, log-likelihood.

### Tolerances
- Deterministic algorithms (OLS, chi-square, t-test): `1e-6` or tighter
- Iterative algorithms (IRLS, EM, Newton): `1e-4` (may need `1e-2` for convergence-sensitive quantities)
- P-values near 0 or 1: use absolute tolerance `1e-4`
- Random/stochastic methods (bootstrap): compare distribution properties, not exact values

### Fixture Format
```json
{
  "method": "welchANOVA",
  "r_code": "oneway.test(y ~ group, data = df, var.equal = FALSE)",
  "data": { "groups": [[1,2,3], [4,5,6], [7,8,9]] },
  "expected": {
    "statistic": 9.9391,
    "df_num": 2,
    "df_den": 9.8506,
    "pValue": 0.004338
  }
}
```

### When to Run
- After implementing any new statistical function
- After modifying the internals of an existing statistical function
- Before marking any stats task as complete

## Visualization Testing Protocol
- You CANNOT see rendered plots. Never assume a plot looks correct.
- After creating or modifying any visualization:
  1. Generate a test HTML page that renders the plot with synthetic data.
  2. Save it to ./tmp/ with a descriptive name.
  3. Ask the user: "Test page saved to ./tmp/[name].html. Want to review it, or skip and continue?"
  4. Do NOT delete test HTML files from ./tmp/ unless the user agrees.
- This applies to ANY visual change: colors, spacing, fonts, annotations, axis labels, brackets, layout.

## Reporting Changes
- After every change, explain in plain English what was done and why.
- Write as if explaining to a colleague who hasn't seen the code.
- Use text-based illustrations:
  - Before/after comparisons showing old vs new behavior
  - ASCII diagrams for data flow or architecture changes
  - Example input → output showing how behavior differs
  - Key code snippets (not entire files)
- Include: what changed, why, how it affects existing behavior, side effects to watch for.
- Example:
  ```
  Added significance brackets to the violin-box plot.

  Before: Plot showed violins + boxes + jitter with stats only in subtitle.
    [violin] [box] [dots]
    subtitle: F(2,87) = 4.21, p = .018

  After: Plot now shows pairwise comparison brackets above groups.
    [violin] [box] [dots]
       |------*------|     ← bracket with p = .012
              |--ns--|     ← bracket with p = .340
    subtitle: F(2,87) = 4.21, p = .018, η² = 0.09

  Files changed:
    - src/viz/plots/violin-box.ts: added bracket rendering logic
    - src/viz/components/brackets.ts: new file, bracket layout algorithm
    - tests/viz/violin-box.test.ts: added bracket position tests
  ```

## Project Learning Log
- Maintain LEARNINGS.md in project root. Read it at the start of every task.
- After completing any task, append discoveries about:
  - jStat/simple-statistics quirks (wrong defaults, missing features, undocumented behavior)
  - D3 rendering gotchas (SVG vs Canvas, text measurement, responsive sizing)
  - Cross-validation discrepancies between TypeScript and R/Python
  - Browser compatibility issues
  - Performance bottlenecks with large datasets
- Format: `### YYYY-MM-DD\n- [topic]: what was learned`

## Session Handoff
- At end of session or when context is long, write HANDOFF.md with:
  - Completed: what was done, which files
  - Current State: what works, what's broken, what's partial
  - Key Decisions: architecture/design choices and reasoning
  - Open Issues: bugs, edge cases, unresolved questions
  - Next Steps: prioritized list with enough detail to start immediately
  - Context: environment details, dependencies, prerequisites
- At start of every session, read HANDOFF.md and LEARNINGS.md first.

## Change Summary
- After multi-file tasks, update CHANGES.md (newest first):
  ```
  ### YYYY-MM-DD — [brief description]
  - file.ts: what changed and why
  - Tests: what was tested, pass/fail
  ```

## Task Completion Checklist
When completing any task:
1. Write the code.
2. Write tests for all new/modified functions.
3. Run all tests and verify they pass.
4. Generate test HTML for any visual changes and present to user.
5. Update LEARNINGS.md with anything discovered.
6. Update CHANGES.md with what was modified.
7. Update HANDOFF.md with current project state.
8. Report changes in plain English with illustrations.
9. Only then report the task as complete.
No shortcuts. These are part of the task definition.

## Workflow
- Complete the entire task before asking for user input or confirmation.
- Do not stop after each small step. Keep going until done.
- Only ask for clarification if genuinely blocked.

## Git
- Do not run any git commands (push, commit, add, merge, rebase, reset) without explicitly asking the user first.
- Only ask about git once the entire task is fully complete.
- Present a summary of all changes before asking whether to commit.
