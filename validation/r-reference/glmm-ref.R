#!/usr/bin/env Rscript
# Generate reference values for logistic GLMM numerical equivalence tests.
# Uses lme4::glmer() with family = binomial (Laplace approximation).
#
# Usage: Rscript validation/r-reference/glmm-ref.R
# Output: tests/fixtures/glmm-ref.json

library(lme4)
library(jsonlite)

set.seed(42)

# ─── Test 1: Basic random intercept model ─────────────────────────────────
# 200 obs, 10 groups, 1 fixed predictor

n <- 200
nG <- 10
group <- rep(1:nG, each = n/nG)

# Use a simple deterministic PRNG approach for reproducibility
# (matches the LCG used in our TypeScript tests)
x <- rnorm(n)
b <- rnorm(nG, sd = 1.0)  # random intercepts
eta <- -0.5 + 0.8 * x + b[group]
y <- rbinom(n, 1, plogis(eta))

df <- data.frame(y = y, x = x, group = factor(group))

mod <- glmer(y ~ x + (1|group), data = df, family = binomial)
s <- summary(mod)

# Extract all quantities
test1 <- list(
  data = list(
    y = as.numeric(y),
    x = as.numeric(x),
    group = as.numeric(group)
  ),
  fixedEffects = list(
    intercept = list(
      estimate = as.numeric(fixef(mod)[1]),
      se = as.numeric(s$coefficients[1, 2]),
      z = as.numeric(s$coefficients[1, 3]),
      p = as.numeric(s$coefficients[1, 4])
    ),
    x = list(
      estimate = as.numeric(fixef(mod)[2]),
      se = as.numeric(s$coefficients[2, 2]),
      z = as.numeric(s$coefficients[2, 3]),
      p = as.numeric(s$coefficients[2, 4])
    )
  ),
  varianceComponents = list(
    intercept = as.numeric(VarCorr(mod)$group[1, 1])
  ),
  logLik = as.numeric(logLik(mod)),
  aic = as.numeric(AIC(mod)),
  bic = as.numeric(BIC(mod)),
  deviance = as.numeric(deviance(mod)),
  nObs = nobs(mod),
  nGroups = nG
)

# ─── Test 2: Multiple predictors ──────────────────────────────────────────
# 200 obs, 10 groups, 2 fixed predictors

x1 <- rnorm(n)
x2 <- rnorm(n)
b2 <- rnorm(nG, sd = 0.8)
eta2 <- 0.3 + 1.0 * x1 - 0.5 * x2 + b2[group]
y2 <- rbinom(n, 1, plogis(eta2))

df2 <- data.frame(y = y2, x1 = x1, x2 = x2, group = factor(group))
mod2 <- glmer(y ~ x1 + x2 + (1|group), data = df2, family = binomial)
s2 <- summary(mod2)

test2 <- list(
  data = list(
    y = as.numeric(y2),
    x1 = as.numeric(x1),
    x2 = as.numeric(x2),
    group = as.numeric(group)
  ),
  fixedEffects = list(
    intercept = list(
      estimate = as.numeric(fixef(mod2)[1]),
      se = as.numeric(s2$coefficients[1, 2]),
      z = as.numeric(s2$coefficients[1, 3]),
      p = as.numeric(s2$coefficients[1, 4])
    ),
    x1 = list(
      estimate = as.numeric(fixef(mod2)[2]),
      se = as.numeric(s2$coefficients[2, 2]),
      z = as.numeric(s2$coefficients[2, 3]),
      p = as.numeric(s2$coefficients[2, 4])
    ),
    x2 = list(
      estimate = as.numeric(fixef(mod2)[3]),
      se = as.numeric(s2$coefficients[3, 2]),
      z = as.numeric(s2$coefficients[3, 3]),
      p = as.numeric(s2$coefficients[3, 4])
    )
  ),
  varianceComponents = list(
    intercept = as.numeric(VarCorr(mod2)$group[1, 1])
  ),
  logLik = as.numeric(logLik(mod2)),
  aic = as.numeric(AIC(mod2)),
  bic = as.numeric(BIC(mod2)),
  deviance = as.numeric(deviance(mod2)),
  nObs = nobs(mod2),
  nGroups = nG
)

# ─── Write output ─────────────────────────────────────────────────────────

result <- list(
  test1 = test1,
  test2 = test2
)

output_path <- "tests/fixtures/glmm-ref.json"
writeLines(toJSON(result, digits = 10, auto_unbox = TRUE), output_path)
cat("Written to", output_path, "\n")
