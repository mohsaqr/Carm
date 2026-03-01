#!/usr/bin/env Rscript
# Generate 100 synthetic GLMM datasets, fit with lme4 and glmmTMB,
# output all reference values to JSON for cross-validation.
#
# Usage: Rscript validation/r-reference/glmm-100-datasets.R
# Output: tests/fixtures/glmm-100-ref.json

suppressWarnings(suppressMessages({
  library(lme4)
  library(glmmTMB)
  library(jsonlite)
}))

set.seed(2026)

results <- list()

for (ds in 1:100) {
  # Vary parameters across datasets
  n_per_group <- sample(10:30, 1)
  n_groups <- sample(5:20, 1)
  n <- n_per_group * n_groups

  beta0 <- runif(1, -2, 2)
  beta1 <- runif(1, -1.5, 1.5)
  sigma_b <- runif(1, 0.3, 2.5)

  # Generate data
  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n)
  b <- rnorm(n_groups, sd = sigma_b)
  eta <- beta0 + beta1 * x + b[group]
  y <- rbinom(n, 1, plogis(eta))

  # Skip degenerate cases (all 0 or all 1)
  if (length(unique(y)) < 2) next

  df <- data.frame(y = y, x = x, group = factor(group))

  # Fit lme4
  m_lme4 <- NULL
  tryCatch(
    m_lme4 <- suppressWarnings(glmer(y ~ x + (1|group), data = df, family = binomial,
                                      control = glmerControl(optimizer = "bobyqa"))),
    error = function(e) {}
  )
  if (is.null(m_lme4)) next

  # Fit glmmTMB
  m_tmb <- NULL
  tryCatch(
    m_tmb <- suppressWarnings(glmmTMB(y ~ x + (1|group), data = df, family = binomial)),
    error = function(e) {}
  )
  if (is.null(m_tmb)) next

  s_lme4 <- summary(m_lme4)
  s_tmb <- summary(m_tmb)

  results[[length(results) + 1]] <- list(
    id = ds,
    params = list(
      n = n,
      nGroups = n_groups,
      nPerGroup = n_per_group,
      trueBeta0 = beta0,
      trueBeta1 = beta1,
      trueSigmaB = sigma_b
    ),
    data = list(
      y = as.numeric(y),
      x = as.numeric(x),
      group = as.numeric(group)
    ),
    lme4 = list(
      intercept = as.numeric(fixef(m_lme4)[1]),
      x = as.numeric(fixef(m_lme4)[2]),
      interceptSE = as.numeric(s_lme4$coefficients[1, 2]),
      xSE = as.numeric(s_lme4$coefficients[2, 2]),
      variance = as.numeric(VarCorr(m_lme4)$group[1, 1]),
      logLik = as.numeric(logLik(m_lme4)),
      aic = as.numeric(AIC(m_lme4))
    ),
    glmmTMB = list(
      intercept = as.numeric(fixef(m_tmb)$cond[1]),
      x = as.numeric(fixef(m_tmb)$cond[2]),
      interceptSE = as.numeric(s_tmb$coefficients$cond[1, 2]),
      xSE = as.numeric(s_tmb$coefficients$cond[2, 2]),
      variance = as.numeric(VarCorr(m_tmb)$cond$group[1, 1]),
      logLik = as.numeric(logLik(m_tmb)),
      aic = as.numeric(AIC(m_tmb))
    )
  )

  if (ds %% 20 == 0) cat("  done", ds, "/", 100, "\n")
}

cat("Generated", length(results), "datasets\n")

writeLines(toJSON(results, digits = 10, auto_unbox = TRUE),
           "tests/fixtures/glmm-100-ref.json")
cat("Written to tests/fixtures/glmm-100-ref.json\n")
