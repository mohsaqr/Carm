#!/usr/bin/env Rscript
# Round 4: Regression Suite — R Reference Generator
# Generates reference values for 500 datasets using Saqrlab::simulate_data("prediction")
# Compares: simple linear, multiple, polynomial, logistic regression, diagnostics (leverage, Cook's D, VIF)
#
# Usage: Rscript validation/r-reference/regression-ref.R

suppressPackageStartupMessages({
  library(Saqrlab)
  library(jsonlite)
  library(car)
})

N <- 500
cat("Generating", N, "regression reference datasets...\n")
results <- vector("list", N)

for (i in seq_len(N)) {
  d <- tryCatch(simulate_data("prediction", seed = i), error = function(e) NULL)
  if (is.null(d)) { results[[i]] <- list(seed = i, error = "simulate_failed"); next }

  y  <- as.numeric(d$y)
  x1 <- as.numeric(d$x1)
  x2 <- as.numeric(d$x2)
  x3 <- as.numeric(d$x3)
  x4 <- as.numeric(d$x4)
  n  <- length(y)

  # ── Simple linear regression: y ~ x1 ─────────────────────────────────
  simple_res <- tryCatch({
    fit <- lm(y ~ x1)
    s <- summary(fit)
    cf <- s$coefficients
    ci <- confint(fit)
    f_stat <- s$fstatistic  # c(value, numdf, dendf)
    f_pval <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)

    list(
      intercept = cf[1, "Estimate"],
      slope = cf[2, "Estimate"],
      interceptSE = cf[1, "Std. Error"],
      slopeSE = cf[2, "Std. Error"],
      interceptT = cf[1, "t value"],
      slopeT = cf[2, "t value"],
      interceptP = cf[1, "Pr(>|t|)"],
      slopeP = cf[2, "Pr(>|t|)"],
      r2 = s$r.squared,
      adjR2 = s$adj.r.squared,
      fStatistic = f_stat[[1]],
      fPvalue = f_pval,
      interceptCILower = ci[1, 1],
      interceptCIUpper = ci[1, 2],
      slopeCILower = ci[2, 1],
      slopeCIUpper = ci[2, 2]
    )
  }, error = function(e) list(error = e$message))

  # ── Multiple regression: y ~ x1 + x2 + x3 + x4 ─────────────────────
  mult_res <- tryCatch({
    fit <- lm(y ~ x1 + x2 + x3 + x4)
    s <- summary(fit)
    cf <- s$coefficients
    f_stat <- s$fstatistic
    f_pval <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)

    list(
      coefficients = as.numeric(cf[, "Estimate"]),
      se = as.numeric(cf[, "Std. Error"]),
      t = as.numeric(cf[, "t value"]),
      p = as.numeric(cf[, "Pr(>|t|)"]),
      r2 = s$r.squared,
      adjR2 = s$adj.r.squared,
      fStatistic = f_stat[[1]],
      fPvalue = f_pval,
      aic = AIC(fit),
      bic = BIC(fit),
      coeffNames = rownames(cf)
    )
  }, error = function(e) list(error = e$message))

  # ── Polynomial regression: y ~ poly(x1, 2, raw=TRUE) ────────────────
  poly_res <- tryCatch({
    fit <- lm(y ~ poly(x1, 2, raw = TRUE))
    s <- summary(fit)
    cf <- s$coefficients

    list(
      coefficients = as.numeric(cf[, "Estimate"]),
      se = as.numeric(cf[, "Std. Error"]),
      t = as.numeric(cf[, "t value"]),
      p = as.numeric(cf[, "Pr(>|t|)"]),
      r2 = s$r.squared,
      adjR2 = s$adj.r.squared,
      coeffNames = rownames(cf)
    )
  }, error = function(e) list(error = e$message))

  # ── Logistic regression: ybin ~ x1 + x2 + x3 + x4 ──────────────────
  logistic_res <- tryCatch({
    ybin <- as.integer(y > median(y))
    fit <- glm(ybin ~ x1 + x2 + x3 + x4, family = binomial)
    s <- summary(fit)
    cf <- s$coefficients

    list(
      coefficients = as.numeric(cf[, "Estimate"]),
      se = as.numeric(cf[, "Std. Error"]),
      z = as.numeric(cf[, "z value"]),
      p = as.numeric(cf[, "Pr(>|z|)"]),
      aic = fit$aic,
      converged = fit$converged,
      coeffNames = rownames(cf)
    )
  }, error = function(e) list(error = e$message))

  # ── Diagnostics (from multiple regression) ───────────────────────────
  diag_res <- tryCatch({
    fit <- lm(y ~ x1 + x2 + x3 + x4)
    lev <- hatvalues(fit)
    cook <- cooks.distance(fit)
    vif_vals <- car::vif(fit)

    list(
      leverage = as.numeric(lev[seq_len(min(10, n))]),
      cooksD = as.numeric(cook[seq_len(min(10, n))]),
      vif = as.numeric(vif_vals),
      vifNames = names(vif_vals)
    )
  }, error = function(e) list(error = e$message))

  results[[i]] <- list(
    seed = i,
    n = n,
    y = y,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    simple = simple_res,
    multiple = mult_res,
    polynomial = poly_res,
    logistic = logistic_res,
    diagnostics = diag_res
  )

  if (i %% 100 == 0) cat("  Completed", i, "/", N, "\n")
}

outPath <- file.path("validation", "data", "regression-ref.json")
dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(results, digits = 12, auto_unbox = TRUE, na = "null"), outPath)
cat("Wrote", outPath, "(", sum(vapply(results, function(r) is.null(r$error), logical(1))), "valid datasets)\n")
