#!/usr/bin/env Rscript
# Generate reference values for distribution functions cross-validation.
# Run: Rscript validation/r-reference/distributions-ref.R
# Output: tests/fixtures/distributions-ref.json

library(jsonlite)

ref <- list()

# ─── PDFs ────────────────────────────────────────────────────────────────────

ref$normalPDF <- list(
  list(x = 0,   mu = 0, sigma = 1,   expected = dnorm(0, 0, 1)),
  list(x = 1,   mu = 0, sigma = 1,   expected = dnorm(1, 0, 1)),
  list(x = -2,  mu = 0, sigma = 1,   expected = dnorm(-2, 0, 1)),
  list(x = 3,   mu = 5, sigma = 2,   expected = dnorm(3, 5, 2)),
  list(x = 0.5, mu = 0, sigma = 0.5, expected = dnorm(0.5, 0, 0.5))
)

ref$tPDF <- list(
  list(x = 0,  df = 1,  expected = dt(0, 1)),
  list(x = 1,  df = 5,  expected = dt(1, 5)),
  list(x = -2, df = 10, expected = dt(-2, 10)),
  list(x = 3,  df = 30, expected = dt(3, 30)),
  list(x = 0,  df = 2,  expected = dt(0, 2))
)

ref$chiSqPDF <- list(
  list(x = 1,  df = 1,  expected = dchisq(1, 1)),
  list(x = 3,  df = 3,  expected = dchisq(3, 3)),
  list(x = 5,  df = 5,  expected = dchisq(5, 5)),
  list(x = 10, df = 10, expected = dchisq(10, 10)),
  list(x = 0,  df = 3,  expected = dchisq(0, 3))
)

ref$fPDF <- list(
  list(x = 1, df1 = 5, df2 = 10, expected = df(1, 5, 10)),
  list(x = 2, df1 = 3, df2 = 20, expected = df(2, 3, 20)),
  list(x = 0.5, df1 = 10, df2 = 5, expected = df(0.5, 10, 5)),
  list(x = 3, df1 = 2, df2 = 50, expected = df(3, 2, 50))
)

ref$exponentialPDF <- list(
  list(x = 0, rate = 1,   expected = dexp(0, 1)),
  list(x = 1, rate = 1,   expected = dexp(1, 1)),
  list(x = 2, rate = 0.5, expected = dexp(2, 0.5)),
  list(x = 0.5, rate = 3, expected = dexp(0.5, 3))
)

ref$uniformPDF <- list(
  list(x = 0.5, min = 0, max = 1, expected = dunif(0.5, 0, 1)),
  list(x = 3, min = 2, max = 5, expected = dunif(3, 2, 5)),
  list(x = 0, min = 0, max = 1, expected = dunif(0, 0, 1)),
  list(x = -1, min = 0, max = 1, expected = dunif(-1, 0, 1))
)

ref$gammaPDF <- list(
  list(x = 1, shape = 1, rate = 1, expected = dgamma(1, 1, 1)),
  list(x = 2, shape = 2, rate = 1, expected = dgamma(2, 2, 1)),
  list(x = 3, shape = 5, rate = 2, expected = dgamma(3, 5, 2)),
  list(x = 0.5, shape = 0.5, rate = 1, expected = dgamma(0.5, 0.5, 1)),
  list(x = 10, shape = 3, rate = 0.5, expected = dgamma(10, 3, 0.5))
)

ref$betaPDF <- list(
  list(x = 0.5, shape1 = 1, shape2 = 1, expected = dbeta(0.5, 1, 1)),
  list(x = 0.3, shape1 = 2, shape2 = 5, expected = dbeta(0.3, 2, 5)),
  list(x = 0.8, shape1 = 5, shape2 = 2, expected = dbeta(0.8, 5, 2)),
  list(x = 0.5, shape1 = 0.5, shape2 = 0.5, expected = dbeta(0.5, 0.5, 0.5)),
  list(x = 0.1, shape1 = 10, shape2 = 3, expected = dbeta(0.1, 10, 3))
)

# ─── CDFs ────────────────────────────────────────────────────────────────────

ref$exponentialCDF <- list(
  list(x = 0, rate = 1, expected = pexp(0, 1)),
  list(x = 1, rate = 1, expected = pexp(1, 1)),
  list(x = 2, rate = 0.5, expected = pexp(2, 0.5)),
  list(x = 0.5, rate = 3, expected = pexp(0.5, 3))
)

ref$uniformCDF <- list(
  list(x = 0.5, min = 0, max = 1, expected = punif(0.5, 0, 1)),
  list(x = 3, min = 2, max = 5, expected = punif(3, 2, 5)),
  list(x = -1, min = 0, max = 1, expected = punif(-1, 0, 1)),
  list(x = 2, min = 0, max = 1, expected = punif(2, 0, 1))
)

ref$gammaCDF <- list(
  list(x = 1, shape = 1, rate = 1, expected = pgamma(1, 1, 1)),
  list(x = 2, shape = 2, rate = 1, expected = pgamma(2, 2, 1)),
  list(x = 5, shape = 5, rate = 2, expected = pgamma(5, 5, 2)),
  list(x = 0.5, shape = 0.5, rate = 1, expected = pgamma(0.5, 0.5, 1))
)

ref$betaCDF <- list(
  list(x = 0.5, shape1 = 1, shape2 = 1, expected = pbeta(0.5, 1, 1)),
  list(x = 0.3, shape1 = 2, shape2 = 5, expected = pbeta(0.3, 2, 5)),
  list(x = 0.8, shape1 = 5, shape2 = 2, expected = pbeta(0.8, 5, 2)),
  list(x = 0.1, shape1 = 10, shape2 = 3, expected = pbeta(0.1, 10, 3))
)

# ─── Quantiles ───────────────────────────────────────────────────────────────

ref$exponentialQuantile <- list(
  list(p = 0.5, rate = 1, expected = qexp(0.5, 1)),
  list(p = 0.95, rate = 2, expected = qexp(0.95, 2)),
  list(p = 0.1, rate = 0.5, expected = qexp(0.1, 0.5))
)

ref$uniformQuantile <- list(
  list(p = 0.5, min = 0, max = 1, expected = qunif(0.5, 0, 1)),
  list(p = 0.25, min = 2, max = 8, expected = qunif(0.25, 2, 8))
)

ref$gammaQuantile <- list(
  list(p = 0.5, shape = 1, rate = 1, expected = qgamma(0.5, 1, 1)),
  list(p = 0.95, shape = 5, rate = 2, expected = qgamma(0.95, 5, 2)),
  list(p = 0.1, shape = 2, rate = 0.5, expected = qgamma(0.1, 2, 0.5))
)

ref$betaQuantile <- list(
  list(p = 0.5, shape1 = 2, shape2 = 5, expected = qbeta(0.5, 2, 5)),
  list(p = 0.95, shape1 = 5, shape2 = 2, expected = qbeta(0.95, 5, 2)),
  list(p = 0.1, shape1 = 10, shape2 = 3, expected = qbeta(0.1, 10, 3))
)

ref$fQuantile <- list(
  list(p = 0.5, df1 = 5, df2 = 10, expected = qf(0.5, 5, 10)),
  list(p = 0.95, df1 = 3, df2 = 20, expected = qf(0.95, 3, 20)),
  list(p = 0.99, df1 = 10, df2 = 50, expected = qf(0.99, 10, 50)),
  list(p = 0.025, df1 = 2, df2 = 30, expected = qf(0.025, 2, 30))
)

# ─── Discrete PMFs ───────────────────────────────────────────────────────────

ref$binomialPMF <- list(
  list(k = 3, n = 10, prob = 0.5, expected = dbinom(3, 10, 0.5)),
  list(k = 0, n = 5, prob = 0.3, expected = dbinom(0, 5, 0.3)),
  list(k = 10, n = 10, prob = 0.9, expected = dbinom(10, 10, 0.9)),
  list(k = 5, n = 20, prob = 0.25, expected = dbinom(5, 20, 0.25))
)

ref$binomialCDF <- list(
  list(k = 3, n = 10, prob = 0.5, expected = pbinom(3, 10, 0.5)),
  list(k = 0, n = 5, prob = 0.3, expected = pbinom(0, 5, 0.3)),
  list(k = 7, n = 10, prob = 0.5, expected = pbinom(7, 10, 0.5)),
  list(k = 5, n = 20, prob = 0.25, expected = pbinom(5, 20, 0.25))
)

ref$binomialQuantile <- list(
  list(p = 0.5, n = 10, prob = 0.5, expected = qbinom(0.5, 10, 0.5)),
  list(p = 0.95, n = 20, prob = 0.3, expected = qbinom(0.95, 20, 0.3)),
  list(p = 0.1, n = 50, prob = 0.5, expected = qbinom(0.1, 50, 0.5))
)

ref$poissonPMF <- list(
  list(k = 3, lambda = 5, expected = dpois(3, 5)),
  list(k = 0, lambda = 1, expected = dpois(0, 1)),
  list(k = 10, lambda = 10, expected = dpois(10, 10)),
  list(k = 5, lambda = 2.5, expected = dpois(5, 2.5))
)

ref$poissonCDF <- list(
  list(k = 3, lambda = 5, expected = ppois(3, 5)),
  list(k = 0, lambda = 1, expected = ppois(0, 1)),
  list(k = 10, lambda = 10, expected = ppois(10, 10)),
  list(k = 5, lambda = 2.5, expected = ppois(5, 2.5))
)

ref$poissonQuantile <- list(
  list(p = 0.5, lambda = 5, expected = qpois(0.5, 5)),
  list(p = 0.95, lambda = 10, expected = qpois(0.95, 10)),
  list(p = 0.1, lambda = 2.5, expected = qpois(0.1, 2.5))
)

# ─── MLE Fitting ─────────────────────────────────────────────────────────────

# Normal
set.seed(42)
norm_data <- rnorm(100, mean = 5, sd = 2)
norm_fit <- MASS::fitdistr(norm_data, "normal")
ref$fitNormal <- list(
  data = norm_data,
  mu = unname(norm_fit$estimate["mean"]),
  sigma = unname(norm_fit$estimate["sd"]),
  logLik = as.numeric(logLik(norm_fit)),
  aic = AIC(norm_fit),
  bic = BIC(norm_fit)
)

# Exponential
set.seed(42)
exp_data <- rexp(100, rate = 2)
exp_fit <- MASS::fitdistr(exp_data, "exponential")
ref$fitExponential <- list(
  data = exp_data,
  rate = unname(exp_fit$estimate["rate"]),
  logLik = as.numeric(logLik(exp_fit)),
  aic = AIC(exp_fit),
  bic = BIC(exp_fit)
)

# Gamma
set.seed(42)
gamma_data <- rgamma(100, shape = 3, rate = 2)
gamma_fit <- MASS::fitdistr(gamma_data, "gamma")
ref$fitGamma <- list(
  data = gamma_data,
  shape = unname(gamma_fit$estimate["shape"]),
  rate = unname(gamma_fit$estimate["rate"]),
  logLik = as.numeric(logLik(gamma_fit)),
  aic = AIC(gamma_fit),
  bic = BIC(gamma_fit)
)

# Beta — MASS::fitdistr doesn't support beta, use manual MLE
set.seed(42)
beta_data <- rbeta(100, shape1 = 2, shape2 = 5)
# Fit via optim
beta_nll <- function(par) {
  -sum(dbeta(beta_data, par[1], par[2], log = TRUE))
}
beta_opt <- optim(c(2, 5), beta_nll, method = "L-BFGS-B", lower = c(0.01, 0.01))
ref$fitBeta <- list(
  data = beta_data,
  shape1 = beta_opt$par[1],
  shape2 = beta_opt$par[2],
  logLik = -beta_opt$value,
  aic = 2 * beta_opt$value + 2 * 2,
  bic = 2 * beta_opt$value + 2 * log(100)
)

# Poisson
set.seed(42)
pois_data <- rpois(100, lambda = 4.5)
# MLE for Poisson: lambda = mean
ref$fitPoisson <- list(
  data = pois_data,
  lambda = mean(pois_data),
  logLik = sum(dpois(pois_data, mean(pois_data), log = TRUE))
)

# ─── Goodness-of-Fit Tests ───────────────────────────────────────────────────

# Anderson-Darling (normal case — uses nortest::ad.test)
set.seed(42)
ad_data <- rnorm(50, mean = 10, sd = 3)
ad_result <- nortest::ad.test(ad_data)
ref$andersonDarling <- list(
  data = ad_data,
  statistic = unname(ad_result$statistic),
  pValue = ad_result$p.value
)

# Anderson-Darling on non-normal data
set.seed(42)
ad_exp_data <- rexp(50, rate = 1)
ad_exp_result <- nortest::ad.test(ad_exp_data)
ref$andersonDarlingNonNormal <- list(
  data = ad_exp_data,
  statistic = unname(ad_exp_result$statistic),
  pValue = ad_exp_result$p.value
)

# Kolmogorov-Smirnov (known params)
set.seed(42)
ks_data <- rnorm(50, mean = 0, sd = 1)
ks_result <- ks.test(ks_data, "pnorm", mean = 0, sd = 1)
ref$kolmogorovSmirnov <- list(
  data = ks_data,
  statistic = unname(ks_result$statistic),
  pValue = ks_result$p.value
)

# KS on exponential data tested against normal
set.seed(42)
ks_exp_data <- rexp(50, rate = 1)
# Test against exponential with known params
ks_exp_result <- ks.test(ks_exp_data, "pexp", rate = 1)
ref$kolmogorovSmirnovExp <- list(
  data = ks_exp_data,
  statistic = unname(ks_exp_result$statistic),
  pValue = ks_exp_result$p.value
)

# Write JSON
output_path <- "tests/fixtures/distributions-ref.json"
writeLines(toJSON(ref, digits = 15, auto_unbox = TRUE, pretty = TRUE), output_path)
cat("Written:", output_path, "\n")
