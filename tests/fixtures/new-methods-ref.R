#!/usr/bin/env Rscript
# Generate R reference values for 14 new statistical methods.
# Run from repo root: Rscript tests/fixtures/new-methods-ref.R
# Output: tests/fixtures/new-methods-ref.json

library(jsonlite)

results <- list()

# ─── 1. Welch's ANOVA ───────────────────────────────────────────────────
cat("1. Welch's ANOVA...\n")
g1 <- c(2, 3, 7, 2, 6)
g2 <- c(10, 8, 12, 6, 11)
g3 <- c(8, 7, 5, 10, 9)
df <- data.frame(
  y = c(g1, g2, g3),
  group = factor(rep(c("A","B","C"), each=5))
)
w <- oneway.test(y ~ group, data = df, var.equal = FALSE)
results$welchANOVA <- list(
  data = list(g1 = g1, g2 = g2, g3 = g3),
  statistic = unname(w$statistic),
  df_num = unname(w$parameter[1]),
  df_den = unname(w$parameter[2]),
  pValue = w$p.value
)

# Also with heterogeneous variances
g1h <- c(1, 2, 3)
g2h <- c(10, 20, 30, 40, 50)
g3h <- c(5, 6, 7, 8)
dfh <- data.frame(
  y = c(g1h, g2h, g3h),
  group = factor(rep(c("A","B","C"), c(3,5,4)))
)
wh <- oneway.test(y ~ group, data = dfh, var.equal = FALSE)
results$welchANOVA_hetero <- list(
  data = list(g1 = g1h, g2 = g2h, g3 = g3h),
  statistic = unname(wh$statistic),
  df_num = unname(wh$parameter[1]),
  df_den = unname(wh$parameter[2]),
  pValue = wh$p.value
)

# ─── 2. Mood's Median Test ──────────────────────────────────────────────
cat("2. Mood's Median Test...\n")
# R doesn't have a direct mood.test for k-sample; compute manually
m_g1 <- c(1, 2, 3, 4, 5)
m_g2 <- c(6, 7, 8, 9, 10)
m_g3 <- c(3, 4, 5, 6, 7)
all_vals <- c(m_g1, m_g2, m_g3)
grand_median <- median(all_vals)
# Count above vs at-or-below
above <- c(sum(m_g1 > grand_median), sum(m_g2 > grand_median), sum(m_g3 > grand_median))
below <- c(sum(m_g1 <= grand_median), sum(m_g2 <= grand_median), sum(m_g3 <= grand_median))
tbl <- rbind(above, below)
ct <- chisq.test(tbl, correct = FALSE)
results$moodsMedian <- list(
  data = list(g1 = m_g1, g2 = m_g2, g3 = m_g3),
  grandMedian = grand_median,
  table_above = above,
  table_below = below,
  statistic = unname(ct$statistic),
  df = unname(ct$parameter),
  pValue = ct$p.value
)

# ─── 3. Cochran's Q Test ───────────────────────────────────────────────
cat("3. Cochran's Q Test...\n")
# Use manual computation matching the formula
cq_data <- matrix(c(
  1, 1, 0,
  1, 0, 0,
  1, 1, 1,
  0, 0, 0,
  1, 1, 0,
  1, 0, 1,
  0, 1, 0,
  1, 1, 0
), nrow = 8, byrow = TRUE)
n_subj <- nrow(cq_data)
k_cond <- ncol(cq_data)
Cj <- colSums(cq_data)
Ri <- rowSums(cq_data)
N <- sum(cq_data)
Q <- (k_cond * (k_cond - 1) * sum((Cj - N/k_cond)^2)) / (k_cond * sum(Ri) - sum(Ri^2))
cq_df <- k_cond - 1
cq_p <- pchisq(Q, df = cq_df, lower.tail = FALSE)
results$cochranQ <- list(
  data = as.list(as.data.frame(t(cq_data))),  # list of rows
  dataMatrix = lapply(1:nrow(cq_data), function(i) as.numeric(cq_data[i,])),
  n = n_subj,
  k = k_cond,
  statistic = Q,
  df = cq_df,
  pValue = cq_p,
  colSums = as.numeric(Cj),
  rowSums = as.numeric(Ri)
)

# ─── 4. McNemar's Test ─────────────────────────────────────────────────
cat("4. McNemar's Test...\n")
# b = 15, c = 5 (off-diagonal of paired 2x2)
mn_b <- 15
mn_c <- 5
mn_mat <- matrix(c(20, mn_b, mn_c, 10), nrow = 2)
mn_nocorr <- mcnemar.test(mn_mat, correct = FALSE)
mn_corr <- mcnemar.test(mn_mat, correct = TRUE)
results$mcnemar <- list(
  b = mn_b,
  c = mn_c,
  noCorrection = list(
    statistic = unname(mn_nocorr$statistic),
    df = unname(mn_nocorr$parameter),
    pValue = mn_nocorr$p.value
  ),
  withCorrection = list(
    statistic = unname(mn_corr$statistic),
    df = unname(mn_corr$parameter),
    pValue = mn_corr$p.value
  ),
  oddsRatio = mn_b / mn_c
)

# ─── 5. Binomial Test ──────────────────────────────────────────────────
cat("5. Binomial Test...\n")
bt1 <- binom.test(7, 10, p = 0.5)
bt2 <- binom.test(3, 20, p = 0.5)
bt3 <- binom.test(15, 15, p = 0.5)  # all successes
bt4 <- binom.test(0, 10, p = 0.3)   # zero successes
results$binomial <- list(
  test1 = list(
    successes = 7, trials = 10, p0 = 0.5,
    pValue = bt1$p.value,
    ci_lower = bt1$conf.int[1],
    ci_upper = bt1$conf.int[2],
    estimate = unname(bt1$estimate)
  ),
  test2 = list(
    successes = 3, trials = 20, p0 = 0.5,
    pValue = bt2$p.value,
    ci_lower = bt2$conf.int[1],
    ci_upper = bt2$conf.int[2],
    estimate = unname(bt2$estimate)
  ),
  test3 = list(
    successes = 15, trials = 15, p0 = 0.5,
    pValue = bt3$p.value,
    ci_lower = bt3$conf.int[1],
    ci_upper = bt3$conf.int[2],
    estimate = unname(bt3$estimate)
  ),
  test4 = list(
    successes = 0, trials = 10, p0 = 0.3,
    pValue = bt4$p.value,
    ci_lower = bt4$conf.int[1],
    ci_upper = bt4$conf.int[2],
    estimate = unname(bt4$estimate)
  )
)

# ─── 6. Proportions Z-Test ─────────────────────────────────────────────
cat("6. Proportions Z-Test...\n")
# 1-sample
pt1 <- prop.test(30, 100, p = 0.5, correct = FALSE)
# 2-sample
pt2 <- prop.test(c(30, 50), c(100, 100), correct = FALSE)
# 2-sample with Yates
pt3 <- prop.test(c(30, 50), c(100, 100), correct = TRUE)
# Note: prop.test returns chi-sq = z^2; z = sqrt(chi-sq) with sign
z1 <- sqrt(unname(pt1$statistic)) * sign(30/100 - 0.5)
z2 <- sqrt(unname(pt2$statistic)) * sign(30/100 - 50/100)
results$proportionsZ <- list(
  oneSample = list(
    x = 30, n = 100, p0 = 0.5,
    chiSq = unname(pt1$statistic),
    z = z1,
    pValue = pt1$p.value,
    ci_lower = pt1$conf.int[1],
    ci_upper = pt1$conf.int[2]
  ),
  twoSample = list(
    x1 = 30, n1 = 100, x2 = 50, n2 = 100,
    chiSq = unname(pt2$statistic),
    z = z2,
    pValue = pt2$p.value,
    ci_lower = pt2$conf.int[1],
    ci_upper = pt2$conf.int[2]
  ),
  twoSampleYates = list(
    x1 = 30, n1 = 100, x2 = 50, n2 = 100,
    chiSq = unname(pt3$statistic),
    pValue = pt3$p.value
  )
)

# ─── 7. Point-Biserial Correlation ─────────────────────────────────────
cat("7. Point-Biserial Correlation...\n")
pb_bin <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
pb_cont <- c(3.1, 2.5, 4.0, 3.8, 2.9, 6.2, 5.5, 7.1, 6.8, 5.9)
pb <- cor.test(pb_bin, pb_cont)
results$pointBiserial <- list(
  binary = pb_bin,
  continuous = pb_cont,
  r = unname(pb$estimate),
  statistic = unname(pb$statistic),
  df = unname(pb$parameter),
  pValue = pb$p.value,
  ci_lower = pb$conf.int[1],
  ci_upper = pb$conf.int[2]
)

# ─── 8. Cramér's V with CI (chi-square part) ──────────────────────────
cat("8. Cramér's V...\n")
cv_obs <- matrix(c(10, 20, 30, 40, 15, 25), nrow = 2)
cv_chi <- chisq.test(cv_obs, correct = FALSE)
cv_n <- sum(cv_obs)
cv_minDim <- min(nrow(cv_obs), ncol(cv_obs)) - 1
cv_v <- sqrt(unname(cv_chi$statistic) / (cv_n * cv_minDim))
results$cramersV <- list(
  observed = lapply(1:nrow(cv_obs), function(i) as.numeric(cv_obs[i,])),
  chiSq = unname(cv_chi$statistic),
  df = unname(cv_chi$parameter),
  pValue = cv_chi$p.value,
  cramersV = cv_v,
  n = cv_n
)

# ─── 9. Two-Way ANOVA ─────────────────────────────────────────────────
cat("9. Two-Way ANOVA...\n")
set.seed(42)
tw_n <- 60
tw_a <- rep(c("Low","Med","High"), each = 20)
tw_b <- rep(c("Male","Female"), times = 30)
# Generate y with main effects + interaction
tw_y <- 10 +
  ifelse(tw_a == "Med", 2, ifelse(tw_a == "High", 5, 0)) +
  ifelse(tw_b == "Female", 1, 0) +
  ifelse(tw_a == "High" & tw_b == "Female", 3, 0) +
  rnorm(tw_n, 0, 2)
tw_df <- data.frame(y = tw_y, A = factor(tw_a), B = factor(tw_b))

# Type II (default car::Anova)
tw_fit <- lm(y ~ A * B, data = tw_df)
if (requireNamespace("car", quietly = TRUE)) {
  tw_type2 <- car::Anova(tw_fit, type = 2)
  tw_type3 <- car::Anova(tw_fit, type = 3)
  results$twoWayANOVA <- list(
    data = list(y = tw_y, factorA = tw_a, factorB = tw_b),
    typeII = list(
      ss_A = tw_type2["A", "Sum Sq"],
      ss_B = tw_type2["B", "Sum Sq"],
      ss_AB = tw_type2["A:B", "Sum Sq"],
      ss_Residuals = tw_type2["Residuals", "Sum Sq"],
      df_A = tw_type2["A", "Df"],
      df_B = tw_type2["B", "Df"],
      df_AB = tw_type2["A:B", "Df"],
      df_Residuals = tw_type2["Residuals", "Df"],
      F_A = tw_type2["A", "F value"],
      F_B = tw_type2["B", "F value"],
      F_AB = tw_type2["A:B", "F value"],
      p_A = tw_type2["A", "Pr(>F)"],
      p_B = tw_type2["B", "Pr(>F)"],
      p_AB = tw_type2["A:B", "Pr(>F)"]
    ),
    typeIII = list(
      ss_A = tw_type3["A", "Sum Sq"],
      ss_B = tw_type3["B", "Sum Sq"],
      ss_AB = tw_type3["A:B", "Sum Sq"],
      ss_Residuals = tw_type3["Residuals", "Sum Sq"],
      F_A = tw_type3["A", "F value"],
      F_B = tw_type3["B", "F value"],
      F_AB = tw_type3["A:B", "F value"],
      p_A = tw_type3["A", "Pr(>F)"],
      p_B = tw_type3["B", "Pr(>F)"],
      p_AB = tw_type3["A:B", "Pr(>F)"]
    )
  )
} else {
  # Fallback: just anova() (Type I)
  tw_aov <- anova(tw_fit)
  results$twoWayANOVA <- list(
    data = list(y = tw_y, factorA = tw_a, factorB = tw_b),
    typeI = list(
      ss_A = tw_aov["A", "Sum Sq"],
      ss_B = tw_aov["B", "Sum Sq"],
      ss_AB = tw_aov["A:B", "Sum Sq"],
      ss_Residuals = tw_aov["Residuals", "Sum Sq"]
    )
  )
}

# ─── 10. ANCOVA ─────────────────────────────────────────────────────────
cat("10. ANCOVA...\n")
set.seed(123)
an_n <- 30
an_group <- rep(c("Control","Treatment"), each = 15)
an_covariate <- rnorm(an_n, mean = 50, sd = 10)
an_y <- 20 + ifelse(an_group == "Treatment", 8, 0) + 0.3 * an_covariate + rnorm(an_n, 0, 3)
an_df <- data.frame(y = an_y, group = factor(an_group), cov = an_covariate)

an_fit <- lm(y ~ group + cov, data = an_df)
if (requireNamespace("car", quietly = TRUE)) {
  an_type3 <- car::Anova(an_fit, type = 3)
  # Adjusted means at grand mean of covariate
  grand_cov <- mean(an_covariate)
  an_coefs <- coef(an_fit)
  adj_control <- an_coefs["(Intercept)"] + an_coefs["cov"] * grand_cov
  adj_treatment <- an_coefs["(Intercept)"] + an_coefs["groupTreatment"] + an_coefs["cov"] * grand_cov
  results$ancova <- list(
    data = list(y = an_y, group = an_group, covariate = an_covariate),
    ss_group = an_type3["group", "Sum Sq"],
    ss_cov = an_type3["cov", "Sum Sq"],
    ss_Residuals = an_type3["Residuals", "Sum Sq"],
    df_group = an_type3["group", "Df"],
    df_cov = an_type3["cov", "Df"],
    df_Residuals = an_type3["Residuals", "Df"],
    F_group = an_type3["group", "F value"],
    F_cov = an_type3["cov", "F value"],
    p_group = an_type3["group", "Pr(>F)"],
    p_cov = an_type3["cov", "Pr(>F)"],
    adjustedMeans = list(
      Control = unname(adj_control),
      Treatment = unname(adj_treatment)
    ),
    grandMeanCov = grand_cov
  )
}

# ─── 11. Quasi-Poisson ─────────────────────────────────────────────────
cat("11. Quasi-Poisson...\n")
set.seed(99)
qp_n <- 50
qp_x <- runif(qp_n, 0, 5)
qp_mu <- exp(0.5 + 0.4 * qp_x)
# Overdispersed counts (negative binomial as overdispersed Poisson)
qp_y <- rnbinom(qp_n, mu = qp_mu, size = 3)
qp_fit <- glm(qp_y ~ qp_x, family = quasipoisson)
qp_s <- summary(qp_fit)
results$quasiPoisson <- list(
  data = list(y = as.numeric(qp_y), x = qp_x),
  coefficients = list(
    intercept = list(
      estimate = unname(coef(qp_fit)[1]),
      se = unname(qp_s$coefficients[1, "Std. Error"]),
      z = unname(qp_s$coefficients[1, "t value"]),
      pValue = unname(qp_s$coefficients[1, "Pr(>|t|)"])
    ),
    x = list(
      estimate = unname(coef(qp_fit)[2]),
      se = unname(qp_s$coefficients[2, "Std. Error"]),
      z = unname(qp_s$coefficients[2, "t value"]),
      pValue = unname(qp_s$coefficients[2, "Pr(>|t|)"])
    )
  ),
  dispersion = unname(qp_s$dispersion),
  deviance = unname(qp_fit$deviance),
  nullDeviance = unname(qp_fit$null.deviance)
)

# ─── 12. Negative Binomial ─────────────────────────────────────────────
cat("12. Negative Binomial...\n")
if (requireNamespace("MASS", quietly = TRUE)) {
  set.seed(77)
  nb_n <- 80
  nb_x <- runif(nb_n, 0, 4)
  nb_mu <- exp(1.0 + 0.5 * nb_x)
  nb_y <- rnbinom(nb_n, mu = nb_mu, size = 2.5)
  nb_fit <- MASS::glm.nb(nb_y ~ nb_x)
  nb_s <- summary(nb_fit)
  results$negativeBinomial <- list(
    data = list(y = as.numeric(nb_y), x = nb_x),
    coefficients = list(
      intercept = list(
        estimate = unname(coef(nb_fit)[1]),
        se = unname(nb_s$coefficients[1, "Std. Error"]),
        z = unname(nb_s$coefficients[1, "z value"]),
        pValue = unname(nb_s$coefficients[1, "Pr(>|z|)"])
      ),
      x = list(
        estimate = unname(coef(nb_fit)[2]),
        se = unname(nb_s$coefficients[2, "Std. Error"]),
        z = unname(nb_s$coefficients[2, "z value"]),
        pValue = unname(nb_s$coefficients[2, "Pr(>|z|)"])
      )
    ),
    theta = nb_fit$theta,
    aic = AIC(nb_fit),
    logLik = as.numeric(logLik(nb_fit)),
    deviance = unname(nb_fit$deviance)
  )
} else {
  cat("  MASS not available, skipping NB\n")
}

# ─── 13. Ordinal Logistic Regression ──────────────────────────────────
cat("13. Ordinal Logistic...\n")
if (requireNamespace("MASS", quietly = TRUE)) {
  set.seed(55)
  ol_n <- 100
  ol_x1 <- rnorm(ol_n)
  ol_x2 <- rnorm(ol_n)
  # Generate ordinal outcome with 4 categories
  ol_latent <- 1.0 * ol_x1 + 0.5 * ol_x2 + rlogis(ol_n)
  ol_y <- cut(ol_latent, breaks = c(-Inf, -1, 0, 1, Inf), labels = FALSE)
  ol_df <- data.frame(y = factor(ol_y), x1 = ol_x1, x2 = ol_x2)
  ol_fit <- MASS::polr(y ~ x1 + x2, data = ol_df, method = "logistic")
  ol_s <- summary(ol_fit)
  # polr coefficients (note: polr uses P(Y<=j) = logistic(zeta_j - eta))
  # where eta = x'beta. Thresholds = zeta (Intercepts in output)
  results$ordinalLogistic <- list(
    data = list(y = as.numeric(ol_y), x1 = ol_x1, x2 = ol_x2),
    coefficients = list(
      x1 = list(
        estimate = unname(ol_fit$coefficients["x1"]),
        se = unname(ol_s$coefficients["x1", "Std. Error"])
      ),
      x2 = list(
        estimate = unname(ol_fit$coefficients["x2"]),
        se = unname(ol_s$coefficients["x2", "Std. Error"])
      )
    ),
    thresholds = list(
      t1_2 = unname(ol_fit$zeta["1|2"]),
      t2_3 = unname(ol_fit$zeta["2|3"]),
      t3_4 = unname(ol_fit$zeta["3|4"])
    ),
    logLik = as.numeric(logLik(ol_fit)),
    aic = AIC(ol_fit),
    nCategories = length(levels(ol_df$y))
  )
} else {
  cat("  MASS not available, skipping ordinal\n")
}

# ─── 14. Bootstrap CI (reference: mean of known data) ─────────────────
cat("14. Bootstrap CI (reference statistics)...\n")
# For bootstrap we just verify the point estimate and that CI contains it
boot_data <- c(2.3, 4.1, 3.7, 5.2, 1.8, 4.5, 3.9, 6.1, 2.7, 3.5)
results$bootstrap <- list(
  data = boot_data,
  trueMean = mean(boot_data),
  trueSD = sd(boot_data),
  trueMedian = median(boot_data)
)

# ─── Write JSON ─────────────────────────────────────────────────────────
cat("\nWriting JSON...\n")
json_out <- toJSON(results, auto_unbox = TRUE, digits = 10, pretty = TRUE)
writeLines(json_out, "tests/fixtures/new-methods-ref.json")
cat("Done! Written to tests/fixtures/new-methods-ref.json\n")
