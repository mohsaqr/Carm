#!/usr/bin/env Rscript
# Round 1: T-Test & Descriptive Suite — R Reference Generator
# Generates reference values for 1000 datasets using Saqrlab::simulate_data("ttest")
# Compares: descriptive stats, t-tests, nonparametric tests, effect sizes, frequency tests

suppressPackageStartupMessages({
  library(Saqrlab)
  library(jsonlite)
  library(effsize)
  library(e1071)
})

N <- 1000
cat("Generating", N, "ttest reference datasets...\n")
results <- vector("list", N)

for (i in seq_len(N)) {
  d <- tryCatch(simulate_data("ttest", seed = i), error = function(e) NULL)
  if (is.null(d)) { results[[i]] <- list(seed = i, error = "simulate_failed"); next }

  g <- split(d$score, d$group)
  x1 <- as.numeric(g[[1]])
  x2 <- as.numeric(g[[2]])
  n1 <- length(x1); n2 <- length(x2)
  nmin <- min(n1, n2)

  # ── Descriptive per group ──────────────────────────────────────────────
  make_desc <- function(x) {
    n <- length(x)
    ci <- if (n >= 2) as.numeric(t.test(x)$conf.int) else c(NA, NA)
    sw <- if (n >= 3 && n <= 5000) shapiro.test(x) else NULL
    list(
      n = n,
      mean = mean(x),
      median = median(x),
      sd = sd(x),
      variance = var(x),
      se = sd(x) / sqrt(n),
      trimmedMean = mean(x, trim = 0.1),
      skewness = if (n >= 3) e1071::skewness(x, type = 2) else NULL,
      kurtosis = if (n >= 4) e1071::kurtosis(x, type = 2) else NULL,
      ciLower = ci[1],
      ciUpper = ci[2],
      shapiroW = if (!is.null(sw)) sw$statistic[[1]] else NULL,
      shapiroP = if (!is.null(sw)) sw$p.value else NULL
    )
  }
  desc1 <- make_desc(x1)
  desc2 <- make_desc(x2)

  # ── T-tests ────────────────────────────────────────────────────────────
  welch <- t.test(x1, x2, var.equal = FALSE)
  student <- t.test(x1, x2, var.equal = TRUE)
  paired_t <- if (nmin >= 2) t.test(x1[1:nmin], x2[1:nmin], paired = TRUE) else NULL
  less_t <- t.test(x1, x2, var.equal = FALSE, alternative = "less")
  greater_t <- t.test(x1, x2, var.equal = FALSE, alternative = "greater")

  # ── Mann-Whitney U (no continuity correction to match Carm) ────────────
  mwu <- wilcox.test(x1, x2, exact = FALSE, correct = FALSE)
  mwu_less <- wilcox.test(x1, x2, alternative = "less", exact = FALSE, correct = FALSE)

  # ── Wilcoxon signed-rank (paired) ──────────────────────────────────────
  # R returns V = W+ (sum of positive ranks)
  # Carm returns min(W+, W-) — we save both for comparison
  wsr <- NULL
  wsr_min <- NULL
  if (nmin >= 2) {
    diffs <- x1[1:nmin] - x2[1:nmin]
    nonzero <- diffs[diffs != 0]
    if (length(nonzero) >= 1) {
      wsr <- tryCatch(
        wilcox.test(x1[1:nmin], x2[1:nmin], paired = TRUE),
        error = function(e) NULL
      )
      if (!is.null(wsr)) {
        V <- wsr$statistic[[1]]  # W+
        nn <- length(nonzero)
        maxW <- nn * (nn + 1) / 2
        wsr_min <- min(V, maxW - V)  # Carm's convention
      }
    }
  }

  # ── Effect sizes ───────────────────────────────────────────────────────
  cd <- effsize::cohen.d(x1, x2)
  hg <- effsize::cohen.d(x1, x2, hedges.correction = TRUE)

  # Rank-biserial from Mann-Whitney U
  U1 <- mwu$statistic[[1]]
  rb <- 1 - 2 * U1 / (n1 * n2)

  # Cohen's d CI — Hedges & Olkin normal approximation (matches Carm)
  d_val <- cd$estimate
  se_d <- sqrt((n1 + n2) / (n1 * n2) + d_val^2 / (2 * (n1 + n2)))
  d_ci <- d_val + c(-1, 1) * qnorm(0.975) * se_d

  # ── Frequency tests ────────────────────────────────────────────────────
  # Goodness-of-fit: test equal proportions across groups
  obs <- as.numeric(table(d$group))
  gof <- chisq.test(obs)

  # 2x2 contingency: group × median_split(score)
  med_score <- median(d$score)
  above <- as.integer(d$score > med_score)
  ct <- table(d$group, above)

  chi_noY <- chi_Y <- fish <- NULL
  ct_a <- ct_b <- ct_c <- ct_d <- NA
  if (nrow(ct) == 2 && ncol(ct) == 2) {
    ct_a <- ct[1, 1]; ct_b <- ct[1, 2]
    ct_c <- ct[2, 1]; ct_d <- ct[2, 2]
    chi_noY <- suppressWarnings(chisq.test(ct, correct = FALSE))
    chi_Y <- suppressWarnings(chisq.test(ct, correct = TRUE))
    fish <- tryCatch(fisher.test(ct), error = function(e) NULL)
  }

  results[[i]] <- list(
    seed = i, n1 = n1, n2 = n2,
    x1 = x1, x2 = x2,
    desc1 = desc1, desc2 = desc2,
    welch = list(
      statistic = welch$statistic[[1]], pValue = welch$p.value,
      df = welch$parameter[[1]], ciLower = welch$conf.int[1], ciUpper = welch$conf.int[2]
    ),
    student = list(
      statistic = student$statistic[[1]], pValue = student$p.value,
      df = student$parameter[[1]], ciLower = student$conf.int[1], ciUpper = student$conf.int[2]
    ),
    paired = if (!is.null(paired_t)) list(
      statistic = paired_t$statistic[[1]], pValue = paired_t$p.value,
      df = paired_t$parameter[[1]], ciLower = paired_t$conf.int[1], ciUpper = paired_t$conf.int[2]
    ) else NULL,
    lessT = list(statistic = less_t$statistic[[1]], pValue = less_t$p.value),
    greaterT = list(statistic = greater_t$statistic[[1]], pValue = greater_t$p.value),
    mannWhitney = list(statistic = mwu$statistic[[1]], pValue = mwu$p.value),
    mannWhitneyLess = list(statistic = mwu_less$statistic[[1]], pValue = mwu_less$p.value),
    wilcoxon = if (!is.null(wsr)) list(
      V = wsr$statistic[[1]],  # R convention: W+
      W = wsr_min,             # Carm convention: min(W+, W-)
      pValue = wsr$p.value
    ) else NULL,
    cohensD = list(value = cd$estimate),
    hedgesG = list(value = hg$estimate),
    rankBiserial = rb,
    cohensDCI = list(lower = d_ci[1], upper = d_ci[2]),
    goodnessOfFit = list(
      statistic = gof$statistic[[1]], pValue = gof$p.value, df = gof$parameter[[1]]
    ),
    contingency = list(a = ct_a, b = ct_b, c = ct_c, d = ct_d),
    chiSquareNoYates = if (!is.null(chi_noY)) list(
      statistic = chi_noY$statistic[[1]], pValue = chi_noY$p.value, df = chi_noY$parameter[[1]]
    ) else NULL,
    chiSquareYates = if (!is.null(chi_Y)) list(
      statistic = chi_Y$statistic[[1]], pValue = chi_Y$p.value, df = chi_Y$parameter[[1]]
    ) else NULL,
    fisherExact = if (!is.null(fish)) list(
      pValue = fish$p.value, oddsRatio = fish$estimate[[1]]
    ) else NULL
  )

  if (i %% 100 == 0) cat("  Completed", i, "/", N, "\n")
}

outPath <- file.path("validation", "data", "ttest-ref.json")
dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(results, digits = 12, auto_unbox = TRUE, na = "null"), outPath)
cat("Wrote", outPath, "(", sum(sapply(results, function(r) is.null(r$error))), "valid datasets)\n")
