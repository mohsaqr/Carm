#!/usr/bin/env Rscript
# Round 3: ANOVA & Post-Hoc Suite — R Reference Generator
# Generates reference values for 500 datasets using Saqrlab::simulate_data("anova")
# Compares: ANOVA, effect sizes, Kruskal-Wallis, Tukey HSD, Games-Howell,
#           Dunn test (3 correction methods), Friedman test, LMM
#
# Usage: Rscript validation/r-reference/anova-ref.R

suppressPackageStartupMessages({
  library(Saqrlab)
  library(jsonlite)
  library(rstatix)
  library(dunn.test)
  library(lme4)
})

N <- 500
cat("Generating", N, "ANOVA reference datasets...\n")
results <- vector("list", N)

for (i in seq_len(N)) {
  d <- tryCatch(simulate_data("anova", seed = i), error = function(e) NULL)
  if (is.null(d)) { results[[i]] <- list(seed = i, error = "simulate_failed"); next }

  # Ensure group is a factor
  d$group <- as.factor(d$group)
  group_levels <- sort(levels(d$group))
  groups <- split(d$score, d$group)
  # Reorder groups alphabetically for consistent output
  groups <- groups[group_levels]

  k <- length(groups)
  n_total <- nrow(d)
  group_ns <- vapply(groups, length, integer(1))

  # ── One-Way ANOVA ────────────────────────────────────────────────────
  aov_res <- tryCatch({
    aov_fit <- aov(score ~ group, data = d)
    aov_sum <- summary(aov_fit)[[1]]
    dfB <- aov_sum[["Df"]][1]
    dfW <- aov_sum[["Df"]][2]
    ssB <- aov_sum[["Sum Sq"]][1]
    ssW <- aov_sum[["Sum Sq"]][2]
    msB <- aov_sum[["Mean Sq"]][1]
    msW <- aov_sum[["Mean Sq"]][2]
    fStat <- aov_sum[["F value"]][1]
    pVal <- aov_sum[["Pr(>F)"]][1]
    ssT <- ssB + ssW

    list(
      fStatistic = fStat,
      pValue = pVal,
      dfBetween = dfB,
      dfWithin = dfW,
      ssBetween = ssB,
      ssWithin = ssW,
      msBetween = msB,
      msWithin = msW,
      ssTotal = ssT,
      aov_fit = aov_fit,
      msW_val = msW,
      dfW_val = dfW,
      ssB_val = ssB,
      ssT_val = ssT,
      dfB_val = dfB
    )
  }, error = function(e) list(error = e$message))

  if (!is.null(aov_res$error)) {
    results[[i]] <- list(seed = i, error = paste("anova_failed:", aov_res$error))
    next
  }

  # ── Effect sizes ──────────────────────────────────────────────────────
  etaSq <- aov_res$ssB_val / aov_res$ssT_val
  omegaSq <- max(0, (aov_res$ssB_val - aov_res$dfB_val * aov_res$msW_val) /
                      (aov_res$ssT_val + aov_res$msW_val))

  # ── Kruskal-Wallis ───────────────────────────────────────────────────
  kw_res <- tryCatch({
    kw <- kruskal.test(score ~ group, data = d)
    H <- kw$statistic[[1]]
    etaSqKW <- (H - k + 1) / (n_total - k)
    list(
      H = H,
      pValue = kw$p.value,
      df = kw$parameter[[1]],
      etaSquaredKW = etaSqKW
    )
  }, error = function(e) list(error = e$message))

  # ── Tukey HSD ────────────────────────────────────────────────────────
  tukey_res <- tryCatch({
    thsd <- TukeyHSD(aov_res$aov_fit)$group
    pair_names <- rownames(thsd)
    # Sort pairs alphabetically for consistent ordering
    ord <- order(pair_names)
    pair_names <- pair_names[ord]
    lapply(seq_along(pair_names), function(j) {
      idx <- ord[j]
      list(
        pair = pair_names[j],
        diff = thsd[idx, "diff"],
        pAdj = thsd[idx, "p adj"],
        ciLower = thsd[idx, "lwr"],
        ciUpper = thsd[idx, "upr"]
      )
    })
  }, error = function(e) list(error = e$message))

  # ── Games-Howell ──────────────────────────────────────────────────────
  gh_res <- tryCatch({
    gh <- rstatix::games_howell_test(d, score ~ group)
    gh <- gh[order(paste0(gh$group1, "-", gh$group2)), ]
    lapply(seq_len(nrow(gh)), function(j) {
      list(
        group1 = as.character(gh$group1[j]),
        group2 = as.character(gh$group2[j]),
        diff = gh$estimate[j],
        pAdj = gh$p.adj[j],
        ciLower = gh$conf.low[j],
        ciUpper = gh$conf.high[j]
      )
    })
  }, error = function(e) list(error = e$message))

  # ── Dunn test (3 correction methods) ─────────────────────────────────
  run_dunn <- function(method_name) {
    tryCatch({
      # dunn.test prints to console; suppress it
      invisible(capture.output(
        dt <- dunn.test::dunn.test(d$score, d$group, method = method_name)
      ))
      # dt$comparisons gives pair labels, dt$Z and dt$P.adjusted give results
      pair_labels <- dt$comparisons
      ord <- order(pair_labels)
      pair_labels <- pair_labels[ord]
      lapply(seq_along(pair_labels), function(j) {
        idx <- ord[j]
        list(
          pair = pair_labels[j],
          Z = dt$Z[idx],
          pAdj = dt$P.adjusted[idx]
        )
      })
    }, error = function(e) list(error = e$message))
  }

  dunn_bonf <- run_dunn("bonferroni")
  dunn_holm <- run_dunn("holm")
  dunn_bh <- run_dunn("bh")

  # ── Friedman test ────────────────────────────────────────────────────
  # Reshape to balanced matrix: take min group size from each group
  friedman_res <- tryCatch({
    min_n <- min(group_ns)
    if (min_n < 2 || k < 2) stop("insufficient data for Friedman")
    # Build matrix: rows = "subjects", columns = groups
    fmat <- do.call(cbind, lapply(groups, function(g) g[seq_len(min_n)]))
    ft <- friedman.test(fmat)
    list(
      chi2 = ft$statistic[[1]],
      pValue = ft$p.value,
      df = ft$parameter[[1]],
      minGroupSize = min_n
    )
  }, error = function(e) list(error = e$message))

  # ── Linear Mixed Model ──────────────────────────────────────────────
  lmm_res <- tryCatch({
    lmm_fit <- lme4::lmer(score ~ 1 + (1 | group), data = d, REML = TRUE)
    vc <- as.data.frame(VarCorr(lmm_fit))
    # vc has columns: grp, var1, var2, vcov, sdcor
    group_var <- vc$vcov[vc$grp == "group"]
    resid_var <- vc$vcov[vc$grp == "Residual"]
    intercept <- fixef(lmm_fit)[[1]]
    icc <- group_var / (group_var + resid_var)
    list(
      intercept = intercept,
      groupVar = group_var,
      residVar = resid_var,
      icc = icc
    )
  }, error = function(e) list(error = e$message))

  # ── Remove internal objects from aov_res before saving ───────────────
  aov_out <- list(
    fStatistic = aov_res$fStatistic,
    pValue = aov_res$pValue,
    dfBetween = aov_res$dfBetween,
    dfWithin = aov_res$dfWithin,
    ssBetween = aov_res$ssBetween,
    ssWithin = aov_res$ssWithin,
    msBetween = aov_res$msBetween,
    msWithin = aov_res$msWithin,
    ssTotal = aov_res$ssTotal
  )

  results[[i]] <- list(
    seed = i,
    nTotal = n_total,
    k = k,
    groupLevels = group_levels,
    groups = lapply(groups, as.numeric),
    anova = aov_out,
    etaSquared = etaSq,
    omegaSquared = omegaSq,
    kruskalWallis = kw_res,
    tukeyHSD = tukey_res,
    gamesHowell = gh_res,
    dunnBonferroni = dunn_bonf,
    dunnHolm = dunn_holm,
    dunnBH = dunn_bh,
    friedman = friedman_res,
    lmm = lmm_res
  )

  if (i %% 100 == 0) cat("  Completed", i, "/", N, "\n")
}

outPath <- file.path("validation", "data", "anova-ref.json")
dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(results, digits = 12, auto_unbox = TRUE, na = "null"), outPath)
cat("Wrote", outPath, "(", sum(vapply(results, function(r) is.null(r$error), logical(1))), "valid datasets)\n")
