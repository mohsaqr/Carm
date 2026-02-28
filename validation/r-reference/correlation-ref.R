#!/usr/bin/env Rscript
# Round 2: Correlation Suite — R Reference Generator
# Generates reference values for 1000 datasets using Saqrlab::simulate_data("correlation")
# Compares: Pearson, Spearman, Kendall, partial correlation, correlation matrices
#
# Usage: Rscript validation/r-reference/correlation-ref.R

suppressPackageStartupMessages({
  library(Saqrlab)
  library(jsonlite)
  library(ppcor)
})

N <- 1000
cat("Generating", N, "correlation reference datasets...\n")
results <- vector("list", N)

for (i in seq_len(N)) {
  d <- tryCatch(simulate_data("correlation", seed = i), error = function(e) NULL)
  if (is.null(d)) { results[[i]] <- list(seed = i, error = "simulate_failed"); next }

  # Extract the first 4 numeric columns (x1, x2, x3, x4)
  num_cols <- names(d)[vapply(d, is.numeric, logical(1))]
  if (length(num_cols) < 4) {
    results[[i]] <- list(seed = i, error = "insufficient_numeric_columns")
    next
  }
  x1 <- as.numeric(d[[num_cols[1]]])
  x2 <- as.numeric(d[[num_cols[2]]])
  x3 <- as.numeric(d[[num_cols[3]]])
  x4 <- as.numeric(d[[num_cols[4]]])
  n <- length(x1)

  # Build matrix of 4 columns for matrix operations
  mat <- cbind(x1, x2, x3, x4)

  # ── Pairwise indices (all 6 pairs of 4 columns, 0-indexed) ────────────
  pairs <- list(
    c(0L, 1L), c(0L, 2L), c(0L, 3L),
    c(1L, 2L), c(1L, 3L), c(2L, 3L)
  )
  cols <- list(x1, x2, x3, x4)

  # ── Pearson correlations (all 6 pairs) ────────────────────────────────
  pearson <- lapply(pairs, function(p) {
    xi <- cols[[p[1] + 1L]]
    xj <- cols[[p[2] + 1L]]
    ct <- tryCatch(cor.test(xi, xj, method = "pearson"), error = function(e) NULL)
    if (is.null(ct)) return(list(i = p[1], j = p[2], error = "cor_test_failed"))
    list(
      i = p[1], j = p[2],
      r = ct$estimate[[1]],
      t = ct$statistic[[1]],
      pValue = ct$p.value,
      ciLower = ct$conf.int[1],
      ciUpper = ct$conf.int[2]
    )
  })

  # ── Spearman correlations (all 6 pairs) ───────────────────────────────
  spearman <- lapply(pairs, function(p) {
    xi <- cols[[p[1] + 1L]]
    xj <- cols[[p[2] + 1L]]
    ct <- tryCatch(
      suppressWarnings(cor.test(xi, xj, method = "spearman")),
      error = function(e) NULL
    )
    if (is.null(ct)) return(list(i = p[1], j = p[2], error = "cor_test_failed"))
    list(
      i = p[1], j = p[2],
      rho = ct$estimate[[1]],
      pValue = ct$p.value
    )
  })

  # ── Kendall correlations (all 6 pairs) ────────────────────────────────
  kendall <- lapply(pairs, function(p) {
    xi <- cols[[p[1] + 1L]]
    xj <- cols[[p[2] + 1L]]
    ct <- tryCatch(cor.test(xi, xj, method = "kendall"), error = function(e) NULL)
    if (is.null(ct)) return(list(i = p[1], j = p[2], error = "cor_test_failed"))
    list(
      i = p[1], j = p[2],
      tau = ct$estimate[[1]],
      pValue = ct$p.value
    )
  })

  # ── Partial correlations ──────────────────────────────────────────────
  # Test 1: x1 vs x2, controlling for x3 and x4
  partial1 <- tryCatch({
    pc <- ppcor::pcor.test(x1, x2, cbind(x3, x4))
    list(
      statistic = pc$statistic,
      pValue = pc$p.value,
      estimate = pc$estimate
    )
  }, error = function(e) list(error = e$message))

  # Test 2: x1 vs x3, controlling for x2 (single control)
  partial2 <- tryCatch({
    pc <- ppcor::pcor.test(x1, x3, cbind(x2))
    list(
      statistic = pc$statistic,
      pValue = pc$p.value,
      estimate = pc$estimate
    )
  }, error = function(e) list(error = e$message))

  # ── Correlation matrices (4x4, saved as row-major lists) ──────────────
  corrMatPearson <- cor(mat, method = "pearson")
  corrMatSpearman <- cor(mat, method = "spearman")
  corrMatKendall <- cor(mat, method = "kendall")

  # Convert to row-major list of lists
  mat_to_list <- function(m) {
    lapply(seq_len(nrow(m)), function(r) as.numeric(m[r, ]))
  }

  results[[i]] <- list(
    seed = i,
    n = n,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    pearson = pearson,
    spearman = spearman,
    kendall = kendall,
    partial1 = partial1,
    partial2 = partial2,
    corrMatrixPearson = mat_to_list(corrMatPearson),
    corrMatrixSpearman = mat_to_list(corrMatSpearman),
    corrMatrixKendall = mat_to_list(corrMatKendall)
  )

  if (i %% 100 == 0) cat("  Completed", i, "/", N, "\n")
}

outPath <- file.path("validation", "data", "correlation-ref.json")
dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(results, digits = 12, auto_unbox = TRUE, na = "null"), outPath)
cat("Wrote", outPath, "(", sum(vapply(results, function(r) is.null(r$error), logical(1))), "valid datasets)\n")
