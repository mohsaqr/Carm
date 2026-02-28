#!/usr/bin/env Rscript
# PCA Cross-Validation: R Reference Generator
# Generates 500 synthetic correlation datasets via Saqrlab::simulate_data("correlation")
# and runs PCA with varimax rotation to produce reference values.
#
# Computes per dataset:
#   - prcomp() with center=TRUE, scale.=TRUE
#   - Eigenvalues, variance explained, cumulative variance
#   - Loadings (rotation matrix), first 5 scores
#   - Varimax rotation on first k components
#
# Usage: Rscript validation/r-reference/pca-ref.R

suppressPackageStartupMessages({
  library(Saqrlab)
  library(jsonlite)
})

N <- 500
cat("Generating", N, "PCA reference datasets (Saqrlab correlation)...\n")
results <- vector("list", N)

for (i in seq_len(N)) {
  d <- tryCatch(simulate_data("correlation", seed = i), error = function(e) NULL)
  if (is.null(d)) {
    results[[i]] <- list(seed = i, error = "simulate_failed")
    next
  }

  tryCatch({
    # Extract numeric columns (x1, x2, x3, x4)
    num_cols <- c("x1", "x2", "x3", "x4")
    d_mat <- as.matrix(d[, num_cols, drop = FALSE])
    n <- nrow(d_mat)
    p <- ncol(d_mat)

    # ── PCA ────────────────────────────────────────────────────────────────
    pca <- prcomp(d_mat, center = TRUE, scale. = TRUE)
    pca_summary <- summary(pca)

    eigenvalues <- as.numeric(pca$sdev^2)
    # NOTE: Do NOT use pca_summary$importance — it rounds to 5 decimal places!
    # Compute full-precision values directly from eigenvalues.
    var_explained <- eigenvalues / sum(eigenvalues)
    cum_var <- cumsum(var_explained)

    # Loadings (rotation matrix): p x p, save as row-major list-of-lists
    loadings_mat <- pca$rotation  # p x p
    loadings_rows <- lapply(seq_len(p), function(r) as.numeric(loadings_mat[r, ]))

    # Scores (first 5 rows): min(5, n) x p, save as row-major list-of-lists
    n_scores <- min(5L, n)
    scores_mat <- pca$x[seq_len(n_scores), , drop = FALSE]
    scores_rows <- lapply(seq_len(n_scores), function(r) as.numeric(scores_mat[r, ]))

    # ── Varimax rotation ───────────────────────────────────────────────────
    k <- min(3L, p)
    v <- varimax(pca$rotation[, seq_len(k), drop = FALSE])
    v_loadings <- v$loadings
    class(v_loadings) <- NULL  # convert from "loadings" to plain matrix
    v_mat <- matrix(as.numeric(v_loadings), nrow = p, ncol = k)
    varimax_rows <- lapply(seq_len(p), function(r) as.numeric(v_mat[r, ]))

    # ── Raw data ───────────────────────────────────────────────────────────
    x1 <- as.numeric(d_mat[, "x1"])
    x2 <- as.numeric(d_mat[, "x2"])
    x3 <- as.numeric(d_mat[, "x3"])
    x4 <- as.numeric(d_mat[, "x4"])

    results[[i]] <- list(
      seed = i,
      n = n,
      p = p,
      data = list(x1 = x1, x2 = x2, x3 = x3, x4 = x4),
      eigenvalues = eigenvalues,
      varianceExplained = var_explained,
      cumulativeVariance = cum_var,
      loadings = loadings_rows,
      scores = scores_rows,
      varimax = list(
        k = k,
        loadings = varimax_rows
      )
    )

    if (i %% 100 == 0) cat("  Completed", i, "/", N, "\n")

  }, error = function(e) {
    cat(sprintf("[%3d] ERROR: %s\n", i, e$message))
    results[[i]] <<- list(seed = i, error = e$message)
  })
}

# Summary
n_ok <- sum(vapply(results, function(r) is.null(r$error), logical(1)))
n_fail <- N - n_ok
cat(sprintf("\nDone: %d passed, %d failed out of %d\n", n_ok, n_fail, N))

outPath <- file.path("validation", "data", "pca-ref.json")
dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(results, digits = 12, auto_unbox = TRUE, na = "null"), outPath)
cat(sprintf("Saved to %s\n", outPath))
