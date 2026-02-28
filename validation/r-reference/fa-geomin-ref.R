#!/usr/bin/env Rscript
# Cross-validate geomin rotation against GPArotation
# Generates reference values for both delta = 0.01 (GPArotation default)
# and delta = 0.001 (lavaan default)
suppressPackageStartupMessages({ library(psych); library(GPArotation); library(jsonlite) })

# Use same synthetic datasets from existing cross-validation
ref <- fromJSON("validation/data/fa-crossval-data.json", simplifyVector = FALSE)

# ── delta = 0.01 (GPArotation default) ──────────────────────────────────────
results_d01 <- list()
cat("=== Geomin cross-validation (delta = 0.01) ===\n")

for (idx in seq_along(ref$datasets)) {
  ds <- ref$datasets[[idx]]
  X <- do.call(rbind, lapply(ds$data, function(row) as.numeric(unlist(row))))
  colnames(X) <- ds$variableNames
  k <- ds$params$k
  n <- ds$params$n

  tryCatch({
    fa_none <- suppressWarnings(fa(X, nfactors = k, fm = "ml", rotate = "none", warnings = FALSE))
    unrotated <- fa_none$loadings

    geo <- geominQ(unrotated, delta = 0.01)

    L_geo <- unname(matrix(as.numeric(geo$loadings), nrow = ncol(X), ncol = k))
    Phi_geo <- unname(as.matrix(geo$Phi))

    results_d01[[length(results_d01) + 1]] <- list(
      id = ds$id,
      n = n,
      p = ncol(X),
      k = k,
      loadings = L_geo,
      phi = Phi_geo,
      converged = geo$convergence
    )

    cat(sprintf("[%3d] n=%d p=%d k=%d converged=%s\n", ds$id, n, ncol(X), k, geo$convergence))

  }, error = function(e) {
    cat(sprintf("[%3d] n=%d p=%d k=%d ERROR: %s\n", ds$id, n, ncol(X), k, e$message))
  })
}

writeLines(toJSON(results_d01, digits = 10, auto_unbox = TRUE), "validation/data/fa-geomin-ref.json")
cat(sprintf("\nSaved %d geomin (delta=0.01) reference results\n", length(results_d01)))

# ── delta = 0.001 (lavaan default) ──────────────────────────────────────────
# NOTE: With delta=0.001, the geomin landscape has more local optima.
# Use randomStarts=100 to match Carm's multi-start behaviour (default 50)
# and lavaan's default (rstarts=100). Without random starts, GPArotation
# often converges to worse local optima than Carm/lavaan.
results_d001 <- list()
cat("\n=== Geomin cross-validation (delta = 0.001, randomStarts=100) ===\n")

for (idx in seq_along(ref$datasets)) {
  ds <- ref$datasets[[idx]]
  X <- do.call(rbind, lapply(ds$data, function(row) as.numeric(unlist(row))))
  colnames(X) <- ds$variableNames
  k <- ds$params$k
  n <- ds$params$n

  tryCatch({
    fa_none <- suppressWarnings(fa(X, nfactors = k, fm = "ml", rotate = "none", warnings = FALSE))
    unrotated <- fa_none$loadings

    geo <- geominQ(unrotated, delta = 0.001, randomStarts = 100)

    L_geo <- unname(matrix(as.numeric(geo$loadings), nrow = ncol(X), ncol = k))
    Phi_geo <- unname(as.matrix(geo$Phi))

    results_d001[[length(results_d001) + 1]] <- list(
      id = ds$id,
      n = n,
      p = ncol(X),
      k = k,
      loadings = L_geo,
      phi = Phi_geo,
      converged = geo$convergence
    )

    cat(sprintf("[%3d] n=%d p=%d k=%d converged=%s\n", ds$id, n, ncol(X), k, geo$convergence))

  }, error = function(e) {
    cat(sprintf("[%3d] n=%d p=%d k=%d ERROR: %s\n", ds$id, n, ncol(X), k, e$message))
  })
}

writeLines(toJSON(results_d001, digits = 10, auto_unbox = TRUE), "validation/data/fa-geomin-d001-ref.json")
cat(sprintf("\nSaved %d geomin (delta=0.001) reference results\n", length(results_d001)))

# ── Real dataset references for both deltas ──────────────────────────────────
cat("\n=== Real dataset geomin ===\n")
X_real <- as.matrix(read.csv("~/Downloads/rraw_dataaw_data.csv"))

for (delta in c(0.01, 0.001)) {
  cat(sprintf("\n--- delta = %s ---\n", delta))
  real_results <- list()
  rs <- if (delta == 0.001) 100 else 0  # random starts for small delta
  for (k in c(3, 4, 5, 6)) {
    fa_none <- suppressWarnings(fa(X_real, nfactors = k, fm = "ml", rotate = "none", warnings = FALSE))
    geo <- geominQ(fa_none$loadings, delta = delta, randomStarts = rs)
    real_results[[as.character(k)]] <- list(
      k = k,
      loadings = unname(matrix(as.numeric(geo$loadings), nrow = ncol(X_real), ncol = k)),
      phi = unname(as.matrix(geo$Phi)),
      converged = geo$convergence
    )
    cat(sprintf("  k=%d converged=%s\n", k, geo$convergence))
  }
  suffix <- if (delta == 0.01) "" else "-d001"
  fname <- sprintf("validation/data/fa-geomin-real%s-ref.json", suffix)
  writeLines(toJSON(real_results, digits = 10, auto_unbox = TRUE), fname)
  cat(sprintf("Saved to %s\n", fname))
}
