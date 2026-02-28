#!/usr/bin/env Rscript
# Cross-validate quartimax and target rotation against GPArotation
# Generates reference values for both quartimax (orthogonal GPForth) and
# target rotation (oblique GPFoblq with targetQ)
suppressPackageStartupMessages({ library(psych); library(GPArotation); library(jsonlite) })

# Use same synthetic datasets from existing cross-validation
ref <- fromJSON("validation/data/fa-crossval-data.json", simplifyVector = FALSE)

# ── Quartimax rotation (orthogonal via GPForth) ──────────────────────────────

results_quartimax <- list()
cat("=== Quartimax cross-validation ===\n")

for (idx in seq_along(ref$datasets)) {
  ds <- ref$datasets[[idx]]
  X <- do.call(rbind, lapply(ds$data, function(row) as.numeric(unlist(row))))
  colnames(X) <- ds$variableNames
  k <- ds$params$k
  n <- ds$params$n

  tryCatch({
    fa_none <- suppressWarnings(fa(X, nfactors = k, fm = "ml", rotate = "none", warnings = FALSE))
    unrotated <- fa_none$loadings

    # Use GPForth for orthogonal quartimax
    qrt <- GPForth(unrotated, method = "quartimax")

    L_qrt <- unname(matrix(as.numeric(qrt$loadings), nrow = ncol(X), ncol = k))
    # T'T should be I for orthogonal
    T_mat <- unname(as.matrix(qrt$Th))
    # f is NULL in GPForth; extract final criterion from Table[last, 2]
    f_val <- qrt$Table[nrow(qrt$Table), 2]

    results_quartimax[[length(results_quartimax) + 1]] <- list(
      id = ds$id,
      n = n,
      p = ncol(X),
      k = k,
      loadings = L_qrt,
      T_matrix = T_mat,
      criterion = f_val,
      converged = qrt$convergence
    )

    cat(sprintf("[%3d] n=%d p=%d k=%d f=%.6f converged=%s\n",
                ds$id, n, ncol(X), k, f_val, qrt$convergence))

  }, error = function(e) {
    cat(sprintf("[%3d] n=%d p=%d k=%d ERROR: %s\n", ds$id, n, ncol(X), k, e$message))
  })
}

writeLines(toJSON(results_quartimax, digits = 10, auto_unbox = TRUE),
           "validation/data/fa-quartimax-ref.json")
cat(sprintf("\nSaved %d quartimax reference results\n\n", length(results_quartimax)))

# ── Target rotation (oblique via GPFoblq with targetQ) ───────────────────────

results_target <- list()
cat("=== Target rotation cross-validation ===\n")

# For target rotation, we generate a target from a rounded version of the
# geomin solution (threshold at 0.3), then rotate toward it.
for (idx in seq_along(ref$datasets)) {
  ds <- ref$datasets[[idx]]
  X <- do.call(rbind, lapply(ds$data, function(row) as.numeric(unlist(row))))
  colnames(X) <- ds$variableNames
  k <- ds$params$k
  n <- ds$params$n

  tryCatch({
    fa_none <- suppressWarnings(fa(X, nfactors = k, fm = "ml", rotate = "none", warnings = FALSE))
    unrotated <- fa_none$loadings

    # First get a geomin solution to derive a target
    geo <- geominQ(unrotated, delta = 0.01)
    L_geo <- as.numeric(geo$loadings)
    dim(L_geo) <- c(ncol(X), k)

    # Create target: keep loadings > 0.3 in absolute value, set others to 0
    Target <- ifelse(abs(L_geo) > 0.3, L_geo, 0)

    # Full target rotation (all weights = 1, equivalent to targetQ)
    tgt <- targetQ(unrotated, Target = Target)

    L_tgt <- unname(matrix(as.numeric(tgt$loadings), nrow = ncol(X), ncol = k))
    Phi_tgt <- unname(as.matrix(tgt$Phi))
    # f may be NULL; extract from Table if needed
    f_val <- if (!is.null(tgt$f)) tgt$f else tgt$Table[nrow(tgt$Table), 2]

    # Save unrotated loadings for sign alignment in cross-validation
    L_unrot <- unname(matrix(as.numeric(unrotated), nrow = ncol(X), ncol = k))

    results_target[[length(results_target) + 1]] <- list(
      id = ds$id,
      n = n,
      p = ncol(X),
      k = k,
      target = unname(Target),
      unrotated = L_unrot,
      loadings = L_tgt,
      phi = Phi_tgt,
      criterion = f_val,
      converged = tgt$convergence
    )

    cat(sprintf("[%3d] n=%d p=%d k=%d f=%.6f converged=%s\n",
                ds$id, n, ncol(X), k, f_val, tgt$convergence))

  }, error = function(e) {
    cat(sprintf("[%3d] n=%d p=%d k=%d ERROR: %s\n", ds$id, n, ncol(X), k, e$message))
  })
}

writeLines(toJSON(results_target, digits = 10, auto_unbox = TRUE),
           "validation/data/fa-target-ref.json")
cat(sprintf("\nSaved %d target reference results\n", length(results_target)))
