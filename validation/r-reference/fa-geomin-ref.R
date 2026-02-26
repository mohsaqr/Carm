#!/usr/bin/env Rscript
# Cross-validate geomin rotation against GPArotation
suppressPackageStartupMessages({ library(psych); library(GPArotation); library(jsonlite) })

# Use same synthetic datasets from existing cross-validation
ref <- fromJSON("validation/data/fa-crossval-data.json", simplifyVector = FALSE)

results <- list()
cat("=== Geomin cross-validation ===\n")

for (idx in seq_along(ref$datasets)) {
  ds <- ref$datasets[[idx]]
  X <- do.call(rbind, lapply(ds$data, function(row) as.numeric(unlist(row))))
  colnames(X) <- ds$variableNames
  k <- ds$params$k
  n <- ds$params$n

  tryCatch({
    # Run ML extraction without rotation first
    fa_none <- suppressWarnings(fa(X, nfactors = k, fm = "ml", rotate = "none", warnings = FALSE))
    unrotated <- fa_none$loadings

    # Apply geomin oblique rotation via GPArotation
    geo <- geominQ(unrotated, delta = 0.01)

    L_geo <- unname(matrix(as.numeric(geo$loadings), nrow = ncol(X), ncol = k))
    Phi_geo <- unname(as.matrix(geo$Phi))

    results[[length(results) + 1]] <- list(
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

writeLines(toJSON(results, digits = 10, auto_unbox = TRUE), "validation/data/fa-geomin-ref.json")
cat(sprintf("\nSaved %d geomin reference results\n", length(results)))

# Also test on real dataset
cat("\n=== Real dataset geomin ===\n")
X_real <- as.matrix(read.csv("~/Downloads/rraw_dataaw_data.csv"))
for (k in c(3, 5)) {
  fa_none <- suppressWarnings(fa(X_real, nfactors = k, fm = "ml", rotate = "none", warnings = FALSE))
  geo <- geominQ(fa_none$loadings, delta = 0.01)
  L <- unname(matrix(as.numeric(geo$loadings), nrow = ncol(X_real), ncol = k))
  cat(sprintf("k=%d converged=%s\n", k, geo$convergence))
  cat("First 6 loadings:\n")
  for (i in 1:6) {
    cat(sprintf("  %s: %s\n", colnames(X_real)[i],
      paste(sprintf("%.6f", L[i,]), collapse = "  ")))
  }
  cat("Phi:\n")
  print(round(geo$Phi, 4))
}

# Save real dataset reference
real_results <- list()
for (k in c(3, 4, 5, 6)) {
  fa_none <- suppressWarnings(fa(X_real, nfactors = k, fm = "ml", rotate = "none", warnings = FALSE))
  geo <- geominQ(fa_none$loadings, delta = 0.01)
  real_results[[as.character(k)]] <- list(
    k = k,
    loadings = unname(matrix(as.numeric(geo$loadings), nrow = ncol(X_real), ncol = k)),
    phi = unname(as.matrix(geo$Phi)),
    converged = geo$convergence
  )
}
writeLines(toJSON(real_results, digits = 10, auto_unbox = TRUE), "validation/data/fa-geomin-real-ref.json")
cat("Saved real dataset geomin reference\n")
