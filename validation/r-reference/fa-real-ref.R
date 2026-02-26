#!/usr/bin/env Rscript
# Cross-validate FA on user's real dataset
suppressPackageStartupMessages({ library(psych); library(jsonlite) })

X <- read.csv("~/Downloads/rraw_dataaw_data.csv")
X <- as.matrix(X)

cat("Dataset:", nrow(X), "x", ncol(X), "\n")

# Test multiple factor counts
results <- list()
for (k in c(3, 4, 5, 6)) {
  cat(sprintf("\n=== k=%d factors ===\n", k))

  tryCatch({
    fa_res <- suppressWarnings(fa(X, nfactors = k, fm = "ml", rotate = "promax", warnings = FALSE))

    # Extract loadings
    L <- unname(matrix(as.numeric(fa_res$loadings), nrow = ncol(X), ncol = k))

    # KMO
    kmo_res <- KMO(X)

    # Bartlett
    n <- nrow(X)
    R <- cor(X, use = "pairwise.complete.obs")
    p <- ncol(X)
    chi_sq <- -((n - 1) - (2 * p + 5) / 6) * log(det(R))
    df <- p * (p - 1) / 2

    cat(sprintf("KMO overall: %.6f\n", kmo_res$MSA))
    cat(sprintf("Bartlett chi²: %.4f, df: %d\n", chi_sq, df))
    cat(sprintf("Uniquenesses: %s\n", paste(sprintf("%.6f", fa_res$uniquenesses), collapse = ", ")))

    # Factor correlations
    if (!is.null(fa_res$Phi)) {
      cat("Factor correlations (Phi):\n")
      print(round(fa_res$Phi, 6))
    }

    results[[as.character(k)]] <- list(
      k = k,
      loadings = L,
      uniquenesses = as.numeric(fa_res$uniquenesses),
      variableNames = colnames(X),
      kmo = kmo_res$MSA,
      bartlett_chi2 = chi_sq,
      bartlett_df = df,
      phi = if (!is.null(fa_res$Phi)) unname(as.matrix(fa_res$Phi)) else NULL
    )

    cat("Loadings (first 6 rows):\n")
    for (i in 1:min(6, nrow(L))) {
      cat(sprintf("  %s:", colnames(X)[i]))
      for (j in 1:k) cat(sprintf(" %9.6f", L[i, j]))
      cat("\n")
    }

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}

# Save reference data
writeLines(toJSON(results, digits = 10, auto_unbox = TRUE), "validation/data/fa-real-crossval-ref.json")
cat("\nSaved reference to validation/data/fa-real-crossval-ref.json\n")
