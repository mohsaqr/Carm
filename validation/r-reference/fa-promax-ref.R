#!/usr/bin/env Rscript
# Factor Analysis Cross-Validation: R Reference Values
# Generates 100 synthetic factor-structured datasets and runs:
#   - psych::fa() with ML extraction and promax rotation
#   - psych::KMO() for sampling adequacy
#   - psych::cortest.bartlett() for sphericity test
#   - Eigenvalue decomposition of correlation matrix
#
# Usage: Rscript validation/r-reference/fa-promax-ref.R

suppressPackageStartupMessages({
  library(psych)
  library(jsonlite)
})

set.seed(42)

N_DATASETS <- 100

# ── Generate dataset configurations ──────────────────────────────────────

configs <- vector("list", N_DATASETS)
for (i in seq_len(N_DATASETS)) {
  n <- sample(c(100, 150, 200, 300, 500), 1)
  k <- sample(2:5, 1)
  p <- k * sample(c(3, 4, 5), 1)
  loading_strength <- round(runif(1, 0.45, 0.85), 2)
  configs[[i]] <- list(n = n, p = p, k = k, ls = loading_strength, seed = 1000L + i)
}

# ── Run FA on each dataset ───────────────────────────────────────────────

results <- vector("list", N_DATASETS)
ok_count <- 0L
fail_count <- 0L

for (i in seq_len(N_DATASETS)) {
  cfg <- configs[[i]]
  set.seed(cfg$seed)
  n <- cfg$n; p <- cfg$p; k <- cfg$k; ls <- cfg$ls

  # Build loading matrix (p x k) with simple structure + small cross-loadings
  items_per <- p %/% k
  rem <- p %% k
  L <- matrix(0, p, k)
  row_idx <- 1L
  for (f in seq_len(k)) {
    ni <- items_per + as.integer(f <= rem)
    for (j in seq_len(ni)) {
      L[row_idx, f] <- ls + runif(1, -0.1, 0.1)
      for (other in seq_len(k)) {
        if (other != f) L[row_idx, other] <- runif(1, -0.05, 0.15)
      }
      row_idx <- row_idx + 1L
    }
  }

  # Generate data: X = F %*% t(L) + E
  Fscores <- matrix(rnorm(n * k), n, k)
  uniq <- runif(p, 0.2, 0.6)
  E <- vapply(seq_len(p), function(j) rnorm(n, 0, sqrt(uniq[j])), numeric(n))
  X <- Fscores %*% t(L) + E
  colnames(X) <- paste0("V", seq_len(p))

  tryCatch({
    fa_res <- suppressWarnings(fa(X, nfactors = k, fm = "ml", rotate = "promax", warnings = FALSE))
    R <- cor(X)
    kmo_res <- KMO(R)
    bart_res <- cortest.bartlett(R, n = n)
    eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

    # CFI computation
    chi_m <- as.numeric(fa_res$STATISTIC)
    df_m  <- as.numeric(fa_res$dof)
    chi_n <- as.numeric(fa_res$null.chisq)
    df_n  <- as.numeric(fa_res$null.dof)
    cfi   <- 1 - max(chi_m - df_m, 0) / max(chi_n - df_n, 0)

    # SRMR from residual matrix
    resid_mat <- fa_res$residual
    srmr <- if (!is.null(resid_mat)) {
      lower_tri <- resid_mat[lower.tri(resid_mat)]
      sqrt(mean(lower_tri^2))
    } else { NA_real_ }

    # AIC (SEM convention: chi2 + 2*q where q = free params)
    q <- p * k + p - k * (k - 1L) / 2  # loadings + uniquenesses - rotation constraints
    aic <- chi_m + 2 * q

    loadings_mat <- unname(matrix(as.numeric(fa_res$loadings), nrow = p, ncol = k))

    results[[i]] <- list(
      id = i,
      params = list(n = n, p = p, k = k, loadingStrength = ls, seed = cfg$seed),
      data = unname(lapply(seq_len(n), function(r) as.numeric(X[r, ]))),
      variableNames = unname(colnames(X)),
      rResults = list(
        eigenvalues = as.numeric(eig),
        loadings = unname(lapply(seq_len(p), function(r) loadings_mat[r, ])),
        communalities = as.numeric(fa_res$communality),
        uniqueness = as.numeric(fa_res$uniquenesses),
        fit = list(
          chiSq     = chi_m,
          df        = df_m,
          pValue    = as.numeric(fa_res$PVAL),
          rmsea     = as.numeric(fa_res$RMSEA[1]),
          rmseaLower = as.numeric(fa_res$RMSEA[2]),
          rmseaUpper = as.numeric(fa_res$RMSEA[3]),
          cfi       = cfi,
          tli       = as.numeric(fa_res$TLI),
          srmr      = srmr,
          aic       = aic,
          bic       = as.numeric(fa_res$BIC)
        ),
        kmo = list(
          overall = as.numeric(kmo_res$MSA),
          perItem = as.numeric(kmo_res$MSAi)
        ),
        bartlett = list(
          chiSq  = as.numeric(bart_res$chisq),
          df     = as.numeric(bart_res$df),
          pValue = as.numeric(bart_res$p.value)
        ),
        nFactors = k
      )
    )
    ok_count <- ok_count + 1L
    cat(sprintf("[%3d] n=%3d p=%2d k=%d ls=%.2f — OK\n", i, n, p, k, ls))
  }, error = function(e) {
    fail_count <<- fail_count + 1L
    cat(sprintf("[%3d] n=%3d p=%2d k=%d ls=%.2f — FAIL: %s\n", i, n, p, k, ls, e$message))
    results[[i]] <<- NULL
  })
}

# Remove NULLs
results <- results[!vapply(results, is.null, logical(1))]

cat(sprintf("\nDone: %d passed, %d failed\n", ok_count, fail_count))
write_json(list(datasets = results), "validation/data/fa-crossval-data.json", auto_unbox = TRUE, digits = 12)
cat(sprintf("Saved %d datasets to validation/data/fa-crossval-data.json\n", length(results)))
