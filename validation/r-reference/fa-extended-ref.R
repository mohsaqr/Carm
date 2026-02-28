#!/usr/bin/env Rscript
# FA Extended Cross-Validation: R Reference Generator
# Generates 200 synthetic factor analysis datasets via Saqrlab::simulate_data("factor_analysis")
# and runs multiple extraction/rotation combos, CFA, and diagnostics.
#
# Computes per dataset:
#   - PAF + oblimin, PAF + quartimin
#   - ML + oblimin (with fit indices), ML + varimax
#   - CFA via lavaan (fit measures + standardized loadings)
#   - Diagnostics: KMO, Bartlett's test
#
# Usage: Rscript validation/r-reference/fa-extended-ref.R

suppressPackageStartupMessages({
  library(Saqrlab)
  library(jsonlite)
  library(psych)
  library(GPArotation)
  library(lavaan)
})

N <- 200
cat("Generating", N, "FA extended reference datasets (Saqrlab factor_analysis)...\n")
results <- vector("list", N)

for (i in seq_len(N)) {
  d <- tryCatch(simulate_data("factor_analysis", seed = i), error = function(e) NULL)
  if (is.null(d)) {
    results[[i]] <- list(seed = i, error = "simulate_failed")
    next
  }

  tryCatch({
    nf <- attr(d, "n_factors")
    if (is.null(nf) || is.na(nf)) {
      results[[i]] <- list(seed = i, error = "no_n_factors_attr")
      next
    }

    # Extract numeric columns only
    num_idx <- vapply(d, is.numeric, logical(1))
    d_mat <- as.matrix(d[, num_idx, drop = FALSE])
    n <- nrow(d_mat)
    p <- ncol(d_mat)
    var_names <- colnames(d_mat)

    # Ensure nf is feasible
    if (nf > p || nf < 1) {
      results[[i]] <- list(seed = i, error = sprintf("nf=%d out of range for p=%d", nf, p))
      next
    }

    # Helper: convert loadings to row-major list-of-lists
    loadings_to_rows <- function(fa_result) {
      L <- matrix(as.numeric(fa_result$loadings), nrow = p, ncol = nf)
      lapply(seq_len(p), function(r) as.numeric(L[r, ]))
    }

    # Helper: extract Phi (factor correlation matrix)
    phi_to_rows <- function(fa_result) {
      if (is.null(fa_result$Phi)) return(NULL)
      Phi <- unname(as.matrix(fa_result$Phi))
      lapply(seq_len(nrow(Phi)), function(r) as.numeric(Phi[r, ]))
    }

    # ── PAF + oblimin ──────────────────────────────────────────────────
    paf_oblimin <- tryCatch({
      fa_res <- suppressWarnings(fa(d_mat, nfactors = nf, fm = "pa", rotate = "oblimin", warnings = FALSE))
      list(
        loadings = loadings_to_rows(fa_res),
        communalities = as.numeric(fa_res$communality),
        phi = phi_to_rows(fa_res)
      )
    }, error = function(e) list(error = e$message))

    # ── PAF + quartimin ────────────────────────────────────────────────
    paf_quartimin <- tryCatch({
      fa_res <- suppressWarnings(fa(d_mat, nfactors = nf, fm = "pa", rotate = "quartimin", warnings = FALSE))
      list(
        loadings = loadings_to_rows(fa_res),
        communalities = as.numeric(fa_res$communality)
      )
    }, error = function(e) list(error = e$message))

    # ── ML + oblimin (with fit indices) ────────────────────────────────
    ml_oblimin <- tryCatch({
      fa_res <- suppressWarnings(fa(d_mat, nfactors = nf, fm = "ml", rotate = "oblimin", warnings = FALSE))

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

      list(
        loadings = loadings_to_rows(fa_res),
        fit = list(
          chiSq  = chi_m,
          df     = df_m,
          pValue = as.numeric(fa_res$PVAL),
          rmsea  = as.numeric(fa_res$RMSEA[1]),
          cfi    = cfi,
          tli    = as.numeric(fa_res$TLI),
          srmr   = srmr
        )
      )
    }, error = function(e) list(error = e$message))

    # ── ML + varimax ──────────────────────────────────────────────────
    ml_varimax <- tryCatch({
      fa_res <- suppressWarnings(fa(d_mat, nfactors = nf, fm = "ml", rotate = "varimax", warnings = FALSE))
      list(
        loadings = loadings_to_rows(fa_res)
      )
    }, error = function(e) list(error = e$message))

    # ── CFA via lavaan ─────────────────────────────────────────────────
    cfa_result <- tryCatch({
      # Build model from the true factor structure attribute
      true_loadings <- attr(d, "loadings")
      if (is.null(true_loadings)) {
        # Fallback: use PAF oblimin result to assign items to factors
        if (!is.null(paf_oblimin$loadings)) {
          # Each item assigned to factor with highest absolute loading
          L_mat <- do.call(rbind, paf_oblimin$loadings)
          assignments <- apply(abs(L_mat), 1, which.max)
        } else {
          stop("No loadings available for CFA model construction")
        }
      } else {
        # true_loadings is a matrix: items x factors
        # Assign each item to factor with highest loading
        if (is.matrix(true_loadings)) {
          assignments <- apply(abs(true_loadings), 1, which.max)
        } else {
          # It might be a vector if single factor
          assignments <- rep(1L, p)
        }
      }

      # Build lavaan model string
      factors <- sort(unique(assignments))
      model_parts <- vapply(factors, function(f) {
        items <- var_names[assignments == f]
        sprintf("F%d =~ %s", f, paste(items, collapse = " + "))
      }, character(1))
      model_string <- paste(model_parts, collapse = "\n")

      fit <- cfa(model_string, data = as.data.frame(d_mat))

      # Fit measures
      fm <- fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli", "srmr"))

      # Standardized loadings
      std_sol <- standardizedSolution(fit)
      loading_rows <- std_sol[std_sol$op == "=~", , drop = FALSE]

      std_loadings <- lapply(seq_len(nrow(loading_rows)), function(r) {
        list(
          factor = loading_rows$lhs[r],
          item = loading_rows$rhs[r],
          est = loading_rows$est.std[r],
          se = loading_rows$se[r],
          pvalue = loading_rows$pvalue[r]
        )
      })

      list(
        model = model_string,
        fitMeasures = list(
          chisq  = as.numeric(fm["chisq"]),
          df     = as.numeric(fm["df"]),
          pvalue = as.numeric(fm["pvalue"]),
          rmsea  = as.numeric(fm["rmsea"]),
          cfi    = as.numeric(fm["cfi"]),
          tli    = as.numeric(fm["tli"]),
          srmr   = as.numeric(fm["srmr"])
        ),
        standardizedLoadings = std_loadings
      )
    }, error = function(e) list(error = e$message))

    # ── Diagnostics: KMO & Bartlett ────────────────────────────────────
    diagnostics <- tryCatch({
      R <- cor(d_mat)
      kmo_res <- KMO(R)
      bart_res <- cortest.bartlett(R, n = n)
      list(
        kmo = as.numeric(kmo_res$MSA),
        bartlett = list(
          chiSq  = as.numeric(bart_res$chisq),
          pValue = as.numeric(bart_res$p.value)
        )
      )
    }, error = function(e) list(error = e$message))

    # ── Raw data (columns as arrays) ──────────────────────────────────
    data_cols <- lapply(var_names, function(v) as.numeric(d_mat[, v]))
    names(data_cols) <- var_names

    results[[i]] <- list(
      seed = i,
      n = n,
      p = p,
      nf = nf,
      variableNames = var_names,
      data = data_cols,
      paf_oblimin = paf_oblimin,
      paf_quartimin = paf_quartimin,
      ml_oblimin = ml_oblimin,
      ml_varimax = ml_varimax,
      cfa = cfa_result,
      diagnostics = diagnostics
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

outPath <- file.path("validation", "data", "fa-extended-ref.json")
dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(results, digits = 12, auto_unbox = TRUE, na = "null"), outPath)
cat(sprintf("Saved to %s\n", outPath))
