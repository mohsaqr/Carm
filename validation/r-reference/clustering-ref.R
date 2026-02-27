#!/usr/bin/env Rscript
# Clustering Cross-Validation: R Reference Generator
# Generates 200 synthetic cluster datasets via Saqrlab::simulate_data("clusters")
# and runs K-Means, GMM (mclust), HAC, and silhouette analysis.
#
# Computes per dataset:
#   - K-Means (Lloyd algorithm, nstart=500, reproducible via set.seed)
#   - GMM with mclust: VVV, VVI, EEI covariance models
#   - HAC: ward.D2, complete, single, average linkages
#   - Silhouette scores from K-Means partition
#
# Usage: Rscript validation/r-reference/clustering-ref.R

suppressPackageStartupMessages({
  library(Saqrlab)
  library(jsonlite)
  library(mclust)
  library(cluster)
})

N <- 200
cat("Generating", N, "clustering reference datasets (Saqrlab clusters)...\n")
results <- vector("list", N)

for (i in seq_len(N)) {
  d <- tryCatch(simulate_data("clusters", seed = i), error = function(e) NULL)
  if (is.null(d)) {
    results[[i]] <- list(seed = i, error = "simulate_failed")
    next
  }

  tryCatch({
    features <- as.matrix(d[, c("x1", "x2")])
    true_cluster <- as.integer(d$true_cluster)
    k <- max(true_cluster)
    n <- nrow(features)

    # ── K-Means (Lloyd, nstart=1, reproducible) ─────────────────────────
    set.seed(i)
    km <- kmeans(features, centers = k, nstart = 500, algorithm = "Lloyd", iter.max = 100)

    # Sort centroids by first column for consistent ordering
    km_ord <- order(km$centers[, 1])
    km_centroids <- unname(km$centers[km_ord, , drop = FALSE])
    km_centroids_rows <- lapply(seq_len(k), function(r) as.numeric(km_centroids[r, ]))

    # Remap labels to match sorted centroid order
    km_label_map <- integer(k)
    km_label_map[km_ord] <- seq_len(k)
    km_labels <- km_label_map[km$cluster]

    km_withinss <- as.numeric(km$withinss[km_ord])

    km_result <- list(
      centroids = km_centroids_rows,
      withinss = km_withinss,
      totss = km$totss,
      tot_withinss = km$tot.withinss,
      labels = km_labels
    )

    # ── GMM (mclust) ───────────────────────────────────────────────────
    run_gmm <- function(modelName) {
      tryCatch({
        m <- Mclust(features, G = k, modelNames = modelName)
        if (is.null(m)) return(NULL)

        # Sort means by first dimension
        means_mat <- t(m$parameters$mean)  # k x 2
        gmm_ord <- order(means_mat[, 1])
        sorted_means <- means_mat[gmm_ord, , drop = FALSE]
        means_rows <- lapply(seq_len(k), function(r) as.numeric(sorted_means[r, ]))

        # Remap labels to match sorted order
        gmm_label_map <- integer(k)
        gmm_label_map[gmm_ord] <- seq_len(k)
        gmm_labels <- gmm_label_map[m$classification]

        list(
          means = means_rows,
          loglik = m$loglik,
          bic = m$bic,
          labels = gmm_labels
        )
      }, error = function(e) {
        list(error = e$message)
      })
    }

    gmm_vvv <- run_gmm("VVV")
    gmm_vvi <- run_gmm("VVI")
    gmm_eei <- run_gmm("EEI")

    # ── HAC ────────────────────────────────────────────────────────────
    dist_mat <- dist(features)

    run_hac <- function(method) {
      tryCatch({
        hc <- hclust(dist_mat, method = method)
        labels <- as.integer(cutree(hc, k = k))
        list(
          heights = sort(hc$height),
          labels = labels
        )
      }, error = function(e) {
        list(error = e$message)
      })
    }

    hac_ward <- run_hac("ward.D2")
    hac_complete <- run_hac("complete")
    hac_single <- run_hac("single")
    hac_average <- run_hac("average")

    # ── Silhouette (from K-Means partition) ────────────────────────────
    sil <- tryCatch({
      s <- silhouette(km$cluster, dist_mat)
      sil_vals <- as.numeric(s[, "sil_width"])
      n_sil <- min(20L, length(sil_vals))
      list(
        mean = mean(sil_vals),
        scores = sil_vals[seq_len(n_sil)]
      )
    }, error = function(e) {
      list(error = e$message)
    })

    # ── Raw data ───────────────────────────────────────────────────────
    x1 <- as.numeric(features[, "x1"])
    x2 <- as.numeric(features[, "x2"])

    results[[i]] <- list(
      seed = i,
      n = n,
      k = k,
      data = list(x1 = x1, x2 = x2),
      true_cluster = true_cluster,
      kmeans = km_result,
      gmm_vvv = gmm_vvv,
      gmm_vvi = gmm_vvi,
      gmm_eei = gmm_eei,
      hac_ward = hac_ward,
      hac_complete = hac_complete,
      hac_single = hac_single,
      hac_average = hac_average,
      silhouette = sil
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

outPath <- file.path("validation", "data", "clustering-ref.json")
dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(results, digits = 12, auto_unbox = TRUE, na = "null"), outPath)
cat(sprintf("Saved to %s\n", outPath))
