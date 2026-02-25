#!/usr/bin/env Rscript
# Generate reference values for DBSCAN cross-validation.
# Output: JSON file with exact numerical results from R.
#
# R DBSCAN convention:
#   - Noise points are labeled 0
#   - Clusters are 1-indexed (1, 2, 3, ...)
# TypeScript convention may differ (e.g. -1 for noise, 0-indexed clusters).
# The test harness must handle this mapping.

library(dbscan)
library(cluster)
suppressMessages(library(jsonlite))

cat("=== R DBSCAN Cross-Validation Reference ===\n")

# ─── 1. SYNTHETIC DATA ─────────────────────────────────────────────────
# 3 well-separated clusters, 10 points each, 2D
set.seed(42)
data <- rbind(
  matrix(rnorm(20, mean = 0, sd = 0.5), ncol = 2),   # 10 points, cluster 1
  matrix(rnorm(20, mean = 5, sd = 0.5), ncol = 2),   # 10 points, cluster 2
  matrix(rnorm(20, mean = 10, sd = 0.5), ncol = 2)   # 10 points, cluster 3
)

cat("Data dimensions:", nrow(data), "x", ncol(data), "\n")
cat("Data head:\n")
print(head(data, 6))

# ─── 2. DBSCAN ─────────────────────────────────────────────────────────
# eps=1.5 should easily capture all 3 clusters with no noise
db <- dbscan::dbscan(data, eps = 1.5, minPts = 3)

cat("\nDBSCAN result:\n")
print(db)
cat("Cluster labels (0=noise):", db$cluster, "\n")
cat("Number of clusters (excluding noise):", max(db$cluster), "\n")
cat("Number of noise points:", sum(db$cluster == 0), "\n")

# Core/border point identification
# A core point has >= minPts neighbors within eps (including itself)
dist_matrix <- as.matrix(dist(data))
n <- nrow(data)
is_core <- vapply(seq_len(n), function(i) {
  sum(dist_matrix[i, ] <= 1.5) >= 3
}, logical(1))

cat("Core points:", sum(is_core), "\n")
cat("Border points:", sum(!is_core & db$cluster > 0), "\n")

# ─── 3. k-NN DISTANCES ─────────────────────────────────────────────────
# k=3 to match minPts=3
# kNNdist(data, k) returns a vector of k-th NN distances (length n)
# kNNdist(data, k, all=TRUE) returns an n x k matrix of all 1..k NN distances
knn_k_only <- dbscan::kNNdist(data, k = 3)             # vector: k-th NN distance per point
knn_all <- dbscan::kNNdist(data, k = 3, all = TRUE)    # n x k matrix: all NN distances

# For the "knee plot", sort the k-th NN distances descending
k_distances <- sort(knn_k_only, decreasing = TRUE)

cat("\nk-NN distances (k=3, sorted descending):\n")
print(round(k_distances, 6))

# ─── 4. SILHOUETTE ─────────────────────────────────────────────────────
# Silhouette only for non-noise points (cluster > 0)
# R's silhouette() will error or give meaningless results if noise (0) is included
non_noise_idx <- which(db$cluster > 0)
cat("\nNon-noise points:", length(non_noise_idx), "of", n, "\n")

if (length(non_noise_idx) == n) {
  # No noise points — compute silhouette on all data
  sil <- cluster::silhouette(db$cluster, dist(data))
} else if (length(non_noise_idx) > 0) {
  # Subset to non-noise points only
  sil <- cluster::silhouette(db$cluster[non_noise_idx], dist(data[non_noise_idx, ]))
} else {
  sil <- NULL
}

if (!is.null(sil)) {
  sil_summary <- summary(sil)
  sil_mean <- mean(sil[, "sil_width"])
  sil_per_cluster <- sil_summary$clus.avg.widths

  cat("Silhouette mean:", sil_mean, "\n")
  cat("Silhouette per cluster:", sil_per_cluster, "\n")
  cat("Silhouette widths (per point):\n")
  print(round(sil[, "sil_width"], 6))
} else {
  sil_mean <- NA
  sil_per_cluster <- NA
  cat("No non-noise points — silhouette not computed.\n")
}

# ─── 5. ADDITIONAL DBSCAN WITH DIFFERENT eps ─────────────────────────
# eps=0.5 should produce more noise points (tighter clusters)
db_tight <- dbscan::dbscan(data, eps = 0.5, minPts = 3)

cat("\nDBSCAN (eps=0.5, minPts=3):\n")
print(db_tight)
cat("Cluster labels:", db_tight$cluster, "\n")
cat("Noise count:", sum(db_tight$cluster == 0), "\n")

# Silhouette for tight eps (only non-noise)
non_noise_tight <- which(db_tight$cluster > 0)
if (length(non_noise_tight) > 1 && length(unique(db_tight$cluster[non_noise_tight])) > 1) {
  sil_tight <- cluster::silhouette(db_tight$cluster[non_noise_tight], dist(data[non_noise_tight, ]))
  sil_tight_mean <- mean(sil_tight[, "sil_width"])
  cat("Silhouette mean (tight):", sil_tight_mean, "\n")
} else {
  sil_tight <- NULL
  sil_tight_mean <- NA
  cat("Cannot compute silhouette for tight eps (insufficient non-noise clusters).\n")
}

# ─── 6. DBSCAN WITH eps WHERE SOME NOISE EXISTS ─────────────────────
# eps=1.0 might produce a few noise points between clusters
db_mid <- dbscan::dbscan(data, eps = 1.0, minPts = 3)

cat("\nDBSCAN (eps=1.0, minPts=3):\n")
print(db_mid)
cat("Cluster labels:", db_mid$cluster, "\n")
cat("Noise count:", sum(db_mid$cluster == 0), "\n")

non_noise_mid <- which(db_mid$cluster > 0)
if (length(non_noise_mid) > 1 && length(unique(db_mid$cluster[non_noise_mid])) > 1) {
  sil_mid <- cluster::silhouette(db_mid$cluster[non_noise_mid], dist(data[non_noise_mid, ]))
  sil_mid_mean <- mean(sil_mid[, "sil_width"])
  sil_mid_per_point <- as.numeric(sil_mid[, "sil_width"])
  cat("Silhouette mean (mid):", sil_mid_mean, "\n")
} else {
  sil_mid <- NULL
  sil_mid_mean <- NA
  sil_mid_per_point <- NA
  cat("Cannot compute silhouette for mid eps.\n")
}

# ─── 7. WRITE OUTPUT ───────────────────────────────────────────────────

output <- list(
  # Raw data for reconstruction
  data = data,

  # Convention note
  convention = "R DBSCAN: 0=noise, 1..K=cluster IDs (1-indexed)",

  # Main run: eps=1.5, minPts=3
  dbscan_eps1_5 = list(
    eps = 1.5,
    minPts = 3,
    clusters = as.numeric(db$cluster),
    nClusters = max(db$cluster),
    nNoise = sum(db$cluster == 0),
    isCore = as.logical(is_core),
    nCore = sum(is_core),
    nBorder = sum(!is_core & db$cluster > 0),
    silhouetteMean = sil_mean,
    silhouettePerCluster = as.numeric(sil_per_cluster),
    silhouetteWidths = if (!is.null(sil)) as.numeric(sil[, "sil_width"]) else NA,
    nonNoiseIndices = as.numeric(non_noise_idx)  # 1-indexed
  ),

  # k-NN distances
  knnDistances = list(
    k = 3,
    sortedDescending = as.numeric(k_distances),
    kthNeighbor = as.numeric(knn_k_only),  # vector: k-th NN distance per point (unsorted)
    allNeighbors = knn_all                  # n x k matrix: distances to 1st, 2nd, 3rd NN
  ),

  # Tight run: eps=0.5, minPts=3
  dbscan_eps0_5 = list(
    eps = 0.5,
    minPts = 3,
    clusters = as.numeric(db_tight$cluster),
    nClusters = max(db_tight$cluster),
    nNoise = sum(db_tight$cluster == 0),
    silhouetteMean = sil_tight_mean,
    nonNoiseIndices = as.numeric(non_noise_tight)  # 1-indexed
  ),

  # Mid run: eps=1.0, minPts=3
  dbscan_eps1_0 = list(
    eps = 1.0,
    minPts = 3,
    clusters = as.numeric(db_mid$cluster),
    nClusters = max(db_mid$cluster),
    nNoise = sum(db_mid$cluster == 0),
    silhouetteMean = sil_mid_mean,
    silhouetteWidths = if (!is.null(sil_mid)) sil_mid_per_point else NA,
    nonNoiseIndices = as.numeric(non_noise_mid)  # 1-indexed
  )
)

outfile <- "tests/r_dbscan_reference.json"
write_json(output, outfile, digits = 10, auto_unbox = TRUE, pretty = TRUE)
cat("\nReference written to", outfile, "\n")
