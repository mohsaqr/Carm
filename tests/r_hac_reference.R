#!/usr/bin/env Rscript
# Generate reference values for Hierarchical Agglomerative Clustering (HAC) cross-validation.
# Output: JSON file with exact numerical results from R.
#
# R hclust convention:
#   - merge matrix: negative indices = original observations (1-indexed),
#     positive indices = previously merged clusters (1-indexed).
#     e.g. merge[i,] = c(-3, -7) means observation 3 and 7 were merged at step i.
#     e.g. merge[i,] = c(-5, 2) means observation 5 was merged with cluster formed at step 2.
#   - order: 1-indexed leaf ordering for the dendrogram.
#   - height: merge distance at each step (length n-1).
#   - cutree: returns 1-indexed cluster labels.

library(cluster)
suppressMessages(library(jsonlite))

cat("=== R HAC Cross-Validation Reference ===\n")

# ─── 1. SYNTHETIC DATA ─────────────────────────────────────────────────
# Same data as DBSCAN reference: 3 well-separated clusters, 10 points each, 2D
set.seed(42)
data <- rbind(
  matrix(rnorm(20, mean = 0, sd = 0.5), ncol = 2),   # 10 points, cluster 1
  matrix(rnorm(20, mean = 5, sd = 0.5), ncol = 2),   # 10 points, cluster 2
  matrix(rnorm(20, mean = 10, sd = 0.5), ncol = 2)   # 10 points, cluster 3
)

cat("Data dimensions:", nrow(data), "x", ncol(data), "\n")
cat("Data head:\n")
print(head(data, 6))

# Distance matrix (used by all methods)
d <- dist(data)

# ─── 2. WARD.D2 ────────────────────────────────────────────────────────
# Ward.D2 uses squared Euclidean distances internally (Lance-Williams formula)
hc_ward <- hclust(d, method = "ward.D2")

cat("\n--- Ward.D2 ---\n")
cat("Merge matrix (first 5 rows):\n")
print(head(hc_ward$merge, 5))
cat("Heights (first 10):\n")
print(round(head(hc_ward$height, 10), 6))
cat("Order:\n")
print(hc_ward$order)

# Cut into k=3
cut_ward_3 <- cutree(hc_ward, k = 3)
cat("Cutree k=3:", cut_ward_3, "\n")

# Cut into k=2 and k=4 for additional validation
cut_ward_2 <- cutree(hc_ward, k = 2)
cut_ward_4 <- cutree(hc_ward, k = 4)
cat("Cutree k=2:", cut_ward_2, "\n")
cat("Cutree k=4:", cut_ward_4, "\n")

# Cophenetic correlation
coph_ward <- cor(cophenetic(hc_ward), d)
cat("Cophenetic correlation:", coph_ward, "\n")

# Silhouette for k=3
sil_ward <- cluster::silhouette(cut_ward_3, d)
sil_ward_mean <- mean(sil_ward[, "sil_width"])
sil_ward_per_cluster <- summary(sil_ward)$clus.avg.widths
cat("Silhouette mean (k=3):", sil_ward_mean, "\n")
cat("Silhouette per cluster:", sil_ward_per_cluster, "\n")

# ─── 3. SINGLE LINKAGE ─────────────────────────────────────────────────
hc_single <- hclust(d, method = "single")

cat("\n--- Single Linkage ---\n")
cat("Merge matrix (first 5 rows):\n")
print(head(hc_single$merge, 5))
cat("Heights (first 10):\n")
print(round(head(hc_single$height, 10), 6))
cat("Order:\n")
print(hc_single$order)

cut_single_3 <- cutree(hc_single, k = 3)
cat("Cutree k=3:", cut_single_3, "\n")

coph_single <- cor(cophenetic(hc_single), d)
cat("Cophenetic correlation:", coph_single, "\n")

sil_single <- cluster::silhouette(cut_single_3, d)
sil_single_mean <- mean(sil_single[, "sil_width"])
cat("Silhouette mean (k=3):", sil_single_mean, "\n")

# ─── 4. COMPLETE LINKAGE ───────────────────────────────────────────────
hc_complete <- hclust(d, method = "complete")

cat("\n--- Complete Linkage ---\n")
cat("Merge matrix (first 5 rows):\n")
print(head(hc_complete$merge, 5))
cat("Heights (first 10):\n")
print(round(head(hc_complete$height, 10), 6))
cat("Order:\n")
print(hc_complete$order)

cut_complete_3 <- cutree(hc_complete, k = 3)
cat("Cutree k=3:", cut_complete_3, "\n")

coph_complete <- cor(cophenetic(hc_complete), d)
cat("Cophenetic correlation:", coph_complete, "\n")

sil_complete <- cluster::silhouette(cut_complete_3, d)
sil_complete_mean <- mean(sil_complete[, "sil_width"])
cat("Silhouette mean (k=3):", sil_complete_mean, "\n")

# ─── 5. AVERAGE LINKAGE (UPGMA) ────────────────────────────────────────
hc_average <- hclust(d, method = "average")

cat("\n--- Average Linkage ---\n")
cat("Merge matrix (first 5 rows):\n")
print(head(hc_average$merge, 5))
cat("Heights (first 10):\n")
print(round(head(hc_average$height, 10), 6))
cat("Order:\n")
print(hc_average$order)

cut_average_3 <- cutree(hc_average, k = 3)
cat("Cutree k=3:", cut_average_3, "\n")

coph_average <- cor(cophenetic(hc_average), d)
cat("Cophenetic correlation:", coph_average, "\n")

sil_average <- cluster::silhouette(cut_average_3, d)
sil_average_mean <- mean(sil_average[, "sil_width"])
cat("Silhouette mean (k=3):", sil_average_mean, "\n")

# ─── 6. EXTRACT HELPER ─────────────────────────────────────────────────

extract_hclust <- function(hc, cut3, cut_others, sil, coph, name) {
  sil_summary <- summary(sil)
  list(
    method = name,
    # Merge matrix: n-1 rows x 2 cols
    # Negative = original obs (1-indexed), Positive = merged cluster (1-indexed)
    merge = hc$merge,
    # Heights at each merge step (length n-1)
    height = as.numeric(hc$height),
    # Leaf order for dendrogram (1-indexed)
    order = as.numeric(hc$order),
    # Labels (if any)
    labels = if (!is.null(hc$labels)) hc$labels else NULL,
    # Cut results
    cutree_k3 = as.numeric(cut3),
    cutree_others = cut_others,
    # Cophenetic correlation
    copheneticCorrelation = coph,
    # Silhouette for k=3
    silhouetteMean = mean(sil[, "sil_width"]),
    silhouettePerCluster = as.numeric(sil_summary$clus.avg.widths),
    silhouetteWidths = as.numeric(sil[, "sil_width"])
  )
}

ward_ref <- extract_hclust(
  hc_ward, cut_ward_3,
  list(k2 = as.numeric(cut_ward_2), k4 = as.numeric(cut_ward_4)),
  sil_ward, coph_ward, "ward.D2"
)

single_ref <- extract_hclust(
  hc_single, cut_single_3,
  list(),
  sil_single, coph_single, "single"
)

complete_ref <- extract_hclust(
  hc_complete, cut_complete_3,
  list(),
  sil_complete, coph_complete, "complete"
)

average_ref <- extract_hclust(
  hc_average, cut_average_3,
  list(),
  sil_average, coph_average, "average"
)

# ─── 7. DISTANCE MATRIX ────────────────────────────────────────────────
# Store the full distance matrix for cross-validation
# dist() returns lower triangle; convert to full matrix
dist_full <- as.matrix(d)

cat("\n--- Distance matrix (first 5x5) ---\n")
print(round(dist_full[1:5, 1:5], 6))

# ─── 8. WRITE OUTPUT ───────────────────────────────────────────────────

output <- list(
  # Raw data for reconstruction
  data = data,

  # Convention notes
  convention = list(
    merge = "Negative indices = original observations (1-indexed). Positive indices = merged clusters (1-indexed, referring to row in merge matrix).",
    order = "1-indexed leaf ordering for dendrogram plot.",
    height = "Merge distance at each step. Length = n-1.",
    cutree = "1-indexed cluster labels.",
    wardD2 = "Ward.D2 uses sqrt of Ward's criterion. Heights are not squared distances but sqrt of the Lance-Williams merge criterion."
  ),

  # Distance matrix (full n x n)
  distanceMatrix = dist_full,

  # Linkage methods
  ward = ward_ref,
  single = single_ref,
  complete = complete_ref,
  average = average_ref
)

outfile <- "tests/r_hac_reference.json"
write_json(output, outfile, digits = 10, auto_unbox = TRUE, pretty = TRUE)
cat("\nReference written to", outfile, "\n")
