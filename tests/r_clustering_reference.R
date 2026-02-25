#!/usr/bin/env Rscript
# Generate reference values for clustering cross-validation.
# Output: JSON file with exact numerical results from R.

library(mclust)
library(poLCA)
library(MASS)
suppressMessages(library(jsonlite))

cat("=== R Clustering Cross-Validation Reference ===\n")

# ─── 1. FIXED GMM DATA ───────────────────────────────────────────────────
# 2D data: 3 clusters with known means, unit covariance
set.seed(42)
n_per <- 30
g1 <- mvrnorm(n_per, mu = c(0, 0), Sigma = diag(2))
g2 <- mvrnorm(n_per, mu = c(6, 0), Sigma = diag(2))
g3 <- mvrnorm(n_per, mu = c(3, 5), Sigma = diag(2))
gmm_data <- rbind(g1, g2, g3)

# Fit GMM with mclust — VVV model, K=3
fit_vvv <- Mclust(gmm_data, G = 3, modelNames = "VVV", verbose = FALSE)
# Fit GMM — EII model, K=3
fit_eii <- Mclust(gmm_data, G = 3, modelNames = "EII", verbose = FALSE)
# Fit GMM — VVI model, K=3
fit_vvi <- Mclust(gmm_data, G = 3, modelNames = "VVI", verbose = FALSE)

extract_gmm <- function(fit, name) {
  probs <- fit$z                        # posterior conditional probs (N × K)
  n <- fit$n
  K <- fit$G

  # Normalized entropy: 1 + sum(probs * log(probs)) / (n * log(K))
  entropy <- 1 + sum(probs * log(pmax(probs, 1e-300))) / (n * log(K))

  # Case-specific entropy contributions: Ei = 1 + rowSums(probs * log(probs)) / log(K)
  Ei <- 1 + rowSums(probs * log(pmax(probs, 1e-300))) / log(K)

  # Average posterior probabilities per component
  clusters <- fit$classification
  avepp <- vapply(seq_len(K), function(k) {
    assigned <- which(clusters == k)
    if (length(assigned) == 0) return(0)
    mean(apply(probs[assigned, , drop = FALSE], 1, max))
  }, numeric(1))

  # ICL = BIC + 2 * raw_entropy (where raw = -sum(probs * log(probs)))
  raw_entropy <- -sum(probs * log(pmax(probs, 1e-300)))
  icl <- (-fit$bic) + 2 * raw_entropy  # standard BIC + 2E

  list(
    model = name,
    k = K,
    weights = as.numeric(fit$parameters$pro),
    means = t(fit$parameters$mean),  # K × D
    logLikelihood = fit$loglik,
    bic = -fit$bic,  # mclust returns negative BIC (higher = better), we want standard BIC (lower = better)
    df = attr(fit$bic, "df"),
    classification = as.numeric(fit$classification),
    entropy = entropy,
    caseEntropy = as.numeric(Ei),
    avepp = as.numeric(avepp),
    icl = icl
  )
}

gmm_vvv <- extract_gmm(fit_vvv, "VVV")
gmm_eii <- extract_gmm(fit_eii, "EII")
gmm_vvi <- extract_gmm(fit_vvi, "VVI")

cat("GMM VVV LL:", gmm_vvv$logLikelihood, "\n")
cat("GMM VVV BIC:", gmm_vvv$bic, "\n")
cat("GMM VVV means:\n"); print(round(t(fit_vvv$parameters$mean), 4))
cat("GMM VVV entropy:", gmm_vvv$entropy, "\n")
cat("GMM VVV AvePP:", gmm_vvv$avepp, "\n")
cat("GMM VVV ICL:", gmm_vvv$icl, "\n")
cat("GMM EII LL:", gmm_eii$logLikelihood, "\n")
cat("GMM EII entropy:", gmm_eii$entropy, "\n")
cat("GMM VVI LL:", gmm_vvi$logLikelihood, "\n")
cat("GMM VVI entropy:", gmm_vvi$entropy, "\n")

# ─── 2. FIXED LCA DATA ───────────────────────────────────────────────────
# Binary data: 2 classes, 5 items
set.seed(42)
n_lca <- 100
# Class 1: items 1-3 high, 4-5 low
c1 <- matrix(0, nrow = n_lca/2, ncol = 5)
for (i in 1:(n_lca/2)) {
  c1[i,] <- c(rbinom(1,1,0.9), rbinom(1,1,0.85), rbinom(1,1,0.8),
              rbinom(1,1,0.15), rbinom(1,1,0.1))
}
# Class 2: items 1-3 low, 4-5 high
c2 <- matrix(0, nrow = n_lca/2, ncol = 5)
for (i in 1:(n_lca/2)) {
  c2[i,] <- c(rbinom(1,1,0.1), rbinom(1,1,0.15), rbinom(1,1,0.2),
              rbinom(1,1,0.85), rbinom(1,1,0.9))
}
lca_data_raw <- rbind(c1, c2)

# poLCA requires data as data.frame with values 1 and 2 (not 0 and 1!)
lca_df <- as.data.frame(lca_data_raw + 1)  # Convert 0/1 to 1/2
colnames(lca_df) <- paste0("V", 1:5)

# Fit LCA with poLCA, K=2, single start
set.seed(42)
f <- cbind(V1, V2, V3, V4, V5) ~ 1
fit_lca <- poLCA(f, lca_df, nclass = 2, nrep = 1, verbose = FALSE)

# Extract — poLCA returns P(x=2|class) in probs, we want P(x=1) = P(original=1)
# probs is a list of K matrices, one per item. Each matrix is K_classes × 2 (prob of value 1, prob of value 2)
rho_r <- matrix(0, nrow = 2, ncol = 5)
for (d in 1:5) {
  for (k in 1:2) {
    rho_r[k, d] <- fit_lca$probs[[d]][k, 2]  # Column 2 = P(value=2) = P(original=1)
  }
}

lca_ref <- list(
  k = 2,
  priorWeights = as.numeric(fit_lca$P),
  rho = rho_r,
  logLikelihood = fit_lca$llik,
  bic = fit_lca$bic,
  aic = fit_lca$aic,
  df = 5 * 2 + 1,  # k*m + (k-1)
  classification = as.numeric(fit_lca$predclass)
)

cat("\nLCA LL:", lca_ref$logLikelihood, "\n")
cat("LCA BIC:", lca_ref$bic, "\n")
cat("LCA AIC:", lca_ref$aic, "\n")
cat("LCA weights:", lca_ref$priorWeights, "\n")
cat("LCA rho:\n"); print(round(rho_r, 4))

# ─── 3. FIXED KMEANS DATA ────────────────────────────────────────────────
# Use same gmm_data, fit with K=3 Lloyd's algorithm
set.seed(42)
# Use fixed starting centers to ensure determinism
# Pick first point from each true cluster
init_centers <- gmm_data[c(1, n_per+1, 2*n_per+1), ]
km <- kmeans(gmm_data, centers = init_centers, algorithm = "Lloyd", iter.max = 300)

kmeans_ref <- list(
  k = 3,
  centroids = km$centers,
  labels = as.numeric(km$cluster),
  inertia = km$tot.withinss,
  iterations = km$iter
)

cat("\nKMeans inertia:", kmeans_ref$inertia, "\n")
cat("KMeans iterations:", kmeans_ref$iterations, "\n")
cat("KMeans centers:\n"); print(round(km$centers, 4))

# ─── 4. WRITE OUTPUT ─────────────────────────────────────────────────────

output <- list(
  gmm_data = gmm_data,
  gmm_vvv = gmm_vvv,
  gmm_eii = gmm_eii,
  gmm_vvi = gmm_vvi,
  lca_data = lca_data_raw,
  lca = lca_ref,
  kmeans = kmeans_ref
)

outfile <- "tests/r_clustering_reference.json"
write_json(output, outfile, digits = 10, auto_unbox = TRUE, pretty = TRUE)
cat("\nReference written to", outfile, "\n")
