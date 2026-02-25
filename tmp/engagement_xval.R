#!/usr/bin/env Rscript
# Cross-validation: mclust VVI G=3 with priorControl on engagement data
library(mclust)
suppressMessages(library(jsonlite))

# Read data (European format: ; separator, , decimal)
data <- read.csv2("tmp/engagement.csv", stringsAsFactors = FALSE)
cat("Dimensions:", nrow(data), "x", ncol(data), "\n")
cat("Column names:\n")
print(names(data))

# Identify the z-scored engagement variables
zcols <- grep("^ZPRE_ENG", names(data), value = TRUE)
cat("\nZ-scored engagement columns:", zcols, "\n")

# Extract and clean (convert comma decimals, remove NAs/spaces)
x <- data[, zcols]
# Some values may have commas as decimals
for (col in names(x)) {
  x[[col]] <- as.numeric(gsub(",", ".", x[[col]]))
}
cat("x dimensions before NA removal:", nrow(x), "x", ncol(x), "\n")
cat("NAs per column:", sapply(x, function(c) sum(is.na(c))), "\n")

# Remove rows with any NA
x <- x[complete.cases(x), ]
cat("x dimensions after NA removal:", nrow(x), "x", ncol(x), "\n")
cat("Summary:\n")
print(summary(x))

# ── Fit mclust: VVI, G=3, with prior ─────────────────────────────────
mod <- Mclust(x, modelNames = "VVI", G = 3, prior = priorControl())

cat("\n=== mclust Results ===\n")
cat("Model:", mod$modelName, "\n")
cat("G:", mod$G, "\n")
cat("Log-likelihood:", mod$loglik, "\n")
cat("BIC (standard, lower=better):", -mod$bic, "\n")
cat("DF:", attr(mod$bic, "df"), "\n")

cat("\nMeans:\n")
print(round(t(mod$parameters$mean), 4))

cat("\nWeights:", round(mod$parameters$pro, 4), "\n")
cat("Classification table:\n")
print(table(mod$classification))

# ── Entropy (mclust convention) ───────────────────────────────────────
probs <- mod$z
n <- mod$n
K <- mod$G

E <- 1 + sum(probs * log(pmax(probs, 1e-300))) / (n * log(K))
Ei <- 1 + rowSums(probs * log(pmax(probs, 1e-300))) / log(K)
cat("\nEntropy (normalized):", E, "\n")
cat("mean(Ei):", mean(Ei), "\n")

# AvePP
clusters <- mod$classification
avepp <- sapply(1:K, function(k) {
  assigned <- which(clusters == k)
  if (length(assigned) == 0) return(0)
  mean(apply(probs[assigned, , drop = FALSE], 1, max))
})
cat("AvePP:", round(avepp, 4), "\n")

# ICL
raw_entropy <- -sum(probs * log(pmax(probs, 1e-300)))
icl <- (-mod$bic) + 2 * raw_entropy
cat("ICL:", icl, "\n")

# ── Variance (diagonal for VVI) ──────────────────────────────────────
cat("\nVariances (VVI diagonal):\n")
for (k in 1:K) {
  cat("  Cluster", k, ":", round(diag(mod$parameters$variance$sigma[,,k]), 4), "\n")
}

# ── Export ────────────────────────────────────────────────────────────
output <- list(
  data = as.matrix(x),
  n = nrow(x),
  d = ncol(x),
  varNames = names(x),
  model = "VVI",
  k = K,
  weights = as.numeric(mod$parameters$pro),
  means = t(mod$parameters$mean),
  logLikelihood = mod$loglik,
  bic = -mod$bic,
  df = as.numeric(attr(mod$bic, "df")),
  classification = as.numeric(mod$classification),
  entropy = E,
  caseEntropy = as.numeric(Ei),
  avepp = as.numeric(avepp),
  icl = icl,
  posteriors = probs
)

outfile <- "tmp/engagement_mclust_ref.json"
write_json(output, outfile, digits = 12, auto_unbox = TRUE, pretty = TRUE)
cat("\nReference written to", outfile, "\n")
