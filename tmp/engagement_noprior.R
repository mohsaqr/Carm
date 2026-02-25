#!/usr/bin/env Rscript
library(mclust)
suppressMessages(library(jsonlite))

data <- read.csv2("tmp/engagement.csv", stringsAsFactors = FALSE)
x <- data[, c("ZPRE_ENG_EMOC", "ZPRE_ENG_COGN", "ZPRE_ENG_COND")]
for (col in names(x)) x[[col]] <- as.numeric(gsub(",", ".", x[[col]]))
x <- x[complete.cases(x), ]

mod <- Mclust(x, modelNames = "VVI", G = 3)

probs <- mod$z
n <- mod$n
K <- mod$G

entropy <- 1 + sum(probs * log(pmax(probs, 1e-300))) / (n * log(K))
Ei <- 1 + rowSums(probs * log(pmax(probs, 1e-300))) / log(K)
clusters <- mod$classification

avepp <- sapply(1:K, function(k) {
  assigned <- which(clusters == k)
  if (length(assigned) == 0) return(0)
  mean(apply(probs[assigned, , drop = FALSE], 1, max))
})

raw_entropy <- -sum(probs * log(pmax(probs, 1e-300)))
icl <- (-mod$bic) + 2 * raw_entropy
df <- nMclustParams(mod$modelName, mod$d, mod$G)
aic <- 2 * df - 2 * mod$loglik

# Per-cluster descriptives
per_cluster <- lapply(1:K, function(k) {
  idx <- which(clusters == k)
  sub <- as.matrix(x[idx, ])
  list(
    n = length(idx),
    means = colMeans(sub),
    sds = apply(sub, 2, sd),
    vars = diag(mod$parameters$variance$sigma[,,k]),
    avepp = avepp[k],
    entropy_mean = mean(Ei[idx]),
    entropy_sd = sd(Ei[idx]),
    entropy_min = min(Ei[idx]),
    entropy_max = max(Ei[idx])
  )
})

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
  aic = aic,
  df = df,
  icl = icl,
  classification = as.numeric(clusters),
  entropy = entropy,
  caseEntropy = as.numeric(Ei),
  avepp = as.numeric(avepp),
  posteriors = probs,
  perCluster = per_cluster
)

write_json(output, "tmp/engagement_noprior_ref.json", digits = 12, auto_unbox = TRUE, pretty = TRUE)

cat("=== mclust VVI K=3 (NO prior) ===\n")
cat("LL:", mod$loglik, "\n")
cat("BIC:", -mod$bic, "\n")
cat("AIC:", aic, "\n")
cat("DF:", df, "\n")
cat("ICL:", icl, "\n")
cat("Entropy:", entropy, "\n")
cat("AvePP:", round(avepp, 4), "\n")
cat("Weights:", round(mod$parameters$pro, 4), "\n")
cat("Sizes:", table(clusters), "\n")
cat("Means:\n"); print(round(t(mod$parameters$mean), 4))
