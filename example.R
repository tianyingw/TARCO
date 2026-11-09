# Run from this directory: Rscript --vanilla example.R
source("R/project_admm.R")
source("R/tarco.R")
source("R/grouping.R")

set.seed(2026)
n <- 100L
p <- 8L
taxa <- paste0("Taxon", seq_len(p))

# A balanced tree: leaves, pairs, and two four-leaf subtrees; omit the root.
leaf_sets <- c(as.list(seq_len(p)), list(1:2, 3:4, 5:6, 7:8, 1:4, 5:8))
A <- vapply(leaf_sets, function(s) as.numeric(seq_len(p) %in% s), numeric(p))
rownames(A) <- taxa

# Multiplicative measurement error becomes additive after the ALR transform.
latent_log <- matrix(rnorm(n * p), n, p)
tau <- 0.2
observed_log <- latent_log + matrix(rnorm(n * p, sd = tau), n, p)
composition <- exp(observed_log)
composition <- composition / rowSums(composition)
X <- sweep(log(composition), 1, log(composition[, p]), "-")
colnames(X) <- taxa
beta_true <- c(rep(0.8, 4), rep(-0.8, 4))
y <- as.numeric(2 + latent_log %*% beta_true + rnorm(n, sd = 0.5))

# Known covariance of the ALR measurement error; the reference column is zero.
Sigma <- matrix(0, p, p)
Sigma[1:(p - 1), 1:(p - 1)] <- tau^2 * (diag(p - 1) + 1)
folds <- sample(rep(1:5, length.out = n))
lambdas <- exp(seq(-4, 1, length.out = 20))
fit <- fit_tarco(X, y, A, Sigma, lambdas, folds)

# Group the fitted beta coefficients by their range within each subtree.
partition <- TreePartition(fit$beta, A, threshold = 0.01)
coefficients <- data.frame(taxon = taxa, true_beta = beta_true,
                           beta = fit$beta, group = partition$labels)
groups <- do.call(rbind, lapply(seq_along(partition$groups), function(k) {
  idx <- partition$groups[[k]]
  data.frame(group = k, taxa = paste(taxa[idx], collapse = ", "),
             mean_beta = mean(fit$beta[idx]),
             min_beta = min(fit$beta[idx]), max_beta = max(fit$beta[idx]))
}))

cat("\nSelected lambda:", fit$lambda, "\nSolver status:", fit$status, "\n")
cat("Recovered intercept:", fit$intercept, "\n\nCoefficients:\n")
print(coefficients, row.names = FALSE, digits = 4)
cat("\nGroups (diameter threshold = 0.01):\n")
print(groups, row.names = FALSE, digits = 4)
dir.create("output", showWarnings = FALSE)
write.csv(coefficients, "output/coefficients.csv", row.names = FALSE)
write.csv(groups, "output/groups.csv", row.names = FALSE)
