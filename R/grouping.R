# Maximal subtrees whose coefficient range is at most the threshold.
TreePartition <- function(beta, A, threshold = 0, tol = 1e-10) {
  beta <- as.numeric(beta)
  stopifnot(length(beta) == nrow(A), all(is.finite(beta)), threshold >= 0)
  # The root belongs to the grouping rule, but not to the fitted dictionary.
  sets <- c(lapply(seq_len(ncol(A)), function(k) which(A[, k] != 0)),
            list(seq_along(beta)))
  keys <- vapply(sets, paste, collapse = ",", FUN.VALUE = character(1))
  sets <- sets[!duplicated(keys)]
  eligible <- vapply(sets, function(g) diff(range(beta[g])) <= threshold + tol,
                     logical(1))
  candidates <- which(eligible)
  candidates <- candidates[order(lengths(sets)[candidates], decreasing = TRUE)]
  labels <- integer(length(beta))
  kept <- list()
  for (k in candidates) {
    g <- sets[[k]]
    if (all(labels[g] == 0L)) {
      kept[[length(kept) + 1L]] <- g
      labels[g] <- length(kept)
    }
  }
  stopifnot(all(labels > 0L))
  list(labels = labels, groups = kept)
}
