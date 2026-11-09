# Corrected quadratic loss with training-sample centering.
suppressPackageStartupMessages(library(CVXR))
cvxr_new <- packageVersion("CVXR") >= "1.8.0"
solve_cvx <- function(problem, ...) {
  if (!cvxr_new) return(solve(problem, ...))
  objective <- psolve(problem, ...)
  list(status = status(problem), value = objective,
       getValue = function(x) as.matrix(value(x)))
}

project_maxnorm <- function(S, eps = 1e-6) {
  S <- (S + t(S)) / 2
  raw_eig <- eigen(S, symmetric = TRUE)
  if (min(raw_eig$values) >= 0)
    return(list(Q = S, status = "already_psd", adjustment = 0, distance = 0))
  # Preserve exact zero rows (including the embedded ALR reference).
  active <- which(rowSums(abs(S)) > 0)
  if (length(active) < nrow(S)) {
    small <- project_maxnorm(S[active, active, drop = FALSE], eps)
    Qfull <- matrix(0, nrow(S), ncol(S))
    Qfull[active, active] <- small$Q
    small$Q <- Qfull
    if(!is.null(small$certificate)) {
      dual_full <- matrix(0,nrow(S),ncol(S))
      dual_full[active,active] <- small$certificate$dual
      small$certificate$dual <- dual_full
    }
    return(small)
  }
  p <- nrow(S)
  Q <- if (cvxr_new) Variable(c(p, p), PSD = TRUE) else Variable(p, p, PSD = TRUE)
  # Scale the epigraph problem for numerical stability.
  magnitude <- max(1, max(abs(S)))
  target <- S / magnitude
  radius <- Variable(nonneg = TRUE)
  problem <- Problem(Minimize(radius), list(Q-target <= radius, target-Q <= radius))
  attempts <- list()
  record_attempt <- function(label) {
    list(label=label, status=result$status, stats=if(cvxr_new) solver_stats(problem) else NULL)
  }
  result <- solve_cvx(problem, solver = "SCS", abstol = eps, reltol = eps, feastol = eps, num_iter = 50000)
  attempts[[length(attempts)+1L]] <- record_attempt("original_50000")
  if (identical(result$status, "optimal_inaccurate")) {
    cat("Retrying projection with a larger iteration budget\n"); flush.console()
    result <- solve_cvx(problem, solver = "SCS", abstol = eps, reltol = eps,
                        feastol = eps, num_iter = 200000, warm_start = TRUE,
                        acceleration_lookback = 0)
    attempts[[length(attempts)+1L]] <- record_attempt("original_200000_no_acceleration")
  }
  # Use ADMM if SCS does not converge.
  if (!identical(result$status, "optimal")) {
    cat("Retrying the same projection with a primal-dual certified ADMM solver\n")
    flush.console()
    rescued <- project_admm(S,eps=eps)
    attempts[[length(attempts)+1L]] <- list(label="ADMM_primal_dual_certificate",
      status=rescued$status,certificate=rescued$certificate[c("upper","lower","gap","iterations","elapsed")])
    rescued$attempts <- attempts
    return(rescued)
  }
  if (!identical(result$status, "optimal"))
    stop("Projection failed strict status check: ", result$status)
  q <- magnitude * result$getValue(Q); q <- (q + t(q)) / 2
  eg <- eigen(q, symmetric = TRUE)
  clean <- eg$vectors %*% (pmax(eg$values, 0) * t(eg$vectors))
  adjustment <- max(abs(clean - q))
  if (adjustment > 20 * eps * max(1, max(abs(S))))
    stop("Projection eigenvalue clipping exceeds tolerance.")
  list(Q = clean, status = result$status, adjustment = adjustment,
       distance = max(abs(clean - S)), attempts=attempts)
}

moments <- function(X, y, Sigma, center = NULL, intercept = TRUE,
                    training = TRUE, corrected = TRUE) {
  n <- nrow(X)
  stopifnot(n >= 2L, length(y) == n, all(dim(Sigma) == c(ncol(X), ncol(X))),
            all(is.finite(X)), all(is.finite(y)), all(is.finite(Sigma)))
  if (!training && intercept && is.null(center))
    stop("Validation must use means from the training fold.")
  if (is.null(center)) center <- if (intercept)
    list(x = colMeans(X), y = mean(y)) else list(x = rep(0, ncol(X)), y = 0)
  V <- sweep(X, 2, center$x); z <- as.numeric(y) - center$y
  # Sample centering removes 1/n of the expected training noise covariance.
  factor <- if (training && intercept) (n - 1) / n else 1
  S <- crossprod(V) / n - if (corrected) factor * Sigma else 0
  pr <- if (corrected) project_maxnorm(S) else
    list(Q = S, status = "uncorrected", adjustment = 0, distance = 0)
  c(pr, list(q = as.numeric(crossprod(V, z)) / n, center = center, covariance_factor = factor, n = n,
             centered = intercept, training = training))
}

solve_path <- function(moment, F, weights, lambdas) {
  eg <- eigen(moment$Q, symmetric = TRUE)
  L <- sqrt(pmax(eg$values, 0)) * t(eg$vectors)
  theta <- Variable(ncol(F)); lam <- Parameter(nonneg = TRUE)
  b <- F %*% theta
  problem <- Problem(Minimize(sum_squares(L %*% b) / 2 - sum(moment$q * b) +
                               lam * sum(abs(weights * theta))),
                     list(sum(b) == 0))
  lapply(lambdas, function(lambda) {
    value(lam) <- lambda
    res <- tryCatch(solve_cvx(problem, solver = "ECOS", abstol = 1e-8,
                          reltol = 1e-8, feastol = 1e-8), error = function(e) NULL)
    # Exclude unbounded candidates.
    if (!is.null(res) && identical(res$status, "unbounded"))
      return(list(beta = NULL, status = "unbounded"))
    if (is.null(res) || !identical(res$status, "optimal"))
      res <- tryCatch(solve_cvx(problem, solver = "SCS", abstol = 1e-6, reltol = 1e-6,
                            feastol = 1e-6, num_iter = 50000), error = function(e) NULL)
    if (is.null(res) || !identical(res$status, "optimal"))
      return(list(beta = NULL, status = if (is.null(res)) "error" else res$status))
    beta <- as.numeric(F %*% res$getValue(theta))
    if (any(!is.finite(beta)) || abs(sum(beta)) > 1e-5 * max(1, sum(abs(beta))))
      return(list(beta = NULL, status = "constraint_failure"))
    list(beta = beta, status = res$status, objective = res$value)
  })
}

prepare_cv <- function(X, y, Sigma, folds, intercept = TRUE, corrected = TRUE) {
  folds_unique <- sort(unique(folds))
  fold_moments <- lapply(folds_unique, function(k) {
    cat("Preparing fold", k, "of", length(folds_unique), "\n"); flush.console()
    train <- folds != k; test <- !train
    tr <- moments(X[train, , drop = FALSE], y[train], Sigma,
                  intercept = intercept, corrected = corrected)
    va <- moments(X[test, , drop = FALSE], y[test], Sigma, center = tr$center,
                  intercept = intercept, training = FALSE, corrected = corrected)
    list(train = tr, validation = va)
  })
  list(folds = fold_moments, full = moments(X, y, Sigma,
             intercept = intercept, corrected = corrected))
}

fit_prepared <- function(pre, F, weights, lambdas) {
  losses <- matrix(Inf, nrow = length(pre$folds), ncol = length(lambdas))
  status <- matrix("", nrow(losses), ncol(losses))
  for (k in seq_along(pre$folds)) {
    cat("Fitting fold", k, "with", length(lambdas), "penalties\n"); flush.console()
    fit <- solve_path(pre$folds[[k]]$train, F, weights, lambdas)
    va <- pre$folds[[k]]$validation
    for (l in seq_along(fit)) {
      status[k, l] <- fit[[l]]$status
      b <- fit[[l]]$beta
      if (!is.null(b)) losses[k, l] <- as.numeric(crossprod(b, va$Q %*% b)) / 2 - sum(va$q * b)
    }
  }
  score <- colMeans(losses)
  if (!any(is.finite(score))) stop("All candidate penalties failed on at least one fold.")
  j <- which.min(score)
  selected <- solve_path(pre$full, F, weights, lambdas[j])[[1]]
  if (is.null(selected$beta)) stop("Selected full-sample fit failed: ", selected$status)
  selected$intercept <- pre$full$center$y - sum(pre$full$center$x * selected$beta)
  selected$lambda <- lambdas[j]
  selected$boundary <- j %in% c(1L, length(lambdas))
  selected$cv <- data.frame(lambda = lambdas, loss = score)
  selected$fold_status <- status
  selected
}


# X contains observed ALR coordinates, with a zero reference column.
fit_tarco <- function(X, y, A, Sigma, lambdas, folds, alpha = 0) {
  X <- as.matrix(X); A <- as.matrix(A); Sigma <- as.matrix(Sigma)
  stopifnot(nrow(A) == ncol(X), all(A %in% c(0, 1)),
            all(colSums(A) > 0), length(folds) == nrow(X),
            !anyNA(folds), length(unique(folds)) >= 2L,
            length(lambdas) > 0L, all(is.finite(lambdas)), all(lambdas > 0),
            length(alpha) == 1L, is.finite(alpha))
  pre <- prepare_cv(X, y, Sigma, folds)
  fit_prepared(pre, A, colSums(A)^alpha, lambdas)
}
