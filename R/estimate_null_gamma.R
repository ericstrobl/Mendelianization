# Estimate the null covariance matrix of multivariate z statistics while
# adaptively trimming variants with the largest multivariate statistics.
estimate_null_gamma <- function(
    Z,
    threshold_grid = seq(0, 0.20, length.out = 20),
    tolerance = 1e-4,
    max_iterations = 15,
    stop_after = 2
) {
  # Number of variants and outcomes, respectively.
  n <- nrow(Z)
  m <- ncol(Z)
  
  # Symmetrize a matrix and floor very small eigenvalues to ensure that it is
  # positive definite and therefore suitable for Cholesky decomposition.
  make_pd <- function(G) {
    e <- eigen((G + t(G)) / 2, symmetric = TRUE)
    e$values <- pmax(e$values, max(e$values) * 1e-8)
    G <- e$vectors %*% diag(e$values) %*% t(e$vectors)
    dimnames(G) <- list(colnames(Z), colnames(Z))
    G
  }
  
  # Robustly estimate each outcome's marginal null variance from the median
  # squared z statistic, whose null median is the median of chi-square(1).
  scale2 <- apply(Z, 2,function(z) median(z^2) / qchisq(0.5, 1))
  D <- diag(sqrt(scale2))
  # Combine robust marginal scales with the genome-wide correlation matrix to
  # obtain the positive-definite starting estimate.
  Gamma0 <- D %*% cov2cor(crossprod(Z) / n) %*% D
  Gamma0 <- make_pd(Gamma0)
  
  # Iteratively estimate Gamma for one candidate trimming proportion.
  fit_gamma <- function(threshold) {
    Gamma <- Gamma0
    # Retain observations inside the central chi-square region. The correction
    # reverses the covariance shrinkage induced by elliptical truncation.
    cutoff <- qchisq(1 - threshold, m)
    correction <- pchisq(cutoff, m + 2) / pchisq(cutoff, m)
    
    for (iteration in seq_len(max_iterations)) {
      # Compute each variant's squared Mahalanobis distance under the current fit.
      Omega <- chol2inv(chol(Gamma))
      Q <- rowSums((Z %*% Omega) * Z)
      keep <- Q <= cutoff
      # Re-estimate Gamma from the retained central observations and apply the
      # analytical truncation correction.
      Gamma_new <- crossprod(Z[keep, , drop = FALSE]) /
        sum(keep) / correction
      Gamma_new <- make_pd(Gamma_new)
      # Stop when the relative Frobenius change in Gamma is sufficiently small.
      change <- sqrt(sum((Gamma_new - Gamma)^2) / sum(Gamma^2))
      Gamma <- Gamma_new
      if (change < tolerance) break
    }
    Gamma
  }
  
  # Store the fitted covariance and conditional KS distance for each threshold.
  fits <- list()
  KS <- numeric()
  increases <- 0
  
  # Evaluate trimming proportions in ascending order.
  for (i in seq_along(threshold_grid)) {
    threshold <- threshold_grid[i]
    message(sprintf("Testing threshold %d/%d: %.4f",i, length(threshold_grid), threshold))
    
    # Recompute multivariate p values using the covariance fitted at this threshold.
    Gamma <- fit_gamma(threshold)
    Omega <- chol2inv(chol(Gamma))
    Q <- rowSums((Z %*% Omega) * Z)
    p <- pchisq(Q, m, lower.tail = FALSE)
    # Exclude p <= 0.05 because they are enriched for genuine associations, then
    # rescale the remaining conditional Uniform(0.05, 1) values to Uniform(0, 1).
    u <- sort((p[p > 0.05] - 0.05) / 0.95)
    k <- length(u)
    
    # Compute the two-sided one-sample KS distance from Uniform(0, 1).
    KS[i] <- if (k < 2) Inf else max(
      abs(seq_len(k) / k - u),
      abs((seq_len(k) - 1) / k - u)
    )
    fits[[i]] <- Gamma
    message(sprintf("Conditional KS: %.5f", KS[i]))
    
    # Stop searching after the KS distance worsens for the requested number of
    # consecutive thresholds.
    increases <- if (i > 1 && KS[i] > KS[i - 1]) increases + 1 else 0
    if (increases >= stop_after) {
      message("KS increased twice consecutively; stopping.")
      break
    }
  }
  
  # Select the evaluated threshold with the smallest conditional KS distance.
  selected <- which.min(KS)
  message(sprintf(
    "Selected threshold %.4f with KS %.5f",
    threshold_grid[selected], KS[selected]
  ))
  
  # Return the selected covariance estimate and its trimming proportion.
  list(
    Gamma = fits[[selected]],
    lambda = threshold_grid[selected]
  )
}