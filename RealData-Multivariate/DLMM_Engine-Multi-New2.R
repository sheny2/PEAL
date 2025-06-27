library(Matrix)
library(mvtnorm)

## 1. Modified Summary Statistics Function for Multivariate Case ##
lmm.get.summary3.multivariate <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(Y))
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  py <- ncol(Y)

  ShXYZ <- list()
  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    wth <- weights[id.site == sh]
    Xh <- X[id.site == sh, ]
    Yh <- Y[id.site == sh, ]
    Zh <- Z[[h]]  # Changed from Z[h][[1]] to Z[[h]]

    # Modified summary statistics for multivariate case
    ShX  <- crossprod((Xh * wth), Xh)  # px × px
    ShXZ <- crossprod((Xh * wth), Zh)  # px × pz
    ShXY <- crossprod((Xh * wth), Yh)  # px × py
    ShZ  <- crossprod((Zh * wth), Zh)  # pz × pz
    ShZY <- crossprod((Zh * wth), Yh)  # pz × py
    ShY  <- crossprod((Yh * sqrt(wth))) # py × py
    Nh <- sum(id.site == sh)

    ShXYZ[[sh]] <- list(ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ = ShZ, ShZY = ShZY, ShY = ShY, Nh = Nh)
  }

  return(ShXYZ)
}

## 2. Function to Estimate Covariance Parameters (Run at Central Site) ##
estimate_covariances <- function(ShXYZ, par.init = NULL, reml = TRUE, verbose = TRUE) {
  K <- length(ShXYZ)
  py <- ncol(ShXYZ[[1]]$ShY)
  px <- ncol(ShXYZ[[1]]$ShX)

  # Initialize parameters: [sigma_u, sigma_v1, ..., sigma_vK, sigma, rho]
  if (is.null(par.init)) {
    par.init <- c(rep(1, K+1), 1, 0.1)
  }

  cat("par.init are ", par.init, "\n")

  # Negative log-likelihood function
  neg_loglik <- function(par) {
    sigma_u <- par[1]
    sigma_v <- par[2:(K+1)]
    sigma <- par[K+2]
    rho <- par[K+3]

    # Validate parameters
    if (sigma_u <= 0 || any(sigma_v <= 0) || sigma <= 0) return(Inf)
    if (rho <= -1/(py-1) || rho >= 1) return(Inf)

    # Create exchangeable Sigma_e
    Sigma_e <- matrix(rho*sigma^2, py, py)
    diag(Sigma_e) <- sigma^2

    # Pre-compute Sigma_e inverse and logdet
    Sigma_e_inv <- solve(Sigma_e)
    logdet_Sigma_e <- as.numeric(determinant(Sigma_e, logarithm = TRUE)$modulus)

    logdet_sum <- 0
    quad_form_sum <- matrix(0, py, py)
    N <- 0

    for (h in 1:K) {
      sh <- names(ShXYZ)[h]
      ShZ <- ShXYZ[[sh]]$ShZ
      ShZY <- ShXYZ[[sh]]$ShZY
      ShY <- ShXYZ[[sh]]$ShY
      Nh <- ShXYZ[[sh]]$Nh

      pzh <- ncol(ShZ)

      # Construct exchangeable D_h matrix for this hospital
      D_h <- diag(c(sigma_u^2, rep(sigma_v[h]^2, pzh-1)), pzh)

      # Compute log|I + Z'Z D| using matrix determinant lemma
      logdet_term <- as.numeric(determinant(diag(1, pzh) + ShZ %*% D_h, logarithm = TRUE)$modulus)
      logdet_sum <- logdet_sum + logdet_term

      # Compute quadratic form Y'(I - Z(D^{-1} + Z'Z)^{-1}Z')Y using Woodbury identity
      D_h_inv <- solve(D_h)
      W_h <- solve(D_h_inv + ShZ)
      quad_term <- ShY - t(ShZY) %*% W_h %*% ShZY
      quad_form_sum <- quad_form_sum + quad_term

      N <- N + Nh
    }

    # Compute log-likelihood
    if (reml) {
      # REML version: includes additional term for fixed effects
      loglik <- -0.5 * (logdet_sum + (N - px) * logdet_Sigma_e +
                          sum(diag(Sigma_e_inv %*% quad_form_sum)))
    } else {
      # ML version
      loglik <- -0.5 * (logdet_sum + N * logdet_Sigma_e +
                          sum(diag(Sigma_e_inv %*% quad_form_sum)))
    }

    return(-loglik)
  }

  # Optimize with constraints: sigma > 0, -1/(py-1) < rho < 1
  opt <- optim(par.init, neg_loglik, method = "L-BFGS-B",
               lower = c(rep(1e-6, K+1), 1e-6, -1/(py-1)+1e-6),
               upper = c(rep(Inf, K+1), Inf, 1-1e-6))

  # Extract parameters
  sigma_u <- opt$par[1]
  sigma_v <- opt$par[2:(K+1)]
  sigma <- opt$par[K+2]
  rho <- opt$par[K+3]

  # Create exchangeable Sigma_e
  Sigma_e <- matrix(rho*sigma^2, py, py)
  diag(Sigma_e) <- sigma^2

  cat("Number of function evaluations used to estimate D and Sigma:", opt$counts[1], '\n')
  cat(sigma_u, '\n', sigma_v, '\n', sigma, rho, '\n')

  list(sigma_u = sigma_u, sigma_v = sigma_v, sigma = sigma, rho = rho,
       Sigma_e = Sigma_e, opt = opt)
}




## 3. Function to Estimate Fixed Effects (Run at Local Sites) ##
estimate_fixed_effects <- function(Y, X, Z, id.site, D, Sigma_e, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(Y))
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  py <- ncol(Y)
  K <- length(id.site.uniq)

  sum_xt_sigmax <- matrix(0, px * py, px * py)
  sum_xt_sigmay <- matrix(0, px * py, py)

  # Pre-compute inverse of Sigma_e once
  Sigma_e_inv <- solve(Sigma_e)

  # Pre-allocate memory for repeated operations
  eye_py <- diag(py)
  eye_py_kron <- kronecker(eye_py, matrix(1, 1, 1)) # Placeholder for efficient kronecker later

  for (h in 1:K) {
    sh <- id.site.uniq[h]
    idx <- id.site == sh
    wth <- weights[idx]
    Xh <- X[idx, , drop = FALSE]
    Yh <- Y[idx, , drop = FALSE]
    Zh <- Z[[h]]

    pzh <- ncol(Zh)
    Nh <- nrow(Zh)

    # More efficient D_local construction
    if (ncol(D) != pzh) {
      D_local <- matrix(0, pzh, pzh)
      diag(D_local) <- c(D[1,1], rep(D[2,2], pzh - 1))
    } else {
      D_local <- D
    }

    # Efficient kronecker products using linear algebra properties
    # Instead of full kronecker, we compute block-wise operations
    Zh_multi <- matrix(0, Nh * py, pzh * py)
    for (j in 1:py) {
      rows <- ((j-1)*Nh + 1):(j*Nh)
      cols <- ((j-1)*pzh + 1):(j*pzh)
      Zh_multi[rows, cols] <- Zh
    }

    # Efficient construction of D_multi
    D_multi <- matrix(0, pzh * py, pzh * py)
    for (i in 1:py) {
      for (j in 1:py) {
        rows <- ((i-1)*pzh + 1):(i*pzh)
        cols <- ((j-1)*pzh + 1):(j*pzh)
        D_multi[rows, cols] <- ifelse(i == j, D_local, 0)
      }
    }

    # Construct Vh using Woodbury identity to avoid large matrix inversion
    ZhD <- Zh_multi %*% D_multi
    Vh <- tcrossprod(ZhD, Zh_multi)
    diag_Vh <- kronecker(diag(Nh), Sigma_e)
    Vh <- Vh + diag_Vh

    # Efficient inversion using Cholesky decomposition
    chol_Vh <- chol(Vh)
    Vh_inv <- chol2inv(chol_Vh)

    # Weighted Vh_inv
    if (!all(wth == 1)) {
      sqrt_wth <- sqrt(wth)
      Vh_inv <- Vh_inv * tcrossprod(sqrt_wth)
    }

    # Efficient construction of Xh_multi
    Xh_multi <- matrix(0, Nh * py, px * py)
    for (j in 1:py) {
      rows <- ((j-1)*Nh + 1):(j*Nh)
      cols <- ((j-1)*px + 1):(j*px)
      Xh_multi[rows, cols] <- Xh
    }

    # Update sum_xt_sigmax efficiently
    sum_xt_sigmax <- sum_xt_sigmax + crossprod(Xh_multi, Vh_inv %*% Xh_multi)

    # Update sum_xt_sigmay efficiently without full Yh vectorization
    for (j in 1:py) {
      Yh_j <- numeric(Nh * py)
      Yh_j[((j-1)*Nh + 1):(j*Nh)] <- Yh[, j]
      sum_xt_sigmay[, j] <- sum_xt_sigmay[, j] + crossprod(Xh_multi, Vh_inv %*% Yh_j)
    }
  }

  # Solve using Cholesky decomposition for stability and speed
  chol_sum <- chol(sum_xt_sigmax)
  B <- backsolve(chol_sum, forwardsolve(t(chol_sum), sum_xt_sigmay))

  B <- cbind(B[1:9,1],B[10:18,2],B[19:27,3])

  # print(dim(B))

  # Compute standard errors using Cholesky
  B_se <- matrix(sqrt(diag(chol2inv(chol_sum))), px, py)

  list(B = B, B_se = B_se)
}



## 4. Main Fitting Function ##
federated_lmm_multivariate <- function(Y, X, Z, id.site, weights = NULL,
                                       max_iter = 2, tol = 1e-8, verbose = TRUE) {
  # Initialize
  conv <- FALSE
  iter <- 0
  prev_loglik <- -Inf
  D <- NULL
  Sigma_e <- NULL
  B <- NULL

  while (!conv && iter < max_iter) {
    iter <- iter + 1
    if (verbose) cat("\nIteration:", iter, "\n")

    # Step 1: Calculate summary statistics (or residuals if iter > 1)
    if (iter == 1) {
      ShXYZ <- lmm.get.summary3.multivariate(Y, X, Z, id.site, weights)
    } else {
      # Calculate residuals using current estimates
      residuals <- Y - X %*% B
      ShXYZ <- lmm.get.summary3.multivariate(residuals, X, Z, id.site, weights)
    }

    # Step 2: Estimate covariance parameters
    cov_est <- estimate_covariances(ShXYZ, reml = TRUE, verbose = verbose)
    D <- diag(c(cov_est$sigma_u^2, cov_est$sigma_v^2))
    Sigma_e <- cov_est$Sigma_e

    # Step 3: Estimate fixed effects
    fe_est <- estimate_fixed_effects(Y, X, Z, id.site, D, Sigma_e, weights)
    B <- fe_est$B

    # Check convergence
    current_loglik <- -cov_est$opt$value
    if (verbose) {
      cat("Log-likelihood:", current_loglik, "\n")
      cat("Fixed effects:\n")
      print(B)
    }

    if (abs(current_loglik - prev_loglik) < tol && iter > 1) {
      conv <- TRUE
      if (verbose) cat("Convergence reached.\n")
    }
    prev_loglik <- current_loglik
  }

  # Final estimates
  list(B = B,
       B_se = fe_est$B_se,
       D = D,
       cov_est = cov_est,
       Sigma_e = Sigma_e,
       iterations = iter,
       loglik = current_loglik)
}

