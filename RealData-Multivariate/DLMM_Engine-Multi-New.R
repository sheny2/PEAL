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
  px <- ncol(ShXYZ[[1]]$ShX)
  py <- ncol(ShXYZ[[1]]$ShY)

  # Initialize parameters if not provided
  if (is.null(par.init)) {
    # Parameters: [sigma_u, sigma_v1, ..., sigma_vK, vech(Sigma_e)]
    par.init <- c(rep(1, K+1), diag(py)[lower.tri(diag(py), diag = TRUE)])
  }

  # Negative log-likelihood function for optimization
  neg_loglik <- function(par) {
    # Extract parameters
    sigma_u <- par[1]
    sigma_v <- par[2:(K+1)]
    Sigma_e_params <- par[(K+2):length(par)]

    # Reconstruct Sigma_e (ensure positive definite)
    L <- matrix(0, py, py)
    L[lower.tri(L, diag = TRUE)] <- Sigma_e_params
    Sigma_e <- tcrossprod(L)

    logdet_sum <- 0
    quad_form_sum <- matrix(0, py, py)
    N <- 0

    for (h in 1:K) {
      ShZ <- ShXYZ[[h]]$ShZ
      ShZY <- ShXYZ[[h]]$ShZY
      ShY <- ShXYZ[[h]]$ShY
      Nh <- ShXYZ[[h]]$Nh

      pzh <- ncol(ShZ)  # Number of random effects at site h (RI model, it's 1+mh)

      # Construct D with correct dimensions (pzh × pzh)
      D <- diag(c(sigma_u^2, rep(sigma_v[h]^2, (pzh - 1))), pzh)

      Vinv <- solve(D)  # Inverse of random effects covariance

      # Compute Wh = (Vinv + ShZ)^{-1}
      Wh <- solve(Vinv + ShZ)

      # Accumulate terms for covariance estimation
      quad_form_sum <- quad_form_sum + (ShY - t(ShZY) %*% Wh %*% ShZY)
      logdet_sum <- logdet_sum + as.numeric(determinant(diag(1, pzh) + ShZ %*% D, logarithm = TRUE)$modulus)
      N <- N + Nh
    }

    if (reml) {
      Sigma_e <- quad_form_sum / (N - px)
      loglik <- -0.5 * (logdet_sum + (N - px) * determinant(Sigma_e, logarithm = TRUE)$modulus +
                          sum(diag(solve(Sigma_e, quad_form_sum))))
    } else {
      Sigma_e <- quad_form_sum / N
      loglik <- -0.5 * (logdet_sum + N * determinant(Sigma_e, logarithm = TRUE)$modulus +
                          sum(diag(solve(Sigma_e, quad_form_sum))))
    }

    return(-loglik)
  }

  # Optimize covariance parameters
  opt <- optim(par.init, neg_loglik, method = "L-BFGS-B", hessian = T,
               lower = c(rep(1e-6, K+1), rep(-Inf, py*(py+1)/2)))

  cat("The number of function evaluations used is", opt$counts[1], '\n')

  # Extract and return estimated parameters
  sigma_u <- opt$par[1]
  sigma_v <- opt$par[2:(K+1)]
  Sigma_e_params <- opt$par[(K+2):length(opt$par)]

  cat(sigma_u, '\n' sigma_v, '\n' Sigma_e)

  L <- matrix(0, py, py)
  L[lower.tri(L, diag = TRUE)] <- Sigma_e_params
  Sigma_e <- tcrossprod(L)

  list(sigma_u = sigma_u, sigma_v = sigma_v, Sigma_e = Sigma_e, opt = opt)
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
  sum_xt_sigmay <- matrix(0, px * py, py)  # Now (px * py) × py

  Sigma_e_inv <- solve(Sigma_e)

  for (h in 1:K) {
    sh <- id.site.uniq[h]
    wth <- weights[id.site == sh]
    Xh <- X[id.site == sh, ]
    Yh <- Y[id.site == sh, ]  # Nh × py
    Zh <- Z[[h]]

    pzh <- ncol(Zh)
    Nh <- nrow(Zh)

    # Ensure D matches Zh dimensions (pzh × pzh)
    if (ncol(D) != pzh) {
      D_local <- diag(c(D[1,1], rep(D[2,2], (pzh - 1))), pzh)
    } else {
      D_local <- D
    }

    # Construct block-diagonal matrices for multivariate case
    Zh_multi <- kronecker(diag(1, py), Zh)  # (Nh * py) × (pzh * py)
    D_multi <- kronecker(D_local, diag(1, py))  # (pzh * py) × (pzh * py)
    Xh_multi <- kronecker(diag(1, py), Xh)  # (Nh * py) × (px * py)

    # Construct full covariance matrix
    Vh <- Zh_multi %*% D_multi %*% t(Zh_multi) + kronecker(diag(1, Nh), Sigma_e)
    Vh_inv <- solve(Vh)

    # Compute contributions to sum_xt_sigmax and sum_xt_sigmay
    sum_xt_sigmax <- sum_xt_sigmax + t(Xh_multi) %*% (Vh_inv * wth) %*% Xh_multi

    # Properly compute sum_xt_sigmay (avoid vectorizing Yh)
    for (j in 1:py) {
      Yh_j <- matrix(0, Nh * py, 1)
      Yh_j[((j-1)*Nh + 1):(j*Nh)] <- Yh[, j]  # Place response j in the correct block
      sum_xt_sigmay[, j] <- sum_xt_sigmay[, j] + t(Xh_multi) %*% (Vh_inv * wth) %*% Yh_j
    }
  }


  # Estimate fixed effects matrix (reshape to px × py)
  # B <- matrix(solve(sum_xt_sigmax) %*% sum_xt_sigmay, px, py)
  B <- solve(sum_xt_sigmax, sum_xt_sigmay)

  # Compute standard errors (reshape accordingly)
  B_se <- matrix(sqrt(diag(solve(sum_xt_sigmax))), px, py)


  list(B = B, B_se = B_se)

}

## 4. Main Fitting Function ##
federated_lmm_multivariate <- function(Y, X, Z, id.site, weights = NULL,
                                       max_iter = 10, tol = 1e-8, verbose = TRUE) {
  # Step 1: Get initial summary statistics from all sites
  ShXYZ <- lmm.get.summary3.multivariate(Y, X, Z, id.site, weights)

  # Initialize parameters
  conv <- FALSE
  iter <- 0
  prev_loglik <- -Inf
  D <- NULL
  Sigma_e <- NULL

  while (!conv && iter < max_iter) {
    iter <- iter + 1
    if (verbose) cat("Iteration:", iter, "\n")

    # Step 2: Central site estimates covariance parameters
    cov_est <- estimate_covariances(ShXYZ, reml = TRUE, verbose = verbose)

    # Construct D matrix from estimated parameters
    K <- length(unique(id.site))
    D <- diag(c(cov_est$sigma_u^2, cov_est$sigma_v^2))

    # Step 3: Send D and Sigma_e back to sites for fixed effects estimation
    fe_est <- estimate_fixed_effects(Y, X, Z, id.site, D, cov_est$Sigma_e, weights)

    # Check convergence (using log-likelihood)
    current_loglik <- -cov_est$opt$value
    if (verbose) cat("Log-likelihood:", current_loglik, "\n")

    if (abs(current_loglik - prev_loglik) < tol) {
      conv <- TRUE
      if (verbose) cat("Convergence reached.\n")
    }
    prev_loglik <- current_loglik
  }

  # Final estimates
  list(B = fe_est$B,
       B_se = fe_est$B_se,
       D = D,
       Sigma_e = cov_est$Sigma_e,
       iterations = iter,
       loglik = current_loglik)
}
