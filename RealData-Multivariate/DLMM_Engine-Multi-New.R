library(Matrix)
library(mvtnorm)
library(MASS)  # For multivariate normal functions

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
    Zh <- Z[[h]]  # Random effects design matrix for site h

    # Compute summary statistics accounting for multivariate outcome
    ShX  <- crossprod((Xh * wth), Xh)  # px × px
    ShXZ <- lapply(1:py, function(j) crossprod((Xh * wth), Zh))  # List of px × pz matrices (one per outcome)
    ShXY <- crossprod((Xh * wth), Yh)  # px × py
    ShZ  <- lapply(1:py, function(j) crossprod((Zh * wth), Zh))  # List of pz × pz matrices
    ShZY <- lapply(1:py, function(j) crossprod((Zh * wth), Yh[,j]))  # List of pz × 1 vectors
    ShY  <- crossprod((Yh * sqrt(wth))) # py × py covariance of outcomes
    Nh <- sum(id.site == sh)

    ShXYZ[[sh]] <- list(ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ = ShZ, ShZY = ShZY, ShY = ShY, Nh = Nh)
  }

  return(ShXYZ)
}



#-----------------------------------------------------------------------
# Federated EM Algorithm for Multivariate Mixed Models
#-----------------------------------------------------------------------
federated_mv_lmm <- function(ShXYZ_list, max_iter = 100, tol = 1e-6) {

  # Initialize parameters
  px <- ncol(ShXYZ_list[[1]]$ShX)  # Number of covariates
  py <- ncol(ShXYZ_list[[1]]$ShY)   # Number of outcomes
  H <- length(ShXYZ_list)           # Number of sites

  # Initial guesses (can be improved with warm starts)
  beta <- matrix(0, nrow = px, ncol = py)      # Fixed effects
  D <- diag(1, py)                             # Site-level covariance
  Sigma_v <- diag(1, py)                       # Patient-level covariance
  Sigma_e <- diag(1, py)                       # Residual covariance (exchangeable)
  rho <- 0.3                                   # Initial correlation

  for (iter in 1:max_iter) {
    #-------------------------------------------
    # E-step: Compute expected sufficient statistics (locally)
    #-------------------------------------------
    sum_D <- matrix(0, py, py)
    sum_Sigma_v <- matrix(0, py, py)
    sum_Sigma_e <- matrix(0, py, py)
    sum_XX <- matrix(0, px, px)
    sum_XY <- matrix(0, px, py)

    for (h in 1:H) {
      Sh <- ShXYZ_list[[h]]

      cat(dim(Sh$ShZ))
      cat(dim(Sigma_v))
      # Construct marginal covariance for site h
      V_h <- Sh$ShZ %*% Sigma_v %*% t(Sh$ShZ) +
        kronecker(diag(Sh$Nh), D) +
        kronecker(diag(Sh$Nh * ncol(Sh$ShZ)/py), Sigma_e)

      # Compute E[u_h] and E[v_hi] (pseudo-code, adjust for actual random effects structure)
      E_uh <- D %*% solve(V_h) %*% (Sh$ShZY - Sh$ShXZ %*% beta)
      E_vhi <- Sigma_v %*% solve(V_h) %*% (Sh$ShZY - Sh$ShXZ %*% beta)

      # Accumulate sufficient statistics
      sum_D <- sum_D + tcrossprod(E_uh)
      sum_Sigma_v <- sum_Sigma_v + tcrossprod(E_vhi)
      sum_Sigma_e <- sum_Sigma_e + (Sh$ShY - Sh$ShXY %*% beta - t(E_uh) %*% Sh$ShZY - t(E_vhi) %*% Sh$ShZY)

      # For fixed effects update
      sum_XX <- sum_XX + Sh$ShX
      sum_XY <- sum_XY + Sh$ShXY - Sh$ShXZ %*% (E_uh + E_vhi)
    }

    #-------------------------------------------
    # M-step: Update parameters (globally)
    #-------------------------------------------
    # Update D (site-level covariance)
    D_new <- sum_D / H

    # Update Sigma_v (patient-level covariance)
    total_patients <- sum(sapply(ShXYZ_list, function(x) ncol(x$ShZ)/py))
    Sigma_v_new <- sum_Sigma_v / total_patients

    # Update Sigma_e (residual covariance, exchangeable)
    total_visits <- sum(sapply(ShXYZ_list, function(x) x$Nh))
    Sigma_e_new <- sum_Sigma_e / total_visits

    # Enforce exchangeable structure
    sigma_e <- mean(diag(Sigma_e_new))
    off_diag <- mean(Sigma_e_new[lower.tri(Sigma_e_new)])
    rho_new <- off_diag / sigma_e
    Sigma_e_new <- sigma_e * ((1 - rho_new) * diag(py) + rho_new * matrix(1, py, py))

    # Update beta (fixed effects)
    beta_new <- solve(sum_XX) %*% sum_XY

    #-------------------------------------------
    # Check convergence
    #-------------------------------------------
    if (max(abs(beta_new - beta)) < tol &&
        max(abs(D_new - D)) < tol &&
        max(abs(Sigma_v_new - Sigma_v)) < tol &&
        abs(rho_new - rho) < tol) {
      break
    }

    # Update parameters for next iteration
    beta <- beta_new
    D <- D_new
    Sigma_v <- Sigma_v_new
    Sigma_e <- Sigma_e_new
    rho <- rho_new
  }

  # Return results
  list(
    beta = beta,
    D = D,
    Sigma_v = Sigma_v,
    Sigma_e = Sigma_e,
    rho = rho,
    iterations = iter
  )
}
