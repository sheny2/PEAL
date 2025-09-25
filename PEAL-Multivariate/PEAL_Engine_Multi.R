## =========================
## Multivariate PEAL (MV-PEAL) RI only
## Exchangeable residual correlation across outcomes
## =========================

library(tidyverse)
library(data.table)
library(lme4)
library(nlme)
library(Matrix)
library(minqa)

# -------------------------
# helpers for Z_hv
# -------------------------
generate_record_count <- function(data) {
  counts <- table(data[, "n_hi"])
  result_matrix <- cbind(as.numeric(names(counts)), as.numeric(counts))
  colnames(result_matrix) <- c("n_hi", "frequency")
  return(result_matrix)
}

generate_Zhv_matrix <- function(data) {
  record_count_matrix <- generate_record_count(data)
  diagonal_blocks <- list()
  for (i in 1:nrow(record_count_matrix)) {
    n_hi <- record_count_matrix[i, "n_hi"]
    frequency <- record_count_matrix[i, "frequency"]
    identity_block <- diag(frequency)
    ones_vector <- matrix(1, nrow = n_hi, ncol = 1)
    diagonal_blocks[[i]] <- kronecker(identity_block, ones_vector)
  }
  big_matrix <- do.call(Matrix::bdiag, diagonal_blocks)
  big_matrix <- cbind(1, big_matrix)
  return(as.matrix(big_matrix))
}

# -------------------------
# Exchangeable residual correlation utilities
# Rcorr = (1 - rho) I_R + rho * 11^T (scale-free)
# We profile sigma^2 separately; all whitening uses Rcorr^{-1}
# -------------------------
.Rcorr_inv_and_logdet <- function(R, rho, eps = 1e-8) {
  # bounds: rho in (-1/(R-1), 1)
  lower <- -1/(R-1) + eps
  upper <- 1 - eps
  if (rho <= lower || rho >= upper) stop("rho out of admissible range for exchangeable R-corr.")
  # eigenvalues of Rcorr: (1-ρ) with multiplicity (R-1); (1-ρ+ρR) once
  lam1 <- (1 - rho)               # multiplicity R-1
  lam2 <- (1 - rho + rho * R)     # multiplicity 1
  logdet_Rcorr <- (R-1)*log(lam1) + log(lam2)

  # closed form inverse: Rcorr^{-1} = (1/(1-ρ)) [ I - (ρ/(1-ρ+ρR)) 11^T ]
  a <- 1/(1 - rho)
  b <- rho / (1 - rho + rho * R)
  Rinve <- a * (diag(R) - b * matrix(1, R, R))
  list(Rinve = Rinve, logdet_Rcorr = logdet_Rcorr, lower = lower, upper = upper)
}

# -------------------------
# Multivariate summary stats at each site
# Y must be an N x R matrix (one row per visit, R outcomes as columns)
# X: N x p_x (shared across outcomes), Z: list of site-specific Z_h
# returns per-site: ShX, ShXZ, ShZ, ShXY (p_x x R), ShZY (p_zh x R), ShYY (R x R), Nh
# -------------------------
peal.get.summary_mv <- function(Y, X, Z, id.site, weights = NULL) {
  if (!is.matrix(Y)) stop("Y must be an N x R matrix for multivariate.")
  R <- ncol(Y)
  if (is.null(weights)) weights <- rep(1, nrow(X))
  X <- as.matrix(X)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)

  ShXYZ <- list()
  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    idx <- which(id.site == sh)
    wth <- weights[idx]
    Xh <- X[idx, , drop = FALSE]
    Yh <- as.matrix(Y[idx, , drop = FALSE])  # N_h x R
    Zh <- Z[h][[1]]                          # N_h x p_zh

    # weight rows
    Xh_w <- Xh * wth
    Zh_w <- Zh * wth
    Yh_w <- Yh * wth

    ShX  <- crossprod(Xh_w, Xh)             # p_x x p_x
    ShXZ <- crossprod(Xh_w, Zh)             # p_x x p_zh
    ShXY <- crossprod(Xh_w, Yh)             # p_x x R

    ShZ  <- crossprod(Zh_w, Zh)             # p_zh x p_zh
    ShZY <- crossprod(Zh_w, Yh)             # p_zh x R

    ShYY <- crossprod(Yh_w, Yh)             # R x R
    Nh   <- length(idx)

    ShXYZ[[sh]] <- list(ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ = ShZ, ShZY = ShZY, ShYY = ShYY,
                        Nh = Nh)
  }
  attr(ShXYZ, "R") <- ncol(Y)
  return(ShXYZ)
}

# -------------------------
# Multivariate profile likelihood (Woodbury), exchangeable residuals
# Parameters 'par' length = (K + 1) [+ 1 if estimate_rho], same as univariate:
#   par[1] = theta_u   (site RI variance scaled by sigma^2)
#   par[1+h] = theta_vh for site h (patient-level RE variance scaled by sigma^2)
# If estimate_rho = TRUE, append par[length(par)] = rho (in (-1/(R-1), 1))
# -------------------------
lmm.profile3_mv <- function(par,
                            Y, X, Z, id.site, ShXYZ,
                            reml = TRUE, pooled = FALSE,
                            estimate_rho = TRUE, rho_fixed = 0,
                            verbose = FALSE) {

  id.site.uniq <- if (pooled) unique(as.character(id.site)) else names(ShXYZ)
  K  <- length(id.site.uniq)
  px <- if (pooled) ncol(X) else ncol(ShXYZ[[1]]$ShX)
  pz <- K + 1                                 # theta_u + theta_vh (per site)
  R  <- attr(ShXYZ, "R"); if (is.null(R)) stop("ShXYZ missing R attribute.")

  # parse parameters
  if (estimate_rho) {
    if (length(par) != (pz + 1)) stop("par must include rho as last element.")
    rho <- par[pz + 1]
  } else {
    rho <- rho_fixed
    if (length(par) != pz) stop("par length must be K+1 when rho is fixed.")
  }

  # exchangeable Rcorr inverse (scale-free)
  rc <- .Rcorr_inv_and_logdet(R, rho)
  Rinve <- rc$Rinve
  logdet_Rcorr <- rc$logdet_Rcorr

  # accumulators
  lpterm1 <- 0      # log|Rcorr| terms + log|I + Z~ Θ Z~^T| via |Θ| + |Θ^{-1} + Z~^T Z~|
  lpterm2 <- 0      # y~'y~ - y~'Z~(...)Z~'y~
  remlterm <- 0
  bterm1 <- matrix(0, R*px, R*px)  # X~'Γ^{-1}X~
  bterm2 <- rep(0, R*px)           # X~'Γ^{-1}y~
  Nsum <- 0

  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    S <- ShXYZ[[sh]]
    ShX  <- S$ShX
    ShXZ <- S$ShXZ
    ShXY <- S$ShXY      # p_x x R
    ShZ  <- S$ShZ
    ShZY <- S$ShZY      # p_zh x R
    ShYY <- S$ShYY      # R x R
    Nh   <- S$Nh

    Nsum <- Nsum + Nh
    pzh <- ncol(ShZ)

    # Θ block (scaled RE covariance) for this site: diag(theta_u, theta_vh,...)
    theta_u  <- par[1]
    theta_vh <- par[1 + h]
    V0 <- diag(c(theta_u, rep(theta_vh, pzh - 1)), pzh)          # p_zh x p_zh (per outcome)
    V0_inv <- diag(1/diag(V0), pzh, pzh)

    # build Kronecker-whitened design cross-products
    # (~) denotes prewhitening by Rcorr^{-1/2} on outcomes (scale-free)
    ShX_tilde  <- kronecker(Rinve, ShX)                          # (R*px) x (R*px)
    ShXZ_tilde <- kronecker(Rinve, ShXZ)                         # (R*px) x (R*pzh)
    ShZ_tilde  <- kronecker(Rinve, ShZ)                          # (R*pzh) x (R*pzh)

    # Whitened X'Y and Z'Y as vec( X' Y Rinve ), vec( Z' Y Rinve )
    ShXY_tilde_vec <- as.vector(ShXY %*% Rinve)                  # length R*px (col-stacked)
    ShZY_tilde_vec <- as.vector(ShZY %*% Rinve)                  # length R*pzh
    YtildeY <- sum(Rinve * ShYY)                                 # scalar: tr(Rinve * Y'Y)

    # Θ_h and Θ_h^{-1} (block-diagonal across outcomes)
    Theta_h_inv <- as.matrix(bdiag(replicate(R, V0_inv, simplify = FALSE))) # (R*pzh)x(R*pzh)
    logdet_Theta_h <- R * as.numeric(determinant(V0, logarithm = TRUE)$modulus)

    # K_h(Θ) = Θ^{-1} + Z~'Z~   (this is the Woodbury block we solve with)
    A_h <- Theta_h_inv + ShZ_tilde

    # log|Γ_h| = log|Rcorr|^{Nh} + log|Θ_h| + log|A_h|
    lpterm1 <- lpterm1 + Nh * logdet_Rcorr +
      logdet_Theta_h + as.numeric(determinant(A_h, logarithm = TRUE)$modulus)

    # Woodbury reductions
    # X~'Γ^{-1}X~ = X~'X~ - X~'Z~ A_h^{-1} Z~'X~
    Ainv_StXZt <- solve(A_h, t(ShXZ_tilde))                      # (R*pzh) x (R*px)
    bterm1 <- bterm1 + (ShX_tilde - ShXZ_tilde %*% Ainv_StXZt)

    # X~'Γ^{-1}y~ = X~'y~ - X~'Z~ A_h^{-1} Z~'y~
    Ainv_StZY  <- solve(A_h, ShZY_tilde_vec)                     # (R*pzh) x 1
    bterm2 <- bterm2 + (ShXY_tilde_vec - as.vector(ShXZ_tilde %*% Ainv_StZY))

    # y~'Γ^{-1}y~ = y~'y~ - y~'Z~ A_h^{-1} Z~'y~
    lpterm2 <- lpterm2 + (YtildeY - drop(t(ShZY_tilde_vec) %*% Ainv_StZY))
  }

  # Solve for beta (length R*px)
  L <- chol(bterm1)
  beta_hat <- backsolve(L, forwardsolve(t(L), bterm2))

  # Quadratic form for sigma^2 profiling
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * beta_hat) + t(beta_hat) %*% bterm1 %*% beta_hat)

  # Degrees of freedom
  NtotR <- Nsum * R
  pxR   <- px * R

  if (reml) {
    s2 <- qterm / (NtotR - pxR)
    remlterm <- determinant(bterm1 / s2, logarithm = TRUE)$modulus
    # profile REML log-likelihood
    lp <- -(lpterm1 + NtotR * log(s2) + qterm / s2 + remlterm) / 2
  } else {
    s2 <- qterm / NtotR
    lp <- -(lpterm1 + (1 + log(qterm * 2 * pi / NtotR)) * NtotR) / 2
  }

  # convenience (for AIC-like)
  lk <- - (lpterm1 + (1 + log(qterm * 2 * pi / NtotR)) * NtotR) / 2

  res <- list(
    lp = lp, b = beta_hat, s2 = s2, lk = lk, rho = rho,
    allterms = list(lpterm1 = lpterm1, lpterm2 = lpterm2,
                    qterm = qterm, remlterm = remlterm,
                    bterm1 = bterm1, bterm2 = bterm2)
  )
  return(res)
}

# -------------------------
# Top-level fit function for MV-PEAL
# Y: N x R matrix
# corstr: 'exchangeable' (default) or 'independence' (sets rho = 0)
# -------------------------
peal.fit.RI_mv <- function(Y, X, Z, id.site, weights = NULL,
                           pooled = FALSE, reml = TRUE,
                           common.s2 = TRUE,
                           ShXYZ = NULL,
                           corstr = c('exchangeable','independence'),
                           mypar.init = NULL,
                           estimate_rho = TRUE, rho_init = 0.1,
                           hessian = TRUE, verbose = TRUE) {

  corstr <- match.arg(corstr)
  if (!is.matrix(Y)) stop("For multivariate fit, supply Y as N x R matrix.")
  R <- ncol(Y)

  id.site.uniq <- if (pooled) unique(as.character(id.site)) else names(ShXYZ %||% list())
  if (pooled) {
    K  <- length(id.site.uniq)
    px <- ncol(X)
    ShXYZ <- peal.get.summary_mv(Y, X, Z, id.site = id.site, weights = weights)
  } else {
    if (is.null(ShXYZ)) ShXYZ <- peal.get.summary_mv(Y, X, Z, id.site = id.site, weights = weights)
    id.site.uniq <- names(ShXYZ)
    K  <- length(id.site.uniq)
    px <- ncol(ShXYZ[[1]]$ShX)
  }
  pz <- K + 1  # theta_u + theta_vh per site

  # Initialize parameters: thetas (scaled by sigma^2), optionally rho
  if (is.null(mypar.init)) {
    mypar.init <- rep(1, pz)
    if (estimate_rho && corstr == "exchangeable") {
      mypar.init <- c(mypar.init, rho_init)
    }
    if (verbose) cat('Default mypar.init (theta_u + theta_v[h]) =', mypar.init, '\n')
  } else {
    if (estimate_rho && corstr == "exchangeable" && length(mypar.init) == pz) {
      mypar.init <- c(mypar.init, rho_init)
    }
  }

  # objective
  fn <- function(parameter) {
    if (corstr == "independence") {
      return(-lmm.profile3_mv(par = parameter[1:pz],
                              Y, X, Z, id.site, ShXYZ,
                              reml = reml, pooled = FALSE,
                              estimate_rho = FALSE, rho_fixed = 0,
                              verbose = FALSE)$lp)
    } else {
      return(-lmm.profile3_mv(par = parameter,
                              Y, X, Z, id.site, ShXYZ,
                              reml = reml, pooled = FALSE,
                              estimate_rho = TRUE,
                              verbose = FALSE)$lp)
    }
  }

  # bounds: theta params >= 1e-6; rho in (-1/(R-1)+eps, 1-eps)
  lower <- rep(1e-6, length(mypar.init))
  upper <- rep(Inf,  length(mypar.init))
  if (corstr == "exchangeable") {
    rc <- .Rcorr_inv_and_logdet(R, rho = 0)  # just to get bounds
    lower[length(lower)] <- rc$lower
    upper[length(upper)] <- rc$upper - 1e-8
  }

  res <- optim(mypar.init, fn, method = "L-BFGS-B",
               hessian = hessian, lower = lower, upper = upper)

  if (verbose) cat(ifelse(all(res$convergence == 0), "Convergence Reached", "Non-convergence!"),
                   "and",
                   ifelse(all(eigen(res$hessian)$value > 0), "Hessian PD", "Hessian not PD"), "\n",
                   "The number of function evaluations used is ", res$counts[1], '\n')

  # recompute profile at optimum
  if (corstr == "independence") {
    res.profile <- lmm.profile3_mv(par = res$par[1:pz],
                                   Y, X, Z, id.site, ShXYZ,
                                   reml = reml, pooled = FALSE,
                                   estimate_rho = FALSE, rho_fixed = 0,
                                   verbose = FALSE)
    rho_hat <- 0
  } else {
    res.profile <- lmm.profile3_mv(par = res$par,
                                   Y, X, Z, id.site, ShXYZ,
                                   reml = reml, pooled = FALSE,
                                   estimate_rho = TRUE,
                                   verbose = FALSE)
    rho_hat <- res.profile$rho
  }

  s2 <- res.profile$s2

  # Inference: Var(vec(beta)) = s2 * (bterm1)^{-1}
  Vbeta <- solve(res.profile$allterms$bterm1) * s2
  se    <- sqrt(diag(Vbeta))
  wald  <- res.profile$b / se
  lb    <- res.profile$b - 1.96 * se
  ub    <- res.profile$b + 1.96 * se

  # Return residual covariance Sigma_e = s2 * Rcorr(rho)
  Rcorr <- (1 - rho_hat) * diag(R) + rho_hat * matrix(1, R, R)
  Sigma_e <- s2 * Rcorr

  out <- list(
    b = res.profile$b,                # length R*px (stacked by outcomes)
    b.sd = se,
    wald = wald, lb = lb, ub = ub,
    theta = res$par[1:pz],            # scaled RE variances (Θ params)
    rho = rho_hat,
    Sigma_e = Sigma_e,
    s2 = s2,
    opt = res,
    res.profile = res.profile
  )
  return(out)
}

# -------------
# %||% helper
# -------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
