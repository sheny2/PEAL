## =========================
## Multivariate PEAL (MV-PEAL) RI only (not profiling out s2)
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





lmm.noprofile3_mv <- function(par,
                              Y, X, Z, id.site, ShXYZ,
                              reml = TRUE, pooled = FALSE,
                              estimate_rho = TRUE, rho_fixed = 0,
                              verbose = FALSE) {
  
  id.site.uniq <- if (pooled) unique(as.character(id.site)) else names(ShXYZ)
  K  <- length(id.site.uniq)
  px <- if (pooled) ncol(X) else ncol(ShXYZ[[1]]$ShX)
  pz <- K + 1                          # theta_u + theta_vh (per site)
  R  <- attr(ShXYZ, "R"); if (is.null(R)) stop("ShXYZ missing R attribute.")
  
  # -----------------------------
  # Parse parameters
  # par structure:
  #   1..(K+1)         : theta_u, theta_v[h]
  #   (optional) next  : rho (if estimate_rho)
  #   last             : log_s2
  # -----------------------------
  have_rho <- as.logical(estimate_rho)
  if (have_rho) {
    stopifnot(length(par) == (pz + 1 + 1))
    rho     <- par[pz + 1]
    log_s2  <- par[pz + 2]
  } else {
    stopifnot(length(par) == (pz + 1))
    rho     <- rho_fixed
    log_s2  <- par[pz + 1]
  }
  s2 <- exp(log_s2)
  
  # exchangeable Rcorr inverse (scale-free)
  rc <- .Rcorr_inv_and_logdet(R, rho)
  Rinve <- rc$Rinve
  logdet_Rcorr <- rc$logdet_Rcorr
  
  # accumulators
  lpterm1 <- 0
  lpterm2 <- 0
  remlterm <- 0
  bterm1 <- matrix(0, R*px, R*px)
  bterm2 <- rep(0, R*px)
  Nsum <- 0
  
  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    S <- ShXYZ[[sh]]
    ShX  <- S$ShX
    ShXZ <- S$ShXZ
    ShXY <- S$ShXY
    ShZ  <- S$ShZ
    ShZY <- S$ShZY
    ShYY <- S$ShYY
    Nh   <- S$Nh
    
    Nsum <- Nsum + Nh
    pzh <- ncol(ShZ)
    
    theta_u  <- par[1]
    theta_vh <- par[1 + h]
    V0 <- diag(c(theta_u, rep(theta_vh, pzh - 1)), pzh)
    V0_inv <- diag(1/diag(V0), pzh, pzh)
    
    ShX_tilde  <- kronecker(Rinve, ShX)
    ShXZ_tilde <- kronecker(Rinve, ShXZ)
    ShZ_tilde  <- kronecker(Rinve, ShZ)
    
    ShXY_tilde_vec <- as.vector(ShXY %*% Rinve)
    ShZY_tilde_vec <- as.vector(ShZY %*% Rinve)
    YtildeY <- sum(Rinve * ShYY)
    
    Theta_h_inv <- as.matrix(Matrix::bdiag(replicate(R, V0_inv, simplify = FALSE)))
    logdet_Theta_h <- R * as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    
    A_h <- Theta_h_inv + ShZ_tilde
    
    # |Γ_h| decomposition still scale-free w.r.t. Rcorr, but σ^2 will enter globally
    lpterm1 <- lpterm1 + Nh * logdet_Rcorr +
      logdet_Theta_h + as.numeric(determinant(A_h, logarithm = TRUE)$modulus)
    
    Ainv_StXZt <- solve(A_h, t(ShXZ_tilde))
    bterm1 <- bterm1 + (ShX_tilde - ShXZ_tilde %*% Ainv_StXZt)
    
    Ainv_StZY <- solve(A_h, ShZY_tilde_vec)
    bterm2 <- bterm2 + (ShXY_tilde_vec - as.vector(ShXZ_tilde %*% Ainv_StZY))
    
    lpterm2 <- lpterm2 + (YtildeY - drop(t(ShZY_tilde_vec) %*% Ainv_StZY))
  }
  
  # Solve beta (same as before; keep small ridge for PD)
  L <- chol(bterm1 + 1e-6 * diag(nrow(bterm1)))
  beta_hat <- backsolve(L, forwardsolve(t(L), bterm2))
  
  # Quadratic form (scale-free piece; σ^2 enters next)
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * beta_hat) +
                        t(beta_hat) %*% bterm1 %*% beta_hat)
  
  NtotR <- Nsum * R
  pxR   <- px * R
  
  # REML: log|X' Γ^{-1} X| contribution.
  # bterm1 is scale-free w.r.t. Rcorr; true Γ^{-1} has factor 1/s2.
  # So log|X'Γ^{-1}X| = log|bterm1| - pxR*log(s2)
  if (reml) {
    remlterm <- as.numeric(determinant(bterm1, logarithm = TRUE)$modulus) - pxR * log(s2)
    # (RE)ML log-likelihood with explicit s2:
    # lp = -1/2 [ lpterm1 + NtotR*log(s2) + qterm/s2 + remlterm ]
    lp <- -0.5 * (lpterm1 + NtotR*log(s2) + qterm/s2 + remlterm)
  } else {
    # ML (no X-part): lp = -1/2 [ lpterm1 + NtotR*log(s2) + qterm/s2 ]
    lp <- -0.5 * (lpterm1 + NtotR*log(s2) + qterm/s2)
  }
  
  # convenience (for AIC-like; same “pure ML” kernel)
  lk <- -0.5 * (lpterm1 + NtotR * (1 + log(qterm * 2 * pi / NtotR)))
  
  list(
    lp = lp, b = beta_hat, s2 = s2, rho = rho,
    allterms = list(lpterm1 = lpterm1, lpterm2 = lpterm2,
                    qterm = qterm, remlterm = remlterm,
                    bterm1 = bterm1, bterm2 = bterm2)
  )
}


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
  
  if (is.null(ShXYZ)) ShXYZ <- peal.get.summary_mv(Y, X, Z, id.site = id.site, weights = weights)
  id.site.uniq <- names(ShXYZ); K <- length(id.site.uniq)
  px <- ncol(ShXYZ[[1]]$ShX)
  pz <- K + 1
  
  # init: theta's (+ rho if exchangeable) + log_s2
  if (is.null(mypar.init)) {
    base <- rep(1, pz)
    if (estimate_rho && corstr == "exchangeable") base <- c(base, rho_init)
    mypar.init <- c(base, log_s2 = 0)   # log_s2 init at 0 -> s2 = 1
    if (verbose) cat('Default mypar.init =', mypar.init, '\n')
  } else {
    if (estimate_rho && corstr == "exchangeable" && length(mypar.init) == pz) {
      mypar.init <- c(mypar.init, rho_init)
    }
    # ensure we have room for log_s2
    if (length(mypar.init) == (pz + as.integer(corstr=="exchangeable"))) {
      mypar.init <- c(mypar.init, log_s2 = 0)
    }
  }
  
  fn <- function(parameter) {
    if (corstr == "independence") {
      # parameter = [thetas..., log_s2]
      return(-lmm.noprofile3_mv(par = parameter[1:(pz+1)],
                                Y, X, Z, id.site, ShXYZ,
                                reml = reml, pooled = FALSE,
                                estimate_rho = FALSE, rho_fixed = 0,
                                verbose = FALSE)$lp)
    } else {
      # parameter = [thetas..., rho, log_s2]
      return(-lmm.noprofile3_mv(par = parameter,
                                Y, X, Z, id.site, ShXYZ,
                                reml = reml, pooled = FALSE,
                                estimate_rho = TRUE,
                                verbose = FALSE)$lp)
    }
  }
  
  # bounds
  lower <- rep(1e-5, length(mypar.init))
  upper <- rep(Inf,  length(mypar.init))
  if (corstr == "exchangeable") {
    rc <- .Rcorr_inv_and_logdet(R, rho = 0)
    lower[length(lower)-1] <- rc$lower + 1e-8     # bound for rho
    upper[length(upper)-1] <- rc$upper - 1e-8
  }
  # last parameter is log_s2: leave unbounded in R, or clamp wide
  lower[length(lower)] <- -Inf
  upper[length(upper)] <-  Inf
  
  res <- optim(mypar.init, fn, method = "L-BFGS-B",
               hessian = hessian, lower = lower, upper = upper)
  
  if (verbose) cat(ifelse(all(res$convergence == 0), "Convergence Reached", "Non-convergence!"),
                   "and",
                   ifelse(all(eigen(res$hessian)$value > 0), "Hessian PD", "Hessian not PD"), "\n",
                   "The number of function evaluations used is ", res$counts[1], '\n')
  
  # Evaluate at optimum (no profiling)
  res.profile <- if (corstr == "independence") {
    lmm.noprofile3_mv(par = res$par[1:(pz+1)],
                      Y, X, Z, id.site, ShXYZ,
                      reml = reml, pooled = FALSE,
                      estimate_rho = FALSE, rho_fixed = 0,
                      verbose = FALSE)
  } else {
    lmm.noprofile3_mv(par = res$par,
                      Y, X, Z, id.site, ShXYZ,
                      reml = reml, pooled = FALSE,
                      estimate_rho = TRUE,
                      verbose = FALSE)
  }
  
  rho_hat <- if (corstr == "independence") 0 else res.profile$rho
  s2 <- res.profile$s2
  
  Vbeta <- solve(res.profile$allterms$bterm1) * s2
  se    <- sqrt(diag(Vbeta))
  wald  <- res.profile$b / se
  lb    <- res.profile$b - 1.96 * se
  ub    <- res.profile$b + 1.96 * se
  
  Rcorr <- (1 - rho_hat) * diag(R) + rho_hat * matrix(1, R, R)
  Sigma_e <- s2 * Rcorr
  
  list(
    b = res.profile$b, b.sd = se, wald = wald, lb = lb, ub = ub,
    theta = res$par[1:pz], rho = rho_hat, Sigma_e = Sigma_e, s2 = s2,
    opt = res, res.profile = res.profile
  )
}
