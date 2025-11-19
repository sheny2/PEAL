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
# General (unstructured) residual correlation utilities
# Parameterization: free lower-triangular (unit diag) -> S = L L^T -> Corr = cov2cor(S)
# par_corr length = R*(R-1)/2 (row-major lower-tri order: (2,1), (3,1),(3,2), ... )
# Returns: list(Corr, Rinve, logdet)
# -------------------------
.Corr_from_cholfree <- function(R, par_corr, jitter = 1e-8) {
  stopifnot(length(par_corr) == R*(R-1)/2)
  L <- diag(1, R)
  idx <- 1
  for (i in 2:R) {
    for (j in 1:(i-1)) {
      L[i, j] <- par_corr[idx]; idx <- idx + 1
    }
  }
  S <- L %*% t(L)                      # SPD by construction
  Corr <- cov2cor(S)                   # scale to correlation (unit diag)
  # stabilize and invert
  Cj <- Corr + diag(jitter, R)
  U  <- chol(Cj)
  Rinve <- chol2inv(U)                 # C^{-1}
  logdet <- 2 * sum(log(diag(U)))      # log|C|
  list(Corr = Corr, Rinve = Rinve, logdet_Rcorr = logdet)
}



lmm.profile3_mv <- function(par,
                            Y, X, Z, id.site, ShXYZ,
                            reml = TRUE, pooled = FALSE,
                            estimate_rho = TRUE, rho_fixed = 0,
                            corstr = c("exchangeable","independence","unstructured"),
                            verbose = FALSE) {
  
  corstr <- match.arg(corstr)
  id.site.uniq <- if (pooled) unique(as.character(id.site)) else names(ShXYZ)
  K  <- length(id.site.uniq)
  px <- if (pooled) ncol(X) else ncol(ShXYZ[[1]]$ShX)
  pz <- K + 1
  R  <- attr(ShXYZ, "R"); if (is.null(R)) stop("ShXYZ missing R attribute.")
  
  # ----- parse parameters -----
  if (corstr == "independence") {
    if (length(par) != pz) stop("par length must be K+1 for independence.")
    rho <- 0
    rc <- list(Rinve = diag(R), logdet_Rcorr = 0)
  } else if (corstr == "exchangeable") {
    if (estimate_rho) {
      if (length(par) != (pz + 1)) stop("par must include rho as last element.")
      rho <- par[pz + 1]
    } else {
      rho <- rho_fixed
      if (length(par) != pz) stop("par length must be K+1 when rho is fixed.")
    }
    rc <- .Rcorr_inv_and_logdet(R, rho)
  } else { # unstructured
    q <- R*(R-1)/2
    if (!estimate_rho && q > 0) stop("For 'unstructured', estimate_rho must be TRUE.")
    if (length(par) != (pz + q)) stop(sprintf("par must be length %d (thetas) + %d (corr).", pz, q))
    par_corr <- par[(pz + 1):(pz + q)]
    rc <- .Corr_from_cholfree(R, par_corr)
    rho <- NA_real_  # not a single rho; we'll return Corr instead
  }
  
  Rinve <- rc$Rinve
  logdet_Rcorr <- rc$logdet_Rcorr
  
  # ----- accumulators -----
  lpterm1 <- 0; lpterm2 <- 0; remlterm <- 0
  bterm1 <- matrix(0, R*px, R*px)
  bterm2 <- rep(0, R*px)
  Nsum <- 0
  
  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    S <- ShXYZ[[sh]]
    ShX  <- S$ShX;  ShXZ <- S$ShXZ; ShXY <- S$ShXY
    ShZ  <- S$ShZ;  ShZY <- S$ShZY; ShYY <- S$ShYY
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
    
    Theta_h_inv <- as.matrix(bdiag(replicate(R, V0_inv, simplify = FALSE)))
    logdet_Theta_h <- R * as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    
    A_h <- Theta_h_inv + ShZ_tilde
    
    lpterm1 <- lpterm1 + Nh * logdet_Rcorr +
      logdet_Theta_h + as.numeric(determinant(A_h, logarithm = TRUE)$modulus)
    
    Ainv_StXZt <- solve(A_h, t(ShXZ_tilde))
    bterm1 <- bterm1 + (ShX_tilde - ShXZ_tilde %*% Ainv_StXZt)
    
    Ainv_StZY  <- solve(A_h, ShZY_tilde_vec)
    bterm2 <- bterm2 + (ShXY_tilde_vec - as.vector(ShXZ_tilde %*% Ainv_StZY))
    
    lpterm2 <- lpterm2 + (YtildeY - drop(t(ShZY_tilde_vec) %*% Ainv_StZY))
  }
  
  L <- chol(bterm1 + 1e-6 * diag(nrow(bterm1)))
  beta_hat <- backsolve(L, forwardsolve(t(L), bterm2))
  
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * beta_hat) + t(beta_hat) %*% bterm1 %*% beta_hat)
  
  NtotR <- Nsum * R
  pxR   <- px * R
  
  if (reml) {
    rB <- Matrix::rankMatrix(bterm1, tol = 1e-10)
    df <- NtotR - as.numeric(rB)
    s2 <- qterm / df
    remlterm <- determinant(bterm1 / s2, logarithm = TRUE)$modulus
    lp <- -(lpterm1 + NtotR * log(s2) + qterm / s2 + remlterm) / 2
  } else {
    s2 <- qterm / NtotR
    lp <- -(lpterm1 + (1 + log(qterm * 2 * pi / NtotR)) * NtotR) / 2
  }
  
  lk <- - (lpterm1 + (1 + log(qterm * 2 * pi / NtotR)) * NtotR) / 2
  
  list(
    lp = lp, b = beta_hat, s2 = s2, lk = lk, rho = rho,
    # expose Corr if unstructured (else NULL)
    Corr = if (corstr == "unstructured") rc$Corr else NULL,
    allterms = list(lpterm1 = lpterm1, lpterm2 = lpterm2,
                    qterm = qterm, remlterm = remlterm,
                    bterm1 = bterm1, bterm2 = bterm2)
  )
}



peal.fit.RI_mv <- function(Y, X, Z, id.site, weights = NULL,
                           pooled = FALSE, reml = TRUE,
                           common.s2 = TRUE,
                           ShXYZ = NULL,
                           corstr = c('exchangeable','independence','unstructured'),
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
  pz <- K + 1
  q  <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  # ----- init -----
  if (is.null(mypar.init)) {
    mypar.init <- rep(1, pz)
    if (corstr == "exchangeable" && estimate_rho) {
      mypar.init <- c(mypar.init, rho_init)
    } else if (corstr == "unstructured") {
      # start near independence: zeros in cholfree => Corr ~ I
      mypar.init <- c(mypar.init, rep(0, q))
    }
    if (verbose) cat('Default mypar.init (theta_u + theta_v[h] [+ corr params]) =', mypar.init, '\n')
  } else {
    # if provided only thetas, append rho/corr init as needed
    if (length(mypar.init) == pz) {
      if (corstr == "exchangeable" && estimate_rho) mypar.init <- c(mypar.init, rho_init)
      if (corstr == "unstructured")               mypar.init <- c(mypar.init, rep(0, q))
    }
  }
  
  # ----- objective -----
  fn <- function(parameter) {
    if (corstr == "independence") {
      -lmm.profile3_mv(par = parameter[1:pz],
                       Y, X, Z, id.site, ShXYZ,
                       reml = reml, pooled = FALSE,
                       estimate_rho = FALSE, rho_fixed = 0,
                       corstr = "independence",
                       verbose = FALSE)$lp
    } else if (corstr == "exchangeable") {
      -lmm.profile3_mv(par = parameter,
                       Y, X, Z, id.site, ShXYZ,
                       reml = reml, pooled = FALSE,
                       estimate_rho = TRUE,
                       corstr = "exchangeable",
                       verbose = FALSE)$lp
    } else { # unstructured
      -lmm.profile3_mv(par = parameter,
                       Y, X, Z, id.site, ShXYZ,
                       reml = reml, pooled = FALSE,
                       estimate_rho = TRUE,
                       corstr = "unstructured",
                       verbose = FALSE)$lp
    }
  }
  
  # ----- bounds -----
  # thetas >= 1e-5; exchangeable rho within admissible; unstructured corr params are free
  lower <- rep(1e-5, pz)
  upper <- rep(Inf,  pz)
  if (corstr == "exchangeable") {
    rc <- .Rcorr_inv_and_logdet(R, rho = 0)
    lower <- c(lower, rc$lower + 1e-8)
    upper <- c(upper, rc$upper - 1e-8)
  } else if (corstr == "unstructured") {
    lower <- c(lower, rep(-Inf, q))
    upper <- c(upper, rep( Inf, q))
  }
  
  res <- optim(mypar.init, fn, method = "L-BFGS-B",
               hessian = hessian, lower = lower, upper = upper)
  
  if (verbose) cat(ifelse(all(res$convergence == 0), "Convergence Reached", "Non-convergence!"),
                   "and",
                   ifelse(all(eigen(res$hessian)$value > 0), "Hessian PD", "Hessian not PD"), "\n",
                   "Function evaluations:", res$counts[1], '\n')
  
  # ----- recompute profile and assemble outputs -----
  if (corstr == "independence") {
    res.profile <- lmm.profile3_mv(par = res$par[1:pz],
                                   Y, X, Z, id.site, ShXYZ,
                                   reml = reml, pooled = FALSE,
                                   estimate_rho = FALSE, rho_fixed = 0,
                                   corstr = "independence",
                                   verbose = FALSE)
    rho_hat <- 0; Corr_hat <- diag(R)
  } else if (corstr == "exchangeable") {
    res.profile <- lmm.profile3_mv(par = res$par,
                                   Y, X, Z, id.site, ShXYZ,
                                   reml = reml, pooled = FALSE,
                                   estimate_rho = TRUE,
                                   corstr = "exchangeable",
                                   verbose = FALSE)
    rho_hat <- res.profile$rho; Corr_hat <- (1 - rho_hat) * diag(R) + rho_hat * matrix(1, R, R)
  } else {
    res.profile <- lmm.profile3_mv(par = res$par,
                                   Y, X, Z, id.site, ShXYZ,
                                   reml = reml, pooled = FALSE,
                                   estimate_rho = TRUE,
                                   corstr = "unstructured",
                                   verbose = FALSE)
    rho_hat <- NA_real_; Corr_hat <- res.profile$Corr
  }
  
  s2 <- res.profile$s2
  Lb <- chol(res.profile$allterms$bterm1)
  Vbeta <- chol2inv(Lb) * s2
  se    <- sqrt(diag(Vbeta))
  wald  <- res.profile$b / se
  lb    <- res.profile$b - 1.96 * se
  ub    <- res.profile$b + 1.96 * se
  
  Sigma_e <- s2 * Corr_hat
  
  list(
    b = res.profile$b,
    b.sd = se, wald = wald, lb = lb, ub = ub,
    theta = res$par[1:pz],
    rho = if (corstr == "exchangeable") rho_hat else NULL,
    Corr = Corr_hat,             # always return the full Corr
    Sigma_e = Sigma_e,
    s2 = s2,
    opt = res,
    res.profile = res.profile
  )
}


# -------------
# %||% helper
# -------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
