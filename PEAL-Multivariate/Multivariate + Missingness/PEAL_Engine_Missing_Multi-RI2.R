# MV-PEAL with missing outcomes, random intercepts only

library(Matrix)

# Encode which outcomes are observed on a row (length-R "01" string, e.g., "110")
row_pattern_key <- function(yrow) {
  paste0(ifelse(!is.na(yrow), "1", "0"), collapse = "")
}

# Turn a key like "101" into outcome indices (1-based)
key_to_idx <- function(key) {
  which(strsplit(key, "", fixed = TRUE)[[1]] == "1")
}

# Build outcome-selection matrix E_pi (R x s) from a key and R
selector_from_key <- function(key, R) {
  idx <- key_to_idx(key)
  s <- length(idx)
  E <- matrix(0, nrow = R, ncol = s)
  if (s > 0) E[cbind(idx, seq_len(s))] <- 1
  E
}

# Robust site-specific Z_h (site RI + patient-within-site RIs), order-safe.
# data: a data.frame containing at least 'site', 'patient' columns, already ordered as intended.
# Returns a list Z_list with one matrix per site, aligned to the rows of that site.
build_Z_list_by_site <- function(data, site_col = "site", patient_col = "patient") {
  sites <- as.character(data[[site_col]])
  pats  <- as.character(data[[patient_col]])
  site_levels <- unique(sites)
  Z_list <- vector("list", length(site_levels))
  names(Z_list) <- site_levels
  for (sh in site_levels) {
    idx <- which(sites == sh)
    # patient within *this* site
    pat_fac <- factor(pats[idx], levels = unique(pats[idx]))
    # Z_patient: one column per patient
    Z_pat <- model.matrix(~ pat_fac + 0L)
    # prepend site intercept column of ones
    Z_h <- cbind(Intercept_site = 1, Z_pat)
    Z_list[[sh]] <- as.matrix(Z_h)
  }
  Z_list
}



.Rcorr_inv_and_logdet <- function(R, rho, eps = 1e-8) {
  lower <- -1/(R-1) + eps
  upper <- 1 - eps
  if (rho <= lower || rho >= upper) stop("rho out of admissible range.")
  lam1 <- (1 - rho)               # multiplicity R-1
  lam2 <- (1 - rho + rho * R)     # multiplicity 1
  logdet_Rcorr <- (R-1)*log(lam1) + log(lam2)
  a <- 1/(1 - rho)
  b <- rho / (1 - rho + rho * R)
  Rinve <- a * (diag(R) - b * matrix(1, R, R))
  list(Rinve = Rinve, logdet_Rcorr = logdet_Rcorr, lower = lower, upper = upper)
}

# Subset version for size s (used per pattern)
.Rcorr_inv_and_logdet_s <- function(s, rho, eps = 1e-8) {
  if (s < 1) stop("Pattern of size 0 not allowed.")
  if (s == 1) {
    # degenerate: correlation is 1x1
    return(list(Rinve = matrix(1,1,1), logdet_Rcorr = 0,
                lower = -Inf, upper = Inf))
  }
  lower <- -1/(s-1) + eps
  upper <- 1 - eps
  if (rho <= lower || rho >= upper) stop(sprintf("rho out of range for s=%d.", s))
  a <- 1/(1 - rho)
  b <- rho / (1 - rho + rho * s)
  Rinve <- a * (diag(s) - b * matrix(1, s, s))
  logdet_Rcorr <- (s-1)*log(1 - rho) + log(1 - rho + rho * s)
  list(Rinve = Rinve, logdet_Rcorr = logdet_Rcorr, lower = lower, upper = upper)
}




# Build lossless one-shot summaries stratified by (site, outcome-missingness pattern).
# Inputs:
#   data: data.frame with columns site, patient, X1..Xpx, Y1..YR
#   X_cols: names of X columns
#   Y_cols: names of Y columns (length R)
#   weights: optional vector length N (defaults to 1)
# Returns:
#   nested list ShXYZ[[site]][[pattern_key]] with S^{..} and metadata.
peal.get.summary_mv_patterns <- function(data, X_cols, Y_cols,
                                         site_col = "site", patient_col = "patient",
                                         weights = NULL) {
  X <- as.matrix(data[, X_cols, drop = FALSE])
  Y <- as.matrix(data[, Y_cols, drop = FALSE])
  R <- ncol(Y)
  if (is.null(weights)) weights <- rep(1, nrow(X))
  sites <- as.character(data[[site_col]])
  
  Z_list <- build_Z_list_by_site(data, site_col = site_col, patient_col = patient_col)
  
  Sh <- list()
  pat_keys <- apply(Y, 1, row_pattern_key)
  valid_row <- (pat_keys != paste0(rep("0", R), collapse = ""))
  
  for (sh in unique(sites)) {
    site_idx_all <- which(sites == sh)        # global indices for this site
    site_valid   <- valid_row[site_idx_all]    # logical mask, same length as this site's rows
    idx_site     <- site_idx_all[site_valid]   # global indices of valid rows
    if (!length(idx_site)) next
    
    Xh <- X[idx_site, , drop = FALSE]
    Yh <- Y[idx_site, , drop = FALSE]
    wh <- weights[idx_site]
    
    # --- key change: do NOT index Z by global indices ---
    Zh_site <- Z_list[[sh]]                    # already nrow == length(site_idx_all)
    Zh      <- Zh_site[site_valid, , drop = FALSE]  # align to filtered rows
    
    keys_h       <- pat_keys[idx_site]
    keys_h_uniq  <- unique(keys_h)
    
    Sh[[sh]] <- list()
    for (key in keys_h_uniq) {
      rows <- which(keys_h == key)             # positions within *filtered* block
      if (!length(rows)) next
      
      Xs <- Xh[rows, , drop = FALSE]
      Zs <- Zh[rows, , drop = FALSE]
      Ys_full <- Yh[rows, , drop = FALSE]
      
      idx_outcomes <- key_to_idx(key)
      s  <- length(idx_outcomes)
      Ys <- Ys_full[, idx_outcomes, drop = FALSE]
      
      ws <- wh[rows]
      Xw <- Xs * ws
      Zw <- Zs * ws
      
      ShX  <- crossprod(Xw, Xs)
      ShXZ <- crossprod(Xw, Zs)
      ShXY <- crossprod(Xw, Ys)
      
      ShZ  <- crossprod(Zw, Zs)
      ShZY <- crossprod(Zw, Ys)
      ShYY <- crossprod(Ys * ws, Ys)
      
      Sh[[sh]][[key]] <- list(
        key = key, s = s, idx_outcomes = idx_outcomes,
        ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
        ShZ = ShZ, ShZY = ShZY, ShYY = ShYY,
        Nh = nrow(Xs), pzh = ncol(ShZ)
      )
    }
  }
  
  attr(Sh, "R")  <- ncol(Y)
  attr(Sh, "px") <- ncol(X)
  Sh
}



# par: c(theta_u, theta_v_1, ..., theta_v_K[, rho])
# ShPat: nested summaries from peal.get.summary_mv_patterns()
lmm.profile3_mv_patterns <- function(par, ShPat,
                                     reml = TRUE,
                                     estimate_rho = TRUE, rho_fixed = 0,
                                     verbose = FALSE) {
  R  <- attr(ShPat, "R"); if (is.null(R)) stop("Sh summaries missing R.")
  px <- attr(ShPat, "px"); if (is.null(px)) stop("Sh summaries missing px.")
  sites <- names(ShPat)
  K <- length(sites)
  pz <- K + 1  # theta_u + theta_vh per site
  
  # parse parameters
  if (estimate_rho) {
    if (length(par) != (pz + 1)) stop("par must include rho as last element.")
    rho <- par[pz + 1]
  } else {
    rho <- rho_fixed
    if (length(par) != pz) stop("par length must be K+1 when rho is fixed.")
  }
  
  # accumulators
  lpterm1 <- 0
  lpterm2 <- 0
  remlterm <- 0
  bterm1 <- matrix(0, R*px, R*px)
  bterm2 <- rep(0, R*px)
  Nsum <- 0
  
  for (h in seq_len(K)) {
    sh <- sites[h]
    S_h <- ShPat[[sh]]
    if (length(S_h) == 0) next
    
    theta_u  <- par[1]
    theta_vh <- par[1 + h]
    pzh <- S_h[[1]]$pzh
    
    V0 <- diag(c(theta_u, rep(theta_vh, pzh - 1)), pzh)
    V0_inv <- diag(1/diag(V0), pzh, pzh)
    Theta_h_inv <- as.matrix(bdiag(replicate(R, V0_inv, simplify = FALSE)))
    logdet_Theta_h <- R * as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    
    ## ---- Aggregates in FULL (R-augmented) space ----
    SxxR_sum <- matrix(0, R * px, R * px)      # sum of kronecker(Rinve_embed, ShX)
    Sxz_sum  <- matrix(0, R * px, R * pzh)     # sum of kronecker(Rinve_embed, ShXZ)
    Szz_sum  <- matrix(0, R * pzh, R * pzh)    # sum of kronecker(Rinve_embed, ShZ)
    
    sxy_sum  <- rep(0, R * px)                 # vec of stacked columns across outcomes
    szy_sum  <- rep(0, R * pzh)                # vec similarly for ZY
    
    ytildeY_sum       <- 0
    logdet_Rcorr_sum  <- 0
    Nh_site           <- 0
    
    keys <- names(S_h)
    for (key in keys) {
      S   <- S_h[[key]]
      s   <- S$s
      Nh  <- S$Nh
      Nh_site <- Nh_site + Nh
      
      rc_s <- .Rcorr_inv_and_logdet_s(s, rho)
      Rinve_s <- rc_s$Rinve
      logdet_Rcorr_sum <- logdet_Rcorr_sum + Nh * rc_s$logdet_Rcorr
      
      Epi <- selector_from_key(key, R)
      Rinve_embed <- Epi %*% Rinve_s %*% t(Epi)  # R x R
      
      ## Whitened/embedded cross-products (pattern -> full R)
      SxxR_sum <- SxxR_sum + kronecker(Rinve_embed, S$ShX)    # (R*px) x (R*px)
      Sxz_sum  <- Sxz_sum  + kronecker(Rinve_embed, S$ShXZ)   # (R*px) x (R*pzh)
      Szz_sum  <- Szz_sum  + kronecker(Rinve_embed, S$ShZ)    # (R*pzh) x (R*pzh)
      
      ## XY, ZY: build full-R matrices then vec
      idx <- key_to_idx(key)
      ShXY_full <- matrix(0, nrow = nrow(S$ShXY), ncol = R)
      ShZY_full <- matrix(0, nrow = nrow(S$ShZY), ncol = R)
      ShXY_full[, idx] <- S$ShXY %*% Rinve_s
      ShZY_full[, idx] <- S$ShZY %*% Rinve_s
      
      sxy_sum <- sxy_sum + as.vector(ShXY_full)  # length R*px
      szy_sum <- szy_sum + as.vector(ShZY_full)  # length R*pzh
      
      ## y~'y~
      ytildeY_sum <- ytildeY_sum + sum(Rinve_s * S$ShYY)
    }
    
    ## Single per-site Woodbury
    A_h <- Theta_h_inv + Szz_sum
    
    ## log|Gamma_h|: sum log|R_{s}| (by pattern) + log|Theta_h| (once) + log|A_h|
    lpterm1 <- lpterm1 + logdet_Rcorr_sum + logdet_Theta_h +
      as.numeric(determinant(A_h, logarithm = TRUE)$modulus)
    
    ## Reductions
    Ainv_SxzT <- solve(A_h, t(Sxz_sum))                     # (R*pzh) x (R*px)
    bterm1    <- bterm1 + (SxxR_sum - Sxz_sum %*% Ainv_SxzT)
    
    Ainv_szy  <- solve(A_h, szy_sum)                        # (R*pzh) x 1
    bterm2    <- bterm2 + (sxy_sum - as.vector(Sxz_sum %*% Ainv_szy))
    
    lpterm2   <- lpterm2 + (ytildeY_sum - drop(t(szy_sum) %*% Ainv_szy))
    Nsum      <- Nsum + Nh_site
  }
  
  
  # Solve for beta and compute qterm
  L <- chol(bterm1 + 1e-6 * diag(nrow(bterm1)))
  beta_hat <- backsolve(L, forwardsolve(t(L), bterm2))
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * beta_hat) + t(beta_hat) %*% bterm1 %*% beta_hat)
  
  NtotR <- Nsum * R
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
  
  list(
    lp = lp, b = beta_hat, s2 = s2, rho = if (estimate_rho) rho else 0,
    allterms = list(bterm1 = bterm1, bterm2 = bterm2, qterm = qterm,
                    lpterm1 = lpterm1, lpterm2 = lpterm2, remlterm = remlterm)
  )
}





peal.fit.RI_mv_patterns <- function(data, X_cols, Y_cols,
                                    site_col = "site", patient_col = "patient",
                                    weights = NULL, reml = TRUE,
                                    corstr = c("exchangeable","independence"),
                                    mypar.init = NULL, rho_init = 0.1,
                                    hessian = TRUE, verbose = TRUE) {
  
  corstr <- match.arg(corstr)
  R <- length(Y_cols)
  
  # one-shot summaries
  ShPat <- peal.get.summary_mv_patterns(data, X_cols, Y_cols,
                                        site_col = site_col, patient_col = patient_col,
                                        weights = weights)
  sites <- names(ShPat); K <- length(sites)
  pz <- K + 1
  px <- attr(ShPat, "px")
  
  # init params
  if (is.null(mypar.init)) {
    mypar.init <- rep(1, pz)  # theta_u + theta_v[h]
    if (corstr == "exchangeable") mypar.init <- c(mypar.init, rho_init)
    if (verbose) cat("Default mypar.init =", mypar.init, "\n")
  } else {
    if (corstr == "exchangeable" && length(mypar.init) == pz)
      mypar.init <- c(mypar.init, rho_init)
  }
  
  # objective
  fn <- function(parameter) {
    if (corstr == "independence") {
      -lmm.profile3_mv_patterns(parameter[1:pz], ShPat,
                                reml = reml, estimate_rho = FALSE, rho_fixed = 0,
                                verbose = FALSE)$lp
    } else {
      -lmm.profile3_mv_patterns(parameter, ShPat,
                                reml = reml, estimate_rho = TRUE,
                                verbose = FALSE)$lp
    }
  }
  
  # bounds
  lower <- rep(1e-5, length(mypar.init))
  upper <- rep(Inf,  length(mypar.init))
  if (corstr == "exchangeable") {
    # Use s=R for bounds reference; true check is handled inside .Rcorr_inv_and_logdet_s
    rc <- .Rcorr_inv_and_logdet(R, rho = 0)
    lower[length(lower)] <- rc$lower + 1e-8
    upper[length(upper)] <- rc$upper - 1e-8
  }
  
  res <- optim(mypar.init, fn, method = "L-BFGS-B",
               hessian = hessian, lower = lower, upper = upper)
  
  if (verbose) cat(ifelse(res$convergence == 0, "Convergence Reached", "Non-convergence!"),
                   "and",
                   ifelse(all(eigen(res$hessian)$value > 0), "Hessian PD", "Hessian not PD"), "\n",
                   "Function evaluations:", res$counts[1], "\n")
  
  # recompute at optimum
  if (corstr == "independence") {
    prof <- lmm.profile3_mv_patterns(res$par[1:pz], ShPat,
                                     reml = reml, estimate_rho = FALSE, rho_fixed = 0)
    rho_hat <- 0
  } else {
    prof <- lmm.profile3_mv_patterns(res$par, ShPat,
                                     reml = reml, estimate_rho = TRUE)
    rho_hat <- prof$rho
  }
  
  s2 <- prof$s2
  # Var(vec(beta)) = s2 * (bterm1)^{-1}; use Cholesky for stability
  B1 <- prof$allterms$bterm1
  Vbeta <- chol2inv(chol(B1)) * s2
  se    <- sqrt(diag(Vbeta))
  wald  <- prof$b / se
  lb    <- prof$b - 1.96 * se
  ub    <- prof$b + 1.96 * se
  
  Rcorr <- (1 - rho_hat) * diag(R) + rho_hat * matrix(1, R, R)
  Sigma_e <- s2 * Rcorr
  
  list(
    b = prof$b, b.sd = se, wald = wald, lb = lb, ub = ub,
    theta = res$par[1:pz], rho = rho_hat, Sigma_e = Sigma_e,
    s2 = s2, opt = res, res.profile = prof
  )
}












#### EM 

# Build Sigma_e from (s2, rho)
.Sigma_from_rho <- function(R, s2, rho) {
  Rcorr <- (1 - rho) * diag(R) + rho * matrix(1, R, R)
  s2 * Rcorr
}

# E-step: compute completed (full-R) sufficient statistics per site:
#   ShX, ShXZ, ShZ, ShXY_full (px x R), ShZY_full (pzh x R), ShYY_full (R x R), Nh
# Uses current (beta_hat, s2_hat, rho_hat). Random intercepts have mean 0, so
# we condition on fixed-effects mean only (this is exact for RI-only with independent outcomes' RIs).
peal.get.summary_mv_completed_EM <- function(data, X_cols, Y_cols,
                                             site_col = "site", patient_col = "patient",
                                             weights = NULL,
                                             beta_hat, s2_hat, rho_hat) {
  X <- as.matrix(data[, X_cols, drop = FALSE])
  Y <- as.matrix(data[, Y_cols, drop = FALSE])
  R <- ncol(Y)
  px <- ncol(X)
  if (is.null(weights)) weights <- rep(1, nrow(X))
  sites <- as.character(data[[site_col]])
  
  # beta_hat is length R*px stacked by outcomes; reshape to p x R
  if (length(beta_hat) != R * px) stop("beta_hat length must be R*px.")
  B <- matrix(beta_hat, nrow = px, ncol = R)
  
  # robust Z per site (aligned with row order within each site)
  Z_list <- build_Z_list_by_site(data, site_col = site_col, patient_col = patient_col)
  
  # residual covariance
  Sigma <- .Sigma_from_rho(R, s2_hat, rho_hat)
  
  Sh <- list()
  for (sh in unique(sites)) {
    idx <- which(sites == sh)
    if (!length(idx)) next
    Xh <- X[idx,, drop = FALSE]
    Yh <- Y[idx,, drop = FALSE]
    Zh <- Z_list[[sh]]
    wh <- weights[idx]
    
    pzh <- ncol(Zh)
    
    # accumulators (full-R)
    ShX  <- matrix(0, px,  px)
    ShXZ <- matrix(0, px,  pzh)
    ShZ  <- matrix(0, pzh, pzh)
    
    ShXY <- matrix(0, px,  R)
    ShZY <- matrix(0, pzh, R)
    ShYY <- matrix(0, R,   R)
    Nh   <- nrow(Xh)
    
    for (j in seq_len(Nh)) {
      xj <- matrix(Xh[j,, drop = TRUE], nrow = 1)   # 1 x p
      zj <- matrix(Zh[j,, drop = TRUE], nrow = 1)   # 1 x pzh
      yj <- as.numeric(Yh[j,, drop = TRUE])
      wj <- wh[j]
      
      obs <- which(!is.na(yj))
      mis <- which(is.na(yj))
      
      # mean under current beta (RI means are zero)
      eta <- as.numeric(xj %*% B)  # length R
      
      if (length(mis) == 0L) {
        Ey  <- yj
        Eyy <- tcrossprod(yj)       # y y^T
      } else if (length(obs) == 0L) {
        # nothing observed: E[y]=eta; Var[y]=Sigma
        Ey  <- eta
        Eyy <- Sigma + tcrossprod(Ey)
      } else {
        # conditional moments for the missing block
        Soo <- Sigma[obs, obs, drop = FALSE]
        Som <- Sigma[obs, mis, drop = FALSE]
        Smm <- Sigma[mis, mis, drop = FALSE]
        
        # E[y_m | y_o]
        mu_m <- eta[mis] + t(Som) %*% solve(Soo, (yj[obs] - eta[obs]))
        mu_m <- as.numeric(mu_m)
        
        # Var[y_m | y_o]
        V_m  <- Smm - t(Som) %*% solve(Soo, Som)
        
        Ey <- eta
        Ey[obs] <- yj[obs]
        Ey[mis] <- mu_m
        
        # E[yy^T | y_o] = Var(y|y_o) + E[y|y_o]E[y|y_o]^T
        Eyy <- tcrossprod(Ey)
        Eyy[mis, mis] <- Eyy[mis, mis] + V_m
      }
      
      # accumulate with weight
      ShX  <- ShX  + wj * crossprod(xj, xj)
      ShXZ <- ShXZ + wj * crossprod(xj, zj)
      ShZ  <- ShZ  + wj * crossprod(zj, zj)
      
      ShXY <- ShXY + wj * (t(xj) %*% matrix(Ey, nrow = 1))    # (p x 1) %*% (1 x R) -> p x R
      ShZY <- ShZY + wj * (t(zj) %*% matrix(Ey, nrow = 1))    # (pzh x 1) %*% (1 x R) -> pzh x R
      ShYY <- ShYY + wj * Eyy
    }
    
    Sh[[sh]] <- list(ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
                     ShZ = ShZ, ShZY = ShZY, ShYY = ShYY,
                     Nh = Nh, pzh = pzh)
  }
  attr(Sh, "R")  <- R
  attr(Sh, "px") <- px
  return(Sh)
}

# Profile likelihood using full-R summaries (no patterns)
# par = (theta_u, theta_v_1, ..., theta_v_K[, rho])
lmm.profile3_mv_full <- function(par, ShFull,
                                 reml = TRUE,
                                 estimate_rho = TRUE, rho_fixed = 0,
                                 verbose = FALSE) {
  R  <- attr(ShFull, "R")
  px <- attr(ShFull, "px")
  sites <- names(ShFull); K <- length(sites)
  pz <- K + 1
  
  if (estimate_rho) {
    stopifnot(length(par) == pz + 1)
    rho <- par[pz + 1]
  } else {
    stopifnot(length(par) == pz)
    rho <- rho_fixed
  }
  
  rc <- .Rcorr_inv_and_logdet(R, rho)
  Rinve <- rc$Rinve
  logdet_Rcorr <- rc$logdet_Rcorr
  
  lpterm1 <- 0
  lpterm2 <- 0
  remlterm <- 0
  bterm1 <- matrix(0, R * px, R * px)
  bterm2 <- rep(0, R * px)
  Nsum <- 0
  
  for (h in seq_len(K)) {
    sh <- sites[h]
    S  <- ShFull[[sh]]
    if (is.null(S)) next
    
    theta_u  <- par[1]
    theta_vh <- par[1 + h]
    pzh <- S$pzh
    
    V0 <- diag(c(theta_u, rep(theta_vh, pzh - 1)), pzh)
    V0_inv <- diag(1/diag(V0), pzh, pzh)
    Theta_h_inv <- as.matrix(bdiag(replicate(R, V0_inv, simplify = FALSE)))
    logdet_Theta_h <- R * as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    
    # whitened cross-products (full R)
    ShX_tilde  <- kronecker(Rinve, S$ShX)
    ShXZ_tilde <- kronecker(Rinve, S$ShXZ)
    ShZ_tilde  <- kronecker(Rinve, S$ShZ)
    
    ShXY_tilde_vec <- as.vector(S$ShXY %*% Rinve)
    ShZY_tilde_vec <- as.vector(S$ShZY %*% Rinve)
    YtildeY        <- sum(Rinve * S$ShYY)
    
    A_h <- Theta_h_inv + ShZ_tilde
    
    lpterm1 <- lpterm1 + S$Nh * logdet_Rcorr + logdet_Theta_h +
      as.numeric(determinant(A_h, logarithm = TRUE)$modulus)
    
    Ainv_StXZt <- solve(A_h, t(ShXZ_tilde))
    bterm1 <- bterm1 + (ShX_tilde - ShXZ_tilde %*% Ainv_StXZt)
    
    Ainv_StZY <- solve(A_h, ShZY_tilde_vec)
    bterm2 <- bterm2 + (ShXY_tilde_vec - as.vector(ShXZ_tilde %*% Ainv_StZY))
    
    lpterm2 <- lpterm2 + (YtildeY - drop(t(ShZY_tilde_vec) %*% Ainv_StZY))
    
    Nsum <- Nsum + S$Nh
  }
  
  L <- chol(bterm1 + 1e-6 * diag(nrow(bterm1)))
  beta_hat <- backsolve(L, forwardsolve(t(L), bterm2))
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * beta_hat) + t(beta_hat) %*% bterm1 %*% beta_hat)
  
  NtotR <- Nsum * R
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
  
  list(lp = lp, b = beta_hat, s2 = s2, rho = if (estimate_rho) rho else 0,
       allterms = list(bterm1 = bterm1, bterm2 = bterm2, qterm = qterm,
                       lpterm1 = lpterm1, lpterm2 = lpterm2, remlterm = remlterm))
}

# Top-level M-step fit on full summaries (optimize theta and rho)
peal.fit.RI_mv_full_from_EMsumm <- function(ShFull, reml = TRUE,
                                            estimate_rho = TRUE, rho_init = 0.1,
                                            mypar.init = NULL, hessian = TRUE, verbose = TRUE) {
  R  <- attr(ShFull, "R")
  sites <- names(ShFull); K <- length(sites)
  pz <- K + 1
  
  if (is.null(mypar.init)) {
    mypar.init <- rep(1, pz)
    if (estimate_rho) mypar.init <- c(mypar.init, rho_init)
  } else if (estimate_rho && length(mypar.init) == pz) {
    mypar.init <- c(mypar.init, rho_init)
  }
  
  fn <- function(parameter) {
    if (estimate_rho) {
      -lmm.profile3_mv_full(parameter, ShFull, reml = reml, estimate_rho = TRUE)$lp
    } else {
      -lmm.profile3_mv_full(parameter[1:pz], ShFull, reml = reml,
                            estimate_rho = FALSE, rho_fixed = 0)$lp
    }
  }
  
  lower <- rep(1e-5, length(mypar.init)); upper <- rep(Inf, length(mypar.init))
  if (estimate_rho) {
    rc <- .Rcorr_inv_and_logdet(R, rho = 0)
    lower[length(lower)] <- rc$lower + 1e-8
    upper[length(upper)] <- rc$upper - 1e-8
  }
  
  res <- optim(mypar.init, fn, method = "L-BFGS-B", hessian = hessian, lower = lower, upper = upper)
  prof <- if (estimate_rho) {
    lmm.profile3_mv_full(res$par, ShFull, reml = reml, estimate_rho = TRUE)
  } else {
    lmm.profile3_mv_full(res$par[1:pz], ShFull, reml = reml, estimate_rho = FALSE, rho_fixed = 0)
  }
  
  s2 <- prof$s2
  B1 <- prof$allterms$bterm1
  Vbeta <- chol2inv(chol(B1)) * s2
  se    <- sqrt(diag(Vbeta))
  wald  <- prof$b / se
  lb    <- prof$b - 1.96 * se
  ub    <- prof$b + 1.96 * se
  
  rho_hat <- prof$rho
  Rcorr <- (1 - rho_hat) * diag(R) + rho_hat * matrix(1, R, R)
  Sigma_e <- s2 * Rcorr
  
  list(b = prof$b, b.sd = se, wald = wald, lb = lb, ub = ub,
       theta = res$par[1:pz], rho = rho_hat, Sigma_e = Sigma_e,
       s2 = s2, opt = res, res.profile = prof)
}


# EM driver:
#   - Initialize with your pattern-aware fit (M0).
#   - E-step: complete per-site full-R summaries using current params.
#   - M-step: refit using full-R summaries (no patterns).
# Repeat em_iter times (usually 1–2 is enough).
peal.em.fit.RI_mv <- function(data, X_cols, Y_cols,
                              site_col = "site", patient_col = "patient",
                              weights = NULL, reml = TRUE,
                              corstr_init = c("exchangeable","independence"),
                              rho_init = 0.1, mypar.init = NULL,
                              em_iter = 1, verbose = TRUE) {
  
  corstr_init <- match.arg(corstr_init)
  R <- length(Y_cols)
  estimate_rho_all <- (corstr_init == "exchangeable")
  
  init_fit <- peal.fit.RI_mv_patterns(
    data, X_cols, Y_cols,
    site_col = site_col, patient_col = patient_col,
    weights = weights, reml = reml,
    corstr = corstr_init,
    mypar.init = mypar.init,
    rho_init = if (estimate_rho_all) rho_init else 0,
    hessian = TRUE, verbose = verbose
  )
  
  hist <- list(init = init_fit)
  cur  <- init_fit
  
  for (it in seq_len(em_iter)) {
    if (verbose) cat(sprintf("\n[EM] Iteration %d\n", it))
    
    ShFull <- peal.get.summary_mv_completed_EM(
      data, X_cols, Y_cols,
      site_col = site_col, patient_col = patient_col, weights = weights,
      beta_hat = cur$b,
      s2_hat   = cur$s2,
      rho_hat  = if (estimate_rho_all) cur$rho else 0
    )
    
    mfit <- peal.fit.RI_mv_full_from_EMsumm(
      ShFull, reml = reml,
      estimate_rho = estimate_rho_all,
      rho_init     = if (estimate_rho_all) cur$rho else 0,
      mypar.init   = c(cur$theta),
      hessian = TRUE, verbose = verbose
    )
    
    hist[[paste0("em", it)]] <- mfit
    cur <- mfit
  }
  
  cur$history <- hist
  cur$em_iter <- em_iter
  
  # ---- Report theta from the pre-EM (pattern-only) initializer ----
  cur$theta_em <- cur$theta                 # keep EM-updated theta (for reference)
  cur$theta    <- hist$init$theta           # <-- overwrite reported theta
  
  class(cur) <- c("mvpeal_em_fit", class(cur))
  return(cur)
}
