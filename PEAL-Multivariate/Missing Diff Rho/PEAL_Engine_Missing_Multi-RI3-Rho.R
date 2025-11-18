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


# -------------------------
# General (unstructured) residual correlation utilities
# Parameterization: free lower-triangular (unit diag) -> S = L L^T -> Corr = cov2cor(S)
# par_corr length = R*(R-1)/2 ordered as (2,1), (3,1),(3,2), ..., (R,1),...,(R,R-1)
# Returns: list(Corr, Rinve, logdet_Rcorr)
# -------------------------
.Corr_from_cholfree <- function(R, par_corr, jitter = 1e-8) {
  stopifnot(length(par_corr) == R*(R-1)/2)
  L <- diag(1, R)
  idx <- 1
  for (i in 2:R) for (j in 1:(i-1)) { L[i, j] <- par_corr[idx]; idx <- idx + 1 }
  S <- L %*% t(L)                 # SPD by construction
  Corr <- cov2cor(S)              # unit-diagonal correlation
  Cj <- Corr + diag(jitter, R)    # stabilize
  U  <- chol(Cj)
  Rinve <- chol2inv(U)            # Corr^{-1}
  logdet <- 2 * sum(log(diag(U))) # log|Corr|
  list(Corr = Corr, Rinve = Rinve, logdet_Rcorr = logdet)
}

# For a subset of outcomes (pattern of size s) given a FULL Corr (R x R)
# returns inverse and logdet of Corr[o,o].
.Corr_subset_inv_logdet <- function(Corr, o, jitter = 1e-8) {
  Coo <- Corr[o, o, drop = FALSE]
  Coo <- Coo + diag(jitter, nrow(Coo))
  Uoo <- chol(Coo)
  inv <- chol2inv(Uoo)
  ld  <- 2 * sum(log(diag(Uoo)))
  list(Rinve = inv, logdet_Rcorr = ld)
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



lmm.profile3_mv_patterns <- function(par, ShPat,
                                     reml = TRUE,
                                     corstr = c("exchangeable","independence","unstructured"),
                                     estimate_rho = TRUE, rho_fixed = 0,
                                     verbose = FALSE) {
  corstr <- match.arg(corstr)
  R  <- attr(ShPat, "R");  if (is.null(R)) stop("Sh summaries missing R.")
  px <- attr(ShPat, "px"); if (is.null(px)) stop("Sh summaries missing px.")
  sites <- names(ShPat); K <- length(sites)
  pz <- K + 1  # theta_u + theta_vh per site
  q  <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  # ----- parse parameters / correlation structure -----
  if (corstr == "independence") {
    stopifnot(length(par) == pz)
    rho <- 0; Corr <- diag(R)
  } else if (corstr == "exchangeable") {
    if (estimate_rho) {
      stopifnot(length(par) == pz + 1)
      rho <- par[pz + 1]; Corr <- (1 - rho) * diag(R) + rho * matrix(1, R, R)
    } else {
      stopifnot(length(par) == pz)
      rho <- rho_fixed;   Corr <- (1 - rho) * diag(R) + rho * matrix(1, R, R)
    }
  } else { # unstructured
    stopifnot(estimate_rho)
    stopifnot(length(par) == pz + q)
    par_corr <- par[(pz + 1):(pz + q)]
    Corr <- .Corr_from_cholfree(R, par_corr)$Corr
    rho <- NA_real_
  }
  
  # ----- accumulators -----
  lpterm1 <- 0; lpterm2 <- 0; remlterm <- 0
  bterm1 <- matrix(0, R*px, R*px)
  bterm2 <- rep(0, R*px)
  Nsum <- 0
  
  for (h in seq_len(K)) {
    sh <- sites[h]; S_h <- ShPat[[sh]]
    if (length(S_h) == 0) next
    
    theta_u  <- par[1]
    theta_vh <- par[1 + h]
    
    # These don't depend on pattern (only site h)
    pzh <- S_h[[1]]$pzh
    V0 <- diag(c(theta_u, rep(theta_vh, pzh - 1)), pzh)
    V0_inv <- diag(1/diag(V0), pzh, pzh)
    Theta_h_inv <- as.matrix(bdiag(replicate(R, V0_inv, simplify = FALSE)))
    logdet_Theta_h <- R * as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    
    # Per-site sums in FULL R-embedded space
    SxxR_sum <- matrix(0, R*px,  R*px)
    Sxz_sum  <- matrix(0, R*px,  R*pzh)
    Szz_sum  <- matrix(0, R*pzh, R*pzh)
    sxy_sum  <- rep(0, R*px)
    szy_sum  <- rep(0, R*pzh)
    ytildeY_sum      <- 0
    logdet_Rcorr_sum <- 0
    Nh_site          <- 0
    
    for (key in names(S_h)) {
      S <- S_h[[key]]
      o <- S$idx_outcomes
      s <- S$s
      Nh_site <- Nh_site + S$Nh
      
      # inverse + logdet for the pattern's observed sub-corr
      if (corstr == "exchangeable") {
        rc_s <- .Rcorr_inv_and_logdet_s(s, if (isTRUE(estimate_rho)) rho else rho_fixed)
      } else if (corstr == "independence") {
        rc_s <- list(Rinve = diag(s), logdet_Rcorr = 0)
      } else {
        rc_s <- .Corr_subset_inv_logdet(Corr, o)
      }
      Rinve_s <- rc_s$Rinve
      logdet_Rcorr_sum <- logdet_Rcorr_sum + S$Nh * rc_s$logdet_Rcorr
      
      # embed to R x R via E_pi
      Epi <- selector_from_key(key, R)
      Rinve_embed <- Epi %*% Rinve_s %*% t(Epi)
      
      # contributions (pattern -> full R)
      SxxR_sum <- SxxR_sum + kronecker(Rinve_embed, S$ShX)
      Sxz_sum  <- Sxz_sum  + kronecker(Rinve_embed, S$ShXZ)
      Szz_sum  <- Szz_sum  + kronecker(Rinve_embed, S$ShZ)
      
      # XY, ZY: build full-R then vec (prewhiten with Rinve_s on columns o)
      ShXY_full <- matrix(0, nrow = nrow(S$ShXY), ncol = R)
      ShZY_full <- matrix(0, nrow = nrow(S$ShZY), ncol = R)
      ShXY_full[, o] <- S$ShXY %*% Rinve_s
      ShZY_full[, o] <- S$ShZY %*% Rinve_s
      sxy_sum <- sxy_sum + as.vector(ShXY_full)
      szy_sum <- szy_sum + as.vector(ShZY_full)
      
      # y~'y~ = tr(Rinve_s * Y_o'Y_o)
      ytildeY_sum <- ytildeY_sum + sum(Rinve_s * S$ShYY)
    }
    
    # Woodbury at the site level
    A_h <- Theta_h_inv + Szz_sum
    
    lpterm1 <- lpterm1 + logdet_Rcorr_sum + logdet_Theta_h +
      as.numeric(determinant(A_h, logarithm = TRUE)$modulus)
    
    Ainv_StXZt <- solve(A_h, t(Sxz_sum))
    bterm1 <- bterm1 + (SxxR_sum - Sxz_sum %*% Ainv_StXZt)
    
    Ainv_szy <- solve(A_h, szy_sum)
    bterm2 <- bterm2 + (sxy_sum - as.vector(Sxz_sum %*% Ainv_szy))
    
    lpterm2 <- lpterm2 + (ytildeY_sum - drop(t(szy_sum) %*% Ainv_szy))
    Nsum    <- Nsum + Nh_site
  }
  
  # beta, qterm, sigma2 profiling
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
    lp = lp, b = beta_hat, s2 = s2,
    rho = if (corstr == "exchangeable") (if (estimate_rho) rho else rho_fixed) else NA_real_,
    Corr = if (corstr == "unstructured") Corr else NULL,
    allterms = list(bterm1 = bterm1, bterm2 = bterm2, qterm = qterm,
                    lpterm1 = lpterm1, lpterm2 = lpterm2, remlterm = remlterm)
  )
}


peal.fit.RI_mv_patterns <- function(data, X_cols, Y_cols,
                                    site_col = "site", patient_col = "patient",
                                    weights = NULL, reml = TRUE,
                                    corstr = c("exchangeable","independence","unstructured"),
                                    mypar.init = NULL, rho_init = 0.1,
                                    hessian = TRUE, verbose = TRUE) {
  
  corstr <- match.arg(corstr)
  R <- length(Y_cols); q <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  # one-shot summaries
  ShPat <- peal.get.summary_mv_patterns(data, X_cols, Y_cols,
                                        site_col = site_col, patient_col = patient_col,
                                        weights = weights)
  sites <- names(ShPat); K <- length(sites)
  pz <- K + 1
  
  # init params
  if (is.null(mypar.init)) {
    mypar.init <- rep(1, pz)  # theta_u + theta_v[h]
    if (corstr == "exchangeable") mypar.init <- c(mypar.init, rho_init)
    if (corstr == "unstructured") mypar.init <- c(mypar.init, rep(0, q))  # near I
    if (verbose) cat("Default mypar.init =", mypar.init, "\n")
  } else if (length(mypar.init) == pz) {
    if (corstr == "exchangeable") mypar.init <- c(mypar.init, rho_init)
    if (corstr == "unstructured") mypar.init <- c(mypar.init, rep(0, q))
  }
  
  # objective
  fn <- function(parameter) {
    -lmm.profile3_mv_patterns(parameter, ShPat,
                              reml = reml, corstr = corstr,
                              estimate_rho = (corstr != "independence"),
                              rho_fixed = if (corstr == "exchangeable") 0 else 0,
                              verbose = FALSE)$lp
  }
  
  # bounds
  lower <- rep(1e-5, pz); upper <- rep(Inf, pz)
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
  
  if (verbose) cat(ifelse(res$convergence == 0, "Convergence Reached", "Non-convergence!"),
                   "and",
                   ifelse(all(eigen(res$hessian)$value > 0), "Hessian PD", "Hessian not PD"), "\n",
                   "Function evaluations:", res$counts[1], "\n")
  
  prof <- lmm.profile3_mv_patterns(res$par, ShPat,
                                   reml = reml, corstr = corstr,
                                   estimate_rho = (corstr != "independence"))
  
  s2 <- prof$s2
  B1 <- prof$allterms$bterm1
  Vbeta <- chol2inv(chol(B1)) * s2
  se    <- sqrt(diag(Vbeta))
  wald  <- prof$b / se
  lb    <- prof$b - 1.96 * se
  ub    <- prof$b + 1.96 * se
  
  Corr_hat <- if (corstr == "unstructured") prof$Corr else
    if (corstr == "exchangeable") ((1 - prof$rho) * diag(R) + prof$rho * matrix(1, R, R)) else diag(R)
  Sigma_e <- s2 * Corr_hat
  
  list(
    b = prof$b, b.sd = se, wald = wald, lb = lb, ub = ub,
    theta = res$par[1:pz],
    rho = if (corstr == "exchangeable") prof$rho else NULL,
    Corr = Corr_hat, Sigma_e = Sigma_e,
    s2 = s2, opt = res, res.profile = prof, ShPat = ShPat, corstr = corstr
  )
}



















#### EM 
# Profile likelihood using full-R summaries (no patterns)
# par = (theta_u, theta_v_1, ..., theta_v_K[, rho])
lmm.profile3_mv_full <- function(par, ShFull,
                                 reml = TRUE,
                                 corstr = c("exchangeable","independence","unstructured"),
                                 estimate_rho = TRUE, rho_fixed = 0,
                                 verbose = FALSE) {
  corstr <- match.arg(corstr)
  R  <- attr(ShFull, "R")
  px <- attr(ShFull, "px")
  sites <- names(ShFull); K <- length(sites)
  pz <- K + 1
  q  <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  # parse correlation
  if (corstr == "independence") {
    stopifnot(length(par) == pz); Corr <- diag(R); rho <- 0
  } else if (corstr == "exchangeable") {
    if (estimate_rho) { stopifnot(length(par) == pz + 1); rho <- par[pz + 1] }
    else              { stopifnot(length(par) == pz);      rho <- rho_fixed }
    Corr <- (1 - rho) * diag(R) + rho * matrix(1, R, R)
  } else {
    stopifnot(estimate_rho, length(par) == pz + q)
    par_corr <- par[(pz + 1):(pz + q)]
    Corr <- .Corr_from_cholfree(R, par_corr)$Corr
    rho <- NA_real_
  }
  
  # common inverse/logdet for FULL Corr (used R-wise on ShX/Z)
  Cj <- Corr + diag(1e-8, R)
  U  <- chol(Cj)
  Rinve <- chol2inv(U)
  logdet_Rcorr <- 2 * sum(log(diag(U)))
  
  lpterm1 <- 0; lpterm2 <- 0; remlterm <- 0
  bterm1 <- matrix(0, R * px, R * px)
  bterm2 <- rep(0, R * px)
  Nsum <- 0
  
  for (h in seq_len(K)) {
    sh <- sites[h]; S <- ShFull[[sh]]
    if (is.null(S)) next
    
    theta_u  <- par[1]
    theta_vh <- par[1 + h]
    pzh <- S$pzh
    
    V0 <- diag(c(theta_u, rep(theta_vh, pzh - 1)), pzh)
    V0_inv <- diag(1/diag(V0), pzh, pzh)
    Theta_h_inv <- as.matrix(bdiag(replicate(R, V0_inv, simplify = FALSE)))
    logdet_Theta_h <- R * as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    
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
  
  list(lp = lp, b = beta_hat, s2 = s2,
       rho = if (corstr == "exchangeable") (if (estimate_rho) rho else rho_fixed) else NA_real_,
       Corr = if (corstr == "unstructured") Corr else NULL,
       allterms = list(bterm1 = bterm1, bterm2 = bterm2, qterm = qterm,
                       lpterm1 = lpterm1, lpterm2 = lpterm2, remlterm = remlterm))
}

peal.fit.RI_mv_full_from_EMsumm <- function(ShFull, reml = TRUE,
                                            corstr = c("exchangeable","independence","unstructured"),
                                            estimate_rho = TRUE, rho_init = 0.1,
                                            mypar.init = NULL, hessian = TRUE, verbose = TRUE) {
  corstr <- match.arg(corstr)
  R  <- attr(ShFull, "R")
  sites <- names(ShFull); K <- length(sites)
  pz <- K + 1
  q  <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  if (is.null(mypar.init)) {
    mypar.init <- rep(1, pz)
    if (corstr == "exchangeable" && estimate_rho) mypar.init <- c(mypar.init, rho_init)
    if (corstr == "unstructured") mypar.init <- c(mypar.init, rep(0, q))
  } else if (length(mypar.init) == pz) {
    if (corstr == "exchangeable" && estimate_rho) mypar.init <- c(mypar.init, rho_init)
    if (corstr == "unstructured") mypar.init <- c(mypar.init, rep(0, q))
  }
  
  fn <- function(parameter) {
    -lmm.profile3_mv_full(parameter, ShFull, reml = reml,
                          corstr = corstr,
                          estimate_rho = (corstr != "independence"),
                          rho_fixed = if (corstr == "exchangeable") 0 else 0)$lp
  }
  
  lower <- rep(1e-5, pz); upper <- rep(Inf, pz)
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
  
  prof <- lmm.profile3_mv_full(res$par, ShFull, reml = reml,
                               corstr = corstr, estimate_rho = (corstr != "independence"))
  
  s2 <- prof$s2
  B1 <- prof$allterms$bterm1
  Vbeta <- chol2inv(chol(B1)) * s2
  se    <- sqrt(diag(Vbeta))
  wald  <- prof$b / se
  lb    <- prof$b - 1.96 * se
  ub    <- prof$b + 1.96 * se
  
  Corr_hat <- if (corstr == "unstructured") prof$Corr else
    if (corstr == "exchangeable") ((1 - prof$rho) * diag(R) + prof$rho * matrix(1, R, R)) else diag(R)
  Sigma_e <- s2 * Corr_hat
  
  list(b = prof$b, b.sd = se, wald = wald, lb = lb, ub = ub,
       theta = res$par[1:pz],
       rho = if (corstr == "exchangeable") prof$rho else NULL,
       Corr = Corr_hat, Sigma_e = Sigma_e,
       s2 = s2, opt = res, res.profile = prof)
}



peal.complete_from_patterns <- function(ShPat, beta_hat, s2_hat,
                                        Corr_hat = NULL, rho_hat = NULL, ridge = 1e-10) {
  R  <- attr(ShPat, "R");  if (is.null(R)) stop("ShPat missing attr R")
  px <- attr(ShPat, "px"); if (is.null(px)) stop("ShPat missing attr px")
  sites <- names(ShPat)
  B <- matrix(beta_hat, nrow = px, ncol = R)
  
  # Accept either Corr_hat or rho_hat for backward compatibility
  if (is.null(Corr_hat)) {
    if (is.null(rho_hat)) stop("Provide Corr_hat (preferred) or rho_hat.")
    Corr_hat <- (1 - rho_hat) * diag(R) + rho_hat * matrix(1, R, R)
  }
  Sigma <- s2_hat * Corr_hat
  
  safe_solve <- function(A, b = NULL, add = ridge) {
    if (is.null(b)) solve(A + add * diag(nrow(A)))
    else            solve(A + add * diag(nrow(A)), b)
  }
  
  ShFull <- vector("list", length = length(sites)); names(ShFull) <- sites
  
  for (sh in sites) {
    S_h <- ShPat[[sh]]
    if (length(S_h) == 0L) { ShFull[[sh]] <- NULL; next }
    
    ShX  <- matrix(0, px, px)
    ShXZ <- NULL; pzh <- NULL
    ShZ  <- NULL
    ShXY_full <- matrix(0, px, R)
    ShZY_full <- NULL
    ShYY_full <- matrix(0, R, R)
    Nh_site   <- 0L
    
    for (key in names(S_h)) {
      S <- S_h[[key]]
      o <- S$idx_outcomes; s <- length(o)
      m <- setdiff(seq_len(R), o)
      Nh_site <- Nh_site + S$Nh
      if (is.null(ShXZ)) { ShXZ <- matrix(0, px,  S$pzh); pzh <- S$pzh }
      if (is.null(ShZ))  { ShZ  <- matrix(0, S$pzh, S$pzh) }
      if (is.null(ShZY_full)) ShZY_full <- matrix(0, pzh, R)
      
      ShX  <- ShX  + S$ShX
      ShXZ <- ShXZ + S$ShXZ
      ShZ  <- ShZ  + S$ShZ
      
      XYo  <- S$ShXY
      ZYo  <- S$ShZY
      YYoo <- S$ShYY
      B_o  <- B[, o, drop = FALSE]
      
      XEo <- S$ShX  %*% B_o
      ZEo <- t(S$ShXZ) %*% B_o
      
      Xres_o <- XYo - XEo
      Zres_o <- ZYo - ZEo
      
      ShXY_full[, o] <- ShXY_full[, o] + XYo
      ShZY_full[, o] <- ShZY_full[, o] + ZYo
      ShYY_full[o, o] <- ShYY_full[o, o] + YYoo
      if (length(m) == 0L) next
      
      # conditional components from general Sigma
      Soo <- Sigma[o, o, drop = FALSE]
      Som <- Sigma[o, m, drop = FALSE]
      Smm <- Sigma[m, m, drop = FALSE]
      Sooinv <- safe_solve(Soo)
      W <- Sooinv %*% Som
      
      B_m  <- B[, m, drop = FALSE]
      XE_m <- S$ShX  %*% B_m
      ZE_m <- t(S$ShXZ) %*% B_m
      
      ShXY_full[, m] <- ShXY_full[, m] + XE_m + Xres_o %*% W
      ShZY_full[, m] <- ShZY_full[, m] + ZE_m + Zres_o %*% W
      
      YoEta_m <- t(XYo) %*% B_m
      YoResT  <- YYoo - t(XYo) %*% B_o
      ShYY_full[o, m] <- ShYY_full[o, m] + YoEta_m + YoResT %*% W
      ShYY_full[m, o] <- t(ShYY_full[o, m])
      
      Vm <- Smm - t(Som) %*% Sooinv %*% Som
      Eta_mEta_mT <- t(B_m) %*% S$ShX %*% B_m
      Res_oRes_oT <- YYoo - t(XYo) %*% B_o - t(B_o) %*% XYo + t(B_o) %*% S$ShX %*% B_o
      Prop_m <- t(W) %*% Res_oRes_oT %*% W
      Cross  <- t(B_m) %*% Xres_o %*% W + t(W) %*% t(Xres_o) %*% B_m
      
      ShYY_full[m, m] <- ShYY_full[m, m] +
        S$Nh * Vm + Eta_mEta_mT + Prop_m + Cross
    }
    
    ShFull[[sh]] <- list(
      ShX = ShX, ShXZ = ShXZ, ShZ = ShZ,
      ShXY = ShXY_full, ShZY = ShZY_full, ShYY = ShYY_full,
      Nh = Nh_site, pzh = pzh
    )
  }
  
  attr(ShFull, "R")  <- R
  attr(ShFull, "px") <- px
  ShFull
}


# EM driver:
peal.em.fit.RI_mv <- function(data, X_cols, Y_cols,
                              site_col = "site", patient_col = "patient",
                              weights = NULL, reml = TRUE,
                              corstr_init = c("exchangeable","independence","unstructured"),
                              rho_init = 0.1, mypar.init = NULL,
                              em_iter = 1, verbose = TRUE) {
  
  corstr_init <- match.arg(corstr_init)
  R <- length(Y_cols)
  estimate_corr <- (corstr_init != "independence")
  
  # One-shot collection (sites send exactly once)
  ShPat <- peal.get.summary_mv_patterns(
    data, X_cols, Y_cols,
    site_col = site_col, patient_col = patient_col, weights = weights
  )
  
  # Initial pattern-aware fit (M0)
  init_fit <- peal.fit.RI_mv_patterns(
    data, X_cols, Y_cols,
    site_col = site_col, patient_col = patient_col, weights = weights,
    reml = reml, corstr = corstr_init,
    mypar.init = mypar.init,
    rho_init = if (corstr_init == "exchangeable") rho_init else 0,
    hessian = TRUE, verbose = verbose
  )
  
  hist <- list(init = init_fit)
  cur  <- init_fit
  
  for (it in seq_len(em_iter)) {
    if (verbose) cat(sprintf("\n[EM] Iteration %d (one-shot, no new site stats)\n", it))
    
    # E-step (lead): complete full-R summaries from pattern summaries
    ShFull <- peal.complete_from_patterns(
      ShPat,
      beta_hat = cur$b,
      s2_hat   = cur$s2,
      Corr_hat = cur$Corr,           # <-- use full Corr if available
      rho_hat  = cur$rho             # fallback for exchangeable
    )
    
    # M-step: optimize theta and correlation params on completed summaries
    mfit <- peal.fit.RI_mv_full_from_EMsumm(
      ShFull, reml = reml,
      corstr = corstr_init,
      estimate_rho = estimate_corr,
      rho_init     = if (!is.null(cur$rho)) cur$rho else 0,
      mypar.init   = c(cur$theta),   # use θ from first step (your request)
      hessian = TRUE, verbose = verbose
    )
    
    hist[[paste0("em", it)]] <- mfit
    cur <- mfit
  }
  
  # report θ from pre-EM initializer; keep EM-updated θ separately
  cur$history <- hist
  cur$em_iter <- em_iter
  cur$theta_em <- cur$theta
  cur$theta    <- hist$init$theta
  class(cur) <- c("mvpeal_em_fit", class(cur))
  cur
}
