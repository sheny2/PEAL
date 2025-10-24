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
  
  # build Z per site (robust)
  Z_list <- build_Z_list_by_site(data, site_col = site_col, patient_col = patient_col)
  
  Sh <- list()
  site_levels <- unique(sites)
  
  # Precompute pattern key per row
  pat_keys <- apply(Y, 1, row_pattern_key)
  print(table(pat_keys))
  valid_row <- (pat_keys != paste0(rep("0", R), collapse = ""))  # at least one outcome observed
  
  for (sh in site_levels) {
    idx_site <- which(sites == sh & valid_row)
    if (!length(idx_site)) next
    Xh <- X[idx_site, , drop = FALSE]
    Yh <- Y[idx_site, , drop = FALSE]
    wh <- weights[idx_site]
    Zh <- Z_list[[sh]]           # aligned to idx_site order by construction in build_Z_list_by_site
    
    # derive pattern keys in this site
    keys_h <- pat_keys[idx_site]
    keys_h_uniq <- unique(keys_h)
    
    Sh[[sh]] <- list()
    for (key in keys_h_uniq) {
      rows <- which(keys_h == key)
      if (!length(rows)) next
      
      # rows within the site's submatrix
      Xs <- Xh[rows, , drop = FALSE]
      Zs <- Zh[rows, , drop = FALSE]
      Ys_full <- Yh[rows, , drop = FALSE]
      
      idx_outcomes <- key_to_idx(key)
      s <- length(idx_outcomes)
      Ys <- Ys_full[, idx_outcomes, drop = FALSE]
      
      ws <- wh[rows]
      Xw <- Xs * ws
      Zw <- Zs * ws
      Yw <- Ys * ws
      
      ShX  <- crossprod(Xw, Xs)         # p_x x p_x
      ShXZ <- crossprod(Xw, Zs)         # p_x x p_zh
      ShXY <- crossprod(Xw, Ys)         # p_x x s
      
      ShZ  <- crossprod(Zw, Zs)         # p_zh x p_zh
      ShZY <- crossprod(Zw, Ys)         # p_zh x s
      
      ShYY <- crossprod(Yw, Ys)         # s x s
      Nhp  <- nrow(Xs)
      
      Sh[[sh]][[key]] <- list(
        key = key, s = s, idx_outcomes = idx_outcomes,
        ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
        ShZ = ShZ, ShZY = ShZY, ShYY = ShYY,
        Nh = Nhp,
        pzh = ncol(ShZ)   # convenience
      )
    }
  }
  attr(Sh, "R") <- length(Y_cols)
  attr(Sh, "px") <- ncol(X)
  return(Sh)
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
    
    # get pzh from any pattern in this site
    pzh <- S_h[[1]]$pzh
    # per-site V0 and its inverse (per outcome block)
    V0 <- diag(c(theta_u, rep(theta_vh, pzh - 1)), pzh)
    V0_inv <- diag(1/diag(V0), pzh, pzh)
    Theta_h_inv <- as.matrix(bdiag(replicate(R, V0_inv, simplify = FALSE))) # (R*pzh)x(R*pzh)
    logdet_Theta_h <- R * as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    
    # loop over patterns in this site
    keys <- names(S_h)
    for (key in keys) {
      S <- S_h[[key]]
      s <- S$s
      Nh <- S$Nh
      idx_out <- S$idx_outcomes
      Nsum <- Nsum + Nh
      
      # R_{rho,s}^{-1} and its logdet
      rc_s <- .Rcorr_inv_and_logdet_s(s, rho)
      Rinve_s <- rc_s$Rinve
      logdet_Rcorr_s <- rc_s$logdet_Rcorr
      
      # Embed to R x R
      Epi <- selector_from_key(key, R)
      Rinve_embed <- Epi %*% Rinve_s %*% t(Epi)   # R x R
      
      # Whitened cross-products
      # Kronecker with embedded Rinve on outcome axis
      ShX_tilde  <- kronecker(Rinve_embed, S$ShX)      # (R*px)x(R*px)
      ShXZ_tilde <- kronecker(Rinve_embed, S$ShXZ)     # (R*px)x(R*pzh)
      ShZ_tilde  <- kronecker(Rinve_embed, S$ShZ)      # (R*pzh)x(R*pzh)
      
      # Whiten X'Y and Z'Y into full R stacks:
      # Build full (px x R) matrix with zeros, place ShXY %*% Rinve_s into outcome cols idx_out
      ShXY_full <- matrix(0, nrow = nrow(S$ShXY), ncol = R)
      ShZY_full <- matrix(0, nrow = nrow(S$ShZY), ncol = R)
      ShXY_full[, idx_out] <- S$ShXY %*% Rinve_s
      ShZY_full[, idx_out] <- S$ShZY %*% Rinve_s
      ShXY_tilde_vec <- as.vector(ShXY_full)         # length R*px
      ShZY_tilde_vec <- as.vector(ShZY_full)         # length R*pzh
      
      # Y~'Y~ = tr(Rinv_s * ShYY)
      YtildeY <- sum(Rinve_s * S$ShYY)
      
      # A_{h,pi} and terms
      A_hpi <- Theta_h_inv + ShZ_tilde
      
      # log|Gamma_{h,pi}| piece
      lpterm1 <- lpterm1 + Nh * logdet_Rcorr_s + logdet_Theta_h +
        as.numeric(determinant(A_hpi, logarithm = TRUE)$modulus)
      
      # Woodbury reductions
      Ainv_StXZt <- solve(A_hpi, t(ShXZ_tilde))   # (R*pzh) x (R*px)
      bterm1 <- bterm1 + (ShX_tilde - ShXZ_tilde %*% Ainv_StXZt)
      
      Ainv_StZY  <- solve(A_hpi, ShZY_tilde_vec)  # (R*pzh) x 1
      bterm2 <- bterm2 + (ShXY_tilde_vec - as.vector(ShXZ_tilde %*% Ainv_StZY))
      
      lpterm2 <- lpterm2 + (YtildeY - drop(t(ShZY_tilde_vec) %*% Ainv_StZY))
    }
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


