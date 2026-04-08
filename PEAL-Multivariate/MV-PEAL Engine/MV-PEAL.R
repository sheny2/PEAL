library(Matrix)
library(Rcpp)
library(RcppArmadillo)

# ----------------------------------------------------------------
# LOAD C++ BACKEND
# ----------------------------------------------------------------
# Ensure peal_core.cpp is in the working directory
if(file.exists("peal_core.cpp")) {
  sourceCpp("peal_core.cpp")
} else {
  warning("peal_core.cpp not found. Rcpp functionality will be unavailable.")
}

# ==============================================================================
# SECTION 1: SHARED UTILITIES (Correlation & Matrix Helpers)
# ==============================================================================

# General (unstructured) residual correlation utility
# Parameterization: free lower-triangular (unit diag) -> S = L L^T -> Corr = cov2cor(S)
.Corr_from_cholfree <- function(R, par_corr, jitter = 1e-8) {
  # par_corr length = R*(R-1)/2
  L <- diag(1, R)
  idx <- 1
  for (i in 2:R) for (j in 1:(i-1)) {
    L[i, j] <- par_corr[idx]
    idx <- idx + 1
  }
  S <- L %*% t(L)                 
  Corr <- cov2cor(S)              
  
  Cj <- Corr + diag(jitter, R)    
  U  <- chol(Cj)
  Rinve <- chol2inv(U)            
  logdet <- 2 * sum(log(diag(U))) 
  list(Corr = Corr, Rinve = Rinve, logdet_Rcorr = logdet)
}

# Exchangeable residual correlation utility
.Rcorr_inv_and_logdet_s <- function(s, rho, eps = 1e-8) {
  if (s < 1) stop("Pattern of size 0 not allowed.")
  if (s == 1) {
    return(list(Rinve = matrix(1,1,1), logdet_Rcorr = 0, lower = -Inf, upper = Inf))
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

.Rcorr_inv_and_logdet <- function(R, rho, eps = 1e-8) {
  .Rcorr_inv_and_logdet_s(R, rho, eps)
}

# ==============================================================================
# SECTION 2: PURE R IMPLEMENTATION
# ==============================================================================

# Helper: Encode which outcomes are observed
row_pattern_key <- function(yrow) {
  paste0(ifelse(!is.na(yrow), "1", "0"), collapse = "")
}

key_to_idx <- function(key) {
  which(strsplit(key, "", fixed = TRUE)[[1]] == "1")
}

selector_from_key <- function(key, R) {
  idx <- key_to_idx(key)
  s <- length(idx)
  E <- matrix(0, nrow = R, ncol = s)
  if (s > 0) E[cbind(idx, seq_len(s))] <- 1
  E
}

build_Z_list_by_site <- function(data, site_col = "site", patient_col = "patient") {
  sites <- as.character(data[[site_col]])
  pats  <- as.character(data[[patient_col]])
  site_levels <- unique(sites)
  Z_list <- vector("list", length(site_levels))
  names(Z_list) <- site_levels
  
  for (sh in site_levels) {
    idx <- which(sites == sh)
    pat_fac <- factor(pats[idx], levels = unique(pats[idx]))
    Z_pat <- model.matrix(~ pat_fac + 0L)
    Z_h <- cbind(Intercept_site = 1, Z_pat)
    Z_list[[sh]] <- as.matrix(Z_h)
  }
  Z_list
}

# For a subset of outcomes given a FULL Corr
.Corr_subset_inv_logdet <- function(Corr, o, jitter = 1e-8) {
  Coo <- Corr[o, o, drop = FALSE]
  Coo <- Coo + diag(jitter, nrow(Coo))
  Uoo <- chol(Coo)
  inv <- chol2inv(Uoo)
  ld  <- 2 * sum(log(diag(Uoo)))
  list(Rinve = inv, logdet_Rcorr = ld)
}

# R-based Summary Generator
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
    site_idx_all <- which(sites == sh)
    site_valid   <- valid_row[site_idx_all]
    idx_site     <- site_idx_all[site_valid]
    if (!length(idx_site)) next
    
    Xh <- X[idx_site, , drop = FALSE]
    Yh <- Y[idx_site, , drop = FALSE]
    wh <- weights[idx_site]
    
    Zh_site <- Z_list[[sh]]
    Zh      <- Zh_site[site_valid, , drop = FALSE]
    
    keys_h       <- pat_keys[idx_site]
    keys_h_uniq  <- unique(keys_h)
    
    Sh[[sh]] <- list()
    for (key in keys_h_uniq) {
      rows <- which(keys_h == key)
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

# R-based Profile Likelihood
lmm.profile3_mv_patterns <- function(par, ShPat,
                                     reml = TRUE,
                                     corstr = c("exchangeable","independence","unstructured"),
                                     estimate_rho = TRUE, rho_fixed = 0,
                                     verbose = FALSE) {
  corstr <- match.arg(corstr)
  R  <- attr(ShPat, "R")
  px <- attr(ShPat, "px")
  sites <- names(ShPat); K <- length(sites)
  pz <- K + 1
  q  <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  if (corstr == "independence") {
    rho <- 0; Corr <- diag(R)
  } else if (corstr == "exchangeable") {
    if (estimate_rho) rho <- par[pz + 1] else rho <- rho_fixed
    Corr <- (1 - rho) * diag(R) + rho * matrix(1, R, R)
  } else { # unstructured
    par_corr <- par[(pz + 1):(pz + q)]
    Corr_list <- .Corr_from_cholfree(R, par_corr)
    Corr <- Corr_list$Corr
    rho  <- NA_real_
  }
  
  lpterm1 <- 0; lpterm2 <- 0; remlterm <- 0
  bterm1 <- matrix(0, R*px, R*px)
  bterm2 <- rep(0, R*px)
  Nsum_rows <- 0; Nobs_total <- 0
  
  for (h in seq_len(K)) {
    sh <- sites[h]; S_h <- ShPat[[sh]]
    if (length(S_h) == 0) next
    
    theta_u  <- par[1]
    theta_vh <- par[1 + h]
    pzh <- S_h[[1]]$pzh
    
    V0 <- diag(c(theta_u, rep(theta_vh, pzh - 1)), pzh)
    V0_inv <- diag(1/diag(V0), pzh, pzh)
    Theta_h_inv <- as.matrix(bdiag(replicate(R, V0_inv, simplify = FALSE)))
    logdet_Theta_h <- R * as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    
    SxxR_sum <- matrix(0, R*px,  R*px)
    Sxz_sum  <- matrix(0, R*px,  R*pzh)
    Szz_sum  <- matrix(0, R*pzh, R*pzh)
    sxy_sum  <- rep(0, R*px)
    szy_sum  <- rep(0, R*pzh)
    ytildeY_sum <- 0; logdet_Rcorr_sum <- 0; Nh_site <- 0
    
    for (key in names(S_h)) {
      S <- S_h[[key]]
      o <- S$idx_outcomes; s <- S$s
      Nh_site <- Nh_site + S$Nh
      Nobs_total <- Nobs_total + S$Nh * s
      
      if (corstr == "exchangeable") {
        rc_s <- .Rcorr_inv_and_logdet_s(s, if(estimate_rho) rho else rho_fixed)
      } else if (corstr == "independence") {
        rc_s <- list(Rinve = diag(s), logdet_Rcorr = 0)
      } else {
        rc_s <- .Corr_subset_inv_logdet(Corr, o)
      }
      Rinve_s <- rc_s$Rinve
      logdet_Rcorr_sum <- logdet_Rcorr_sum + S$Nh * rc_s$logdet_Rcorr
      
      Epi <- selector_from_key(key, R)
      Rinve_embed <- Epi %*% Rinve_s %*% t(Epi)
      
      SxxR_sum <- SxxR_sum + kronecker(Rinve_embed, S$ShX)
      Sxz_sum  <- Sxz_sum  + kronecker(Rinve_embed, S$ShXZ)
      Szz_sum  <- Szz_sum  + kronecker(Rinve_embed, S$ShZ)
      
      ShXY_full <- matrix(0, nrow = nrow(S$ShXY), ncol = R)
      ShZY_full <- matrix(0, nrow = nrow(S$ShZY), ncol = R)
      ShXY_full[, o] <- S$ShXY %*% Rinve_s
      ShZY_full[, o] <- S$ShZY %*% Rinve_s
      sxy_sum <- sxy_sum + as.vector(ShXY_full)
      szy_sum <- szy_sum + as.vector(ShZY_full)
      
      ytildeY_sum <- ytildeY_sum + sum(Rinve_s * S$ShYY)
    }
    
    A_h <- Theta_h_inv + Szz_sum
    lpterm1 <- lpterm1 + logdet_Rcorr_sum + logdet_Theta_h +
      as.numeric(determinant(A_h, logarithm = TRUE)$modulus)
    
    Ainv_StXZt <- solve(A_h, t(Sxz_sum))
    bterm1 <- bterm1 + (SxxR_sum - Sxz_sum %*% Ainv_StXZt)
    Ainv_szy <- solve(A_h, szy_sum)
    bterm2 <- bterm2 + (sxy_sum - as.vector(Sxz_sum %*% Ainv_szy))
    lpterm2 <- lpterm2 + (ytildeY_sum - drop(t(szy_sum) %*% Ainv_szy))
    Nsum_rows <- Nsum_rows + Nh_site
  }
  
  L <- chol(bterm1 + 1e-6 * diag(nrow(bterm1)))
  beta_hat <- backsolve(L, forwardsolve(t(L), bterm2))
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * beta_hat) + t(beta_hat) %*% bterm1 %*% beta_hat)
  Ntot <- Nobs_total
  
  if (reml) {
    rB <- Matrix::rankMatrix(bterm1, tol = 1e-10)
    df <- Ntot - as.numeric(rB)
    s2 <- qterm / df
    remlterm <- determinant(bterm1 / s2, logarithm = TRUE)$modulus
    lp <- -(lpterm1 + Ntot * log(s2) + qterm / s2 + remlterm) / 2
  } else {
    s2 <- qterm / Ntot
    lp <- -(lpterm1 + (1 + log(qterm * 2 * pi / Ntot)) * Ntot) / 2
  }
  
  list(lp = lp, b = beta_hat, s2 = s2,
       rho = if (corstr == "exchangeable") (if (estimate_rho) rho else rho_fixed) else NA_real_,
       Corr = if (corstr == "unstructured") Corr else NULL)
}

# Main Fit Function
peal.fit.R <- function(data, X_cols, Y_cols,
                       site_col = "site", patient_col = "patient",
                       weights = NULL, reml = TRUE,
                       corstr = c("exchangeable","independence","unstructured"),
                       mypar.init = NULL, rho_init = 0.1,
                       hessian = TRUE, verbose = TRUE) {
  
  corstr <- match.arg(corstr)
  R <- length(Y_cols)
  q <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  ShPat <- peal.get.summary_mv_patterns(data, X_cols, Y_cols, site_col, patient_col, weights)
  sites <- names(ShPat); K <- length(sites)
  pz <- K + 1
  
  if (is.null(mypar.init)) {
    mypar.init <- rep(1, pz)
    if (corstr == "exchangeable") mypar.init <- c(mypar.init, rho_init)
    if (corstr == "unstructured") mypar.init <- c(mypar.init, rep(0, q))
  }
  
  fn <- function(parameter) {
    -lmm.profile3_mv_patterns(
      parameter, ShPat, reml = reml, corstr = corstr,
      estimate_rho = (corstr != "independence"),
      rho_fixed = 0
    )$lp
  }
  
  lower <- rep(1e-5, pz)
  upper <- rep(Inf,  pz)
  if (corstr == "exchangeable") {
    rc <- .Rcorr_inv_and_logdet(R, rho = 0)
    lower <- c(lower, rc$lower + 1e-8); upper <- c(upper, rc$upper - 1e-8)
  } else if (corstr == "unstructured") {
    lower <- c(lower, rep(-Inf, q)); upper <- c(upper, rep( Inf, q))
  }
  
  res <- optim(mypar.init, fn, method = "L-BFGS-B", hessian = hessian, lower = lower, upper = upper)
  
  prof <- lmm.profile3_mv_patterns(res$par, ShPat, reml = reml, corstr = corstr, estimate_rho = (corstr != "independence"))
  
  list(b = prof$b, s2 = prof$s2, theta = res$par[1:pz],
       rho = prof$rho, Corr = prof$Corr, opt = res)
}


# ==============================================================================
# SECTION 3: RCPP / HYBRID IMPLEMENTATION
# ==============================================================================

peal.get.summary_fast_wrapper <- function(data, X_cols, Y_cols, site_col, patient_col) {
  X <- as.matrix(data[, X_cols, drop = FALSE])
  Y <- as.matrix(data[, Y_cols, drop = FALSE])
  sites <- as.character(data[[site_col]])
  patients <- as.character(data[[patient_col]])
  
  # Call C++
  Sh <- get_summaries_fast_cpp(X, Y, sites, patients)
  
  attr(Sh, "R") <- ncol(Y)
  attr(Sh, "px") <- ncol(X)
  return(Sh)
}

peal.profile_fast <- function(par, ShPat, reml = TRUE, corstr = "exchangeable", 
                              estimate_rho = TRUE, rho_fixed = 0) {
  R  <- attr(ShPat, "R"); px <- attr(ShPat, "px")
  K  <- length(ShPat); pz <- K + 1
  q  <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  Corr_full <- diag(R); rho <- 0
  
  if (corstr == "unstructured") {
    par_corr <- par[(pz + 1):(pz + q)]
    Corr_full <- .Corr_from_cholfree(R, par_corr)$Corr
    rho <- NA
  } else if (corstr == "exchangeable" && estimate_rho) {
    rho <- par[pz+1]
  } else if (corstr == "exchangeable") {
    rho <- rho_fixed
  }
  
  res_cpp <- compute_site_stats_cpp(ShPat, par, R, px, corstr, rho_fixed, Corr_full, estimate_rho)
  bterm1 <- res_cpp$bterm1; bterm2 <- res_cpp$bterm2
  
  L <- chol(bterm1 + 1e-6 * diag(nrow(bterm1)))
  beta_hat <- backsolve(L, forwardsolve(t(L), bterm2))
  
  qterm <- as.numeric(res_cpp$lpterm2 - 2 * sum(bterm2 * beta_hat) + t(beta_hat) %*% bterm1 %*% beta_hat)
  Ntot <- res_cpp$Nobs_total
  
  if (reml) {
    rB <- Matrix::rankMatrix(bterm1, tol = 1e-10)
    df <- Ntot - as.numeric(rB)
    s2 <- qterm / df
    remlterm <- determinant(bterm1/s2, logarithm=TRUE)$modulus
    lp <- -(res_cpp$lpterm1 + Ntot * log(s2) + qterm / s2 + remlterm) / 2
  } else {
    s2 <- qterm / Ntot
    lp <- -(res_cpp$lpterm1 + (1 + log(qterm * 2 * pi / Ntot)) * Ntot) / 2
  }
  list(lp = lp, b = beta_hat, s2 = s2, rho = rho, Corr = Corr_full)
}

peal.gradient.hybrid <- function(par, ShPat, reml, corstr, estimate_rho, rho_fixed) {
  prof <- peal.profile_fast(par, ShPat, reml, corstr, estimate_rho, rho_fixed)
  R  <- attr(ShPat, "R"); px <- attr(ShPat, "px"); K  <- length(ShPat); pz <- K + 1
  
  # Analytic Gradient for Variance Components (Thetas)
  grad_theta <- compute_gradient_theta_cpp(
    ShPat, par, R, px, corstr, rho_fixed, prof$Corr, estimate_rho, 
    prof$b, prof$s2
  )
  
  # Numeric Gradient for Correlation Parameters (Rho)
  grad_rho <- numeric(0)
  if (estimate_rho) {
    eps <- 1e-5; base_lp <- prof$lp
    if (corstr == "exchangeable") {
      par_up <- par; idx <- length(par)
      par_up[idx] <- par_up[idx] + eps
      lp_up <- peal.profile_fast(par_up, ShPat, reml, corstr, estimate_rho, rho_fixed)$lp
      grad_rho <- (lp_up - base_lp) / eps
    } else if (corstr == "unstructured") {
      q <- R * (R - 1) / 2; grad_rho <- numeric(q)
      start_idx <- length(par) - q + 1
      for (i in 1:q) {
        par_up <- par; idx <- start_idx + (i - 1)
        par_up[idx] <- par_up[idx] + eps
        lp_up <- peal.profile_fast(par_up, ShPat, reml, corstr, estimate_rho, rho_fixed)$lp
        grad_rho[i] <- (lp_up - base_lp) / eps
      }
    }
  }
  
  grad_full <- c(grad_theta, grad_rho)
  grad_full[is.na(grad_full)] <- 0
  return(-grad_full) # Negative for minimization
}

# Main Rcpp Fit Function
peal.fit.cpp <- function(data, X_cols, Y_cols, 
                         site_col, patient_col,
                         corstr, reml, verbose) {
  
  if(verbose) cat("Calculating local summaries (C++)...\n")
  ShPat <- peal.get.summary_fast_wrapper(data, X_cols, Y_cols, site_col, patient_col)
  
  R <- length(Y_cols); K <- length(ShPat); pz <- K + 1
  q <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  init_par <- c(1.0, rep(1.0, K)) 
  if (corstr == "exchangeable") init_par <- c(init_par, 0.1)
  if (corstr == "unstructured") init_par <- c(init_par, rep(0, q))
  
  fn <- function(par) {
    if(any(par[1:pz] <= 0)) return(1e10)
    if(corstr == "exchangeable") {
      if(par[pz+1] <= -0.5 || par[pz+1] >= 0.99) return(1e10)
    }
    prof <- peal.profile_fast(par, ShPat, reml=reml, corstr=corstr)
    return(-prof$lp)
  }
  
  gr <- function(par) {
    if(any(par[1:pz] <= 0)) return(rep(0, length(par))) 
    peal.gradient.hybrid(par, ShPat, reml, corstr, estimate_rho = (corstr != "independence"), rho_fixed = 0)
  }
  
  if(verbose) cat("Optimizing with Analytic Gradients...\n")
  opt <- optim(init_par, fn, gr = gr, method = "L-BFGS-B", 
               lower = c(rep(1e-4, pz), if(corstr=="exchangeable") -0.5 else rep(-Inf, q)),
               upper = c(rep(Inf, pz),  if(corstr=="exchangeable") 0.99 else rep(Inf, q)),
               hessian = TRUE)
  
  final_prof <- peal.profile_fast(opt$par, ShPat, reml=reml, corstr=corstr)
  res_cpp <- compute_site_stats_cpp(ShPat, opt$par, attr(ShPat,"R"), attr(ShPat,"px"), 
                                    corstr, 0, final_prof$Corr, TRUE)
  XtVinvX <- res_cpp$bterm1
  VarBeta <- solve(XtVinvX) * final_prof$s2
  se_beta <- sqrt(diag(VarBeta))
  
  list(b = final_prof$b, b.sd = se_beta, s2 = final_prof$s2,
       theta = opt$par[1:pz], rho = final_prof$rho, Corr = final_prof$Corr,
       opt = opt)
}

# ==============================================================================
# SECTION 4: MAIN INTERFACE
# ==============================================================================

#' MV-PEAL Fit Function
#' 
#' Main entry point for fitting the MV-PEAL model.
peal.fit <- function(data, X_cols, Y_cols, 
                     site_col = "site", patient_col = "patient",
                     weights = NULL, 
                     corstr = c("exchangeable", "independence", "unstructured"), 
                     reml = TRUE, 
                     verbose = TRUE,
                     use_rcpp = TRUE) {
  
  corstr <- match.arg(corstr)
  t0 <- Sys.time()
  
  if (use_rcpp) {
    if (!is.null(weights)) {
      warning("Weights are not supported in the C++ implementation. Ignoring weights or switch use_rcpp=FALSE.")
    }
    # Check if compiled
    if (!exists("get_summaries_fast_cpp")) {
      stop("C++ backend not loaded. Ensure 'peal_core.cpp' is present and compiled, or set use_rcpp=FALSE.")
    }
    
    result <- peal.fit.cpp(data, X_cols, Y_cols, site_col, patient_col, corstr, reml, verbose)
    result$method <- "Rcpp (Analytic Gradient)"
    
  } else {
    # R Implementation
    if (verbose) cat("Using Pure R implementation...\n")
    result <- peal.fit.R(data, X_cols, Y_cols, site_col, patient_col, weights, reml, corstr, verbose = verbose)
    
    # Standardize output names slightly to match C++ return
    result$b.sd <- sqrt(diag(solve(result$opt$hessian))) # Approx if not calc in profile
    # Note: peal.fit.R actually re-calculates Vbeta if we used the full logic, 
    # but here we simply return the object structure.
    result$method <- "R (no Rcpp)"
  }
  
  
  if (verbose) cat(ifelse(all(result$opt$convergence == 0, eigen(result$opt$hessian)$value > 0),
                          "Convergence Reached", "Non-convergence!"), 'and',
                   ifelse(all(eigen(result$opt$hessian)$value > 0),
                          "Hessian PD", "Hessian not PD"), '\n',
                   "The number of function evaluations used is", result$opt$counts[1], '\n')
  
  result$time <- Sys.time() - t0
  return(result)
}