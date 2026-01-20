library(Matrix)
library(Rcpp)
library(RcppArmadillo)

# ----------------------------------------------------------------
# LOAD C++ BACKEND
# ----------------------------------------------------------------
sourceCpp("peal_core.cpp") 

# ----------------------------------------------------------------
# 1. Data Preparation
# ----------------------------------------------------------------

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

# ----------------------------------------------------------------
# 2. Optimized Profile Likelihood Function (Wrapper)
# ----------------------------------------------------------------

.Corr_from_cholfree <- function(R, par_corr) {
  L <- diag(1, R); idx <- 1
  for (i in 2:R) for (j in 1:(i-1)) { L[i, j] <- par_corr[idx]; idx <- idx + 1 }
  S <- L %*% t(L); cov2cor(S)
}

peal.profile_fast <- function(par, ShPat, reml = TRUE, corstr = "exchangeable", 
                              estimate_rho = TRUE, rho_fixed = 0) {
  
  R  <- attr(ShPat, "R")
  px <- attr(ShPat, "px")
  # Get K from list length
  K  <- length(ShPat)
  pz <- K + 1
  q  <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  Corr_full <- diag(R)
  rho <- 0
  
  if (corstr == "unstructured") {
    par_corr <- par[(pz + 1):(pz + q)]
    Corr_full <- .Corr_from_cholfree(R, par_corr)
    rho <- NA
  } else if (corstr == "exchangeable" && estimate_rho) {
    rho <- par[pz+1]
  } else if (corstr == "exchangeable") {
    rho <- rho_fixed
  }
  
  # Call C++
  res_cpp <- compute_site_stats_cpp(ShPat, par, R, px, corstr, rho_fixed, Corr_full, estimate_rho)
  
  bterm1 <- res_cpp$bterm1
  bterm2 <- res_cpp$bterm2
  
  L <- chol(bterm1 + 1e-6 * diag(nrow(bterm1)))
  beta_hat <- backsolve(L, forwardsolve(t(L), bterm2))
  
  qterm <- as.numeric(
    res_cpp$lpterm2 - 2 * sum(bterm2 * beta_hat) + t(beta_hat) %*% bterm1 %*% beta_hat
  )
  
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

# ----------------------------------------------------------------
# 3. Main Fit Function
# ----------------------------------------------------------------

peal.fit.fast <- function(data, X_cols, Y_cols, 
                          site_col = "site", patient_col = "patient",
                          corstr = "exchangeable", reml = TRUE, verbose = TRUE) {
  
  # 1. Get Summaries (Now calls C++)
  t0 <- Sys.time()
  if(verbose) cat("Calculating local summaries (C++)...\n")
  ShPat <- peal.get.summary_fast_wrapper(data, X_cols, Y_cols, site_col, patient_col)
  
  R <- length(Y_cols)
  K <- length(ShPat)
  pz <- K + 1
  q <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  # 2. Init
  init_par <- c(1.0, rep(1.0, K)) 
  if (corstr == "exchangeable") init_par <- c(init_par, 0.1)
  if (corstr == "unstructured") init_par <- c(init_par, rep(0, q))
  
  # 3. Optimize
  fn <- function(par) {
    if(any(par[1:pz] <= 0)) return(1e10)
    if(corstr == "exchangeable") {
      if(par[pz+1] <= -0.5 || par[pz+1] >= 0.99) return(1e10)
    }
    prof <- peal.profile_fast(par, ShPat, reml=reml, corstr=corstr)
    return(-prof$lp)
  }
  
  if(verbose) cat("Optimizing...\n")
  opt <- optim(init_par, fn, method = "L-BFGS-B", 
               lower = c(rep(1e-4, pz), if(corstr=="exchangeable") -0.5 else rep(-Inf, q)),
               upper = c(rep(Inf, pz),  if(corstr=="exchangeable") 0.99 else rep(Inf, q)),
               hessian = TRUE)
  
  # 4. Final Estimates
  final_prof <- peal.profile_fast(opt$par, ShPat, reml=reml, corstr=corstr)
  
  res_cpp <- compute_site_stats_cpp(ShPat, opt$par, attr(ShPat,"R"), attr(ShPat,"px"), 
                                    corstr, 0, final_prof$Corr, TRUE)
  
  XtVinvX <- res_cpp$bterm1
  VarBeta <- solve(XtVinvX) * final_prof$s2
  se_beta <- sqrt(diag(VarBeta))
  
  list(
    beta = final_prof$b,
    se_beta = se_beta,
    sigma2 = final_prof$s2,
    theta = opt$par[1:pz],
    rho = final_prof$rho,
    Corr = final_prof$Corr,
    opt = opt,
    time = Sys.time() - t0
  )
}



peal.gradient.hybrid <- function(par, ShPat, reml, corstr, estimate_rho, rho_fixed) {
  
  # 1. Get current Profile stats (needed for s2 and beta to compute analytic gradient)
  prof <- peal.profile_fast(par, ShPat, reml, corstr, estimate_rho, rho_fixed)
  
  R  <- attr(ShPat, "R")
  px <- attr(ShPat, "px")
  K  <- length(ShPat)
  pz <- K + 1 # Number of variance parameters (Thetas)
  
  # 2. Analytic Gradient for Variance Components (Thetas) [Length = pz]
  # Returns d(LogLik) / d(Theta)
  grad_theta <- compute_gradient_theta_cpp(
    ShPat, par, R, px, corstr, rho_fixed, prof$Corr, estimate_rho, 
    prof$b, prof$s2
  )
  
  # 3. Numeric Gradient for Correlation Parameters (Rho)
  grad_rho <- numeric(0)
  
  if (estimate_rho) {
    eps <- 1e-5
    base_lp <- prof$lp
    
    if (corstr == "exchangeable") {
      # Single correlation parameter at the end
      par_up <- par
      idx <- length(par)
      par_up[idx] <- par_up[idx] + eps
      
      lp_up <- peal.profile_fast(par_up, ShPat, reml, corstr, estimate_rho, rho_fixed)$lp
      grad_rho <- (lp_up - base_lp) / eps
      
    } else if (corstr == "unstructured") {
      # Multiple correlation parameters at the end
      q <- R * (R - 1) / 2
      grad_rho <- numeric(q)
      
      # The correlation params are the last q elements of par
      start_idx <- length(par) - q + 1
      
      for (i in 1:q) {
        par_up <- par
        idx <- start_idx + (i - 1)
        par_up[idx] <- par_up[idx] + eps
        
        lp_up <- peal.profile_fast(par_up, ShPat, reml, corstr, estimate_rho, rho_fixed)$lp
        grad_rho[i] <- (lp_up - base_lp) / eps
      }
    }
  }
  
  # 4. Combine and Return
  grad_full <- c(grad_theta, grad_rho)
  
  # Safety: If params are out of bounds, finite diff might produce NA
  grad_full[is.na(grad_full)] <- 0
  
  # Verify length matches 'par' to prevent optim error
  if(length(grad_full) != length(par)) {
    stop(paste("Gradient length mismatch. Expected", length(par), "got", length(grad_full)))
  }
  
  # Return NEGATIVE because optim minimizes Negative Log Likelihood
  return(-grad_full) 
}


# ----------------------------------------------------------------
# 5. UPDATED FIT FUNCTION
# ----------------------------------------------------------------
peal.fit.fast <- function(data, X_cols, Y_cols, 
                          site_col = "site", patient_col = "patient",
                          corstr = "exchangeable", reml = TRUE, verbose = TRUE) {
  
  t0 <- Sys.time()
  if(verbose) cat("Calculating local summaries (C++)...\n")
  ShPat <- peal.get.summary_fast_wrapper(data, X_cols, Y_cols, site_col, patient_col)
  
  R <- length(Y_cols); K <- length(ShPat); pz <- K + 1
  q <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  init_par <- c(1.0, rep(1.0, K)) 
  if (corstr == "exchangeable") init_par <- c(init_par, 0.1)
  if (corstr == "unstructured") init_par <- c(init_par, rep(0, q))
  
  # OBJECTIVE FUNCTION
  fn <- function(par) {
    if(any(par[1:pz] <= 0)) return(1e10)
    if(corstr == "exchangeable") {
      if(par[pz+1] <= -0.5 || par[pz+1] >= 0.99) return(1e10)
    }
    prof <- peal.profile_fast(par, ShPat, reml=reml, corstr=corstr)
    return(-prof$lp)
  }
  
  # GRADIENT FUNCTION
  gr <- function(par) {
    # Simple bound check before calling expensive gradient
    # Variances (first pz params) must be positive
    if(any(par[1:pz] <= 0)) return(rep(0, length(par))) 
    
    peal.gradient.hybrid(
      par, ShPat, reml, 
      corstr, 
      estimate_rho = (corstr != "independence"), 
      rho_fixed = 0
    )
  }
  
  if(verbose) cat("Optimizing with Analytic Gradients...\n")
  
  # PASS 'gr' TO OPTIM
  opt <- optim(init_par, fn, 
               gr = gr,
               method = "L-BFGS-B", 
               lower = c(rep(1e-4, pz), if(corstr=="exchangeable") -0.5 else rep(-Inf, q)),
               upper = c(rep(Inf, pz),  if(corstr=="exchangeable") 0.99 else rep(Inf, q)),
               hessian = TRUE) # Hessian is cheap at the very end
  
  final_prof <- peal.profile_fast(opt$par, ShPat, reml=reml, corstr=corstr)
  res_cpp <- compute_site_stats_cpp(ShPat, opt$par, attr(ShPat,"R"), attr(ShPat,"px"), 
                                    corstr, 0, final_prof$Corr, TRUE)
  XtVinvX <- res_cpp$bterm1
  VarBeta <- solve(XtVinvX) * final_prof$s2
  se_beta <- sqrt(diag(VarBeta))
  
  list(b = final_prof$b, b.sd = se_beta, s2 = final_prof$s2,
       theta = opt$par[1:pz], rho = final_prof$rho, Corr = final_prof$Corr,
       opt = opt, time = Sys.time() - t0)
}