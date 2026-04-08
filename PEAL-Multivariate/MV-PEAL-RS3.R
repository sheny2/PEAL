library(Matrix); library(Rcpp); library(RcppArmadillo)
sourceCpp("peal_rs_core2.cpp")

.Corr_from_cholfree <- function(R, par_corr) {
  L <- diag(1, R); idx <- 1
  for (i in 2:R) for (j in 1:(i-1)) { L[i, j] <- par_corr[idx]; idx <- idx + 1 }
  return(cov2cor(L %*% t(L)))
}

peal.fit.rs <- function(data, X_cols, Y_cols, slope_col, site_col="site", patient_col="patient", 
                        corstr=c("exchangeable", "unstructured"), reml=TRUE, verbose=TRUE) {
  corstr <- match.arg(corstr)
  t0 <- Sys.time()
  X <- as.matrix(data[, X_cols]); Y <- as.matrix(data[, Y_cols]); X_slp <- as.numeric(data[[slope_col]])
  ShPat <- get_summaries_rs_cpp(X, Y, X_slp, as.character(data[[site_col]]), as.character(data[[patient_col]]))
  
  R <- ncol(Y); px <- ncol(X); K <- length(ShPat)
  pz <- 1 + K + K 
  q  <- if(corstr == "unstructured") R*(R-1)/2 else 1
  attr(ShPat, "R") <- R; attr(ShPat, "px") <- px
  
  # OPTIMIZATION 1: Initialize variance components on the LOG scale
  init_par <- c(log(1.0), rep(log(0.5), K), rep(log(0.2), K), rep(0.1, q))
  
  # Bounds are now completely open (-Inf) for the variance components!
  lower_bnds <- c(rep(-Inf, pz), if(corstr=="unstructured") rep(-Inf, q) else -0.5)
  upper_bnds <- c(rep(Inf, pz), if(corstr=="unstructured") rep(Inf, q) else 0.99)
  
  fn <- function(par) {
    # Transform variance components back to linear scale for C++ math
    par_orig <- par
    par_orig[1:pz] <- exp(par[1:pz]) 
    
    Cf <- if(corstr=="unstructured") .Corr_from_cholfree(R, par_orig[(pz+1):(pz+q)]) else (1-par_orig[pz+1])*diag(R) + par_orig[pz+1]*matrix(1,R,R)
    res <- compute_site_stats_rs_cpp(ShPat, par_orig, R, px, Cf)
    
    # OPTIMIZATION 3: Safeguard the Cholesky decomposition
    lp <- tryCatch({
      L <- chol(res$b1 + 1e-7*diag(nrow(res$b1)))
      b <- backsolve(L, forwardsolve(t(L), res$b2))
      qt <- as.numeric(res$lp2 - 2*sum(res$b2*b) + t(b)%*%res$b1%*%b)
      s2 <- qt / (if(reml) res$Ntot - px*R else res$Ntot)
      log_det_b1_s2 <- if(reml) 2 * sum(log(diag(L))) - nrow(res$b1) * log(s2) else 0
      
      -(res$lp1 + res$Ntot*log(s2) + qt/s2 + log_det_b1_s2) / 2
    }, error = function(e) {
      return(-1e10) # Return massive penalty if matrix is computationally singular
    })
    
    return(-lp)
  }
  
  gr <- function(par) {
    par_orig <- par
    par_orig[1:pz] <- exp(par[1:pz])
    
    Cf <- if(corstr=="unstructured") .Corr_from_cholfree(R, par_orig[(pz+1):(pz+q)]) else (1-par_orig[pz+1])*diag(R) + par_orig[pz+1]*matrix(1,R,R)
    res <- compute_site_stats_rs_cpp(ShPat, par_orig, R, px, Cf)
    
    L <- chol(res$b1 + 1e-7*diag(nrow(res$b1)))
    b <- backsolve(L, forwardsolve(t(L), res$b2))
    qt <- as.numeric(res$lp2 - 2*sum(res$b2*b) + t(b)%*%res$b1%*%b)
    s2 <- qt / (if(reml) res$Ntot - px*R else res$Ntot)
    
    # Get the raw gradient from C++
    g_theta_raw <- compute_gradient_rs_cpp(ShPat, par_orig, R, px, Cf, b, s2)
    
    # OPTIMIZATION 1 (Chain Rule): Adjust gradient for the log-transformation
    g_theta_log <- g_theta_raw[1:pz] * par_orig[1:pz]
    
    g_corr <- numeric(q)
    eps <- 1e-4 # Slightly larger epsilon for central difference stability
    
    # OPTIMIZATION 2: Central finite differences for the correlation parameters
    for(i in 1:q) { 
      idx <- pz + i
      p_up <- par; p_up[idx] <- min(p_up[idx] + eps, upper_bnds[idx])
      p_down <- par; p_down[idx] <- max(p_down[idx] - eps, lower_bnds[idx])
      
      # Central difference: (f(x+h) - f(x-h)) / 2h
      g_corr[i] <- (fn(p_up) - fn(p_down)) / (p_up[idx] - p_down[idx])
    }
    
    return(c(-g_theta_log, g_corr)) # Ensure we negate both since we are minimizing
  }
  
  opt <- nlminb(start = init_par, 
                objective = fn, 
                gradient = gr, 
                lower = lower_bnds, 
                upper = upper_bnds,
                control = list(eval.max = 2000, iter.max = 1500, rel.tol = 1e-8))
  
  # Re-transform the final optimized variance components back to original scale
  final_par_orig <- opt$par
  final_par_orig[1:pz] <- exp(opt$par[1:pz])
  
  final_Cf <- if(corstr=="unstructured") .Corr_from_cholfree(R, final_par_orig[(pz+1):(pz+q)]) else (1-final_par_orig[pz+1])*diag(R) + final_par_orig[pz+1]*matrix(1,R,R)
  stats <- compute_site_stats_rs_cpp(ShPat, final_par_orig, R, px, final_Cf)
  beta_f <- solve(stats$b1) %*% stats$b2
  qt_f <- as.numeric(stats$lp2 - 2*sum(stats$b2*beta_f) + t(beta_f)%*%stats$b1%*%beta_f)
  s2_f <- qt_f / (if(reml) stats$Ntot - px*R else stats$Ntot) 
  
  return(list(
    beta = beta_f,
    sigma2_residual = s2_f,
    sigma_site_int = sqrt(final_par_orig[1] * s2_f), 
    sigma_pat_int = sqrt(final_par_orig[2:(K+1)] * s2_f),
    sigma_pat_slope = sqrt(final_par_orig[(K+2):pz] * s2_f),
    Corr = final_Cf,
    opt = opt, 
    time = Sys.time() - t0
  ))
}