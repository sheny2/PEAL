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
  
  init_par <- c(1.0, rep(0.5, K), rep(0.2, K), rep(0.1, q))
  
  # Set bounds up here so the gradient function can see them
  lower_bnds <- c(rep(1e-4, pz), if(corstr=="unstructured") rep(-Inf, q) else -0.5)
  upper_bnds <- c(rep(Inf, pz), if(corstr=="unstructured") rep(Inf, q) else 0.99)
  
  fn <- function(par) {
    # REMOVED: if(any(par[1:pz] <= 0)) return(1e10)
    # nlminb will natively respect the 1e-4 lower bound, so we don't need the manual cliff.
    
    Cf <- if(corstr=="unstructured") .Corr_from_cholfree(R, par[(pz+1):(pz+q)]) else (1-par[pz+1])*diag(R) + par[pz+1]*matrix(1,R,R)
    res <- compute_site_stats_rs_cpp(ShPat, par, R, px, Cf)
    
    L <- chol(res$b1 + 1e-7*diag(nrow(res$b1)))
    b <- backsolve(L, forwardsolve(t(L), res$b2))
    qt <- as.numeric(res$lp2 - 2*sum(res$b2*b) + t(b)%*%res$b1%*%b)
    s2 <- qt / (if(reml) res$Ntot - px*R else res$Ntot)
    lp <- -(res$lp1 + res$Ntot*log(s2) + qt/s2 + (if(reml) determinant(res$b1/s2)$modulus else 0)) / 2
    return(-lp)
  }
  
  gr <- function(par) {
    Cf <- if(corstr=="unstructured") .Corr_from_cholfree(R, par[(pz+1):(pz+q)]) else (1-par[pz+1])*diag(R) + par[pz+1]*matrix(1,R,R)
    res <- compute_site_stats_rs_cpp(ShPat, par, R, px, Cf)
    L <- chol(res$b1 + 1e-7*diag(nrow(res$b1))); b <- backsolve(L, forwardsolve(t(L), res$b2))
    qt <- as.numeric(res$lp2 - 2*sum(res$b2*b) + t(b)%*%res$b1%*%b)
    s2 <- qt / (if(reml) res$Ntot - px*R else res$Ntot)
    
    g_theta <- compute_gradient_rs_cpp(ShPat, par, R, px, Cf, b, s2)
    
    # Boundary-safe finite differences for the correlation parameters
    g_corr <- numeric(q)
    eps <- 1e-5
    base_lp <- -fn(par)
    
    for(i in 1:q) { 
      idx <- pz + i
      if (par[idx] + eps <= upper_bnds[idx]) {
        # Forward difference (normal)
        p_up <- par; p_up[idx] <- p_up[idx] + eps
        g_corr[i] <- (-fn(p_up) - base_lp) / eps
      } else {
        # Backward difference (if we are hugging the upper bound)
        p_down <- par; p_down[idx] <- p_down[idx] - eps
        g_corr[i] <- (base_lp - (-fn(p_down))) / eps 
      }
    }
    
    return(-c(g_theta[1:pz], g_corr))
  }
  
  # Added specific control arguments for noisy/complex matrices
  opt <- nlminb(start = init_par, 
                objective = fn, 
                gradient = gr, 
                lower = lower_bnds, 
                upper = upper_bnds,
                control = list(
                  eval.max = 2000, 
                  iter.max = 1500, 
                  rel.tol = 1e-6,  # Loosen rel.tol slightly to ignore deep-valley numerical noise
                  x.tol = 1e-6
                ))
  
  final_Cf <- if(corstr=="unstructured") .Corr_from_cholfree(R, opt$par[(pz+1):(pz+q)]) else (1-opt$par[pz+1])*diag(R) + opt$par[pz+1]*matrix(1,R,R)
  stats <- compute_site_stats_rs_cpp(ShPat, opt$par, R, px, final_Cf)
  beta_f <- solve(stats$b1) %*% stats$b2
  qt_f <- as.numeric(stats$lp2 - 2*sum(stats$b2*beta_f) + t(beta_f)%*%stats$b1%*%beta_f)
  s2_f <- qt_f / (if(reml) stats$Ntot - px*R else stats$Ntot) 
  
  return(list(
    beta = beta_f,
    sigma2_residual = s2_f,
    sigma_site_int = sqrt(opt$par[1] * s2_f), 
    sigma_pat_int = sqrt(opt$par[2:(K+1)] * s2_f),
    sigma_pat_slope = sqrt(opt$par[(K+2):pz] * s2_f),
    Corr = final_Cf,
    opt = opt, 
    time = Sys.time() - t0
  ))
}