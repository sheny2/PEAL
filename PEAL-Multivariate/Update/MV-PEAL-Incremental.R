library(Matrix)
library(Rcpp)
library(RcppArmadillo)

# LOAD C++
sourceCpp("peal_core_incremental.cpp") 

# ----------------------------------------------------------------
# UTILS
# ----------------------------------------------------------------
.Corr_from_cholfree <- function(R, par_corr) {
  L <- diag(1, R); idx <- 1
  for (i in 2:R) for (j in 1:(i-1)) { L[i, j] <- par_corr[idx]; idx <- idx + 1 }
  S <- L %*% t(L); cov2cor(S)
}

# ----------------------------------------------------------------
# 1. GENERATE REGISTRY
#    Creates a list mapping Site -> [PatientIDs]
# ----------------------------------------------------------------
peal.make.registry <- function(data, site_col, patient_col) {
  sites <- as.character(data[[site_col]])
  pats  <- as.character(data[[patient_col]])
  
  reg <- split(pats, sites)
  reg <- lapply(reg, unique)
  return(reg)
}

# ----------------------------------------------------------------
# 2. FIT FUNCTION (Initial)
# ----------------------------------------------------------------
peal.fit.fast <- function(data, X_cols, Y_cols, site_col, patient_col, 
                          corstr = "exchangeable", reml = TRUE, verbose = TRUE) {
  
  # 1. Build Registry & Summaries
  registry <- peal.make.registry(data, site_col, patient_col)
  
  X <- as.matrix(data[, X_cols, drop=FALSE])
  Y <- as.matrix(data[, Y_cols, drop=FALSE])
  sites <- as.character(data[[site_col]])
  pats <- as.character(data[[patient_col]])
  
  if(verbose) cat("Generating Initial Summaries...\n")
  # Use the controlled generator even for initial fit to ensure consistency
  ShPat <- get_summaries_controlled_cpp(X, Y, sites, pats, registry)
  
  attr(ShPat, "R") <- ncol(Y)
  attr(ShPat, "px") <- ncol(X)
  
  # 2. Optimize
  res <- peal.optimize.internal(ShPat, corstr, reml, verbose)
  
  res$ShPat <- ShPat
  res$registry <- registry
  res$settings <- list(reml=reml, corstr=corstr)
  return(res)
}

# ----------------------------------------------------------------
# 3. UPDATE FUNCTION (Incremental)
# ----------------------------------------------------------------
peal.update <- function(old_model, new_data, X_cols, Y_cols, site_col, patient_col, verbose=TRUE) {
  
  t0 <- Sys.time()
  
  # 1. Update Registry
  if(verbose) cat("Updating Patient Registry...\n")
  new_registry <- peal.make.registry(new_data, site_col, patient_col)
  old_registry <- old_model$registry
  full_registry <- old_registry
  
  # Merge registries
  all_sites <- unique(c(names(old_registry), names(new_registry)))
  
  # Calculate target sizes (pzh = 1 + n_patients)
  target_sizes <- list()
  
  for(site in all_sites) {
    p_old <- if(!is.null(old_registry[[site]])) old_registry[[site]] else character(0)
    p_new <- if(!is.null(new_registry[[site]])) new_registry[[site]] else character(0)
    
    p_full <- unique(c(p_old, p_new))
    full_registry[[site]] <- p_full
    target_sizes[[site]] <- length(p_full) + 1
  }
  
  # 2. Pad Old Summaries
  if(verbose) cat("Padding Old Summaries...\n")
  ShPat_Old <- old_model$ShPat
  ShPat_Old_Padded <- pad_summaries_cpp(ShPat_Old, target_sizes)
  
  # 3. Generate New Summaries (Controlled by Full Registry)
  if(verbose) cat("Generating New Summaries...\n")
  X <- as.matrix(new_data[, X_cols, drop=FALSE])
  Y <- as.matrix(new_data[, Y_cols, drop=FALSE])
  sites <- as.character(new_data[[site_col]])
  pats <- as.character(new_data[[patient_col]])
  
  ShPat_New <- get_summaries_controlled_cpp(X, Y, sites, pats, full_registry)
  
  # 4. Merge Lists
  if(verbose) cat("Merging Summaries...\n")
  ShPat_Combined <- ShPat_Old_Padded
  
  for(site in names(ShPat_New)) {
    if(is.null(ShPat_Combined[[site]])) {
      ShPat_Combined[[site]] <- ShPat_New[[site]]
    } else {
      # Append the list of pattern matrices
      # Since we rely on the likelihood sum loop, simple appending works
      # provided dimensions match (which padding ensured)
      ShPat_Combined[[site]] <- c(ShPat_Combined[[site]], ShPat_New[[site]])
    }
  }
  
  attr(ShPat_Combined, "R") <- attr(old_model$ShPat, "R")
  attr(ShPat_Combined, "px") <- attr(old_model$ShPat, "px")
  
  # 5. Optimize (Warm Start)
  res <- peal.optimize.internal(ShPat_Combined, old_model$settings$corstr, old_model$settings$reml, verbose, old_model$theta)
  
  res$ShPat <- ShPat_Combined
  res$registry <- full_registry
  res$settings <- old_model$settings
  res$time <- Sys.time() - t0
  return(res)
}

# ----------------------------------------------------------------
# 4. INTERNAL OPTIMIZER (Shared)
# ----------------------------------------------------------------
peal.optimize.internal <- function(ShPat, corstr, reml, verbose, init_theta=NULL) {
  
  R <- attr(ShPat, "R"); K <- length(ShPat); pz <- K + 1
  q <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  # Init Par
  if(is.null(init_theta)) {
    init_par <- c(1.0, rep(1.0, K)) 
    if (corstr == "exchangeable") init_par <- c(init_par, 0.1)
    if (corstr == "unstructured") init_par <- c(init_par, rep(0, q))
  } else {
    # If K grew, we need to handle new theta_v
    old_pz <- length(init_theta) - (if(corstr=="exchangeable") 1 else if(corstr=="unstructured") q else 0)
    old_K <- old_pz - 1
    
    if(K > old_K) {
      theta_u <- init_theta[1]
      theta_v <- init_theta[2:old_pz]
      theta_v_new <- rep(mean(theta_v), K - old_K)
      tail_par <- init_theta[(old_pz+1):length(init_theta)]
      init_par <- c(theta_u, theta_v, theta_v_new, tail_par)
    } else {
      init_par <- init_theta
    }
  }
  
  fn <- function(par) {
    if(any(par[1:pz] <= 0)) return(1e10)
    if(corstr == "exchangeable") if(par[pz+1] <= -0.5 || par[pz+1] >= 0.99) return(1e10)
    prof <- peal.profile_fast(par, ShPat, reml, corstr)
    return(-prof$lp)
  }
  
  gr <- function(par) {
    if(any(par[1:pz] <= 0)) return(rep(0, length(par))) 
    peal.gradient.hybrid(par, ShPat, reml, corstr, (corstr!="independence"), 0)
  }
  
  if(verbose) cat("Optimizing...\n")
  opt <- optim(init_par, fn, gr=gr, method = "L-BFGS-B", 
               lower = c(rep(1e-4, pz), if(corstr=="exchangeable") -0.5 else rep(-Inf, q)),
               upper = c(rep(Inf, pz),  if(corstr=="exchangeable") 0.99 else rep(Inf, q)),
               hessian = TRUE)
  
  final_prof <- peal.profile_fast(opt$par, ShPat, reml, corstr)
  
  # SE Calc
  res_cpp <- compute_site_stats_cpp(ShPat, opt$par, R, attr(ShPat,"px"), corstr, 0, final_prof$Corr, TRUE)
  XtVinvX <- res_cpp$bterm1
  VarBeta <- solve(XtVinvX) * final_prof$s2
  se_beta <- sqrt(diag(VarBeta))
  
  list(b = final_prof$b, b.sd = se_beta, s2 = final_prof$s2,
       theta = opt$par[1:pz], rho = final_prof$rho, Corr = final_prof$Corr, opt = opt)
}

# --- Profile/Gradient Wrappers ---
peal.profile_fast <- function(par, ShPat, reml, corstr, estimate_rho=TRUE, rho_fixed=0) {
  R <- attr(ShPat, "R"); px <- attr(ShPat, "px"); K <- length(ShPat); pz <- K + 1
  q <- if (corstr == "unstructured") R*(R-1)/2 else 0
  
  Corr_full <- diag(R); rho <- 0
  if (corstr == "unstructured") {
    par_corr <- par[(pz + 1):(pz + q)]
    Corr_full <- .Corr_from_cholfree(R, par_corr)
    rho <- NA
  } else if (corstr == "exchangeable" && estimate_rho) rho <- par[pz+1]
  
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
  R <- attr(ShPat, "R"); px <- attr(ShPat, "px")
  grad_theta <- compute_gradient_theta_cpp(ShPat, par, R, px, corstr, rho_fixed, prof$Corr, estimate_rho, prof$b, prof$s2)
  
  grad_rho <- numeric(0)
  if (estimate_rho) {
    eps <- 1e-5; base_lp <- prof$lp
    if (corstr == "exchangeable") {
      par_up <- par; idx <- length(par); par_up[idx] <- par_up[idx] + eps
      lp_up <- peal.profile_fast(par_up, ShPat, reml, corstr, estimate_rho, rho_fixed)$lp
      grad_rho <- (lp_up - base_lp) / eps
    } else if (corstr == "unstructured") {
      q <- R * (R - 1) / 2; grad_rho <- numeric(q); start_idx <- length(par) - q + 1
      for (i in 1:q) {
        par_up <- par; idx <- start_idx + (i - 1); par_up[idx] <- par_up[idx] + eps
        lp_up <- peal.profile_fast(par_up, ShPat, reml, corstr, estimate_rho, rho_fixed)$lp
        grad_rho[i] <- (lp_up - base_lp) / eps
      }
    }
  }
  grad_full <- c(grad_theta, grad_rho)
  grad_full[is.na(grad_full)] <- 0
  return(-grad_full) 
}