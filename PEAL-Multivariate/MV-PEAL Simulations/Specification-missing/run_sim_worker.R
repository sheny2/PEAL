# ==============================================================================
# File: run_sim_worker.R
# Usage: Rscript run_sim_worker.R <SLURM_ARRAY_TASK_ID>
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Task ID must be provided.")
task_id <- as.integer(args[1])

library(parallel)
library(doParallel)
library(foreach)
library(parallelly)

# ----------------------------------------------------------------
# 1. PARALLEL & BATCH CONFIGURATION
# ----------------------------------------------------------------
# Detect available cores allocated by SLURM
n_cores <- parallelly::availableCores()
cat(sprintf("Task %d: Detected %d available cores. Initializing cluster...\n", task_id, n_cores))

# 1. Create Cluster Explicitly
cl <- makeCluster(n_cores)

# 2. Initialize Workers (Load Libs + Source Code ONCE)
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  library(MASS)
  library(Matrix)
  library(Rcpp)
  library(RcppArmadillo)
  
  # Ensure file exists on worker node before sourcing
  if(!file.exists("MV-PEAL.R")) stop("MV-PEAL.R not found on worker node.")
  source("MV-PEAL.R")
})

# 3. Register Backend
registerDoParallel(cl)

# Batch Logic
N_scenarios     <- 8
N_reps_per_scen <- 1000
Total_Sims      <- N_scenarios * N_reps_per_scen 

Array_Size      <- 80
Sims_Per_Job    <- Total_Sims / Array_Size 

start_idx <- (task_id - 1) * Sims_Per_Job + 1
end_idx   <- task_id * Sims_Per_Job

cat(sprintf("Task %d running global simulations %d to %d in parallel\n", task_id, start_idx, end_idx))

# --- SCENARIO CONFIGURATION: VARYING MISSINGNESS ---
# We define 8 levels of gamma_0 intercepts.
start_vec <- c(-2.8, -2.4, -2.0)
step <- 0.6

gamma_grid <- list()
for(i in 1:8) {
  gamma_grid[[i]] <- start_vec + (i-1)*step
}

scen_labels <- paste0("Miss_Lev_", 1:8)

# --- GLOBAL PARAMS (FIXED RHO) ---
H <- 10; px <- 2; py <- 3; p_bin <- 1; p_cont <- 1
beta_true <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)
sigma_u <- 1
sigma_v_hosp <- seq(3, 5, length.out = H)
sigma_e <- 3

# FIXED CORRELATION SETTING
fixed_rho_vec <- c(-0.3, 0.3, 0.8) # Rho12, Rho13, Rho23

# Build Fixed Covariance Matrix once
R_mat <- diag(py);
R_mat[1,2]=R_mat[2,1]=fixed_rho_vec[1]
R_mat[1,3]=R_mat[3,1]=fixed_rho_vec[2]
R_mat[2,3]=R_mat[3,2]=fixed_rho_vec[3]
Sigma_e_mat <- diag(sigma_e, py) %*% R_mat %*% diag(sigma_e, py)

# ----------------------------------------------------------------
# 2. PARALLEL LOOP OVER ASSIGNED SIMULATIONS
# ----------------------------------------------------------------
final_batch_df <- foreach(global_sim_id = start_idx:end_idx,
                          .combine = rbind,
                          # Safely allowing automatic export of global variables
                          .packages = c("data.table", "dplyr", "MASS", "Matrix")) %dopar% {
                            
  # --- SETUP PARAMS ---
  scen_idx <- ceiling(global_sim_id / N_reps_per_scen)
  if(scen_idx > length(gamma_grid)) return(NULL)
  
  curr_gamma <- gamma_grid[[scen_idx]]
  curr_label <- scen_labels[scen_idx]
  seed_val   <- global_sim_id * 321
  
  # --- DATA GEN ---
  set.seed(seed_val)
  
  m_hosp <- sample(50:100, H, replace=TRUE)
  id_hosp <- rep(1:H, m_hosp)
  N_pat <- sum(m_hosp)
  n_visits <- sample(5:30, N_pat, replace=TRUE)
  
  obs_site <- rep(id_hosp, n_visits)
  obs_pat  <- rep(1:N_pat, n_visits)
  N_obs <- sum(n_visits)
  
  u_mat <- matrix(0, H, py)
  for(j in 1:py) u_mat[,j] <- rnorm(H, 0, sigma_u)
  u_exp <- u_mat[obs_site, ]
  
  v_mat <- matrix(0, N_pat, py)
  pat_sds <- rep(sigma_v_hosp, m_hosp)
  for(j in 1:py) v_mat[,j] <- rnorm(N_pat, 0, pat_sds)
  v_exp <- v_mat[obs_pat, ]
  
  X_bin  <- matrix(rbinom(N_obs*p_bin, 1, 0.3), N_obs, p_bin)
  X_cont <- matrix(rnorm(N_obs*p_cont), N_obs, p_cont)
  X_hij  <- cbind(X_bin, X_cont)
  eps    <- mvrnorm(N_obs, rep(0,py), Sigma_e_mat)
  
  Y_hij <- matrix(0, N_obs, py)
  for(k in 1:py) {
    Y_hij[,k] <- X_hij %*% beta_true[,k] + u_exp[,k] + v_exp[,k] + eps[,k]
  }
  
  dat <- data.frame(site=factor(obs_site), patient=factor(obs_pat), X_hij, Y_hij)
  X_cols <- paste0("X",1:px); Y_cols <- paste0("Y",1:py)
  names(dat)[3:ncol(dat)] <- c(X_cols, Y_cols)
  
  # --- APPLY MISSINGNESS (VARYING GAMMA) ---
  dat_miss <- dat
  lp_base <- 0.5*dat$X1 + 0.4*dat$X2
  
  for(j in 1:py) {
    # P(Missing) = logit_inv(gamma0[j] + lp)
    pmiss <- plogis(curr_gamma[j] + lp_base)
    is_na <- rbinom(N_obs, 1, pmiss)
    dat_miss[is_na==1, Y_cols[j]] <- NA
  }
  
  # --- FITTING ---
  models <- c("independence", "exchangeable", "unstructured")
  sim_results_list <- list()
  
  for(m in models) {
    fit <- tryCatch(
      peal.fit(dat_miss, X_cols, Y_cols, "site", "patient", corstr=m, reml=TRUE, verbose=FALSE),
      error=function(e) return(NULL)
    )
    
    if(!is.null(fit) && !is.null(fit$b)) {
      r_est <- rep(NA, 3)
      if(m=="exchangeable") r_est[] <- fit$rho
      if(m=="unstructured") { C<-fit$Corr; r_est<-c(C[2,1],C[3,1],C[3,2]) }
      if(m=="independence") r_est[] <- 0
      
      b_vec <- as.vector(fit$b)
      # Safe handling if b.sd is missing in some failed convergence cases
      se_vec <- if(!is.null(fit$b.sd)) as.vector(fit$b.sd) else rep(NA, length(b_vec))
      
      beta_names <- paste0("Beta", 1:length(b_vec))
      se_names   <- paste0("SE", 1:length(se_vec))
      
      df_row <- data.frame(
        SimID = global_sim_id,
        Scenario = curr_label, # e.g., "Miss_Lev_1"
        Model = m,
        sigma_u = sqrt(fit$theta[1]*fit$s2),
        sigma_v = sqrt((fit$theta[2:(H+1)])*fit$s2)[1], 
        sigma_e = sqrt(fit$s2),
        rho12 = r_est[1], rho13 = r_est[2], rho23 = r_est[3]
      )
      df_row[beta_names] <- b_vec
      df_row[se_names] <- se_vec
      
      sim_results_list[[length(sim_results_list)+1]] <- df_row
    }
  }
  
  if(length(sim_results_list) > 0) {
    return(do.call(rbind, sim_results_list))
  } else {
    return(NULL)
  }
}

stopCluster(cl)

# ----------------------------------------------------------------
# 3. SAVE BATCH
# ----------------------------------------------------------------
if(!is.null(final_batch_df) && nrow(final_batch_df) > 0) {
  out_dir <- "results"
  if(!dir.exists(out_dir)) dir.create(out_dir)
  out_file <- file.path(out_dir, paste0("res_batch_", task_id, ".rds"))
  
  saveRDS(final_batch_df, out_file)
  cat("Saved batch of", nrow(final_batch_df)/3, "simulations to:", out_file, "\n")
} else {
  cat("No results generated.\n")
}