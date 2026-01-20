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
n_cores <- parallelly::availableCores()
cat(sprintf("Task %d: Detected %d available cores. Initializing cluster...\n", task_id, n_cores))

# 1. Create the cluster explicitly
cl <- makeCluster(n_cores)

# 2. Export libraries and Source the functions ONCE on every worker
#    This prevents file contention and repeated recompilation
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  library(MASS)
  library(Matrix)
  library(Rcpp)
  library(RcppArmadillo)
  
  # Ensure the worker is in the right directory or provide full path
  if(!file.exists("MV-PEAL.R")) stop("MV-PEAL.R not found on worker")
  source("MV-PEAL.R") 
})

# 3. Register the cluster
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

# Scenarios
rho_grid <- list(
  c(0.05, 0.10, 0.15), c(0.05, 0.15, 0.25),
  c(0.05, 0.20, 0.35), c(0.05, 0.25, 0.45),
  c(0.05, 0.30, 0.55), c(0.05, 0.35, 0.65),
  c(0.05, 0.45, 0.75), c(0.05, 0.50, 0.85)
)
scen_labels <- sapply(rho_grid, function(r) paste0("rho_", paste(sprintf("%0.2f", r), collapse="_")))

# Global Params
H <- 10; px <- 2; py <- 3; p_bin <- 1; p_cont <- 1
beta_true <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)
sigma_u <- 1
sigma_v_hosp <- seq(3, 5, length.out = H)
sigma_e <- 3

# ----------------------------------------------------------------
# 2. PARALLEL LOOP
# ----------------------------------------------------------------

# Note: We do NOT need to source("MV-PEAL.R") inside here anymore.
# We do NOT need .export="peal.fit" because it exists in the worker's global env via clusterEvalQ.

final_batch_df <- foreach(global_sim_id = start_idx:end_idx,
                          .combine = rbind,
                          # Only export data objects needed for generation
                          .export = c("rho_grid", "scen_labels", "N_reps_per_scen", 
                                      "H", "px", "py", "p_bin", "p_cont", 
                                      "beta_true", "sigma_u", "sigma_v_hosp", "sigma_e"),
                          .packages = c("data.table", "dplyr", "MASS", "Matrix")) %dopar% {
                            
  # --- SETUP PARAMS ---
  scen_idx <- ceiling(global_sim_id / N_reps_per_scen)
  if(scen_idx > length(rho_grid)) return(NULL)

  curr_rho   <- rho_grid[[scen_idx]]
  curr_label <- scen_labels[scen_idx]
  seed_val   <- global_sim_id * 321

  # Build Covariance
  R_mat <- diag(py);
  R_mat[1,2]=R_mat[2,1]=curr_rho[1]; R_mat[1,3]=R_mat[3,1]=curr_rho[2]; R_mat[2,3]=R_mat[3,2]=curr_rho[3]
  
  # Check Positive Definite (Prevent C++ Crash)
  eigen_vals <- eigen(R_mat, only.values=TRUE)$values
  if(min(eigen_vals) <= 1e-6) return(NULL) # Skip invalid correlation matrices
  
  Sigma_e_mat <- diag(sigma_e, py) %*% R_mat %*% diag(sigma_e, py)

  # --- DATA GEN ---
  set.seed(seed_val)

  m_hosp <- sample(50:100, H, replace=TRUE)
  id_hosp <- rep(1:H, m_hosp)
  N_pat <- sum(m_hosp)
  n_visits <- sample(5:20, N_pat, replace=TRUE)

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

  dat_miss <- dat
  # MAR Settings
  gamma_1 <- 0.5  
  gamma_2 <- 0.4  
  gamma_0_vec <- c(-2.5, -1.5, -0.5) 

  for(j in seq_along(Y_cols)){
    y_col <- Y_cols[j]
    current_gamma_0 <- gamma_0_vec[j]
    linear_pred <- current_gamma_0 + gamma_1 * dat$X1 + gamma_2 * dat$X2 # Assuming px >= 2
    prob_miss <- 1 / (1 + exp(-linear_pred))
    is_missing <- rbinom(n = nrow(dat), size = 1, prob = prob_miss)
    dat_miss[is_missing == 1, y_col] <- NA
  }

  # --- FITTING ---
  models <- c("independence", "exchangeable", "unstructured")
  sim_results_list <- list() 

  for(m in models) {
    fit <- tryCatch(
      peal.fit(dat_miss, X_cols, Y_cols, "site", "patient", corstr=m, reml=TRUE, verbose=FALSE),
      error=function(e) return(NULL) # Catch R errors
    )

    if(!is.null(fit) && !is.null(fit$b)) {
      r_est <- rep(NA, 3)
      if(m=="exchangeable") r_est[] <- fit$rho
      if(m=="unstructured") { C<-fit$Corr; r_est<-c(C[2,1],C[3,1],C[3,2]) }
      if(m=="independence") r_est[] <- 0

      b_vec <- as.vector(fit$b)
      # Handle cases where b.sd might be null
      se_vec <- if(!is.null(fit$b.sd)) as.vector(fit$b.sd) else rep(NA, length(b_vec))
      
      beta_names <- paste0("Beta", 1:length(b_vec))
      se_names   <- paste0("SE", 1:length(se_vec))

      df_row <- data.frame(
        SimID = global_sim_id,
        Scenario = curr_label,
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

# 3. Stop Cluster
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