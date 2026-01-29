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
  library(mice) # Added for imputation
  
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
start_vec <- c(-3, -2, -1)
step <- 0.5

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

# FIXED CORRELATION SETTING (Used for generation)
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
                          .packages = c("data.table", "dplyr", "MASS", "Matrix", "mice")) %dopar% {
                            
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
  
  # Full Dataset 
  dat_full <- data.frame(site=factor(obs_site), patient=factor(obs_pat), X_hij, Y_hij)
  X_cols <- paste0("X",1:px); Y_cols <- paste0("Y",1:py)
  names(dat_full)[3:ncol(dat_full)] <- c(X_cols, Y_cols)
  
  # ----------------------------------------------------------------
  # REPLACEMENT: SEQUENTIAL MISSINGNESS 
  # ----------------------------------------------------------------
  dat_miss <- dat_full
  lp_base <- 0.2 * dat_full$X1 + 0.1 * dat_full$X2 
  
  # Interaction Strength (Log-Odds shift)
  theta <- 2
  
  # 1. Y1 Missingness: Depends only on X (and Gamma)
  pmiss_1 <- plogis(curr_gamma[1] + lp_base)
  is_na_1 <- rbinom(N_obs, 1, pmiss_1)
  
  # 2. Y2 Missingness: Depends on X, Gamma, AND R1
  pmiss_2 <- plogis(curr_gamma[2] + lp_base + theta * is_na_1)
  is_na_2 <- rbinom(N_obs, 1, pmiss_2)
  
  # 3. Y3 Missingness: Depends on X, Gamma, AND R1, R2
  pmiss_3 <- plogis(curr_gamma[3] + lp_base - theta * (is_na_1 + is_na_2))
  is_na_3 <- rbinom(N_obs, 1, pmiss_3)
  
  all_missing <- (is_na_1 == 1 & is_na_2 == 1 & is_na_3 == 1)
  
  # If all are missing, reset Y1 observed
  is_na_1[all_missing] <- 0
  
  # Apply NAs to Dataset
  dat_miss$Y1[is_na_1 == 1] <- NA
  dat_miss$Y2[is_na_2 == 1] <- NA
  dat_miss$Y3[is_na_3 == 1] <- NA
  
  # -------------------------------------------------------------
  # DEFINE THE FOUR APPROACHES
  # -------------------------------------------------------------
  
  # 1. Full Data (Ideal)
  data_list <- list()
  data_list[["Full"]] <- dat_full
  
  # 2. Observed Data (Proposed Method - handles missing internally)
  data_list[["Observed"]] <- dat_miss
  
  # 3. Complete Case (CC)
  data_list[["CC"]] <- dat_miss[complete.cases(dat_miss[, Y_cols]), ]
  
  # 4. MICE (Imputation on Full Data)
  # -------------------------------------------------------------
  dat_mice <- tryCatch({
    # We impute the full dataset at once. 
    # Note: We exclude 'patient' ID from predictors to avoid rank deficiency/slowness 
    # with high-cardinality factors, but include 'site' (fixed effects).
    pred_mat <- make.predictorMatrix(dat_miss)
    if("patient" %in% rownames(pred_mat)) {
      pred_mat[,"patient"] <- 0
      pred_mat["patient",] <- 0
    }
    
    mice_out <- mice(dat_miss, m=1, method='pmm', 
                     predictorMatrix = pred_mat,
                     maxit=5, printFlag=FALSE)
    complete(mice_out, 1)
  }, error = function(e) return(NULL))
  
  if(!is.null(dat_mice)) {
    data_list[["MICE"]] <- dat_mice
  }
  
  # -------------------------------------------------------------
  # FIT MODELS
  # -------------------------------------------------------------
  
  sim_results_list <- list()
  
  for(approach_name in names(data_list)) {
    
    current_data <- data_list[[approach_name]]
    
    # Skip empty datasets (e.g., if CC removes everything)
    if(nrow(current_data) < 10) next 
    
    # Fit MV-PEAL with Unstructured Correlation for ALL approaches
    fit <- tryCatch(
      peal.fit(
        data        = current_data, 
        X_cols      = X_cols, 
        Y_cols      = Y_cols, 
        site_col    = "site", 
        patient_col = "patient", 
        corstr      = "unstructured", # ALWAYS Unstructured
        reml        = TRUE, 
        verbose     = FALSE
      ),
      error=function(e) return(NULL)
    )
    
    if(!is.null(fit) && !is.null(fit$b)) {
      # Extract Unstructured Correlation parameters
      C <- fit$Corr
      r_est <- c(C[2,1], C[3,1], C[3,2])
      
      b_vec <- as.vector(fit$b)
      # Safe handling if b.sd is missing 
      se_vec <- if(!is.null(fit$b.sd)) as.vector(fit$b.sd) else rep(NA, length(b_vec))
      
      beta_names <- paste0("Beta", 1:length(b_vec))
      se_names   <- paste0("SE", 1:length(se_vec))
      
      df_row <- data.frame(
        SimID    = global_sim_id,
        Scenario = curr_label, 
        Model    = approach_name, # "Full", "Observed", "CC", "MICE"
        sigma_u  = sqrt(fit$theta[1]*fit$s2),
        sigma_v  = sqrt((fit$theta[2:(H+1)])*fit$s2)[1], 
        sigma_e  = sqrt(fit$s2),
        rho12    = r_est[1], 
        rho13    = r_est[2], 
        rho23    = r_est[3]
      )
      df_row[beta_names] <- b_vec
      df_row[se_names]   <- se_vec
      
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
  cat("Saved batch of", nrow(final_batch_df)/length(unique(final_batch_df$Model)), "simulations to:", out_file, "\n")
} else {
  cat("No results generated.\n")
}