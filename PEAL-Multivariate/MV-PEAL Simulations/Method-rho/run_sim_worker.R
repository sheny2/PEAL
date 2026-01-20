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

cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  library(MASS)
  library(Matrix)
  library(Rcpp)
  library(RcppArmadillo)
  library(mice)
  
  if(!file.exists("MV-PEAL.R")) stop("MV-PEAL.R not found on worker node.")
  source("MV-PEAL.R")
})

registerDoParallel(cl)

# Batch Logic
N_scenarios     <- 8
N_reps_per_scen <- 500
Total_Sims      <- N_scenarios * N_reps_per_scen 

Array_Size      <- 80
Sims_Per_Job    <- Total_Sims / Array_Size 

start_idx <- (task_id - 1) * Sims_Per_Job + 1
end_idx   <- task_id * Sims_Per_Job

cat(sprintf("Task %d running global simulations %d to %d in parallel\n", task_id, start_idx, end_idx))

# --- SCENARIO CONFIGURATION: VARYING CORRELATION ---
rho_grid <- list(
  c(0.05, 0.10, 0.15), 
  c(0.05, 0.15, 0.25),
  c(0.05, 0.20, 0.35), 
  c(0.05, 0.25, 0.45),
  c(0.05, 0.30, 0.55), 
  c(0.05, 0.35, 0.65),
  c(0.05, 0.45, 0.75), 
  c(0.05, 0.50, 0.85)
)

scen_labels <- sapply(rho_grid, function(r) {
  paste0("Rho_", paste(sprintf("%.2f", r), collapse = "_"))
})

# --- GLOBAL PARAMS (FIXED MISSINGNESS) ---
H <- 10; px <- 2; py <- 3; p_bin <- 1; p_cont <- 1
beta_true <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)
sigma_u <- 1
sigma_v_hosp <- seq(3, 5, length.out = H)
sigma_e <- 3

# Fixed Missingness Parameters (Medium-Low)
gamma_0_fixed <- c(-2.5, -2.0, -1.5) # Increasing intercept -> slightly more missing in Y2, Y3
gamma_1 <- 0.5
gamma_2 <- 0.4

# ----------------------------------------------------------------
# 2. PARALLEL LOOP
# ----------------------------------------------------------------
final_batch_df <- foreach(global_sim_id = start_idx:end_idx,
                          .combine = rbind,
                          .packages = c("data.table", "dplyr", "MASS", "Matrix", "mice")) %dopar% {
                            
                            # --- SETUP PARAMS ---
                            scen_idx <- ceiling(global_sim_id / N_reps_per_scen)
                            if(scen_idx > length(rho_grid)) return(NULL)
                            
                            curr_rho   <- rho_grid[[scen_idx]]
                            curr_label <- scen_labels[scen_idx]
                            seed_val   <- global_sim_id * 555
                            
                            # Build Covariance Matrix for this specific scenario
                            R_mat <- diag(py)
                            R_mat[1,2] <- R_mat[2,1] <- curr_rho[1]
                            R_mat[1,3] <- R_mat[3,1] <- curr_rho[2]
                            R_mat[2,3] <- R_mat[3,2] <- curr_rho[3]
                            
                            # Ensure Positive Definiteness (safeguard)
                            eig <- eigen(R_mat)
                            if(any(eig$values <= 0)) return(NULL) # Skip invalid
                            
                            Sigma_e_mat <- diag(sigma_e, py) %*% R_mat %*% diag(sigma_e, py)
                            
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
                            
                            # Full Data
                            dat_full <- data.frame(site=factor(obs_site), patient=factor(obs_pat), X_hij, Y_hij)
                            X_cols <- paste0("X",1:px); Y_cols <- paste0("Y",1:py)
                            names(dat_full)[3:ncol(dat_full)] <- c(X_cols, Y_cols)
                            
                            # --- APPLY MISSINGNESS (FIXED GAMMA) ---
                            dat_miss <- dat_full
                            lp_base <- gamma_1*dat_full$X1 + gamma_2*dat_full$X2
                            
                            for(j in 1:py) {
                              pmiss <- plogis(gamma_0_fixed[j] + lp_base)
                              is_na <- rbinom(N_obs, 1, pmiss)
                              dat_miss[is_na==1, Y_cols[j]] <- NA
                            }
                            
                            # -------------------------------------------------------------
                            # DEFINE THE FOUR APPROACHES
                            # -------------------------------------------------------------
                            data_list <- list()
                            data_list[["Full"]]     <- dat_full
                            data_list[["Observed"]] <- dat_miss
                            data_list[["CC"]]       <- dat_miss[complete.cases(dat_miss[, Y_cols]), ]
                            
                            # MICE Imputation
                            dat_mice <- tryCatch({
                              imp_parts <- list()
                              site_indices <- split(1:nrow(dat_miss), dat_miss$site)
                              for(s_id in names(site_indices)) {
                                idx <- site_indices[[s_id]]
                                sub_d <- dat_miss[idx, ]
                                mice_out <- mice(sub_d, m=1, method='pmm', maxit=5, printFlag=FALSE)
                                imp_parts[[s_id]] <- complete(mice_out, 1)
                              }
                              do.call(rbind, imp_parts)
                            }, error = function(e) NULL)
                            
                            if(!is.null(dat_mice)) data_list[["MICE"]] <- dat_mice
                            
                            # -------------------------------------------------------------
                            # FIT MODELS
                            # -------------------------------------------------------------
                            sim_results_list <- list()
                            
                            for(approach_name in names(data_list)) {
                              current_data <- data_list[[approach_name]]
                              if(nrow(current_data) < 10) next 
                              
                              fit <- tryCatch(
                                peal.fit(
                                  data        = current_data, 
                                  X_cols      = X_cols, 
                                  Y_cols      = Y_cols, 
                                  site_col    = "site", 
                                  patient_col = "patient", 
                                  corstr      = "unstructured", 
                                  reml        = TRUE, 
                                  verbose     = FALSE
                                ), error=function(e) return(NULL)
                              )
                              
                              if(!is.null(fit) && !is.null(fit$b)) {
                                C <- fit$Corr
                                r_est <- c(C[2,1], C[3,1], C[3,2])
                                
                                b_vec <- as.vector(fit$b)
                                se_vec <- if(!is.null(fit$b.sd)) as.vector(fit$b.sd) else rep(NA, length(b_vec))
                                
                                beta_names <- paste0("Beta", 1:length(b_vec))
                                se_names   <- paste0("SE", 1:length(se_vec))
                                
                                df_row <- data.frame(
                                  SimID    = global_sim_id,
                                  Scenario = curr_label, 
                                  Model    = approach_name,
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
                            
                            if(length(sim_results_list) > 0) return(do.call(rbind, sim_results_list))
                            else return(NULL)
                          }

stopCluster(cl)

if(!is.null(final_batch_df) && nrow(final_batch_df) > 0) {
  out_dir <- "results"
  if(!dir.exists(out_dir)) dir.create(out_dir)
  out_file <- file.path(out_dir, paste0("res_batch_", task_id, ".rds"))
  saveRDS(final_batch_df, out_file)
}