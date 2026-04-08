# ==============================================================================
# File: mice_diagnostic_sim.R
#
# Purpose:
#   Evaluate MICE's ability to recover distributional properties of three
#   correlated outcomes in a 3-level (site -> patient -> observation)
#   hierarchical data structure, across 8 missingness severity scenarios.
#
# No model fitting (MV-PEAL or otherwise) — pure imputation diagnostics.
#
# Summary measures tracked per rep:
#   (a) POINTWISE (on masked obs only, where truth is known):
#       - Bias   : mean(Y_imp - Y_true) per outcome
#       - RMSE   : sqrt(mean((Y_imp - Y_true)^2)) per outcome
#       - Coverage: % of masked obs where Y_true falls in 95% imputation interval
#                  (across m imputations)
#
#   (b) DISTRIBUTIONAL (full data, imputed vs true):
#       - Mean bias   : mean(Y_imp_full) - mean(Y_true_full)
#       - Variance ratio: var(Y_imp_full) / var(Y_true_full)
#       - Correlation recovery: rho_imp vs rho_true for Y1Y2, Y1Y3, Y2Y3
#
#   (c) MISSINGNESS RATES: actual % missing per outcome per rep
#
# Two MICE strategies compared:
#   - "SiteLocal"  : impute within each site separately (mirrors run_sim_worker.R)
#   - "Global"     : single MICE call on full dataset, site entered as predictor
#
# ==============================================================================
rm(list = ls())
library(parallel)
library(doParallel)
library(foreach)
library(parallelly)
library(data.table)
library(dplyr)
library(MASS)
library(Matrix)
library(mice)

# ----------------------------------------------------------------
# 0.  CONFIGURATION
# ----------------------------------------------------------------
N_REPS_PER_SCEN <- 500     # reps per scenario (increase to 500+ for final run)
N_SCENARIOS     <- 8
M_IMPS          <- 3       # MICE imputations per rep
N_CORES         <- max(1, parallelly::availableCores() - 1)

# ----------------------------------------------------------------
# 1.  DATA-GENERATING PARAMETERS  (identical to run_sim_worker.R)
# ----------------------------------------------------------------
H            <- 10
px           <- 2; py <- 3
beta_true    <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)
sigma_u      <- 1
sigma_v_hosp <- seq(3, 5, length.out = H)
sigma_e      <- 3

fixed_rho_vec <- c(-0.3, 0.3, 0.8)
R_mat         <- diag(py)
R_mat[1,2] <- R_mat[2,1] <- fixed_rho_vec[1]
R_mat[1,3] <- R_mat[3,1] <- fixed_rho_vec[2]
R_mat[2,3] <- R_mat[3,2] <- fixed_rho_vec[3]
Sigma_e_mat <- diag(sigma_e, py) %*% R_mat %*% diag(sigma_e, py)

TRUE_CORR <- list(rho12 = fixed_rho_vec[1],
                  rho13 = fixed_rho_vec[2],
                  rho23 = fixed_rho_vec[3])

# Missingness scenario grid (identical to run_sim_worker.R)
start_vec  <- c(-3, -2, -1)
step       <- 0.5
gamma_grid <- lapply(1:N_SCENARIOS, function(i) start_vec + (i - 1) * step)
scen_labels <- paste0("Miss_Lev_", 1:N_SCENARIOS)

X_cols <- c("X1", "X2")
Y_cols <- c("Y1", "Y2", "Y3")

# ----------------------------------------------------------------
# 2.  HELPER FUNCTIONS
# ----------------------------------------------------------------

#' Simulate one full dataset (before missingness)
simulate_data <- function(seed) {
  set.seed(seed)
  
  m_hosp    <- sample(30:50, H, replace = TRUE)
  id_hosp   <- rep(1:H, m_hosp)
  N_pat     <- sum(m_hosp)
  n_visits  <- sample(5:20, N_pat, replace = TRUE)
  obs_site  <- rep(id_hosp, n_visits)
  obs_pat   <- rep(1:N_pat,  n_visits)
  N_obs     <- sum(n_visits)
  
  u_mat <- matrix(rnorm(H * py, 0, sigma_u), H, py)
  u_exp <- u_mat[obs_site, ]
  
  v_mat   <- matrix(0, N_pat, py)
  pat_sds <- rep(sigma_v_hosp, m_hosp)
  for (j in 1:py) v_mat[, j] <- rnorm(N_pat, 0, pat_sds)
  v_exp <- v_mat[obs_pat, ]
  
  X_hij <- cbind(rbinom(N_obs, 1, 0.3), rnorm(N_obs))
  eps   <- mvrnorm(N_obs, rep(0, py), Sigma_e_mat)
  Y_hij <- matrix(0, N_obs, py)
  for (k in 1:py) {
    Y_hij[, k] <- X_hij %*% beta_true[, k] + u_exp[, k] + v_exp[, k] + eps[, k]
  }
  
  dat <- data.frame(site    = factor(obs_site),
                    patient = factor(obs_pat),
                    X_hij, Y_hij)
  names(dat)[3:7] <- c(X_cols, Y_cols)
  dat
}

#' Mask outcomes according to scenario gamma
mask_data <- function(dat_full, gamma) {
  N_obs    <- nrow(dat_full)
  lp_base  <- 0.2 * dat_full$X1 + 0.1 * dat_full$X2
  theta    <- 2
  
  is_na_1 <- rbinom(N_obs, 1, plogis(gamma[1] + lp_base))
  is_na_2 <- rbinom(N_obs, 1, plogis(gamma[2] + lp_base + theta * is_na_1))
  is_na_3 <- rbinom(N_obs, 1, plogis(gamma[3] + lp_base - theta * (is_na_1 + is_na_2)))
  
  # Prevent all three being simultaneously missing
  all_miss        <- (is_na_1 + is_na_2 + is_na_3) == 3
  is_na_1[all_miss] <- 0
  
  mask <- list(Y1 = is_na_1 == 1,
               Y2 = is_na_2 == 1,
               Y3 = is_na_3 == 1)
  
  dat_miss       <- dat_full
  dat_miss$Y1[mask$Y1] <- NA
  dat_miss$Y2[mask$Y2] <- NA
  dat_miss$Y3[mask$Y3] <- NA
  
  list(dat_miss = dat_miss, mask = mask)
}

#' Pool m completed datasets: return list of pooled Y columns
pool_imputations <- function(imp_list) {
  # For each variable, average across imputations at observation level
  pooled <- Reduce("+", lapply(imp_list, function(d) d[, Y_cols])) / length(imp_list)
  as.data.frame(pooled)
}

#' Compute per-outcome pointwise diagnostics on masked observations
pointwise_diag <- function(imp_list, dat_full, mask) {
  results <- list()
  for (y in Y_cols) {
    idx <- which(mask[[y]])
    if (length(idx) == 0) next
    
    true_vals <- dat_full[[y]][idx]
    
    # Pointwise values across imputations
    imp_mat   <- do.call(cbind, lapply(imp_list, function(d) d[[y]][idx]))
    imp_mean  <- rowMeans(imp_mat)
    
    bias <- mean(imp_mean - true_vals)
    rmse <- sqrt(mean((imp_mean - true_vals)^2))
    
    # 95% interval coverage (min/max across m imputations as proxy interval)
    lo  <- apply(imp_mat, 1, min)
    hi  <- apply(imp_mat, 1, max)
    cov <- mean(true_vals >= lo & true_vals <= hi)
    
    results[[y]] <- c(bias = bias, rmse = rmse, coverage = cov,
                      n_miss = length(idx))
  }
  results
}

#' Compute distributional diagnostics on full (N_obs) dataset
distributional_diag <- function(pooled_Y, dat_full) {
  diag_list <- list()
  
  # Per-outcome mean bias and variance ratio
  for (y in Y_cols) {
    true_mean <- mean(dat_full[[y]])
    true_var  <- var(dat_full[[y]])
    imp_mean  <- mean(pooled_Y[[y]])
    imp_var   <- var(pooled_Y[[y]])
    diag_list[[paste0("mean_bias_", y)]]  <- imp_mean - true_mean
    diag_list[[paste0("var_ratio_", y)]]  <- imp_var  / true_var
  }
  
  # Pairwise correlations
  true_cor  <- cor(dat_full[, Y_cols])
  imp_cor   <- cor(pooled_Y)
  pairs     <- list(c(1,2), c(1,3), c(2,3))
  pair_names <- c("rho12", "rho13", "rho23")
  for (i in seq_along(pairs)) {
    r <- pairs[[i]]
    diag_list[[paste0("cor_bias_", pair_names[i])]] <-
      imp_cor[r[1], r[2]] - true_cor[r[1], r[2]]
    diag_list[[paste0("cor_imp_",  pair_names[i])]] <- imp_cor[r[1], r[2]]
    diag_list[[paste0("cor_true_", pair_names[i])]] <- true_cor[r[1], r[2]]
  }
  
  diag_list
}

# ----------------------------------------------------------------
# 3.  SINGLE-REP FUNCTION
# ----------------------------------------------------------------
run_one_rep <- function(rep_id, scen_idx) {
  seed     <- rep_id * 321 + scen_idx * 7
  gamma    <- gamma_grid[[scen_idx]]
  scen_lbl <- scen_labels[scen_idx]
  
  dat_full <- simulate_data(seed)
  masked   <- mask_data(dat_full, gamma)
  dat_miss <- masked$dat_miss
  mask     <- masked$mask
  
  miss_rates <- sapply(Y_cols, function(y) mean(is.na(dat_miss[[y]])))
  
  # ---------- Strategy A: Site-local MICE ----------
  run_mice_strategy <- function(strategy) {
    tryCatch({
      imp_list <- vector("list", M_IMPS)
      
      if (strategy == "SiteLocal") {
        site_idx <- split(seq_len(nrow(dat_miss)), dat_miss$site)
        for (m in seq_len(M_IMPS)) {
          imp_list[[m]] <- do.call(rbind, lapply(site_idx, function(idx) {
            complete(mice(dat_miss[idx, ], m = 1, maxit = 5,
                          printFlag = FALSE, method = "pmm"), 1)
          }))
        }
      } else {  # Global: site as predictor
        predictor_matrix <- quickpred(dat_miss,
                                      include = c(X_cols, "site"),
                                      mincor  = 0)
        mids_obj <- mice(dat_miss, m = M_IMPS, maxit = 5,
                         predictorMatrix = predictor_matrix,
                         printFlag = FALSE, method = "pmm")
        for (m in seq_len(M_IMPS)) imp_list[[m]] <- complete(mids_obj, m)
      }
      
      pooled_Y <- pool_imputations(imp_list)
      pt_diag  <- pointwise_diag(imp_list, dat_full, mask)
      di_diag  <- distributional_diag(pooled_Y, dat_full)
      
      # Flatten to one row
      row <- data.frame(
        RepID    = rep_id,
        ScenIdx  = scen_idx,
        Scenario = scen_lbl,
        Strategy = strategy,
        stringsAsFactors = FALSE
      )
      
      # Missingness rates
      for (y in Y_cols) row[[paste0("miss_rate_", y)]] <- miss_rates[y]
      
      # Pointwise metrics per outcome
      for (y in Y_cols) {
        if (!is.null(pt_diag[[y]])) {
          row[[paste0("pw_bias_",     y)]] <- pt_diag[[y]]["bias"]
          row[[paste0("pw_rmse_",     y)]] <- pt_diag[[y]]["rmse"]
          row[[paste0("pw_coverage_", y)]] <- pt_diag[[y]]["coverage"]
        } else {
          row[[paste0("pw_bias_",     y)]] <- NA
          row[[paste0("pw_rmse_",     y)]] <- NA
          row[[paste0("pw_coverage_", y)]] <- NA
        }
      }
      
      # Distributional metrics
      for (nm in names(di_diag)) row[[nm]] <- di_diag[[nm]]
      
      row
    }, error = function(e) {
      data.frame(RepID = rep_id, ScenIdx = scen_idx, Scenario = scen_lbl,
                 Strategy = strategy, Error = conditionMessage(e),
                 stringsAsFactors = FALSE)
    })
  }
  
  rbind(run_mice_strategy("SiteLocal"),
        run_mice_strategy("Global"))
}

# ----------------------------------------------------------------
# 4.  PARALLEL EXECUTION
# ----------------------------------------------------------------
cat(sprintf("Launching parallel simulation: %d scenarios x %d reps = %d total reps\n",
            N_SCENARIOS, N_REPS_PER_SCEN, N_SCENARIOS * N_REPS_PER_SCEN))
cat(sprintf("Using %d cores | %d MICE imputations per rep\n\n", N_CORES, M_IMPS))

cl <- makeCluster(N_CORES)
clusterEvalQ(cl, {
  library(data.table); library(dplyr); library(MASS)
  library(Matrix); library(mice)
})
clusterExport(cl, varlist = c(
  "H", "px", "py", "beta_true", "sigma_u", "sigma_v_hosp",
  "sigma_e", "Sigma_e_mat", "gamma_grid", "scen_labels",
  "X_cols", "Y_cols", "M_IMPS",
  "simulate_data", "mask_data", "pool_imputations",
  "pointwise_diag", "distributional_diag", "run_one_rep"
))
registerDoParallel(cl)

# Build task list: all (rep, scenario) combos
tasks <- expand.grid(rep_id   = seq_len(N_REPS_PER_SCEN),
                     scen_idx = seq_len(N_SCENARIOS))

all_results <- foreach(
  i         = seq_len(nrow(tasks)),
  .combine  = rbind,
  .packages = c("dplyr", "MASS", "mice")
) %dopar% {
  run_one_rep(rep_id   = tasks$rep_id[i],
              scen_idx = tasks$scen_idx[i])
}

stopCluster(cl)

# ----------------------------------------------------------------
# 5.  AGGREGATE SUMMARY STATISTICS
# ----------------------------------------------------------------
cat("\nSimulation complete. Computing summaries...\n")

# Keep only successful rows (no Error column content)
# clean_results <- all_results[is.na(all_results$Error) | !("Error" %in% names(all_results)), ]
clean_results <- all_results

numeric_cols <- c(
  paste0("miss_rate_", Y_cols),
  paste0("pw_bias_",     Y_cols),
  paste0("pw_rmse_",     Y_cols),
  paste0("pw_coverage_", Y_cols),
  paste0("mean_bias_",   Y_cols),
  paste0("var_ratio_",   Y_cols),
  paste0("cor_bias_",  c("rho12","rho13","rho23")),
  paste0("cor_imp_",   c("rho12","rho13","rho23")),
  paste0("cor_true_",  c("rho12","rho13","rho23"))
)
numeric_cols <- intersect(numeric_cols, names(clean_results))

summary_tbl <- clean_results %>%
  group_by(Scenario, Strategy) %>%
  summarise(
    N_reps = n(),
    across(
      all_of(numeric_cols),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x,   na.rm = TRUE),
        p25  = ~quantile(.x, 0.25, na.rm = TRUE),
        p75  = ~quantile(.x, 0.75, na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  )

# ----------------------------------------------------------------
# 6.  PRINT KEY METRICS TO CONSOLE
# ----------------------------------------------------------------

key_metrics <- c(
  "miss_rate_Y1__mean", "miss_rate_Y2__mean", "miss_rate_Y3__mean",
  "pw_bias_Y1__mean",   "pw_bias_Y2__mean",   "pw_bias_Y3__mean",
  "pw_rmse_Y1__mean",   "pw_rmse_Y2__mean",   "pw_rmse_Y3__mean",
  "pw_coverage_Y1__mean","pw_coverage_Y2__mean","pw_coverage_Y3__mean",
  "mean_bias_Y1__mean", "mean_bias_Y2__mean", "mean_bias_Y3__mean",
  "var_ratio_Y1__mean", "var_ratio_Y2__mean", "var_ratio_Y3__mean",
  "cor_bias_rho12__mean","cor_bias_rho13__mean","cor_bias_rho23__mean"
)
key_metrics <- intersect(key_metrics, names(summary_tbl))

cat("\n===== MICE DIAGNOSTIC SUMMARY (key metrics) =====\n")
print(summary_tbl %>% dplyr::select(Scenario, Strategy, N_reps, all_of(key_metrics)),
      n = Inf, width = Inf)

# ----------------------------------------------------------------
# 7.  SAVE OUTPUTS
# ----------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")
saveRDS(all_results,  "results/mice_diag_raw.rds")
saveRDS(summary_tbl, "results/mice_diag_summary.rds")
write.csv(summary_tbl, "results/mice_diag_summary.csv", row.names = FALSE)

source("MICE-Results.R")