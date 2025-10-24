library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(MASS)
library(reshape2)
library(ggplot2)

source("PEAL_Engine_Missing_Multi-RI.R")

# -------------------------
# Parallel setup
# -------------------------
N <- 30
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# -------------------------
# DGP Parameters
# -------------------------
H <- 5
m_hosp <- sample(10:30, H, replace = TRUE)
px <- 6
p_bin <- 3
p_cont <- px - p_bin
py <- 3

# Fixed effects
beta <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)

# Variance components (SDs)
sigma_u <- 0.3
sigma_v_hosp <- c(0.61, 0.63, 0.65, 0.67, 0.69)

# sigma_u <- 0.3 * 10
# sigma_v_hosp <- c(0.21, 0.23, 0.25, 0.27, 0.29) * 10

# Residual SD and exchangeable correlation
sigma_e <- 0.4
rho <- 0.5
rho_mat <- matrix(rho, nrow = py, ncol = py); diag(rho_mat) <- 1

# True values
true_beta  <- c(beta)
true_sigma <- c(sigma_u, sigma_v_hosp, sigma_e, rho)

# -------------------------
# Pre-allocate result matrices (keep N columns; fill NA on failure)
# -------------------------
par_names_beta  <- paste0("Beta", 1:(px*py))
par_names_sigma <- c("sigma_u", paste0("sigma_v_", 1:H), "sigma_e", "rho")

result_beta_exch  <- matrix(NA_real_, nrow = length(par_names_beta),  ncol = N,
                            dimnames = list(par_names_beta, paste0("sim", 1:N)))
result_sigma_exch <- matrix(NA_real_, nrow = length(par_names_sigma), ncol = N,
                            dimnames = list(par_names_sigma, paste0("sim", 1:N)))

result_beta_indep  <- matrix(NA_real_, nrow = length(par_names_beta),  ncol = N,
                             dimnames = list(par_names_beta, paste0("sim", 1:N)))
result_sigma_indep <- matrix(NA_real_, nrow = length(par_names_sigma), ncol = N,
                             dimnames = list(par_names_sigma, paste0("sim", 1:N)))

# -------------------------
# Simulation loop: fit BOTH models per replicate
# -------------------------
results <- foreach(k = 1:N, .packages = c("data.table","dplyr","MASS")) %dopar% {
  
  source("PEAL_Engine_Missing_Multi-RI.R")
  set.seed(k)
  
  # ---- IDs and visits
  nn <- rep(m_hosp, times = 1)
  id.hosp <- rep(1:H, times = m_hosp)
  id.pat <- sequence(nn)
  n_visits <- sample(1:20, sum(nn), replace = TRUE)
  
  id.visit <- sequence(n_visits)
  id.hosp.expanded <- rep(id.hosp, times = n_visits)
  id.pat.expanded  <- rep(id.pat,  times = n_visits)
  
  # ---- Random effects (RIs)
  u_h  <- rnorm(H, mean = 0, sd = sigma_u)                                # site
  v_hi <- rnorm(sum(nn), mean = 0, sd = rep(sigma_v_hosp, times = m_hosp)) # patient
  
  # expand to visits
  u_h_patient  <- rep(u_h, times = m_hosp)
  u_h_expanded <- rep(u_h_patient, times = n_visits)
  v_hi_expanded <- rep(v_hi, times = n_visits)
  
  # ---- Covariates
  Nobs <- sum(n_visits)
  X_bin  <- matrix(rbinom(Nobs * p_bin, size = 1, prob = 0.3), nrow = Nobs, ncol = p_bin)
  X_cont <- matrix(rnorm(Nobs * p_cont, mean = 0, sd = 1), nrow = Nobs, ncol = p_cont)
  X_hij  <- cbind(X_bin, X_cont)
  
  # ---- Residuals (exchangeable across outcomes)
  Sigma_eps <- diag(sigma_e, py) %*% rho_mat %*% diag(sigma_e, py)
  epsilon_hij <- MASS::mvrnorm(Nobs, mu = rep(0, py), Sigma = Sigma_eps)
  
  # ---- Outcomes
  Y_hij <- matrix(0, nrow = Nobs, ncol = py)
  for (j in 1:py) {
    Y_hij[, j] <- X_hij %*% beta[, j] + u_h_expanded + v_hi_expanded + epsilon_hij[, j]
  }
  
  # ---- Data table
  three_lvl_dat <- data.table(
    site = id.hosp.expanded,
    patient = id.pat.expanded,
    visit = id.visit,
    X_hij, Y_hij
  ) |> data.frame()
  
  setnames(three_lvl_dat, c("site","patient","visit", paste0("X",1:px), paste0("Y",1:py)))
  
  # ---- Preprocess and order
  visit_count <- three_lvl_dat |>
    dplyr::group_by(site, patient) |>
    dplyr::summarise(total_visits = n(), .groups = "drop")
  
  rearranged_data <- merge(three_lvl_dat, visit_count, by = c("site","patient")) |>
    arrange(site, total_visits, patient) |>
    mutate(site = factor(site))
  
  # Introduce substantial missingness in Y_1 to Y_R
  # Option 1: missingness rate (e.g., 40% missing overall)
  missing_prob <- 0.4
  
  # Introduce missing values completely at random (MCAR)
  rearranged_data_with_missing <- rearranged_data
  
  # First pass: introduce missingness for each Y variable
  for(y_col in paste0("Y", 1:py)) {
    missing_indices <- sample(1:nrow(rearranged_data), 
                              size = round(nrow(rearranged_data) * missing_prob))
    rearranged_data_with_missing[missing_indices, y_col] <- NA
  }
  
  # Ensure no cases where all Ys are missing
  all_y_missing <- apply(rearranged_data_with_missing[, paste0("Y", 1:py), drop = FALSE], 1, function(x) all(is.na(x)))
  num_all_missing <- sum(all_y_missing)
  
  if (num_all_missing > 0) {
    cat("Found", num_all_missing, "cases where all Ys are missing. Fixing...\n")
    
    # For each case where all Ys are missing, randomly select one Y to be observed
    for (i in which(all_y_missing)) {
      y_to_keep <- sample(paste0("Y", 1:py), 1)
      rearranged_data_with_missing[i, y_to_keep] <- rearranged_data[i, y_to_keep]
    }
  }
  
  
  # ---- Helper to pack outputs (into SDs + rho)
  pack_fit <- function(fit, cor_model) {
    if (is.null(fit)) {
      return(list(beta = rep(NA_real_, px*py),
                  sigma = rep(NA_real_, H+3)))
    }
    # Optional Hessian sign check (similar to your previous guard)
    if (!is.null(fit$opt$hessian)) {
      ev <- tryCatch(eigen(fit$opt$hessian)$value, error = function(e) NA)
      if (all(is.finite(ev)) && all(ev < 0)) {
        return(list(beta = rep(NA_real_, px*py),
                    sigma = rep(NA_real_, H+3)))
      }
    }
    
    s2 <- fit$s2
    thetas <- fit$theta
    sigU  <- sqrt(thetas[1] * s2)
    sigVh <- sqrt(thetas[1 + seq_len(H)] * s2)
    sigE  <- sqrt(s2)
    rh    <- if (!is.null(fit$rho) && cor_model == "exchangeable") fit$rho else 0
    
    list(beta = as.numeric(fit$b),
         sigma = c(sigU, sigVh, sigE, rh))
  }
  

  
  # ---- Fit Exchangeable
  fit_exch <- tryCatch(
    peal.fit.RI_mv_patterns(
      data     = rearranged_data_with_missing,
      X_cols   = paste0("X", 1:px),
      Y_cols   = paste0("Y", 1:py),
      site_col = "site",
      patient_col = "patient",
      reml     = TRUE,
      corstr   = "exchangeable",   
      rho_init = 0.1,
      verbose  = TRUE
    ),
    error = function(e) { NULL },
    warning = function(w) { invokeRestart("muffleWarning") }
  )
  
  out_exch <- pack_fit(fit_exch, "exchangeable")
  
  # ---- Fit Independence (rho fixed 0)
  fit_indep <- tryCatch(
    peal.fit.RI_mv_patterns(
      data     = rearranged_data_with_missing,
      X_cols   = paste0("X", 1:px),
      Y_cols   = paste0("Y", 1:py),
      site_col = "site",
      patient_col = "patient",
      reml     = TRUE,
      corstr   = "independence",  
      rho_init = 0.1,
      verbose  = TRUE
    ),
    error = function(e) { NULL },
    warning = function(w) { invokeRestart("muffleWarning") }
  )
  
  out_indep <- pack_fit(fit_indep, "independence")
  
  list(beta_exch  = out_exch$beta,
       sigma_exch = out_exch$sigma,
       beta_indep  = out_indep$beta,
       sigma_indep = out_indep$sigma)
}

# -------------------------
# Store into pre-allocated matrices
# -------------------------
for (k in seq_along(results)) {
  result_beta_exch[,  k] <- results[[k]]$beta_exch
  result_sigma_exch[, k] <- results[[k]]$sigma_exch
  result_beta_indep[,  k] <- results[[k]]$beta_indep
  result_sigma_indep[, k] <- results[[k]]$sigma_indep
}

stopCluster(cl)

# -------------------------
# Long-format data with Model labels
# -------------------------
mk_long <- function(mat, par_names, model_label) {
  df <- reshape2::melt(as.data.frame(mat), value.name = "Estimate")
  colnames(df) <- c("Simulation","Estimate")
  df$Parameter <- rep(par_names, ncol(mat))
  df$Model <- model_label
  df
}

beta_exch_df  <- mk_long(result_beta_exch,  par_names_beta,  "Exchangeable")
beta_indep_df <- mk_long(result_beta_indep, par_names_beta,  "Independence")
beta_df <- rbind(beta_exch_df, beta_indep_df)
beta_df$True_Value <- rep(true_beta, times = ncol(result_beta_exch) + ncol(result_beta_indep))

sigma_exch_df  <- mk_long(result_sigma_exch,  par_names_sigma, "Exchangeable")
sigma_indep_df <- mk_long(result_sigma_indep, par_names_sigma, "Independence")
sigma_df <- rbind(sigma_exch_df, sigma_indep_df)
sigma_df$True_Value <- rep(true_sigma, times = ncol(result_sigma_exch) + ncol(result_sigma_indep))

# -------------------------
# Bias / Variance / MSE summaries (by Model)
# -------------------------
beta_summary <- beta_df %>%
  group_by(Model, Parameter) %>%
  summarise(
    True_Value = first(True_Value),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias       = Mean_Est - True_Value,
    Variance   = var(Estimate, na.rm = TRUE),
    MSE        = Bias^2 + Variance,
    .groups = "drop"
  ) %>% arrange(Parameter, Model)

sigma_summary <- sigma_df %>%
  group_by(Model, Parameter) %>%
  summarise(
    True_Value = first(True_Value),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias       = Mean_Est - True_Value,
    Variance   = var(Estimate, na.rm = TRUE),
    MSE        = Bias^2 + Variance,
    .groups = "drop"
  ) %>% arrange(Parameter, Model)

# -------------------------
# Plots: bias comparisons faceted by model
# -------------------------
p_beta_bias <- beta_df %>%
  mutate(Bias = Estimate - True_Value) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1, width = 0.15, height = 0) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~ Model, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Fixed Effects Bias: Exchangeable vs Independence")

p_sigma_bias <- sigma_df %>%
  mutate(Bias = Estimate - True_Value) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1, width = 0.15, height = 0) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~ Model, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Variance Components Bias: Exchangeable vs Independence")

print(p_beta_bias)
print(p_sigma_bias)

# Inspect numeric summaries in console
beta_summary
sigma_summary
