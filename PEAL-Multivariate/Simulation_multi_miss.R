rm(list = ls())
library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(MASS)
library(ggplot2)
library(Matrix)

# Ensure the consolidated source file is present
if(!file.exists("MV-PEAL.R")) stop("MV-PEAL.R file not found in working directory.")

# Define number of cores
num_cores <- max(1, parallelly::availableCores())
cl <- makeCluster(num_cores)
registerDoParallel(cl)

N_sim = 100  # Number of simulations

# ---------------------------------------------------------
# Global Parameters
# ---------------------------------------------------------
H <- 5  # Number of sites
m_hosp <- sample(30:50, H, replace = TRUE)  # Number of patients per site
px <- 2  # Number of covariates
p_bin <- 1  # Number of binary covariates
p_cont <- px - p_bin  # Number of continuous covariates
py <- 3 # Number of outcomes (multivariate)

# Fixed effects
beta <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)

# Variance components (Standard Deviations)
sigma_u <- 1 
sigma_v_hosp <-  seq(0.5, 1, length.out = H)
sigma_e <- 0.3

# UNSTRUCTURED Correlation Setup
rho_true_vec <- c(0.2, 0.3, 0.7)
Rho_mat <- diag(py)
Rho_mat[2,1] <- Rho_mat[1,2] <- rho_true_vec[1]
Rho_mat[3,1] <- Rho_mat[1,3] <- rho_true_vec[2]
Rho_mat[3,2] <- Rho_mat[2,3] <- rho_true_vec[3]

if(min(eigen(Rho_mat)$values) <= 0) stop("Chosen correlation matrix is not positive definite.")

# Container Definitions
result_beta = matrix(nrow = (px*py), ncol = N_sim)
rownames(result_beta) <- paste0("Beta", 1:(px*py))

result_sigma = matrix(nrow = (H + 1 + 1 + length(rho_true_vec)), ncol = N_sim)
rownames(result_sigma) <- c("sigma_u", paste0("sigma_v_", 1:H), "sigma_e", 
                            "rho_12", "rho_13", "rho_23")

# ---------------------------------------------------------
# Simulation Loop
# ---------------------------------------------------------
results <- foreach(k = 1:N_sim, 
                   .packages = c("data.table", "dplyr", "MASS", "Matrix", "Rcpp", "RcppArmadillo"),
                   .noexport = c("peal.fit")) %dopar% {
                     
   source("MV-PEAL.R")
   set.seed(k)
   
   # --- 1. Data Generation ---
   nn <- rep(m_hosp, times = 1) 
   id.hosp <- rep(1:H, times = m_hosp)
   id.pat <- sequence(nn)
   n_visits <- sample(1:20, sum(nn), replace = TRUE) 
   
   id.visit <- sequence(n_visits)
   id.hosp.expanded <- rep(id.hosp, times = n_visits)
   id.pat.expanded <- rep(id.pat, times = n_visits)
   
   u_h <- rnorm(H, mean = 0, sd = sigma_u) 
   u_h_expanded <- rep(rep(u_h, times = m_hosp), times = n_visits)
   
   # v_hi <- rnorm(sum(nn), mean = 0, sd = rep(sigma_v_hosp, times = m_hosp))
   # v_hi_expanded <- rep(v_hi, times = n_visits)
   
   v_hi_mat <- matrix(0, nrow = sum(nn), ncol = py)
   for(j in 1:py) {
     v_hi_mat[, j] <- rnorm(sum(nn), mean = 0, sd = rep(sigma_v_hosp, times = m_hosp))
   }
   v_hi_expanded_mat <- v_hi_mat[rep(1:nrow(v_hi_mat), times = n_visits), ]
   
   
   X_bin <- matrix(rbinom(sum(n_visits) * p_bin, size = 1, prob = 0.3), nrow = sum(n_visits), ncol = p_bin)
   X_cont <- matrix(rnorm(sum(n_visits) * p_cont, mean = 0, sd = 1), nrow = sum(n_visits), ncol = p_cont)
   X_hij <- cbind(X_bin, X_cont) 
   
   # Multivariate Errors with Unstructured Correlation
   Sigma_err <- diag(sigma_e, py) %*% Rho_mat %*% diag(sigma_e, py)
   epsilon_hij <- mvrnorm(sum(n_visits), mu = rep(0, py), Sigma = Sigma_err)
   
   Y_hij <- matrix(0, nrow = sum(n_visits), ncol = py)
   for (j in 1:py) {
     # Y_hij[, j] <- X_hij %*% beta[, j] + u_h_expanded + v_hi_expanded + epsilon_hij[, j]
     Y_hij[, j] <- X_hij %*% beta[, j] + u_h_expanded + v_hi_expanded_mat[, j] + epsilon_hij[, j]
   }
   
   three_lvl_dat <- data.frame(
     site = factor(id.hosp.expanded),
     patient = factor(id.pat.expanded),
     visit = id.visit,
     X_hij, Y_hij
   )
   
   visit_count <- three_lvl_dat %>%
     group_by(site, patient) %>%
     summarise(total_visits = n(), .groups = "drop")
   
   rearranged_data <- merge(three_lvl_dat, visit_count, by = c("site","patient")) %>%
     arrange(site, total_visits, patient) %>%
     mutate(site = factor(site))
   
   
   X_names <- paste0("X", 1:px)
   Y_names <- paste0("Y", 1:py)
   colnames(rearranged_data)[4:(4+px-1)] <- X_names
   colnames(rearranged_data)[(4+px):(4+px+py-1)] <- Y_names
   
   # -----------------------------------------------------
   # 1.5 Generate MAR Missingness (Variable Levels per Y)
   # -----------------------------------------------------
   # Missingness depends on observed covariates X1 and X2.
   rearranged_data_miss <- rearranged_data
   
   # Define slopes (effect of covariates on missingness)
   # These can remain constant or change if you want different dependencies
   gamma_1 <- 0.5  # Effect of X1
   gamma_2 <- 0.4  # Effect of X2
   
   # Define different intercepts for Y1, Y2, Y3 to control their specific missingness levels
   # Lower value (e.g., -3)  -> Low Missingness
   # Higher value (e.g., -1) -> High Missingness
   gamma_0_vec <- c(-2.5, -1.5, -0.5) 
   
   for(j in seq_along(Y_names)){
     y_col <- Y_names[j]
     current_gamma_0 <- gamma_0_vec[j]
     
     # Calculate probability specifically for this column
     linear_pred <- current_gamma_0 + gamma_1 * rearranged_data$X1 + gamma_2 * rearranged_data$X2
     prob_miss <- 1 / (1 + exp(-linear_pred))
     
     # Apply missingness
     is_missing <- rbinom(n = nrow(rearranged_data), size = 1, prob = prob_miss)
     rearranged_data_miss[is_missing == 1, y_col] <- NA
   }
   
   # -----------------------------------------------------
   # 2. Fit Model using 'rearranged_data_miss'
   # -----------------------------------------------------
   fit_mv <- tryCatch(
     {
       peal.fit(
         data = rearranged_data_miss, 
         X_cols = X_names,
         Y_cols = Y_names,
         site_col = "site",
         patient_col = "patient",
         corstr = "unstructured",   # independent, exchangeable, unstructured
         reml = TRUE,
         verbose = FALSE,
         use_rcpp = TRUE 
       )
     },
     error = function(e) return(NULL)
   )
   
   if (is.null(fit_mv)) return(NULL)
   
   hessian_eigen <- eigen(fit_mv$opt$hessian)$values
   if (any(hessian_eigen <= 1e-6)) return(NULL) 
   
   sigma_re_est <- sqrt(fit_mv$theta * fit_mv$s2)
   sigma_e_est  <- sqrt(fit_mv$s2)
   
   C_est <- fit_mv$Corr
   rhos_est <- c(C_est[2,1], C_est[3,1], C_est[3,2]) 
   
   list(
     beta_res = as.vector(fit_mv$b),
     sigma_res = c(sigma_re_est, sigma_e_est, rhos_est)
   )
}

stopCluster(cl)

# Process Results
results <- results[!sapply(results, is.null)]
cat(sprintf("Successful simulations: %d / %d\n", length(results), N_sim))

if (length(results) > 0) {
  
  result_beta_clean <- matrix(nrow = (px*py), ncol = length(results))
  rownames(result_beta_clean) <- rownames(result_beta)
  
  result_sigma_clean <- matrix(nrow = nrow(result_sigma), ncol = length(results))
  rownames(result_sigma_clean) <- rownames(result_sigma)
  
  for (k in 1:length(results)) {
    result_beta_clean[, k] <- results[[k]]$beta_res
    result_sigma_clean[, k] <- results[[k]]$sigma_res
  }
  
  true_beta_vec <- as.vector(beta)
  true_sigma_vec <- c(sigma_u, sigma_v_hosp, sigma_e, rho_true_vec)
  
  # --- Visualization ---
  
  beta_df <- data.frame(
    Simulation = rep(1:length(results), each=nrow(result_beta_clean)),
    Parameter = rep(rownames(result_beta_clean), times=length(results)),
    Estimate = as.vector(result_beta_clean)
  )
  beta_df$True_Value <- rep(true_beta_vec, times = length(results))
  beta_df$Parameter <- factor(beta_df$Parameter, levels = paste0("Beta", 1:(px*py)))
  
  sigma_df <- data.frame(
    Simulation = rep(1:length(results), each=nrow(result_sigma_clean)),
    Parameter = rep(rownames(result_sigma_clean), times=length(results)),
    Estimate = as.vector(result_sigma_clean)
  )
  sigma_df$True_Value <- rep(true_sigma_vec, times = length(results))
  
  p1 <- beta_df %>% 
    mutate(Bias = Estimate - True_Value) %>%
    ggplot(aes(x = Parameter, y = Bias)) +
    geom_boxplot(fill = "lightblue", alpha = 0.6) +
    geom_hline(yintercept = 0, linetype="dashed", color="red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Bias in Beta (Unstructured Corr) with MAR Missingness", y = "Bias")
  
  print(p1)
  
  p2 <- sigma_df %>% 
    mutate(Bias = Estimate - True_Value) %>%
    ggplot(aes(x = Parameter, y = Bias)) +
    geom_boxplot(fill = "lightgreen", alpha = 0.6) +
    geom_hline(yintercept = 0, linetype="dashed", color="red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Bias in Variance/Corr Components (MAR Missingness)", y = "Bias")
  
  print(p2)
  
  # --- Metrics ---
  print(sigma_df %>%
          group_by(Parameter) %>%
          summarise(
            True = first(True_Value),
            Mean = mean(Estimate, na.rm=TRUE),
            Bias = Mean - True,
            MSE = mean((Estimate - True)^2, na.rm=TRUE)
          ))
}