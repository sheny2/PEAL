library(doParallel)
library(foreach)
library(data.table)
library(dplyr)


source("DLMM-Engine-Multi2.R")

# Define number of cores for parallel execution
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)


N = 30


# Parameters
H <- 5  # Number of sites
m_hosp <- sample(50:100, H, replace = TRUE)  # Number of patients per site
px <- 9  # Number of covariates
p_bin <- 5  # Number of binary covariates
p_cont <- px - p_bin  # Number of continuous covariates
py <- 3 # Number of outcomes (multivariate)

# Fixed effects
beta0 <- rnorm(py, 0, 1)  # Intercept for each outcome
beta <- matrix(runif(px*py, 1, 10), nrow = px, ncol = py)

# Variance components
sigma_u <- 3  # Site-level variance
sigma_v_hosp <- runif(H, min = 1, max = 5)  # Varying sigma_v by hospital

# Exchangeable correlation structure
sigma_e <- 3  # Error variance
rho <- 0.5  # Correlation between outcomes
rho_mat <- matrix(rho, nrow = py, ncol = py)  # Exchangeable correlation matrix
diag(rho_mat) <- 1  # Diagonal elements are 1


result_beta = matrix(nrow = (px+1)*py, ncol = N)
rownames(result_beta) <- c(paste0("beta1-", 0:px),
                           paste0("beta2-", 0:px),
                           paste0("beta3-", 0:px))

result_sigma = matrix(nrow = (H+2), ncol = N)
rownames(result_sigma) <- c("sigma_u", paste0("sigma_v_", 1:H), "sigma_e")

result_rho = c()

# Run simulations in parallel using foreach
results <- foreach(k = 1:N, .packages = c("data.table", "dplyr")) %dopar% {

  source("DLMM-Engine-Multi2.R")

  # Generate data
  nn <- rep(m_hosp, times = 1)  # Number of patients per hospital
  id.hosp <- rep(1:H, times = m_hosp)  # Hospital ID
  id.pat <- sequence(nn)  # Patient ID
  n_visits <- sample(1:30, sum(nn), replace = TRUE)  # Number of visits per patient

  # Expand hospital and patient IDs for visits
  id.visit <- sequence(n_visits)
  id.hosp.expanded <- rep(id.hosp, times = n_visits)
  id.pat.expanded <- rep(id.pat, times = n_visits)

  # Random effects
  u_h <- rnorm(H, mean = 0, sd = sigma_u)  # Hospital effects
  v_hi <- rnorm(sum(nn), mean = 0, sd = rep(sigma_v_hosp, times = m_hosp))  # Patient effects varying by hospital

  # Expansion of hospital effects
  u_h_patient <- rep(u_h, times = m_hosp)
  u_h_expanded <- rep(u_h_patient, times = n_visits)

  # Expansion of patient effects
  v_hi_expanded <- rep(v_hi, times = n_visits)

  # Covariates
  X_bin <- matrix(rbinom(sum(n_visits) * p_bin, size = 1, prob = 0.3), nrow = sum(n_visits), ncol = p_bin)
  X_cont <- matrix(rnorm(sum(n_visits) * p_cont, mean = 0, sd = 1), nrow = sum(n_visits), ncol = p_cont)
  X_hij <- cbind(X_bin, X_cont)  # Combine binary & continuous covariates

  # Generate multivariate errors with exchangeable correlation
  epsilon_hij <- MASS::mvrnorm(sum(n_visits), mu = rep(0, py),
                               Sigma = diag(sigma_e, py) %*% rho_mat %*% diag(sigma_e, py))

  # Compute outcomes
  Y_hij <- matrix(0, nrow = sum(n_visits), ncol = py)
  for (k in 1:py) {
    Y_hij[, k] <- beta0[k] + X_hij %*% beta[, k] + u_h_expanded + v_hi_expanded + epsilon_hij[, k]
  }

  # Create data table
  three_lvl_dat <- data.table(
    site = id.hosp.expanded,
    patient = id.pat.expanded,
    visit = id.visit,
    X_hij, Y_hij
  ) %>% data.frame()

  setnames(three_lvl_dat, c("site", "patient", "visit", paste0("X", 1:px), paste0("Y", 1:py)))

  # Preprocessing
  visit_count <- three_lvl_dat %>%
    dplyr::group_by(site, patient) %>%
    dplyr::summarise(total_visits = n(), .groups = "drop")

  # Reorder data
  rearranged_data <- merge(three_lvl_dat, visit_count, by = c("site", "patient")) %>%
    arrange(site, total_visits, patient) %>%
    mutate(site = factor(site))

  # XYZ
  Y <- as.matrix(rearranged_data[, paste0("Y", 1:py)])
  X <- as.matrix(rearranged_data[, paste0("X", 1:px)])
  X <- cbind(1, X)
  Z <- list()

  for(i in 1:H){
    count_mat = rearranged_data %>%
      filter(site == i) %>%
      group_by(site, patient) %>%
      dplyr::summarise(n_hi = n(), .groups = 'drop')

    Z[[i]] <- (generate_Zhv_matrix(count_mat))
  }

  id.site <- rearranged_data$site

  ShXYZ <- lmm.get.summary3(Y, X, Z, id.site)

  # DLMM
  fit03.dlmm = lmm.fit3(Y = NULL, X = NULL, Z = NULL,
                        id.site = NULL, weights = NULL,
                        pooled = F, reml = T,
                        common.s2 = T,
                        ShXYZ = ShXYZ,  # only need summary stats
                        # corstr = 'independence',
                        corstr = 'exchangeable',
                        mypar.init = NULL)


  # Return the estimated parameters
  list(beta_res = fit03.dlmm$b,
       rho_res = fit03.dlmm$rho,
       sigma_res = c(sqrt(fit03.dlmm$V), sqrt(fit03.dlmm$s2)))
}

# Store results into matrices
for (k in 1:N) {
  result_beta[, k] <- c(results[[k]]$beta_res)
  result_sigma[, k] <- results[[k]]$sigma_res
  result_rho[k] <-  results[[k]]$rho_res
}


stopCluster(cl)


plot(result_rho)


# Sample true parameter values
true_beta <- rbind(beta0, beta)
true_sigma <- c(sigma_u, sigma_v_hosp, sigma_e)


beta_df <- reshape2::melt(as.data.frame(result_beta))
colnames(beta_df) <- c("Simulation", "Estimate")
beta_df$Parameter <- rep(rownames(result_beta), ncol(result_beta))

sigma_df <- reshape2::melt(as.data.frame(result_sigma))
colnames(sigma_df) <- c("Simulation", "Estimate")
sigma_df$Parameter <- rep(rownames(result_sigma), ncol(result_sigma))

beta_df$True_Value <- rep(true_beta, times = ncol(result_beta))
sigma_df$True_Value <- rep(true_sigma, times = ncol(result_sigma))


saveRDS(beta_df, file = "beta_df_EXG.rds")
saveRDS(sigma_df, file = "sigma_df_EXG.rds")
beta_df <- readRDS("beta_df_EXG.rds")
sigma_df <- readRDS("sigma_df_EXG.rds")


beta_df %>% mutate(Bias = Estimate - True_Value) %>%
  # filter(! Parameter %in%  c("beta1-0", "beta2-0", "beta3-0")) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(fill = "lightblue", alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

beta_df %>% mutate(Bias = Estimate - True_Value) %>%
  filter(! Parameter %in%  c("beta1-0", "beta2-0", "beta3-0")) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(fill = "lightblue", alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma_df %>% mutate(Bias = Estimate - True_Value) %>%
  # filter(Parameter == "sigma_e") %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(fill = "lightblue", alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


