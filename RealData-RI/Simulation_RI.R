library(doParallel)
library(foreach)
library(data.table)
library(dplyr)

# source("DLMM_Engine_RI_BFGS-Backup.R")
source("DLMM_Engine_RI_BFGS.R")

# Define number of cores for parallel execution
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

N = 50


# Parameters
H <- 5  # of sites
m_hosp <- sample(50:70, H) # of patients (1k to 3k later)

px <- 9  # of covariates
p_bin <- 5  # of binary X
p_cont <- px - p_bin  # of continuous X


beta0 = 5 # intercept
beta <- c(2, 4, 6, 8, 10, 3, 5, 7, 9)  # Fixed effects for covariates


sigma_e <- 5  # error variance
sigma_u <- 5 # site-level variance
sigma_v_hosp <- runif(H, min = 1, max = 5)  # Varying sigma_v by hospital


result_beta = matrix(nrow = (px+1), ncol = N)
rownames(result_beta) <- paste0("Beta", 0:px)
result_sigma = matrix(nrow = (H+2), ncol = N)
rownames(result_sigma) <- c("sigma_u", paste0("sigma_v_", 1:H), "sigma_e")



# Run simulations in parallel using foreach
results <- foreach(k = 1:N, .packages = c("data.table", "dplyr")) %dopar% {
# for(k in 1:N) {

  # source("DLMM_Engine_RI_BFGS-Backup.R")
  source("DLMM_Engine_RI_BFGS.R")

  # Generate data
  nn <- rep(m_hosp, times = 1)  # Number of patients per hospital
  id.hosp <- rep(1:H, times = m_hosp)   # Hospital ID
  id.pat <- sequence(nn)                # Patient ID
  n_visits <- sample(1:30, sum(nn), replace = TRUE)  # number of visits of patients

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

  # covariates
  X_bin <- matrix(rbinom(sum(n_visits) * p_bin, size = 1, prob = 0.5), nrow = sum(n_visits), ncol = p_bin)
  X_cont <- matrix(rnorm(sum(n_visits) * p_cont, mean = 0, sd = 1), nrow = sum(n_visits), ncol = p_cont)
  X_hij <- cbind(X_bin, X_cont)  # Combine binary & continuous covariates

  epsilon_hij <- rnorm(sum(n_visits), mean = 0, sd = sigma_e)

  # Compute outcome
  y_hij <- beta0 + X_hij %*% beta + u_h_expanded + v_hi_expanded + epsilon_hij


  three_lvl_dat <- data.table(
    site = id.hosp.expanded,
    patient = id.pat.expanded,
    visit = id.visit,
    X_hij,
    Y = y_hij
  ) %>% data.frame()

  # Assign column names for covariates
  setnames(three_lvl_dat, c("site", "patient", "visit", paste0("X", 1:px), "Y"))

  # Calculate the number of visits per patient per site
  visit_count <- three_lvl_dat %>%
    dplyr::group_by(site, patient) %>%
    dplyr::summarise(total_visits = n(), .groups = "drop")

  # Reorder data (as a preprocess probably)
  rearranged_data <- merge(three_lvl_dat, visit_count, by = c("site", "patient")) %>%
    arrange(site, total_visits, patient) %>%
    mutate(site = factor(site))

  # XYZ
  Y <- rearranged_data$Y
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
                        corstr = 'independence',
                        mypar.init = NULL)


  # Return the estimated parameters
  list(beta_res = fit03.dlmm$b, sigma_res = c(sqrt(fit03.dlmm$V), sqrt(fit03.dlmm$s2)))
}

# Store results into matrices
for (k in 1:N) {
  result_beta[, k] <- results[[k]]$beta_res
  result_sigma[, k] <- results[[k]]$sigma_res
}


stopCluster(cl)





# Sample true parameter values
true_beta <- c(beta0,beta)
true_sigma <- c(sigma_u, sigma_v_hosp, sigma_e)


beta_df <- reshape2::melt(as.data.frame(result_beta))
colnames(beta_df) <- c("Simulation", "Estimate")
beta_df$Parameter <- rep(rownames(result_beta), ncol(result_beta))

sigma_df <- reshape2::melt(as.data.frame(result_sigma))
colnames(sigma_df) <- c("Simulation", "Estimate")
sigma_df$Parameter <- rep(rownames(result_sigma), ncol(result_sigma))

beta_df$True_Value <- rep(true_beta, times = ncol(result_beta))
sigma_df$True_Value <- rep(true_sigma, times = ncol(result_sigma))



# saveRDS(beta_df, file = "beta_df_large.rds")
# saveRDS(sigma_df, file = "sigma_df_large.rds")
# beta_df <- readRDS("beta_df_large.rds")
# sigma_df <- readRDS("sigma_df_large.rds")
# saveRDS(beta_df, "beta_df_wrong.rds")
# saveRDS(sigma_df, "sigma_df_wrong.rds")

# beta_df <- readRDS("beta_df.rds") # 70-100
# sigma_df <- readRDS("sigma_df.rds") # 70-100
#
# beta_df <- readRDS("beta_df_wrong.rds") # 200 to 300, old code
# sigma_df <- readRDS("sigma_df_wrong.rds") # 200 to 300, old code
#
# beta_df <- readRDS("beta_df_large.rds")
# sigma_df <- readRDS("sigma_df_large.rds")


beta_df %>% mutate(Bias = Estimate - True_Value) %>%
ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(fill = "lightblue", alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


beta_df %>% mutate(Bias = Estimate - True_Value) %>% filter(Parameter!= "Beta0") %>%
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

