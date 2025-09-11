library(tidyverse)
library(lme4)
library(MASS)
library(data.table)
library(Matrix)
source("PEAL_Engine_Multi.R")

# Parameters
H <- 3  # Number of sites
m_hosp <- sample(50:100, H, replace = TRUE)  # Number of patients per site
px <- 6  # Number of covariates
p_bin <- 3  # Number of binary covariates
p_cont <- px - p_bin  # Number of continuous covariates
py <- 3 # Number of outcomes (multivariate)

# Fixed effects
beta <- matrix(runif(px*py, 1, 10), nrow = px, ncol = py)  # Fixed effects for covariates

# Variance components
sigma_u <- 0.2  # Site-level variance
sigma_v_hosp <- runif(H, min = 0.5, max = 0.8)  # Varying sigma_v by hospital

# Exchangeable correlation structure
sigma_e <- 0.3  # Error variance
rho <- 0.6 # Correlation between outcomes

# Simulation function
run_simulation <- function(iter) {
  # Generate data
  nn <- rep(m_hosp, times = 1)  # Number of patients per hospital
  id.hosp <- rep(1:H, times = m_hosp)  # Hospital ID
  id.pat <- sequence(nn)  # Patient ID
  n_visits <- sample(1:20, sum(nn), replace = TRUE)  # Number of visits per patient

  # Expand hospital and patient IDs for visits
  id.visit <- sequence(n_visits)
  id.hosp.expanded <- rep(id.hosp, times = n_visits)
  id.pat.expanded <- rep(id.pat, times = n_visits)

  # Random effects
  u_h <- rnorm(H, mean = 0, sd = sigma_u)  # Hospital effects
  v_hi <- rnorm(sum(nn), mean = 0, sd = rep(sigma_v_hosp, times = m_hosp))  # Patient effects varying by hospital

  # Expansion of effects
  u_h_patient <- rep(u_h, times = m_hosp)
  u_h_expanded <- rep(u_h_patient, times = n_visits)
  v_hi_expanded <- rep(v_hi, times = n_visits)

  # Covariates
  X_bin <- matrix(rbinom(sum(n_visits) * p_bin, size = 1, prob = 0.3), nrow = sum(n_visits), ncol = p_bin)
  X_cont <- matrix(rnorm(sum(n_visits) * p_cont, mean = 0, sd = 1), nrow = sum(n_visits), ncol = p_cont)
  X_hij <- cbind(X_bin, X_cont)

  # Generate multivariate errors
  rho_mat <- matrix(rho, nrow = py, ncol = py)
  diag(rho_mat) <- 1
  epsilon_hij <- mvrnorm(sum(n_visits), mu = rep(0, py),
                         Sigma = diag(sigma_e, py) %*% rho_mat %*% diag(sigma_e, py))

  # Compute outcomes
  Y_hij <- matrix(0, nrow = sum(n_visits), ncol = py)
  for (k in 1:py) {
    Y_hij[, k] <- X_hij %*% beta[, k] + u_h_expanded + v_hi_expanded + epsilon_hij[, k]
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
  Z <- list()

  for(i in 1:H){
    count_mat = rearranged_data %>%
      filter(site == i) %>%
      group_by(site, patient) %>%
      dplyr::summarise(n_hi = n(), .groups = 'drop')

    Z[[i]] <- (generate_Zhv_matrix(count_mat))
  }

  id.site <- rearranged_data$site

  # Run model
  res <- try(federated_lmm_two_stage(Y, X, Z, id.site))

  if(inherits(res, "try-error")) {
    return(NULL)
  } else {
    return(list(
      B = res$B,
      sigma_u = res$sigma_u,
      sigma_v = res$sigma_v,
      sigma = res$sigma,
      rho = res$rho
    ))
  }
}

# Parallel simulation
library(parallel)
library(foreach)
library(doParallel)

n_sim <- 100
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

sim_results <- foreach(i = 1:n_sim, .packages = c("MASS", "data.table", "dplyr", "Matrix"),
                       .combine = 'rbind') %dopar% {
                         source("PEAL_Engine_Multi.R")
                         run_simulation(i)
                       }

stopCluster(cl)


# Suppose your sim_results is a named list: result.1, result.2, ...
# True values
beta_true     <- beta           # 6x3 matrix
sigma_u_true  <- sigma_u        # scalar
sigma_v_true  <- sigma_v_hosp   # length 3
sigma_true    <- sigma_e        # scalar
rho_true      <- rho            # scalar

# Storage
beta_est <- list()
other_est <- data.frame(
  sigma_u = numeric(),
  sigma_v1 = numeric(),
  sigma_v2 = numeric(),
  sigma_v3 = numeric(),
  sigma = numeric(),
  rho = numeric()
)



# Loop over results
for (nm in 1:dim(sim_results)[1]) {
  res <- sim_results[nm,]

  # Flatten B and store
  beta_est[[nm]] <- as.vector(res$B) - as.vector(beta_true)

  # Store other parameters
  other_est <- rbind(other_est, data.frame(
    sigma_u = res$sigma_u - sigma_u_true,
    sigma_v1 = res$sigma_v[1] - sigma_v_true[1],
    sigma_v2 = res$sigma_v[2] - sigma_v_true[2],
    sigma_v3 = res$sigma_v[3] - sigma_v_true[3],
    sigma   = res$sigma   - sigma_true,
    rho     = res$rho     - rho_true
  ))
}

# Convert beta_est to data frame for ggplot
beta_df <- data.frame(
  bias = unlist(beta_est),
  beta_param = rep(paste0("beta", 1:(length(beta_true))), length(beta_est))
)



# Boxplot for beta bias
p1 <- ggplot(beta_df %>% mutate(beta_param =factor(beta_param, levels = paste0("beta", 1:(px*py)))),
             aes(x = beta_param, y = bias)) +
  geom_boxplot(fill = "#8DA0CB") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Beta parameters", y = "Bias", title = "Bias of Beta Estimates") +
  theme_minimal()

# Boxplot for other parameters bias
other_df <- melt(other_est, variable.name = "param", value.name = "bias")

p2 <- ggplot(other_df, aes(x = param, y = bias)) +
  geom_boxplot(fill = "#FC8D62") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Parameter", y = "Bias", title = "Bias of Variance/Correlation Parameters") +
  theme_minimal()

p1
p2

