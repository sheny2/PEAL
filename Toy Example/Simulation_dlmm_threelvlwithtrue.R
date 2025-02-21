# Load required libraries
library(lme4)
library(dplyr)
library(Matrix)
library(tidyverse)

source('dlmm_threelvlengine.R')

# Define simulation function
simulate_comparison <- function(N = 100) {

  results <- list()

  for (sim in 1:N) {
    # Parameters
    H <- 20  # Number of hospitals
    m <- 50  # Patients per hospital
    beta0 <- 5
    p_x <- 5
    beta <- c(5, 3, -2, 4, 1)  # Fixed effect coefficients for 5 covariates
    sigma_u <- 2  # Hospital-level variance
    sigma_v <- 1.5  # Patient-level variance
    sigma_e <- 1  # Error variance

    three_lvl_dat <- data.frame()

    # Generate data
    for (h in 1:H) {
      u_h <- rnorm(1, mean = 0, sd = sigma_u)  # Random hospital effect
      for (i in 1:m) {
        v_hi <- rnorm(1, mean = 0, sd = sigma_v)  # Random patient effect
        n_i <- sample(1:10, 1)  # Random number of visits for each patient
        for (j in 1:n_i) {
          X_hij <- rbinom(p_x, size = 1, prob = 0.3)  # Generate binary covariates
          epsilon_hij <- rnorm(1, mean = 0, sd = sigma_e)  # Random error
          y_hij <- beta0 + sum(X_hij * beta) + u_h + v_hi + epsilon_hij  # Outcome

          three_lvl_dat <- rbind(three_lvl_dat, c(h, i, j, X_hij, y_hij))
        }
      }
    }

    # Add column names
    colnames(three_lvl_dat) <- c("site", "patient", "visit", paste0("X", 1:p_x), "Y")

    # Preprocess data
    visit_count <- three_lvl_dat %>%
      group_by(site, patient) %>%
      summarise(total_visits = n(), .groups = "drop")

    rearranged_data <- merge(three_lvl_dat, visit_count, by = c("site", "patient")) %>%
      arrange(site, total_visits, patient)

    # Fit LMM
    fit03.lmer <- lmer(Y ~ X1 + X2 + X3 + X4 + X5 + (1 | site / patient), data = rearranged_data, REML = T,
                       control = lmerControl(optimizer="bobyqa"))

    # Prepare data for DLMM
    Y <- rearranged_data$Y
    X <- as.matrix(rearranged_data[, paste0("X", 1:p_x)])
    X <- cbind(1, X)
    count_mat <- rearranged_data %>%
      group_by(site, patient) %>%
      summarise(n_hi = n(), .groups = 'drop')
    Z <- generate_Z_matrix(count_mat, H)
    id.site <- rearranged_data$site

    # Fit DLMM
    ShXYZ <- lmm.get.summary3(Y, X, Z, id.site)
    fit03.dlmm <- lmm.fit3(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                           pooled = F, reml = T,
                           common.s2 = T,
                           ShXYZ = ShXYZ,
                           corstr = 'independence',
                           mypar.init = NULL)

    beta_diff <- c(fit03.dlmm$b - c(beta0,beta))

    dlmm_sigma_u <- sqrt(diag(fit03.dlmm$V))[1]
    dlmm_sigma_v <- sqrt(diag(fit03.dlmm$V))[2]
    dlmm_sigma_e <- sqrt(fit03.dlmm$s2)

    # c(dlmm_sigma_u, dlmm_sigma_v, dlmm_sigma_e)

    var_comp <- as.data.frame(VarCorr(fit03.lmer))
    lmm_sigma_v <- var_comp$sdcor[1]
    lmm_sigma_u <- var_comp$sdcor[2]
    lmm_sigma_e <- var_comp$sdcor[3]

    # c(lmm_sigma_u, lmm_sigma_v, lmm_sigma_e)

    sigma_diff <- c(c(dlmm_sigma_u, dlmm_sigma_v, dlmm_sigma_e) -
                      c(sigma_u, sigma_v, sigma_e))

    results[[sim]] <- list(beta_diff = beta_diff,
                         sigma_diff = sigma_diff)
  }

  return(results)
}

# Run the simulation

simulation_results <- simulate_comparison(N = 10)

p_x <- 5

# Convert results into a data frame for easier analysis
beta_diff_df <- do.call(rbind, lapply(simulation_results, function(x) x$beta_diff)) %>%
  as.data.frame() %>%
  setNames(paste0("beta_diff", 0:p_x))

sigma_diff_df <- do.call(rbind, lapply(simulation_results, function(x) x$sigma_diff)) %>%
  as.data.frame() %>%
  setNames(c("s_u_diff", "s_v_diff", "s_e_diff"))

# Combine results into a single data frame
final_results <- cbind(beta_diff_df, sigma_diff_df)
final_results


final_results %>%
  pivot_longer(1:6, names_to = "Diff", values_to = "Value") %>%
  ggplot(aes(x = Diff, y = Value)) + geom_boxplot()

final_results %>%
  pivot_longer(7:9, names_to = "Diff", values_to = "Value") %>%
  ggplot(aes(x = Diff, y = Value)) + geom_boxplot()
