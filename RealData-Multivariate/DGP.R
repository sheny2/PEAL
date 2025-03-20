library(data.table)
library(dplyr)
library(MASS)  # For multivariate normal distribution


##############################
# Try multivariate outcome with exchangeable correlation structure

# Parameters
H <- 5  # Number of sites
m_hosp <- sample(100:200, H, replace = TRUE)  # Number of patients per site
px <- 9  # Number of covariates
p_bin <- 5  # Number of binary covariates
p_cont <- px - p_bin  # Number of continuous covariates
py <- 3  # Number of outcomes (multivariate)

# Fixed effects
beta0 <- rep(5, py)  # Intercept for each outcome
beta <- matrix(c(2, 4, 6, 8, 10, 3, 5, 7, 9), nrow = px, ncol = py)  # Fixed effects for covariates

# Variance components
sigma_e <- 1  # Error variance
sigma_u <- 3  # Site-level variance
sigma_v_hosp <- runif(H, min = 1, max = 5)  # Varying sigma_v by hospital

# Exchangeable correlation structure
rho <- 0.5  # Correlation between outcomes
rho_mat <- matrix(rho, nrow = py, ncol = py)  # Exchangeable correlation matrix
diag(rho_mat) <- 1  # Diagonal elements are 1
sigma <- 1

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
epsilon_hij <- mvrnorm(sum(n_visits), mu = rep(0, py), Sigma = rho_mat * sigma^2)

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
  X_hij,
  Y1 = Y_hij[, 1],
  Y2 = Y_hij[, 2],
  Y3 = Y_hij[, 3]
) %>% data.frame()

setnames(three_lvl_dat, c("site", "patient", "visit", paste0("X", 1:px), "Y1", "Y2", "Y3"))

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




# Load necessary library for multivariate normal simulation
library(MASS)

# Set seed for reproducibility
set.seed(123)

# Simulation settings
n <- 100    # number of subjects
q <- 2      # number of outcomes per subject

# True fixed effects (for intercept-only model, each outcome has its own intercept)
beta <- c(1, 2)

# True covariance matrices
G <- matrix(c(0.5, 0.2,
              0.2, 0.5), nrow = q)   # Covariance for random effects (subject-level)
R <- matrix(c(1, 0.3,
              0.3, 1), nrow = q)     # Residual covariance

# Simulate subject-level random effects and residual errors
b <- mvrnorm(n, mu = rep(0, q), Sigma = G)
epsilon <- mvrnorm(n, mu = rep(0, q), Sigma = R)

# Generate the response matrix Y (each row is a subject's q-dimensional outcome)
Y <- matrix(0, n, q)
for(i in 1:n) {
  Y[i, ] <- beta + b[i, ] + epsilon[i, ]
}

# Define a function to calculate the log likelihood of the multivariate LMM
logLik_MLMM <- function(Y, beta, G, R) {
  n <- nrow(Y)      # number of subjects
  q <- ncol(Y)      # number of outcomes per subject
  V <- G + R        # marginal covariance for each subject
  # Calculate the log-determinant of V
  logdetV <- as.numeric(determinant(V, logarithm = TRUE)$modulus)
  # Inverse of V
  invV <- solve(V)

  ll <- 0
  # Loop over subjects and sum the log likelihood contributions
  for(i in 1:n) {
    resid <- as.vector(Y[i, ] - beta)
    ll <- ll - 0.5 * ( q * log(2*pi) + logdetV + t(resid) %*% invV %*% resid )
  }
  return(ll)
}

# Compute the log likelihood using the true parameters
ll_value <- logLik_MLMM(Y, beta, G, R)
cat("Log-likelihood:", ll_value, "\n")
