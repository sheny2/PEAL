# Parameters
H <- 5  # of sites
m_hosp <- sample(100:200, H, replace = T) # of patients (1k to 3k later)


px <- 9  # of covariates
p_bin <- 5  # of binary X
p_cont <- px - p_bin  # of continuous X

beta <- c(2, 4, 6, 8, 10, 3, 5, 7, 9)  # Fixed effects for covariates

sigma_e <- 2  # error variance
sigma_u <- 3 # site-level variance
sigma_v_hosp <- runif(H, min = 1, max = 5)  # Varying sigma_v by hospital


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
X_bin <- matrix(rbinom(sum(n_visits) * p_bin, size = 1, prob = 0.3), nrow = sum(n_visits), ncol = p_bin)
X_cont <- matrix(rnorm(sum(n_visits) * p_cont, mean = 0, sd = 5), nrow = sum(n_visits), ncol = p_cont)
X_hij <- cbind(X_bin, X_cont)  # Combine binary & continuous covariates

epsilon_hij <- rnorm(sum(n_visits), mean = 0, sd = sigma_e)

# Compute outcome
y_hij <- X_hij %*% beta + u_h_expanded + v_hi_expanded + epsilon_hij

three_lvl_dat <- data.table(
  site = id.hosp.expanded,
  patient = id.pat.expanded,
  visit = id.visit,
  X_hij,
  Y = y_hij
) %>% data.frame()

setnames(three_lvl_dat, c("site", "patient", "visit", paste0("X", 1:px), "Y"))


library(Matrix)
library(data.table)
library(lme4) # For comparison only

three_level_lmm <- function(data, max_iter = 100, tol = 1e-6) {
  # Prepare data
  y <- data$Y
  X <- as.matrix(data[, paste0("X", 1:px)])
  site <- as.integer(factor(data$site))
  patient <- as.integer(factor(data$patient))

  H <- length(unique(site))  # Number of sites
  n_patients <- length(unique(interaction(site, patient)))  # Total patients
  n_obs <- length(y)  # Total observations

  # Create sparse matrices for random effects
  Z_site <- sparse.model.matrix(~ 0 + factor(site))
  Z_patient <- sparse.model.matrix(~ 0 + interaction(factor(site), factor(patient)))

  # Negative log-likelihood function for optim
  neg_loglik <- function(params) {
    beta <- params[1:px]
    sigma_u <- exp(params[px+1])  # Using exp to ensure positive
    sigma_v_hosp <- exp(params[(px+2):(px+1+H)])
    sigma_e <- exp(params[px+2+H])

    # Construct covariance matrices
    G_site <- sigma_u^2 * Diagonal(H)
    G_patient <- Diagonal(n_patients, rep(sigma_v_hosp^2, as.numeric(table(site)[unique(site)])))
    R <- sigma_e^2 * Diagonal(n_obs)
    Z <- cbind(Z_site, Z_patient)
    G <- bdiag(G_site, G_patient)

    V <- Z %*% G %*% t(Z) + R
    V_chol <- chol(V)

    # Compute log-likelihood
    residuals <- y - X %*% beta
    logdet <- 2 * sum(log(diag(V_chol)))
    quad_form <- crossprod(solve(V_chol, residuals))

    0.5 * (logdet + quad_form + n_obs * log(2*pi))
  }

  # Initial parameter estimates
  init_beta <- solve(crossprod(X), crossprod(X, y))
  init_sigma_u <- sd(tapply(y, site, mean))
  init_sigma_v_hosp <- rep(1, H)
  init_sigma_e <- sd(y - X %*% init_beta)

  # Use optim for estimation
  opt <- optim(
    par = c(init_beta, log(init_sigma_u), log(init_sigma_v_hosp), log(init_sigma_e)),
    fn = neg_loglik,
    method = "L-BFGS-B",
    control = list(maxit = max_iter, factr = tol)
  )

  # Extract parameters
  beta <- opt$par[1:px]
  sigma_u <- exp(opt$par[px+1])
  sigma_v_hosp <- exp(opt$par[(px+2):(px+1+H)])
  sigma_e <- exp(opt$par[px+2+H])

  # Prepare results
  list(
    fixed_effects = beta,
    sigma_u = sigma_u,
    sigma_v_hosp = sigma_v_hosp,
    sigma_e = sigma_e,
    logLik = -opt$value,
    convergence = opt$convergence,
    message = opt$message
  )
}

# Run the model
result <- three_level_lmm(three_lvl_dat)

# Compare with lme4 (for validation)
lme4_fit <- lmer(Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                   (1 | site) + (1 | site:patient),
                 data = three_lvl_dat)

# Print results
cat("Custom Implementation Results:\n")
print(result[1:4])

cat("\nlme4 Results:\n")
print(summary(lme4_fit)$varcor)
print(fixef(lme4_fit))
