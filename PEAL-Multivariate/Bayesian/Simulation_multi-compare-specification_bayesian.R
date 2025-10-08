library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(MASS)
library(reshape2)
library(doFuture)
library(progressr)

source("PEAL_Engine_Multi-RI.R")

quiet_stan <- function(...) {
  con <- file(nullfile(), open = "wt")
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    while (sink.number() > 0) try(sink(), silent = TRUE)
    try(close(con), silent = TRUE)
  }, add = TRUE)
  
  sink(con); sink(con, type = "message")  # silence stdout & messages
  suppressWarnings(                      # silence warnings
    suppressMessages(
      rstan::stan(..., refresh = 0, open_progress = FALSE, verbose = FALSE)
    )
  )
}

# -------------------------
# Parallel setup
# -------------------------
N <- 10
num_cores <- detectCores() / 2
cl <- makeCluster(num_cores)
# registerDoParallel(cl)
registerDoFuture()
plan(cluster, workers = cl)   # reuse your existing cluster 'cl

# -------------------------
# DGP Parameters
# -------------------------
H <- 3
m_hosp <- sample(10:30, H, replace = TRUE)
px <- 6
p_bin <- 3
p_cont <- px - p_bin
py <- 3

# Fixed effects
beta <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)

# Variance components (SDs)
sigma_u <- 0.3
sigma_v_hosp <- c(0.61, 0.65, 0.69)


# Residual SD and exchangeable correlation
sigma_e <- 0.5
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


result_beta_exch_bayeisan  <- matrix(NA_real_, nrow = length(par_names_beta),  ncol = N,
                            dimnames = list(par_names_beta, paste0("sim", 1:N)))
result_sigma_exch_bayeisan <- matrix(NA_real_, nrow = length(par_names_sigma), ncol = N,
                            dimnames = list(par_names_sigma, paste0("sim", 1:N)))

result_beta_indep_bayeisan <- matrix(NA_real_, nrow = length(par_names_beta),  ncol = N,
                             dimnames = list(par_names_beta, paste0("sim", 1:N)))
result_sigma_indep_bayeisan <- matrix(NA_real_, nrow = length(par_names_sigma), ncol = N,
                             dimnames = list(par_names_sigma, paste0("sim", 1:N)))

# -------------------------
# Simulation loop: fit BOTH models per replicate
# -------------------------
results <- progressr::with_progress({
  p <- progressr::progressor(steps = N)
  foreach(k = 1:N,
          .packages = c("data.table","dplyr","MASS","rstan"),
          .options.future = list(seed = TRUE)) %dopar% {
            
            p(sprintf("iteration %d/%d", k, N))  # <-- progress line
  
  source("PEAL_Engine_Multi-RI.R")
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
  
  # ---- Build Y, X, Z
  Y <- as.matrix(rearranged_data[, paste0("Y", 1:py)])
  X <- as.matrix(rearranged_data[, paste0("X", 1:px)])
  X <- cbind(X)
  Z <- vector("list", H)
  
  for (i in 1:H) {
    count_mat <- rearranged_data |>
      filter(site == i) |>
      group_by(site, patient) |>
      dplyr::summarise(n_hi = n(), .groups = 'drop')
    Z[[i]] <- generate_Zhv_matrix(count_mat)
  }
  
  id.site <- rearranged_data$site
  
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
    peal.fit.RI_mv(
      Y = Y, X = X, Z = Z, id.site = id.site,
      pooled = FALSE, weights = NULL, reml = TRUE,
      corstr = "exchangeable",
      estimate_rho = TRUE, rho_init = 0.1,
      verbose = FALSE
    ),
    error = function(e) { NULL },
    warning = function(w) { invokeRestart("muffleWarning") }
  )
  
  out_exch <- pack_fit(fit_exch, "exchangeable")
  
  # ---- Fit Independence (rho fixed 0)
  fit_indep <- tryCatch(
    peal.fit.RI_mv(
      Y = Y, X = X, Z = Z, id.site = id.site,
      pooled = FALSE, weights = NULL, reml = TRUE,
      corstr = "independence",
      estimate_rho = FALSE,
      verbose = FALSE
    ),
    error = function(e) { NULL },
    warning = function(w) { invokeRestart("muffleWarning") }
  )
  
  out_indep <- pack_fit(fit_indep, "independence")
  
  
  
  # Baeysian
  dat <- rearranged_data %>%
    mutate(site = as.integer(factor(site)),
           pat  = as.integer(interaction(site, patient, drop = TRUE)))
  
  N   <- nrow(dat)
  H   <- length(unique(dat$site))
  P   <- length(unique(dat$pat))
  R   <- py
  
  site_id <- dat$site                            # length N, in {1,...,H}
  pat_id  <- dat$pat                             # length N, in {1,...,P}
  
  # Map each patient to its site (needed to allow site-specific patient SDs)
  pat_site <- dat %>%
    distinct(pat = pat_id, site = site_id) %>%
    arrange(pat) %>%
    pull(site)
  
  stan_data <- list(
    N = N, H = H, P = P, R = R, p = px,
    y = Y, X = X, site_id = site_id, pat_id = pat_id, pat_site = pat_site
  )
  
  
  options(mc.cores = parallelly::availableCores())
  rstan_options(auto_write = TRUE)
  
  fit_cs <- quiet_stan(
    file = "mv_lmm_cs.stan",
    data = stan_data,
    chains = 4, iter = 500, seed = 123,
    refresh = 0
  )
  
  # Extract draws and produce tidy summaries:
  post <- rstan::extract(fit_cs)
  
  out_exch_bayes <- list(
    beta_mean = as.vector(apply(post$beta, R, colMeans)),
    sigma_u_hat = mean(post$sigma_u),
    sigma_vh_hat = apply(post$sigma_v_h, 2, mean),
    sigma_e_hat = mean(post$sigma_e),
    rho_hat = mean(post$rho)
  )
  
  
  
  # bayesian indep
  fit_cs <- quiet_stan(
    file = "mv_lmm_cs_indep.stan",
    data = stan_data,
    chains = 4, iter = 500, seed = 123,
    refresh = 0
  )
  
  # Extract draws and produce tidy summaries:
  post <- rstan::extract(fit_cs)
  
  out_indep_bayes <- list(
    beta_mean = as.vector(apply(post$beta, R, colMeans)),
    sigma_u_hat = mean(post$sigma_u),
    sigma_vh_hat = apply(post$sigma_v_h, 2, mean),
    sigma_e_hat = mean(post$sigma_e),
    rho_hat = 0
  )

  
  list(beta_exch  = out_exch$beta,
       sigma_exch = out_exch$sigma,
       beta_indep  = out_indep$beta,
       sigma_indep = out_indep$sigma,
       beta_exch_bayes = out_exch_bayes$beta_mean,
       sigma_exch_bayes = c(out_exch_bayes$sigma_u_hat, out_exch_bayes$sigma_vh_hat, 
                            out_exch_bayes$sigma_e_hat, out_exch_bayes$rho_hat),
       beta_indep_bayes = out_indep_bayes$beta_mean,
       sigma_indep_bayes = c(out_indep_bayes$sigma_u_hat, out_indep_bayes$sigma_vh_hat, 
                            out_indep_bayes$sigma_e_hat, out_indep_bayes$rho_hat)
       )
          }
})



# -------------------------
# Store into pre-allocated matrices
# -------------------------
for (k in seq_along(results)) {
  result_beta_exch[,  k] <- results[[k]]$beta_exch
  result_sigma_exch[, k] <- results[[k]]$sigma_exch
  result_beta_indep[,  k] <- results[[k]]$beta_indep
  result_sigma_indep[, k] <- results[[k]]$sigma_indep
  result_beta_exch_bayeisan[,  k] <- results[[k]]$beta_exch_bayes
  result_sigma_exch_bayeisan[, k] <- results[[k]]$sigma_exch_bayes
  result_beta_indep_bayeisan[,  k] <- results[[k]]$beta_indep_bayes
  result_sigma_indep_bayeisan[, k] <- results[[k]]$sigma_indep_bayes
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
beta_exch_bayes_df <- mk_long(result_beta_exch_bayeisan, par_names_beta,  "Exchangeable-B")
beta_indep_bayes_df <- mk_long(result_beta_indep_bayeisan, par_names_beta,  "Independence-B")
beta_df <- rbind(beta_exch_df, beta_indep_df, beta_exch_bayes_df, beta_indep_bayes_df)
beta_df$True_Value <- rep(true_beta, times = ncol(result_beta_exch) + ncol(result_beta_indep) + 
                            ncol(result_beta_exch_bayeisan) + ncol(result_beta_indep_bayeisan))

sigma_exch_df  <- mk_long(result_sigma_exch,  par_names_sigma, "Exchangeable")
sigma_indep_df <- mk_long(result_sigma_indep, par_names_sigma, "Independence")
sigma_exch_bayes_df  <- mk_long(result_sigma_exch_bayeisan,  par_names_sigma, "Exchangeable-B")
sigma_indep_bayes_df  <- mk_long(result_sigma_indep_bayeisan,  par_names_sigma, "Independence-B")
sigma_df <- rbind(sigma_exch_df, sigma_indep_df, sigma_exch_bayes_df, sigma_indep_bayes_df)
sigma_df$True_Value <- rep(true_sigma, times = ncol(result_sigma_exch) + ncol(result_sigma_indep) + 
                             ncol(result_sigma_exch_bayeisan) + ncol(result_sigma_indep_bayeisan))

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



ggplot(beta_summary, aes(x = Parameter, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(sigma_summary, aes(x = Parameter, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Bias Comparison by Model and Parameter",
       x = "Parameter",
       y = "Bias")

