suppressPackageStartupMessages({
  library(doParallel)
  library(doRNG)
  library(foreach)
  library(data.table)
  library(dplyr)
  library(MASS)
  library(reshape2)
  library(ggplot2)
  library(Matrix)  # for nearPD if needed
})

## ------------------------------------------------
## 1. Global settings and true parameters
## ------------------------------------------------


source("PEAL_Engine_Multi-RI_Rho.R")
source("PEAL_Engine_Missing_Multi_RI_Rho_Pattern2.R")

H      <- 5    # number of sites
px     <- 6
p_bin  <- 3
p_cont <- px - p_bin
py     <- 3   # number of outcomes

X_cols <- paste0("X", 1:px)
Y_cols <- paste0("Y", 1:py)

set.seed(1)
beta_true <- matrix(runif(px * py, -3, 3), nrow = px, ncol = py)

sigma_u_true      <- 0.5
sigma_v_hosp_true <- runif(H, min = 0.5, max = 1)   # fixed across sims
sigma_e_true      <- 0.3

# Non-exchangeable correlation
rho12 <- -0.2
rho13 <- 0.35
rho23 <- 0.9
corr_mat_true <- matrix(c(
  1,     rho12, rho13,
  rho12, 1,     rho23,
  rho13, rho23, 1
), nrow = py, byrow = TRUE)

# Ensure SPD
if (min(eigen(corr_mat_true, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
  corr_mat_true <- as.matrix(nearPD(corr_mat_true, corr = TRUE)$mat)
}

sigma_e_vec <- rep(sigma_e_true, py)
De          <- diag(sigma_e_vec, py)
Sigma_e     <- De %*% corr_mat_true %*% De

# True parameter vectors for later comparison
true_beta_vec       <- c(beta_true)
true_sigmaRE_vec    <- c(sigma_u_true, sigma_v_hosp_true)  # length H+1
true_sigma_e        <- sigma_e_true
true_corr_vec       <- as.vector(corr_mat_true)
par_names_beta      <- paste0("Beta", 1:(px * py))
par_names_sigmaRE   <- c("sigma_u", paste0("sigma_v_", 1:H))
par_name_sigma_e    <- "sigma_e"
corr_param_names    <- paste0("rho_", rep(1:py, each = py), "_", rep(1:py, times = py))

## ------------------------------------------------
## 2. Helper: build Z list (site-wise random intercepts)
## ------------------------------------------------
build_Z_list <- function(dat, H) {
  Z <- vector("list", H)
  for (i in 1:H) {
    count_mat <- dat %>%
      dplyr::filter(site == i) %>%
      dplyr::group_by(site, patient) %>%
      dplyr::summarise(n_hi = n(), .groups = "drop")
    Z[[i]] <- generate_Zhv_matrix(count_mat)
  }
  Z
}

## ------------------------------------------------
## 3. Single replicate DGP: full & MAR-masked data
## ------------------------------------------------
simulate_one <- function(seed) {
  set.seed(seed)
  
  # --- DGP ---
  m_hosp <- sample(50:100, H, replace = TRUE)
  nn     <- rep(m_hosp, times = 1)
  id.hosp <- rep(seq_len(H), times = m_hosp)
  id.pat  <- sequence(nn)
  n_visits <- sample(1:30, sum(nn), replace = TRUE)
  
  id.visit         <- sequence(n_visits)
  id.hosp.expanded <- rep(id.hosp, times = n_visits)
  id.pat.expanded  <- rep(id.pat,  times = n_visits)
  
  # Random effects
  u_h  <- rnorm(H, mean = 0, sd = sigma_u_true)
  v_hi <- rnorm(sum(nn), mean = 0, sd = rep(sigma_v_hosp_true, times = m_hosp))
  
  u_h_patient   <- rep(u_h, times = m_hosp)
  u_h_expanded  <- rep(u_h_patient, times = n_visits)
  v_hi_expanded <- rep(v_hi, times = n_visits)
  
  # Covariates
  Nobs   <- sum(n_visits)
  X_bin  <- matrix(rbinom(Nobs * p_bin, 1, 0.3), nrow = Nobs, ncol = p_bin)
  X_cont <- matrix(rnorm(Nobs * p_cont, 0, 1),   nrow = Nobs, ncol = p_cont)
  X_hij  <- cbind(X_bin, X_cont)
  
  # Residuals (non-exchangeable correlation)
  epsilon_hij <- MASS::mvrnorm(n = Nobs, mu = rep(0, py), Sigma = Sigma_e)
  
  # Outcomes
  Y_hij <- matrix(0, nrow = Nobs, ncol = py)
  for (k in 1:py) {
    Y_hij[, k] <- X_hij %*% beta_true[, k] + u_h_expanded + v_hi_expanded + epsilon_hij[, k]
  }
  
  three_lvl_dat <- data.table(
    site    = id.hosp.expanded,
    patient = id.pat.expanded,
    visit   = id.visit,
    X_hij, Y_hij
  ) %>% as.data.frame()
  
  data.table::setnames(three_lvl_dat,
                       c("site","patient","visit", paste0("X",1:px), paste0("Y",1:py)))
  
  visit_count <- three_lvl_dat %>%
    group_by(site, patient) %>%
    summarise(total_visits = n(), .groups = "drop")
  
  rearranged_data <- merge(three_lvl_dat, visit_count, by = c("site","patient")) %>%
    arrange(site, total_visits, patient) %>%
    mutate(site = factor(site))
  
  ## --- MAR missingness ---
  rearranged_data_with_missing <- rearranged_data
  n <- nrow(rearranged_data)
  
  p1_base <- 0.2
  p2_base <- 0.3
  p3_base <- 0.6
  
  # Y3 depends on Y1,Y2 observed
  Y3_miss_prob <- ifelse(!is.na(rearranged_data$Y1) & !is.na(rearranged_data$Y2),
                         p3_base + 0.2,
                         p3_base - 0.2)
  Y3_miss_prob <- pmin(pmax(Y3_miss_prob, 0), 1)
  Y3_missing   <- rbinom(n, 1, Y3_miss_prob) == 1
  rearranged_data_with_missing$Y3[Y3_missing] <- NA
  
  # Y1, Y2 depend on Y3 observed
  Y1_miss_prob <- ifelse(!is.na(rearranged_data_with_missing$Y3),
                         p1_base + 0.2,
                         p1_base - 0.1)
  Y2_miss_prob <- ifelse(!is.na(rearranged_data_with_missing$Y3),
                         p2_base + 0.2,
                         p2_base - 0.1)
  Y1_miss_prob <- pmin(pmax(Y1_miss_prob, 0), 1)
  Y2_miss_prob <- pmin(pmax(Y2_miss_prob, 0), 1)
  
  Y1_missing <- rbinom(n, 1, Y1_miss_prob) == 1
  Y2_missing <- rbinom(n, 1, Y2_miss_prob) == 1
  
  rearranged_data_with_missing$Y1[Y1_missing] <- NA
  rearranged_data_with_missing$Y2[Y2_missing] <- NA
  
  # Ensure no row with all Y's missing
  all_y_missing <- apply(
    rearranged_data_with_missing[, paste0("Y",1:py), drop = FALSE],
    1, function(x) all(is.na(x))
  )
  if (any(all_y_missing)) {
    for (i in which(all_y_missing)) {
      y_to_keep <- sample(paste0("Y",1:py), 1)
      rearranged_data_with_missing[i, y_to_keep] <-
        rearranged_data[i, y_to_keep]
    }
  }
  
  list(
    full   = rearranged_data,
    miss   = rearranged_data_with_missing
  )
}

## ------------------------------------------------
## 4. Fit helpers and packers
## ------------------------------------------------

# assumes peal.fit.RI_mv and peal.fit.RI_mv_patterns are defined

pack_fit <- function(fit, corstr, py, H) {
  if (is.null(fit)) {
    return(list(
      b        = rep(NA_real_, px * py),
      sigmaRE  = rep(NA_real_, H + 1),
      sigma_e  = NA_real_,
      corr_vec = rep(NA_real_, py * py)
    ))
  }
  b_hat      <- as.numeric(fit$b)
  sigmaRE_hat <- sqrt(fit$theta * fit$s2)  # H+1 random-effect SDs
  sigmae_hat  <- sqrt(fit$s2)
  
  if (!is.null(fit$Corr)) {
    Corr_hat <- fit$Corr
    corr_vec_hat <- as.vector(Corr_hat)
  } else {
    # For independence, Corr may be NULL or identity; we treat as NA
    corr_vec_hat <- rep(NA_real_, py * py)
  }
  
  list(
    b        = b_hat,
    sigmaRE  = sigmaRE_hat,
    sigma_e  = sigmae_hat,
    corr_vec = corr_vec_hat
  )
}

fit_full_models <- function(dat_full) {
  Y <- as.matrix(dat_full[, Y_cols])
  X <- as.matrix(dat_full[, X_cols])
  Z <- build_Z_list(dat_full, H)
  id.site <- dat_full$site
  
  # independence
  fit_in <- tryCatch(
    peal.fit.RI_mv(
      Y = Y, X = X, Z = Z, id.site = id.site,
      weights = NULL, reml = TRUE,
      corstr = "independence",
      estimate_rho = FALSE,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  res_in <- pack_fit(fit_in, "independence", py, H)
  
  # exchangeable
  fit_ex <- tryCatch(
    peal.fit.RI_mv(
      Y = Y, X = X, Z = Z, id.site = id.site,
      weights = NULL, reml = TRUE,
      corstr = "exchangeable",
      estimate_rho = TRUE,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  res_ex <- pack_fit(fit_ex, "exchangeable", py, H)
  
  # unstructured
  fit_un <- tryCatch(
    peal.fit.RI_mv(
      Y = Y, X = X, Z = Z, id.site = id.site,
      weights = NULL, reml = TRUE,
      corstr = "unstructured",
      estimate_rho = TRUE,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  res_un <- pack_fit(fit_un, "unstructured", py, H)
  
  list(
    FULL_IN = res_in,
    FULL_EX = res_ex,
    FULL_UN = res_un
  )
}

fit_pattern_models <- function(dat_miss) {
  # IMPORTANT: use rearranged_data_with_missing here
  fit_in <- tryCatch(
    peal.fit.RI_mv_patterns(
      data         = dat_miss,
      X_cols       = X_cols,
      Y_cols       = Y_cols,
      site_col     = "site",
      patient_col  = "patient",
      corstr       = "independence",
      reml         = TRUE,
      verbose      = FALSE
    ),
    error = function(e) NULL
  )
  res_in <- pack_fit(fit_in, "independence", py, H)
  
  fit_ex <- tryCatch(
    peal.fit.RI_mv_patterns(
      data         = dat_miss,
      X_cols       = X_cols,
      Y_cols       = Y_cols,
      site_col     = "site",
      patient_col  = "patient",
      corstr       = "exchangeable",
      reml         = TRUE,
      verbose      = FALSE
    ),
    error = function(e) NULL
  )
  res_ex <- pack_fit(fit_ex, "exchangeable", py, H)
  
  fit_un <- tryCatch(
    peal.fit.RI_mv_patterns(
      data         = dat_miss,
      X_cols       = X_cols,
      Y_cols       = Y_cols,
      site_col     = "site",
      patient_col  = "patient",
      corstr       = "unstructured",
      reml         = TRUE,
      verbose      = FALSE
    ),
    error = function(e) NULL
  )
  res_un <- pack_fit(fit_un, "unstructured", py, H)
  
  list(
    OBS_IN = res_in,
    OBS_EX = res_ex,
    OBS_UN = res_un
  )
}

## ------------------------------------------------
## 5. One replicate wrapper
## ------------------------------------------------
run_one_rep <- function(seed) {
  sim <- simulate_one(seed)
  full_res <- fit_full_models(sim$full)
  patt_res <- fit_pattern_models(sim$miss)
  
  c(full_res, patt_res)  # named list of six models
}

## ------------------------------------------------
## 6. Parallel simulation
## ------------------------------------------------
N <- 100  # number of replicates

model_names <- c("FULL_IN","FULL_EX","FULL_UN",
                 "OBS_IN","OBS_EX","OBS_UN")

num_cores <- max(1, parallelly::availableCores())
cl <- makeCluster(num_cores)
registerDoParallel(cl)
registerDoRNG(12345)

results <- foreach(k = 1:N, .packages = c("data.table","dplyr","MASS")) %dopar% {

  source("PEAL_Engine_Multi-RI_Rho.R")
  source("PEAL_Engine_Missing_Multi_RI_Rho_Pattern2.R")
  
  run_one_rep(sample(1:1e6, 1))
}

stopCluster(cl)


# saveRDS(results, file = "Simulation6Case_Results.rds")
# results <- readRDS("Simulation6Case_Results.rds")



## ------------------------------------------------
## 7. Collect results
## ------------------------------------------------
init_mat <- function(nr) {
  matrix(NA_real_, nrow = nr, ncol = N,
         dimnames = list(NULL, paste0("sim", 1:N)))
}

store <- lapply(model_names, function(nm) {
  list(
    beta   = init_mat(px * py),
    sigmaRE = init_mat(H + 1),
    sigma_e = init_mat(1),
    corr   = init_mat(py * py)
  )
})
names(store) <- model_names

for (k in seq_along(results)) {
  resk <- results[[k]]
  for (nm in model_names) {
    store[[nm]]$beta[,   k] <- resk[[nm]]$b
    store[[nm]]$sigmaRE[,k] <- resk[[nm]]$sigmaRE
    store[[nm]]$sigma_e[,k] <- resk[[nm]]$sigma_e
    store[[nm]]$corr[,  k]  <- resk[[nm]]$corr_vec
  }
}

## ------------------------------------------------
## 8. Bias summaries
## ------------------------------------------------
mk_long <- function(mat, par_names, model_label, value_nm = "Estimate") {
  df <- reshape2::melt(as.data.frame(mat), value.name = value_nm)
  names(df) <- c("Simulation", value_nm)
  df$Parameter <- rep(par_names, ncol(mat))
  df$Model <- model_label
  df
}

### 8.1 Fixed effects beta
beta_long <- do.call(rbind, lapply(names(store), function(nm) {
  mk_long(store[[nm]]$beta, par_names_beta, nm, "Estimate")
}))
beta_long$True <- rep(true_beta_vec, times = ncol(store[[1]]$beta) * length(store))

beta_summary <- beta_long %>%
  group_by(Model, Parameter) %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias       = Mean_Est - True_Value,
    Variance   = var(Estimate, na.rm = TRUE),
    MSE        = Bias^2 + Variance,
    .groups    = "drop"
  ) %>%
  arrange(Parameter, Model)

### 8.2 Random effects SDs
sigmaRE_long <- do.call(rbind, lapply(names(store), function(nm) {
  mk_long(store[[nm]]$sigmaRE, par_names_sigmaRE, nm, "Estimate")
}))
sigmaRE_long$True <- rep(true_sigmaRE_vec,
                         times = ncol(store[[1]]$sigmaRE) * length(store))

sigmaRE_summary <- sigmaRE_long %>%
  group_by(Model, Parameter) %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias       = Mean_Est - True_Value,
    Variance   = var(Estimate, na.rm = TRUE),
    MSE        = Bias^2 + Variance,
    .groups    = "drop"
  ) %>%
  arrange(Parameter, Model)

### 8.3 Residual SD
sigmae_long <- do.call(rbind, lapply(names(store), function(nm) {
  mk_long(store[[nm]]$sigma_e, par_name_sigma_e, nm, "Estimate")
}))
sigmae_long$True <- rep(true_sigma_e,
                        times = ncol(store[[1]]$sigma_e) * length(store))

sigmae_summary <- sigmae_long %>%
  group_by(Model, Parameter) %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias       = Mean_Est - True_Value,
    Variance   = var(Estimate, na.rm = TRUE),
    MSE        = Bias^2 + Variance,
    .groups    = "drop"
  ) %>%
  arrange(Parameter, Model)

### 8.4 Correlation parameters (full 3x3; you can later filter off-diagonals)
corr_long <- do.call(rbind, lapply(names(store), function(nm) {
  mk_long(store[[nm]]$corr, corr_param_names, nm, "Estimate")
}))
corr_long$True <- rep(true_corr_vec,
                      times = ncol(store[[1]]$corr) * length(store))

corr_summary <- corr_long %>%
  group_by(Model, Parameter) %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias       = Mean_Est - True_Value,
    Variance   = var(Estimate, na.rm = TRUE),
    MSE        = Bias^2 + Variance,
    .groups    = "drop"
  ) %>%
  arrange(Parameter, Model)

## ------------------------------------------------
## 9. Example plots (fixed effects bias, variance comp bias)
## ------------------------------------------------

# Fixed-effects bias (jitter + boxplot)
beta_long$Parameter <- factor(beta_long$Parameter,
                              levels = paste0("Beta", 1:(px*py)))

p_beta_bias <- beta_long %>%
  dplyr::filter(Model %in% c("FULL_UN","OBS_UN")) %>%
  mutate(Bias = Estimate - True) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1, width = 0.15, height = 0, size = 0.3) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~ Model, ncol = 3, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.7) +
  ggtitle("Fixed Effects Bias")
print(p_beta_bias)


# Random-effect SD bias
p_sigmaRE_bias <- sigmaRE_long %>%
  dplyr::filter(Model %in% c("FULL_UN","OBS_UN")) %>%
  mutate(Bias = Estimate - True) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1, width = 0.15, height = 0, size = 0.3) +
  geom_boxplot(fill = "lightgreen", alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~ Model, ncol = 3, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.7) +
  ggtitle("Random-effect SD Bias")
print(p_sigmaRE_bias)

# Residual SD bias
p_sigmae_bias <- sigmae_long %>%
  # dplyr::filter(Model %in% c("FULL_UN","OBS_UN")) %>%
  mutate(Bias = Estimate - True) %>%
  ggplot(aes(x = Model, y = Bias)) +
  geom_boxplot(fill = "orange", alpha = 0.6, outlier.shape = NA) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.7) +
  ggtitle("Residual SD Bias")
print(p_sigmae_bias)

# Correlation bias: filter off-diagonal if desired
corr_long_offdiag <- corr_long %>%
  mutate(i = as.integer(substr(Parameter, 5, 5)),
         j = as.integer(substr(Parameter, 7, 7))) %>%
  filter(i != j)
corr_long_offdiag$Parameter <- factor(corr_long_offdiag$Parameter,
                                      levels = unique(corr_long_offdiag$Parameter))


p_corr_bias <- corr_long_offdiag %>%
  filter(!Model %in% c("FULL_IN", "OBS_IN"))%>% 
  filter(!Model %in% c("FULL_EX", "OBS_EX"))%>%
  dplyr::filter(Model %in% c("FULL_UN","OBS_UN")) %>%
  filter(Parameter %in% unique(corr_long_offdiag$Parameter)[c(1,2,4)]) %>%
  mutate(Bias = Estimate - True) %>%
  ggplot(aes(x = Parameter, y = Bias, fill = Model)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.7) +
  ggtitle("Correlation Parameter Bias (off-diagonals)")
print(p_corr_bias)

gridExtra::grid.arrange(p_beta_bias, p_sigmae_bias, p_sigmaRE_bias, p_corr_bias)

# Variance & Relative Efficiency comparisons
beta_long %>%
  filter(Model %in% c("OBS_EX", "OBS_IN", "OBS_UN")) %>%
  group_by(Model, Parameter) %>%
  summarise(
    Variance = var(Estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Model,
    values_from = Variance
  ) %>%
  mutate(
    RE_EX_vs_UN = OBS_EX / OBS_UN,
    RE_IN_vs_UN = OBS_IN / OBS_UN
  ) %>%
  arrange(Parameter)

sigmaRE_long %>%
  filter(Model %in% c("OBS_EX", "OBS_IN", "OBS_UN")) %>%
  group_by(Model, Parameter) %>%
  summarise(
    Variance = var(Estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Model,
    values_from = Variance
  ) %>%
  mutate(
    RE_EX_vs_UN = OBS_EX / OBS_UN,
    RE_IN_vs_UN = OBS_IN / OBS_UN
  ) %>%
  arrange(Parameter)

sigmae_long %>%
  filter(Model %in% c("OBS_EX", "OBS_IN", "OBS_UN")) %>%
  group_by(Model, Parameter) %>%
  summarise(
    Variance = var(Estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Model,
    values_from = Variance
  ) %>%
  mutate(
    RE_EX_vs_UN = OBS_EX / OBS_UN,
    RE_IN_vs_UN = OBS_IN / OBS_UN
  ) %>%
  arrange(Parameter)


# Numeric summaries:
beta_summary
sigmaRE_summary
sigmae_summary
corr_summary
