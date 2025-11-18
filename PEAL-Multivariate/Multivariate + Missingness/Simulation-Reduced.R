# ============================
# Parallel simulation (FULL + two MICE strategies)
# ============================

suppressPackageStartupMessages({
  library(doParallel)
  library(doRNG)
  library(foreach)
  library(data.table)
  library(dplyr)
  library(MASS)
  library(reshape2)
  library(ggplot2)
  library(mice)
  library(parallel)
})

# ---- Engines (peal.fit.RI_mv + generate_Zhv_matrix)
source("PEAL_Engine_Multi-RI.R")

# -------------------------
# Core settings 
# -------------------------
N      <- 2                # number of replicates
H      <- 5                 # sites
px     <- 6                 # covariates
p_bin  <- 3
p_cont <- px - p_bin
py     <- 3                 # outcomes
m_hosp_range <- 200:300       # patients per site
nvisit_range <- 1:20        # visits per patient

# True parameters (fixed across replicates for clear bias eval)
beta_true <- matrix(runif(px * py, -3, 3), nrow = px, ncol = py)

sigma_u_true    <- 0.5
sigma_v_hosp_tr <- seq(0.4, 1, length.out = H)
sigma_e_true    <- 0.5
rho_true        <- 0.5

rho_mat <- matrix(rho_true, nrow = py, ncol = py); diag(rho_mat) <- 1

true_beta_vec  <- c(beta_true)
true_sigma_vec <- c(sigma_u_true, sigma_v_hosp_tr, sigma_e_true, rho_true)
par_names_beta  <- paste0("Beta", 1:(px * py))
par_names_sigma <- c("sigma_u", paste0("sigma_v_", seq_len(H)), "sigma_e", "rho")

# -------------------------
# DGP generator (one replicate)
# Returns: full (no missing), miss (MAR),
#          mice_site (site-wise MICE), mice_full (single MICE on full data)
# -------------------------
simulate_one <- function(seed) {
  set.seed(seed)
  
  # Sizes
  m_hosp <- sample(m_hosp_range, H, replace = TRUE)
  nn     <- rep(m_hosp, times = 1)
  id.hosp <- rep(seq_len(H), times = m_hosp)
  id.pat  <- sequence(nn)
  n_visits <- sample(nvisit_range, sum(nn), replace = TRUE)
  
  # Expand ID to visits
  id.visit         <- sequence(n_visits)
  id.hosp.expanded <- rep(id.hosp, times = n_visits)
  id.pat.expanded  <- rep(id.pat,  times = n_visits)
  
  # Random effects
  u_h  <- rnorm(H, 0, sigma_u_true)
  v_hi <- rnorm(sum(nn), 0, rep(sigma_v_hosp_tr, times = m_hosp))
  
  # Expand REs
  u_h_patient   <- rep(u_h, times = m_hosp)
  u_h_expanded  <- rep(u_h_patient, times = n_visits)
  v_hi_expanded <- rep(v_hi, times = n_visits)
  
  # Covariates
  Nobs  <- sum(n_visits)
  X_bin  <- matrix(rbinom(Nobs * p_bin, 1, 0.3), nrow = Nobs, ncol = p_bin)
  X_cont <- matrix(rnorm(Nobs * p_cont, 0, 1),   nrow = Nobs, ncol = p_cont)
  X_hij  <- cbind(X_bin, X_cont)
  
  # Residuals (exchangeable across outcomes)
  Sigma_eps <- diag(sigma_e_true, py) %*% rho_mat %*% diag(sigma_e_true, py)
  eps <- MASS::mvrnorm(Nobs, mu = rep(0, py), Sigma = Sigma_eps)
  
  # Outcomes
  Y_hij <- matrix(0, nrow = Nobs, ncol = py)
  for (k in 1:py) {
    Y_hij[, k] <- X_hij %*% beta_true[, k] + u_h_expanded + v_hi_expanded + eps[, k]
  }
  
  # Data
  dat <- data.table(
    site    = id.hosp.expanded,
    patient = id.pat.expanded,
    visit   = id.visit,
    X_hij, Y_hij
  ) |> as.data.frame()
  data.table::setnames(dat, c("site", "patient", "visit",
                              paste0("X", 1:px), paste0("Y", 1:py)))
  
  # Preprocess/order
  visit_count <- dat |>
    dplyr::group_by(site, patient) |>
    dplyr::summarise(total_visits = n(), .groups = "drop")
  
  dat_ord <- merge(dat, visit_count, by = c("site", "patient")) |>
    arrange(site, total_visits, patient) |>
    mutate(site = factor(site))
  
  # --------- Introduce MAR (same as your prior logic) ----------
  n <- nrow(dat_ord)
  dat_miss <- dat_ord
  
  p1_base <- 0.2   # baseline for Y1
  p2_base <- 0.3   # baseline for Y2
  p3_base <- 0.5   # baseline for Y3
  
  # Y3 depends on Y1, Y2 observed
  Y3_miss_prob <- ifelse(!is.na(dat_ord$Y1) & !is.na(dat_ord$Y2),
                         p3_base + 0.25,
                         p3_base - 0.25)
  Y3_miss_prob <- pmin(pmax(Y3_miss_prob, 0), 1)
  Y3_missing <- rbinom(n, 1, Y3_miss_prob) == 1
  dat_miss$Y3[Y3_missing] <- NA
  
  # Y1, Y2 depend on Y3 observed
  Y1_miss_prob <- ifelse(!is.na(dat_miss$Y3), p1_base + 0.2, p1_base - 0.1)
  Y2_miss_prob <- ifelse(!is.na(dat_miss$Y3), p2_base + 0.2, p2_base - 0.1)
  Y1_miss_prob <- pmin(pmax(Y1_miss_prob, 0), 1)
  Y2_miss_prob <- pmin(pmax(Y2_miss_prob, 0), 1)
  
  Y1_missing <- rbinom(n, 1, Y1_miss_prob) == 1
  Y2_missing <- rbinom(n, 1, Y2_miss_prob) == 1
  dat_miss$Y1[Y1_missing] <- NA
  dat_miss$Y2[Y2_missing] <- NA
  
  # Ensure not all Y's missing in a row
  all_y_missing <- apply(dat_miss[, paste0("Y", 1:py), drop = FALSE], 1, function(x) all(is.na(x)))
  if (any(all_y_missing)) {
    for (i in which(all_y_missing)) {
      y_keep <- sample(paste0("Y", 1:py), 1)
      dat_miss[i, y_keep] <- dat_ord[i, y_keep]
    }
  }
  
  # --------- MICE (site-wise 1 imputed dataset per site) ----------
  sites <- unique(dat_miss$site)
  imp_list <- vector("list", length(sites)); names(imp_list) <- as.character(sites)
  for (sid in sites) {
    sd <- dat_miss[dat_miss$site == sid, , drop = FALSE]
    imp_site <- mice::mice(sd, m = 1, printFlag = FALSE)
    imp_list[[as.character(sid)]] <- mice::complete(imp_site)
  }
  dat_mice_site <- do.call(rbind, imp_list)
  
  # Harmonize factor levels for site (site-wise)
  site_lvls <- levels(dat_ord$site)
  dat_mice_site$site <- factor(dat_mice_site$site, levels = site_lvls)
  
  # --------- MICE (imputation from full dataset) ----------
  imp_full <- mice::mice(dat_miss, m = 1, printFlag = FALSE)
  dat_mice_full <- mice::complete(imp_full)
  dat_mice_full$site <- factor(dat_mice_full$site, levels = site_lvls)
  
  list(
    full      = dat_ord,       # complete/no missing
    miss      = dat_miss,      # MAR with missingness
    mice_site = dat_mice_site, # site-wise MICE
    mice_full = dat_mice_full  # full-data MICE
  )
}

# -------------------------
# Build Z block lists for peal.fit.RI_mv (full/mice)
# -------------------------
build_Z_list <- function(dat_one, H) {
  Z <- vector("list", H)
  # make sure site is factor with H levels
  if (!is.factor(dat_one$site)) dat_one$site <- factor(dat_one$site)
  for (i in 1:H) {
    count_mat <- dat_one |>
      filter(site == levels(dat_one$site)[i]) |>
      group_by(site, patient) |>
      summarise(n_hi = n(), .groups = "drop")
    Z[[i]] <- generate_Zhv_matrix(count_mat)
  }
  Z
}

# -------------------------
# Pack fit objects into comparable vectors
# -------------------------
pack_mv <- function(fit, cor_model) {
  if (is.null(fit)) {
    return(list(b = rep(NA_real_, px * py),
                sigma = rep(NA_real_, H + 3)))  # u + v1..vH + e + rho
  }
  s2 <- fit$s2
  theta <- fit$theta
  sigU  <- sqrt(theta[1] * s2)
  sigVh <- sqrt(theta[1 + seq_len(H)] * s2)
  sigE  <- sqrt(s2)
  rh    <- if (!is.null(fit$rho) && cor_model == "exchangeable") fit$rho else NA_real_
  list(b = as.numeric(fit$b),
       sigma = c(sigU, sigVh, sigE, rh))
}

# -------------------------
# Single replicate runner: fit 6 models
# FULL_EX / FULL_IN
# MICE_SITE_EX / MICE_SITE_IN
# MICE_FULL_EX / MICE_FULL_IN
# -------------------------
fit_six <- function(dsets) {
  X_cols <- paste0("X", 1:px)
  Y_cols <- paste0("Y", 1:py)
  
  # -- FULL
  Y_full <- as.matrix(dsets$full[, Y_cols])
  X_full <- as.matrix(dsets$full[, X_cols])
  Z_full <- build_Z_list(dsets$full, H)
  id_full <- dsets$full$site
  
  fit_full_ex <- tryCatch(
    peal.fit.RI_mv(Y = Y_full, X = X_full, Z = Z_full, id.site = id_full,
                   weights = NULL, reml = TRUE,
                   corstr = "exchangeable", estimate_rho = TRUE, rho_init = 0.1,
                   verbose = FALSE),
    error = function(e) NULL
  )
  fit_full_in <- tryCatch(
    peal.fit.RI_mv(Y = Y_full, X = X_full, Z = Z_full, id.site = id_full,
                   weights = NULL, reml = TRUE,
                   corstr = "independence", estimate_rho = FALSE,
                   verbose = FALSE),
    error = function(e) NULL
  )
  
  # -- MICE SITE
  Y_ms <- as.matrix(dsets$mice_site[, Y_cols])
  X_ms <- as.matrix(dsets$mice_site[, X_cols])
  Z_ms <- build_Z_list(dsets$mice_site, H)
  id_ms <- dsets$mice_site$site
  
  fit_ms_ex <- tryCatch(
    peal.fit.RI_mv(Y = Y_ms, X = X_ms, Z = Z_ms, id.site = id_ms,
                   weights = NULL, reml = TRUE,
                   corstr = "exchangeable", estimate_rho = TRUE, rho_init = 0.1,
                   verbose = FALSE),
    error = function(e) NULL
  )
  fit_ms_in <- tryCatch(
    peal.fit.RI_mv(Y = Y_ms, X = X_ms, Z = Z_ms, id.site = id_ms,
                   weights = NULL, reml = TRUE,
                   corstr = "independence", estimate_rho = FALSE,
                   verbose = FALSE),
    error = function(e) NULL
  )
  
  # -- MICE FULL
  Y_mf <- as.matrix(dsets$mice_full[, Y_cols])
  X_mf <- as.matrix(dsets$mice_full[, X_cols])
  Z_mf <- build_Z_list(dsets$mice_full, H)
  id_mf <- dsets$mice_full$site
  
  fit_mf_ex <- tryCatch(
    peal.fit.RI_mv(Y = Y_mf, X = X_mf, Z = Z_mf, id.site = id_mf,
                   weights = NULL, reml = TRUE,
                   corstr = "exchangeable", estimate_rho = TRUE, rho_init = 0.1,
                   verbose = FALSE),
    error = function(e) NULL
  )
  fit_mf_in <- tryCatch(
    peal.fit.RI_mv(Y = Y_mf, X = X_mf, Z = Z_mf, id.site = id_mf,
                   weights = NULL, reml = TRUE,
                   corstr = "independence", estimate_rho = FALSE,
                   verbose = FALSE),
    error = function(e) NULL
  )
  
  list(
    full_ex       = pack_mv(fit_full_ex, "exchangeable"),
    full_in       = pack_mv(fit_full_in, "independence"),
    mice_site_ex  = pack_mv(fit_ms_ex,   "exchangeable"),
    mice_site_in  = pack_mv(fit_ms_in,   "independence"),
    mice_full_ex  = pack_mv(fit_mf_ex,   "exchangeable"),
    mice_full_in  = pack_mv(fit_mf_in,   "independence")
  )
}

# -------------------------
# Parallel run
# -------------------------
num_cores <- max(1, parallelly::availableCores() - 1)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
registerDoRNG(12345)

model_names <- c("FULL_EX","FULL_IN",
                 "MICE_SITE_EX","MICE_SITE_IN",
                 "MICE_FULL_EX","MICE_FULL_IN")

results <- foreach(k = 1:N, .packages = c("data.table","dplyr","MASS","mice")) %dopar% {
  source("PEAL_Engine_Multi-RI.R")
  dsets <- simulate_one(sample(1:1e6, 1))
  fit_six(dsets)
}

stopCluster(cl)

# Optionally save:
# saveRDS(results, file = "result_full_micesix.rds")
# results <- readRDS("result_full_micesix.rds")

# -------------------------
# Collect results into matrices
# -------------------------
init_mat <- function(nr) matrix(NA_real_, nrow = nr, ncol = N,
                                dimnames = list(NULL, paste0("sim", 1:N)))

store <- lapply(model_names, function(nm) {
  list(beta = init_mat(px * py), sigma = init_mat(H + 3))
})
names(store) <- model_names

for (k in seq_along(results)) {
  store$FULL_EX$beta[, k]       <- results[[k]]$full_ex$b
  store$FULL_EX$sigma[, k]      <- results[[k]]$full_ex$sigma
  
  store$FULL_IN$beta[, k]       <- results[[k]]$full_in$b
  store$FULL_IN$sigma[, k]      <- results[[k]]$full_in$sigma
  
  store$MICE_SITE_EX$beta[, k]  <- results[[k]]$mice_site_ex$b
  store$MICE_SITE_EX$sigma[, k] <- results[[k]]$mice_site_ex$sigma
  
  store$MICE_SITE_IN$beta[, k]  <- results[[k]]$mice_site_in$b
  store$MICE_SITE_IN$sigma[, k] <- results[[k]]$mice_site_in$sigma
  
  store$MICE_FULL_EX$beta[, k]  <- results[[k]]$mice_full_ex$b
  store$MICE_FULL_EX$sigma[, k] <- results[[k]]$mice_full_ex$sigma
  
  store$MICE_FULL_IN$beta[, k]  <- results[[k]]$mice_full_in$b
  store$MICE_FULL_IN$sigma[, k] <- results[[k]]$mice_full_in$sigma
}

# -------------------------
# Long frames & summaries
# -------------------------
mk_long <- function(mat, par_names, model_label, value_nm = "Estimate") {
  df <- reshape2::melt(as.data.frame(mat), value.name = value_nm)
  names(df) <- c("Simulation", value_nm)
  df$Parameter <- rep(par_names, ncol(mat))
  df$Model <- model_label
  df
}

# Beta
beta_long <- do.call(rbind, lapply(names(store), function(nm) {
  mk_long(store[[nm]]$beta, par_names_beta, nm, "Estimate")
}))
beta_long$True <- rep(true_beta_vec,
                      times = ncol(store[[1]]$beta) * length(store))

beta_summary <- beta_long %>%
  group_by(Model, Parameter) %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias       = Mean_Est - True_Value,
    Variance   = var(Estimate, na.rm = TRUE),
    MSE        = Bias^2 + Variance,
    .groups = "drop"
  ) %>% arrange(Parameter, Model)

# Sigma (order: sigma_u, sigma_v_1..H, sigma_e, rho)
sigma_long <- do.call(rbind, lapply(names(store), function(nm) {
  mk_long(store[[nm]]$sigma, par_names_sigma, nm, "Estimate")
}))
sigma_long$True <- rep(true_sigma_vec,
                       times = ncol(store[[1]]$sigma) * length(store))

sigma_summary <- sigma_long %>%
  group_by(Model, Parameter) %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias       = Mean_Est - True_Value,
    Variance   = var(Estimate, na.rm = TRUE),
    MSE        = Bias^2 + Variance,
    .groups = "drop"
  ) %>% arrange(Parameter, Model)

# -------------------------
# Quick plots (6 models)
# -------------------------
beta_long$Parameter <- factor(beta_long$Parameter,
                              levels = paste0("Beta",1:(py*px)))

p_beta_bias <- beta_long %>%
  mutate(Bias = Estimate - True) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1, width = 0.15, height = 0, size = 0.3) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~ Model, ncol = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, alpha = 0.4, color = "red") +
  ggtitle("Fixed Effects Bias across 6 models")

p_sigma_bias <- sigma_long %>%
  mutate(Bias = Estimate - True) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1, width = 0.15, height = 0, size = 0.3) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~ Model, ncol = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, alpha = 0.4, color = "red") + 
  ggtitle("Variance Components Bias across 6 models")

print(p_beta_bias)
print(p_sigma_bias)

# Inspect numeric summaries
beta_summary
sigma_summary

# Bar visuals
ggplot(beta_summary %>%
         mutate(Parameter = factor(Parameter,
                                   levels = paste0("Beta", 1:(px*py)))),
       aes(x = Parameter, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(sigma_summary,
       aes(x = Parameter, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(beta_summary %>%
         mutate(Parameter = factor(Parameter,
                                   levels = paste0("Beta", 1:(px*py)))), 
       aes(x = Parameter, y = Variance, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(sigma_summary,
       aes(x = Parameter, y = Variance, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(beta_summary %>%
         mutate(Parameter = factor(Parameter,
                                   levels = paste0("Beta", 1:(px*py)))),
       aes(x = Parameter, y = MSE, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(sigma_summary %>% filter(Parameter != "rho"),
       aes(x = Parameter, y = MSE, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
