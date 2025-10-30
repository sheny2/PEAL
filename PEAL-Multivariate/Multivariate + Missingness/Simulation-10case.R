# ============================
# Parallel simulation: 8 fits/replicate
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
})

source("PEAL_Engine_Multi-RI.R")
source("PEAL_Engine_Missing_Multi-RI3.R")

# -------------------------
# Helper: complete-case filter over Y's
# -------------------------
cc <- function(df, Y_cols = grep("^Y", names(df), value = TRUE)) {
  keep <- stats::complete.cases(df[, Y_cols, drop = FALSE])
  df[keep, , drop = FALSE]
}

# -------------------------
# Core settings 
# -------------------------
N      <- 200                # number of replicates
H      <- 5                  # sites
px     <- 6                  # covariates
p_bin  <- 3
p_cont <- px - p_bin
py     <- 3                  # outcomes
m_hosp_range <- 50:70        # patients per site
nvisit_range <- 1:20         # visits per patient

# True parameters (fixed across replicates for clear bias eval)
set.seed(1)
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
  u_h_patient  <- rep(u_h, times = m_hosp)
  u_h_expanded <- rep(u_h_patient, times = n_visits)
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
  for (k in 1:py) Y_hij[, k] <- X_hij %*% beta_true[, k] + u_h_expanded + v_hi_expanded + eps[, k]
  
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
  
  # Introduce MCAR 
  missing_prob <- 0.4
  dat_miss <- dat_ord
  # for (ycol in paste0("Y", 1:py)) {
  #   idx <- sample.int(nrow(dat_miss), size = round(nrow(dat_miss) * missing_prob))
  #   dat_miss[idx, ycol] <- NA
  # }

  # Introduce MCAR 2.
  # dat_miss$Y1[sample(1:nrow(dat_ord), round(nrow(dat_ord)*0.2))] <- NA
  # dat_miss$Y2[sample(1:nrow(dat_ord), round(nrow(dat_ord)*0.3))] <- NA
  # dat_miss$Y3[sample(1:nrow(dat_ord), round(nrow(dat_ord)*0.6))] <- NA

  # Introduce MAR
  # Step 1. Start with baseline missingness probabilities
  n <- nrow(dat_ord)
  p1_base <- 0.2   # baseline for Y1
  p2_base <- 0.3   # baseline for Y2
  p3_base <- 0.5   # baseline for Y3

  # Step 2. For Y3: missingness depends on whether Y1 and Y2 are observed
  # If both Y1 and Y2 are observed → Y3 more likely missing
  Y3_miss_prob <- ifelse(!is.na(dat_ord$Y1) & !is.na(dat_ord$Y2),
                         p3_base + 0.25,  # boost missingness
                         p3_base - 0.25)  # lower missingness
  Y3_miss_prob <- pmin(pmax(Y3_miss_prob, 0), 1)
  Y3_missing <- rbinom(n, 1, Y3_miss_prob) == 1
  dat_miss$Y3[Y3_missing] <- NA

  # Step 3. For Y1 and Y2: missingness depends on whether Y3 is observed
  # If Y3 is observed → Y1/Y2 more likely missing
  Y1_miss_prob <- ifelse(!is.na(dat_miss$Y3),
                         p1_base + 0.2,
                         p1_base - 0.1)
  Y2_miss_prob <- ifelse(!is.na(dat_miss$Y3),
                         p2_base + 0.2,
                         p2_base - 0.1)
  Y1_miss_prob <- pmin(pmax(Y1_miss_prob, 0), 1)
  Y2_miss_prob <- pmin(pmax(Y2_miss_prob, 0), 1)

  Y1_missing <- rbinom(n, 1, Y1_miss_prob) == 1
  Y2_missing <- rbinom(n, 1, Y2_miss_prob) == 1

  dat_miss$Y1[Y1_missing] <- NA
  dat_miss$Y2[Y2_missing] <- NA
  
  
  all_y_missing <- apply(dat_miss[, paste0("Y", 1:py), drop = FALSE], 1, function(x) all(is.na(x)))
  if (any(all_y_missing)) {
    for (i in which(all_y_missing)) {
      y_keep <- sample(paste0("Y", 1:py), 1)
      dat_miss[i, y_keep] <- dat_ord[i, y_keep]
    }
  }
  
  # CC data
  dat_cc <- cc(dat_miss, paste0("Y", 1:py))
  
  # MICE (site-wise 1 imputed dataset per site)
  sites <- unique(dat_miss$site)
  imp_list <- vector("list", length(sites)); names(imp_list) <- as.character(sites)
  for (sid in sites) {
    sd <- dat_miss[dat_miss$site == sid, , drop = FALSE]
    imp <- mice::mice(sd, m = 1, printFlag = FALSE)
    imp_list[[as.character(sid)]] <- mice::complete(imp)
  }
  dat_mice <- do.call(rbind, imp_list)
  
  list(
    full  = dat_ord,
    miss  = dat_miss,
    cc    = dat_cc,
    mice  = dat_mice
  )
}

# -------------------------
# Build Z block lists for peal.fit.RI_mv (full/cc/mice)
# -------------------------
build_Z_list <- function(dat_one, H) {
  Z <- vector("list", H)
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

pack_patterns <- function(fit, cor_model) {
  if (is.null(fit)) {
    return(list(b = rep(NA_real_, px * py),
                sigma = rep(NA_real_, H + 3)))
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
# Single replicate runner: fit 8 models
# -------------------------
fit_eight <- function(dsets) {
  source("PEAL_Engine_Multi-RI.R")
  source("PEAL_Engine_Missing_Multi-RI3.R")
  
  # Common pieces
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
  
  # -- CC
  Y_cc <- as.matrix(dsets$cc[, Y_cols])
  X_cc <- as.matrix(dsets$cc[, X_cols])
  Z_cc <- build_Z_list(dsets$cc, H)
  id_cc <- dsets$cc$site
  
  fit_cc_ex <- tryCatch(
    peal.fit.RI_mv(Y = Y_cc, X = X_cc, Z = Z_cc, id.site = id_cc,
                   weights = NULL, reml = TRUE,
                   corstr = "exchangeable", estimate_rho = TRUE, rho_init = 0.1,
                   verbose = FALSE),
    error = function(e) NULL
  )
  fit_cc_in <- tryCatch(
    peal.fit.RI_mv(Y = Y_cc, X = X_cc, Z = Z_cc, id.site = id_cc,
                   weights = NULL, reml = TRUE,
                   corstr = "independence", estimate_rho = FALSE,
                   verbose = FALSE),
    error = function(e) NULL
  )
  
  # -- MICE
  Y_mi <- as.matrix(dsets$mice[, Y_cols])
  X_mi <- as.matrix(dsets$mice[, X_cols])
  Z_mi <- build_Z_list(dsets$mice, H)
  id_mi <- dsets$mice$site
  
  fit_mi_ex <- tryCatch(
    peal.fit.RI_mv(Y = Y_mi, X = X_mi, Z = Z_mi, id.site = id_mi,
                   weights = NULL, reml = TRUE,
                   corstr = "exchangeable", estimate_rho = TRUE, rho_init = 0.1,
                   verbose = FALSE),
    error = function(e) NULL
  )
  fit_mi_in <- tryCatch(
    peal.fit.RI_mv(Y = Y_mi, X = X_mi, Z = Z_mi, id.site = id_mi,
                   weights = NULL, reml = TRUE,
                   corstr = "independence", estimate_rho = FALSE,
                   verbose = FALSE),
    error = function(e) NULL
  )
  
  # -- Observed-pattern summaries
  fit_pat_ex <- tryCatch(
    peal.fit.RI_mv_patterns(
      data = dsets$miss, X_cols = X_cols, Y_cols = Y_cols,
      site_col = "site", patient_col = "patient",
      reml = TRUE, corstr = "exchangeable", rho_init = 0.1, verbose = FALSE
    ),
    error = function(e) NULL
  )
  fit_pat_in <- tryCatch(
    peal.fit.RI_mv_patterns(
      data = dsets$miss, X_cols = X_cols, Y_cols = Y_cols,
      site_col = "site", patient_col = "patient",
      reml = TRUE, corstr = "independence", rho_init = 0.1, verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  # -- Observed-pattern summaries with EM
  fit_pat_ex_em <- tryCatch(
    peal.em.fit.RI_mv(
      data     = dsets$miss,
      X_cols   = X_cols, Y_cols   = Y_cols,
      site_col = "site", patient_col = "patient",
      corstr_init = "exchangeable",
      em_iter = 1, reml = TRUE, rho_init = 0.1, verbose = FALSE
    ),
    error = function(e) NULL
  )
  fit_pat_in_em <- tryCatch(
    peal.em.fit.RI_mv(
      data     = dsets$miss,
      X_cols   = X_cols, Y_cols   = Y_cols,
      site_col = "site", patient_col = "patient",
      corstr_init = "independence",
      em_iter = 1, reml = TRUE, rho_init = 0.1, verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  list(
    full_ex  = pack_mv(fit_full_ex, "exchangeable"),
    full_in  = pack_mv(fit_full_in, "independence"),
    cc_ex    = pack_mv(fit_cc_ex,   "exchangeable"),
    cc_in    = pack_mv(fit_cc_in,   "independence"),
    mice_ex  = pack_mv(fit_mi_ex,   "exchangeable"),
    mice_in  = pack_mv(fit_mi_in,   "independence"),
    patt_ex  = pack_patterns(fit_pat_ex, "exchangeable"),
    patt_in  = pack_patterns(fit_pat_in, "independence"),
    patt_em_ex  = pack_patterns(fit_pat_ex_em, "exchangeable"),
    patt_em_in  = pack_patterns(fit_pat_in_em, "independence")
  )
}

# -------------------------
# Parallel run
# -------------------------
num_cores <- max(1, parallelly::availableCores() - 1)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
registerDoRNG(12345)

model_names <- c("FULL_EX","FULL_IN","CC_EX","CC_IN","MICE_EX","MICE_IN","OBS_EX","OBS_IN","OBS_EM_EX","OBS_EM_IN")

results <- foreach(k = 1:N, .packages = c("data.table","dplyr","MASS","mice")) %dopar% {
  dsets <- simulate_one(1000 + k)
  fit_eight(dsets)
}

stopCluster(cl)


# saveRDS(results, file = "result10_difflvl.rds")
results <- readRDS("result10_difflvl.rds")
# results <- readRDS("result10.rds")


# -------------------------
# Collect results into matrices
# -------------------------
# For each model, build matrices: beta (px*py × N) and sigma (H+3 × N)
init_mat <- function(nr) matrix(NA_real_, nrow = nr, ncol = N,
                                dimnames = list(NULL, paste0("sim", 1:N)))

store <- lapply(model_names, function(nm) {
  list(beta = init_mat(px * py), sigma = init_mat(H + 3))
})
names(store) <- model_names

for (k in seq_along(results)) {
  resk <- results[[k]]
  for (nm in names(resk)) {
    store[[toupper(gsub("_", "", nm))]] <- store[[toupper(gsub("_", "", nm))]] # no-op to ensure exists
  }
  # map names exactly
  store$FULL_EX$beta[, k] <- results[[k]]$full_ex$b
  store$FULL_EX$sigma[,k] <- results[[k]]$full_ex$sigma
  
  store$FULL_IN$beta[, k] <- results[[k]]$full_in$b
  store$FULL_IN$sigma[,k] <- results[[k]]$full_in$sigma
  
  store$CC_EX$beta[, k]   <- results[[k]]$cc_ex$b
  store$CC_EX$sigma[,k]   <- results[[k]]$cc_ex$sigma
  
  store$CC_IN$beta[, k]   <- results[[k]]$cc_in$b
  store$CC_IN$sigma[,k]   <- results[[k]]$cc_in$sigma
  
  store$MICE_EX$beta[, k] <- results[[k]]$mice_ex$b
  store$MICE_EX$sigma[,k] <- results[[k]]$mice_ex$sigma
  
  store$MICE_IN$beta[, k] <- results[[k]]$mice_in$b
  store$MICE_IN$sigma[,k] <- results[[k]]$mice_in$sigma
  
  store$OBS_EX$beta[, k]  <- results[[k]]$patt_ex$b
  store$OBS_EX$sigma[,k]  <- results[[k]]$patt_ex$sigma
  
  store$OBS_IN$beta[, k]  <- results[[k]]$patt_in$b
  store$OBS_IN$sigma[,k]  <- results[[k]]$patt_in$sigma
  
  store$OBS_EM_EX$beta[, k]  <- results[[k]]$patt_em_ex$b
  store$OBS_EM_EX$sigma[,k]  <- c(results[[k]]$patt_ex$sigma[1:(H+1)], results[[k]]$patt_em_ex$sigma[(H+2):(H+3)])
  
  store$OBS_EM_IN$beta[, k]  <- results[[k]]$patt_em_in$b
  store$OBS_EM_IN$sigma[,k]  <- c(results[[k]]$patt_in$sigma[1:(H+1)], results[[k]]$patt_em_in$sigma[(H+2):(H+3)])
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
beta_long$True <- rep(true_beta_vec, times = ncol(store[[1]]$beta) * length(store))

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
sigma_long$True <- rep(true_sigma_vec, times = ncol(store[[1]]$sigma) * length(store))

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


# Quick plots 
beta_long$Parameter <- factor(beta_long$Parameter, levels = paste0("Beta",1:(py*px)))
p_beta_bias <- beta_long %>%
  mutate(Bias = Estimate - True) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1, width = 0.15, height = 0, size = 0.3) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) +
  # geom_violin(alpha = 1, trim = F, color = F,  fill = "lightblue") +
  facet_wrap(~ Model, ncol = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1, alpha = 0.4) +
  ggtitle("Fixed Effects Bias across 10 models")

p_sigma_bias <- sigma_long %>%
  mutate(Bias = Estimate - True) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1, width = 0.15, height = 0, size = 0.3) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) +
  # geom_violin(alpha = 1, trim = F, color = F,  fill = "lightblue") +
  facet_wrap(~ Model, ncol = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1, alpha = 0.4) + 
  ggtitle("Variance Components Bias across 10 models")

print(p_beta_bias)
print(p_sigma_bias)

# Inspect numeric summaries
beta_summary
sigma_summary



# Additional Visual
ggplot(beta_summary %>% mutate(Parameter = factor(Parameter, levels = paste0("Beta", 1:(px*py)))), aes(x = Parameter, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(sigma_summary, aes(x = Parameter, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



ggplot(beta_summary %>% mutate(Parameter = factor(Parameter, levels = paste0("Beta", 1:(px*py)))), aes(x = Parameter, y = Variance, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(sigma_summary, aes(x = Parameter, y = Variance, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



ggplot(beta_summary %>% mutate(Parameter = factor(Parameter, levels = paste0("Beta", 1:(px*py)))), aes(x = Parameter, y = MSE, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(sigma_summary %>% filter(Parameter != "rho"), aes(x = Parameter, y = MSE, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 