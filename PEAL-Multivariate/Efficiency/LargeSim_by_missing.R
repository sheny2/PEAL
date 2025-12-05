suppressPackageStartupMessages({
  library(doParallel)
  library(doRNG)
  library(foreach)
  library(data.table)
  library(dplyr)
  library(MASS)
  library(reshape2)
  library(ggplot2)
  library(Matrix)    # for nearPD
  library(gridExtra) # for grid.arrange
})

## ------------------------------------------------
## 1. Global settings and true parameters
## ------------------------------------------------

source("MV-PEAL_Pattern.R")

H      <- 5    # number of sites
px     <- 4
p_bin  <- 2
p_cont <- px - p_bin
py     <- 3   # number of outcomes

X_cols <- paste0("X", 1:px)
Y_cols <- paste0("Y", 1:py)

set.seed(1)
beta_true <- matrix(runif(px * py, -3, 3), nrow = px, ncol = py)

sigma_u_true      <- 0.5
sigma_v_hosp_true <- runif(H, min = 0.5, max = 1)   # fixed across sims
sigma_e_true      <- 0.3

sigma_e_vec <- rep(sigma_e_true, py)
De          <- diag(sigma_e_vec, py)

# True parameter vectors (corr is rho-specific; we fix one rho setting below)
true_beta_vec     <- c(beta_true)
true_sigmaRE_vec  <- c(sigma_u_true, sigma_v_hosp_true)  # length H+1
true_sigma_e      <- sigma_e_true

par_names_beta    <- paste0("Beta", 1:(px * py))
par_names_sigmaRE <- c("sigma_u", paste0("sigma_v_", 1:H))
par_name_sigma_e  <- "sigma_e"
corr_param_names  <- paste0("rho_", rep(1:py, each = py), "_", rep(1:py, times = py))

## ------------------------------------------------
## 2. Fix one correlation setting globally
## ------------------------------------------------
rho12 <- -0.25
rho13 <- -0.30
rho23 <-  0.35

corr_mat_true <- matrix(c(
  1,      rho12,  rho13,
  rho12,  1,      rho23,
  rho13,  rho23,  1
), nrow = py, byrow = TRUE)

# Ensure SPD
if (min(eigen(corr_mat_true, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
  corr_mat_true <- as.matrix(nearPD(corr_mat_true, corr = TRUE)$mat)
}

Sigma_e_global <- De %*% corr_mat_true %*% De
true_corr_vec  <- as.vector(corr_mat_true)

## ------------------------------------------------
## 3. Helper: build Z list (site-wise random intercepts)
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
## 4. Single replicate DGP: MAR-masked data
##    Takes Sigma_e and baseline missingness p1_base, p2_base, p3_base
## ------------------------------------------------
simulate_one <- function(seed, Sigma_e,
                         p1_base, p2_base, p3_base) {
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
  
  # Residuals with given Sigma_e
  epsilon_hij <- MASS::mvrnorm(n = Nobs, mu = rep(0, py), Sigma = Sigma_e)
  
  # Outcomes
  Y_hij <- matrix(0, nrow = Nobs, ncol = py)
  for (k in 1:py) {
    Y_hij[, k] <- X_hij %*% beta_true[, k] +
      u_h_expanded + v_hi_expanded + epsilon_hij[, k]
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
  
  ## --- MAR missingness (parametrized by p*_base) ---
  rearranged_data_with_missing <- rearranged_data
  n <- nrow(rearranged_data)
  
  # Y3 depends on Y1,Y2 observed
  Y3_miss_prob <- ifelse(!is.na(rearranged_data$Y1) & !is.na(rearranged_data$Y2),
                         p3_base + 0.15,
                         p3_base - 0.15)
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
## 5. Fit helpers and packers
## ------------------------------------------------

pack_fit <- function(fit, corstr, py, H) {
  if (is.null(fit)) {
    return(list(
      b        = rep(NA_real_, px * py),
      sigmaRE  = rep(NA_real_, H + 1),
      sigma_e  = NA_real_,
      corr_vec = rep(NA_real_, py * py)
    ))
  }
  b_hat       <- as.numeric(fit$b)
  sigmaRE_hat <- sqrt(fit$theta * fit$s2)  # H+1 random-effect SDs
  sigmae_hat  <- sqrt(fit$s2)
  
  if (!is.null(fit$Corr)) {
    Corr_hat      <- fit$Corr
    corr_vec_hat  <- as.vector(Corr_hat)
  } else {
    corr_vec_hat  <- rep(NA_real_, py * py)
  }
  
  list(
    b        = b_hat,
    sigmaRE  = sigmaRE_hat,
    sigma_e  = sigmae_hat,
    corr_vec = corr_vec_hat
  )
}

fit_pattern_models <- function(dat_miss) {
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
## 6. One replicate wrapper
## ------------------------------------------------
run_one_rep <- function(seed, Sigma_e,
                        p1_base, p2_base, p3_base) {
  sim <- simulate_one(seed, Sigma_e,
                      p1_base = p1_base,
                      p2_base = p2_base,
                      p3_base = p3_base)
  patt_res <- fit_pattern_models(sim$miss)
  c(patt_res)
}

## ------------------------------------------------
## 7. Helper to init matrices and long-format helper
## ------------------------------------------------
init_mat <- function(nr, N) {
  matrix(NA_real_, nrow = nr, ncol = N,
         dimnames = list(NULL, paste0("sim", 1:N)))
}

mk_long <- function(mat, par_names, model_label, value_nm = "Estimate") {
  df <- reshape2::melt(as.data.frame(mat), value.name = value_nm)
  names(df) <- c("Simulation", value_nm)
  df$Parameter <- rep(par_names, ncol(mat))
  df$Model <- model_label
  df
}

## ------------------------------------------------
## 8. BIG DRIVER: run whole pipeline for one (p1_base,p2_base,p3_base)
## ------------------------------------------------
run_sim_for_missing <- function(p1_base, p2_base, p3_base,
                                N = 100,
                                label = NULL) {
  if (is.null(label)) {
    label <- paste0("miss_",
                    paste(sprintf("%0.2f", c(p1_base, p2_base, p3_base)),
                          collapse = "_"))
  }
  cat("=== Running simulation for", label, 
      "with (p1, p2, p3) =",
      sprintf("(%.2f, %.2f, %.2f)\n", p1_base, p2_base, p3_base))
  
  Sigma_e     <- Sigma_e_global
  model_names <- c("OBS_IN","OBS_EX","OBS_UN")
  
  # Parallel simulation
  num_cores <- max(1, parallelly::availableCores())
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  registerDoRNG(12345)
  
  results <- foreach(
    k = 1:N,
    .packages = c("data.table","dplyr","MASS","Matrix","reshape2"),
    .export   = c(
      "run_one_rep",
      "simulate_one",
      "fit_pattern_models",
      "pack_fit",
      "H", "px", "p_bin", "p_cont", "py",
      "X_cols", "Y_cols",
      "sigma_u_true", "sigma_v_hosp_true",
      "beta_true", "De", "Sigma_e_global",
      "true_beta_vec", "true_sigmaRE_vec", "true_sigma_e",
      "par_names_beta", "par_names_sigmaRE", "par_name_sigma_e",
      "corr_param_names", "true_corr_vec"
    )
  ) %dopar% {
    source("MV-PEAL_Pattern.R")
    run_one_rep(
      seed    = sample(1:1e6, 1),
      Sigma_e = Sigma_e,
      p1_base = p1_base,
      p2_base = p2_base,
      p3_base = p3_base
    )
  }
  
  stopCluster(cl)
  
  # Save raw results
  saveRDS(results, file = paste0("Simulation_Results_", label, ".rds"))
  
  ## ------------------------------------------------
  ## Collect results
  ## ------------------------------------------------
  store <- lapply(model_names, function(nm) {
    list(
      beta    = init_mat(px * py, N),
      sigmaRE = init_mat(H + 1, N),
      sigma_e = init_mat(1, N),
      corr    = init_mat(py * py, N)
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
  ## Bias summaries
  ## ------------------------------------------------
  beta_long <- do.call(rbind, lapply(names(store), function(nm) {
    mk_long(store[[nm]]$beta, par_names_beta, nm, "Estimate")
  }))
  beta_long$True <- rep(true_beta_vec,
                        times = ncol(store[[1]]$beta) * length(store))
  
  sigmaRE_long <- do.call(rbind, lapply(names(store), function(nm) {
    mk_long(store[[nm]]$sigmaRE, par_names_sigmaRE, nm, "Estimate")
  }))
  sigmaRE_long$True <- rep(true_sigmaRE_vec,
                           times = ncol(store[[1]]$sigmaRE) * length(store))
  
  sigmae_long <- do.call(rbind, lapply(names(store), function(nm) {
    mk_long(store[[nm]]$sigma_e, par_name_sigma_e, nm, "Estimate")
  }))
  sigmae_long$True <- rep(true_sigma_e,
                          times = ncol(store[[1]]$sigma_e) * length(store))
  
  corr_long <- do.call(rbind, lapply(names(store), function(nm) {
    mk_long(store[[nm]]$corr, corr_param_names, nm, "Estimate")
  }))
  corr_long$True <- rep(true_corr_vec,
                        times = ncol(store[[1]]$corr) * length(store))
  
  ## ------------------------------------------------
  ## Plots
  ## ------------------------------------------------
  beta_long$Parameter <- factor(beta_long$Parameter,
                                levels = paste0("Beta", 1:(px*py)))
  
  p_beta_bias <- beta_long %>%
    mutate(Bias = Estimate - True) %>%
    ggplot(aes(x = Parameter, y = Bias)) +
    geom_jitter(alpha = 0.1, width = 0.15, height = 0, size = 0.3) +
    geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) +
    facet_wrap(~ Model, ncol = 3, scales = "free_y") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.7) +
    ggtitle(paste("Fixed Effects Bias -", label))
  
  p_sigmaRE_bias <- sigmaRE_long %>%
    mutate(Bias = Estimate - True) %>%
    ggplot(aes(x = Parameter, y = Bias)) +
    geom_jitter(alpha = 0.1, width = 0.15, height = 0, size = 0.3) +
    geom_boxplot(fill = "lightgreen", alpha = 0.6, outlier.shape = NA) +
    facet_wrap(~ Model, ncol = 3, scales = "free_y") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.7) +
    ggtitle(paste("Random-effect SD Bias -", label))
  
  p_sigmae_bias <- sigmae_long %>%
    mutate(Bias = Estimate - True) %>%
    ggplot(aes(x = Model, y = Bias)) +
    geom_boxplot(fill = "orange", alpha = 0.6, outlier.shape = NA) +
    theme_bw() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.7) +
    ggtitle(paste("Residual SD Bias -", label))
  
  corr_long_offdiag <- corr_long %>%
    mutate(i = as.integer(substr(Parameter, 5, 5)),
           j = as.integer(substr(Parameter, 7, 7))) %>%
    filter(i != j)
  corr_long_offdiag$Parameter <- factor(corr_long_offdiag$Parameter,
                                        levels = unique(corr_long_offdiag$Parameter))
  
  p_corr_bias <- corr_long_offdiag %>%
    filter(Parameter %in% unique(corr_long_offdiag$Parameter)[c(1,2,4)]) %>%
    mutate(Bias = Estimate - True) %>%
    ggplot(aes(x = Parameter, y = Bias, fill = Model)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    theme_bw() +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.7) +
    ggtitle(paste("Correlation Parameter Bias (off-diagonals) -", label))
  
  grid_plot <- gridExtra::grid.arrange(
    p_beta_bias, p_sigmae_bias, p_sigmaRE_bias, p_corr_bias,
    ncol = 2
  )
  
  ggsave(filename = paste0("BiasGrid_", label, ".png"),
         plot = grid_plot, width = 10, height = 8)
  
  ## ------------------------------------------------
  ## Variance & Relative Efficiency (beta)
  ## ------------------------------------------------
  RE_table <- beta_long %>%
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
  
  print(RE_table)
  saveRDS(RE_table, file = paste0("RE_table_beta_", label, ".rds"))
  
  invisible(list(
    results      = results,
    RE_beta      = RE_table,
    beta_long    = beta_long,
    sigmaRE_long = sigmaRE_long,
    sigmae_long  = sigmae_long,
    corr_long    = corr_long
  ))
}

## ------------------------------------------------
## 9. Define 8 missingness settings and run
## ------------------------------------------------

# p1_base, p2_base are modest and fixed;
# p3_base increases from low to high (more missingness in Y3).
miss_grid <- data.frame(
  setting = paste0("miss_", 1:8),
  p1_base = rep(0.10, 8),
  p2_base = rep(0.20, 8),
  p3_base = seq(0.1, 0.8, length.out = 8)
)

N_reps <- 100  # change if you want more/less Monte Carlo replicates

all_out_miss <- lapply(seq_len(nrow(miss_grid)), function(i) {
  row <- miss_grid[i, ]
  run_sim_for_missing(
    p1_base = row$p1_base,
    p2_base = row$p2_base,
    p3_base = row$p3_base,
    N       = N_reps,
    label   = row$setting
  )
})
names(all_out_miss) <- miss_grid$setting


data.frame(
  setting  = paste0("miss_", 1:8),
  #p1_base  = rep(0.20, 8),
  #p1_base  = rep(0.20, 8),
  p1_base  = seq(0.1, 0.45, length.out = 8),
  p2_base  = seq(0.1, 0.65, length.out = 8),
  p3_base  = seq(0.1, 0.85, length.out = 8)
) %>% mutate(missing_base_avg = round((p1_base+ p2_base+ p3_base)/3,2))
