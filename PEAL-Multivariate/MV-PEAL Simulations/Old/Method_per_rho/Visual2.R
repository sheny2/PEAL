suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(tidyr)
  library(stringr)
})

## ------------------------------------------------
## 0. Global info (match your simulation)
## ------------------------------------------------
H  <- 5
px <- 4
py <- 3

par_names_beta    <- paste0("Beta", 1:(px * py))
par_names_sigmaRE <- c("sigma_u", paste0("sigma_v_", 1:H))
par_name_sigma_e  <- "sigma_e"
corr_param_names  <- paste0("rho_", rep(1:py, each = py), "_", rep(1:py, times = py))

## Helpers (matching sim script)
init_mat <- function(nr, N) {
  matrix(
    NA_real_,
    nrow = nr,
    ncol = N,
    dimnames = list(NULL, paste0("sim", 1:N))
  )
}

mk_long <- function(mat, par_names, model_label, value_nm = "Estimate") {
  df <- reshape2::melt(as.data.frame(mat), value.name = value_nm)
  names(df) <- c("Simulation", value_nm)
  df$Parameter <- rep(par_names, ncol(mat))
  df$Model     <- model_label
  df
}

## A nicer, consistent plotting theme
theme_sim <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(face = "bold"),
      legend.position  = "top",
      legend.title     = element_blank(),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0.5)
    )
}

## ------------------------------------------------
## 1. Recreate rho_grid / rho_meta to attach rho info
## ------------------------------------------------
rho_grid <- list(
  c(0.05, 0.10, 0.15),
  c(0.05, 0.15, 0.25),
  c(0.05, 0.20, 0.35),
  c(0.05, 0.25, 0.45),
  c(0.05, 0.30, 0.55),
  c(0.05, 0.35, 0.65),
  c(0.05, 0.45, 0.75),
  c(0.05, 0.50, 0.85)
)

rho_labels <- sapply(rho_grid, function(r) {
  paste0("rho_", paste(sprintf("%0.2f", r), collapse = "_"))
})

rho_meta <- data.frame(
  setting = rho_labels,
  rho12   = sapply(rho_grid, `[`, 1),
  rho13   = sapply(rho_grid, `[`, 2),
  rho23   = sapply(rho_grid, `[`, 3)
) %>%
  mutate(
    ## Frobenius norm of the three unique off-diagonal elements
    rho_frob = sqrt(rho12^2 + rho13^2 + rho23^2)
  )

## ------------------------------------------------
## 2. Locate all simulation summary files
## ------------------------------------------------
files <- list.files(pattern = "^AllOut_rho_.*\\.rds$")
if (length(files) == 0L) stop("No AllOut_rho_*.rds files found in working directory.")

## Containers for all settings
beta_summary_all    <- list()
beta_long_all       <- list()
sigmaRE_summary_all <- list()
sigmae_summary_all  <- list()
corr_summary_all    <- list()

## ------------------------------------------------
## 3. Loop over each rho setting, compute summaries
## ------------------------------------------------
for (f in files) {
  cat("Processing", f, "...\n")
  obj <- readRDS(f)
  
  ## setting label, e.g. "AllOut_rho_0.05_0.10_0.15.rds" -> "rho_0.05_0.10_0.15"
  lab <- gsub("^AllOut_|\\.rds$", "", f)
  
  ## match to rho_meta
  rg <- rho_meta %>% filter(setting == lab)
  if (nrow(rg) != 1) stop("Cannot match setting label in rho_meta for file: ", f)
  
  rho12    <- rg$rho12
  rho13    <- rg$rho13
  rho23    <- rg$rho23
  rho_frob <- rg$rho_frob
  
  ## Extract long dfs from simulation output
  beta_long <- obj$beta_long %>%
    mutate(
      setting  = lab,
      rho12    = rho12,
      rho13    = rho13,
      rho23    = rho23,
      rho_frob = rho_frob
    )
  
  sigmaRE_long <- obj$sigmaRE_long %>%
    mutate(
      setting  = lab,
      rho12    = rho12,
      rho13    = rho13,
      rho23    = rho23,
      rho_frob = rho_frob
    )
  
  sigmae_long <- obj$sigmae_long %>%
    mutate(
      setting  = lab,
      rho12    = rho12,
      rho13    = rho13,
      rho23    = rho23,
      rho_frob = rho_frob
    )
  
  corr_long <- obj$corr_long %>%
    mutate(
      setting  = lab,
      rho12    = rho12,
      rho13    = rho13,
      rho23    = rho23,
      rho_frob = rho_frob
    )
  
  beta_long_all[[lab]] <- beta_long
  
  ## -----------------------------
  ## 3a. Coverage for betas using b_sd
  ## -----------------------------
  results     <- obj$results
  N           <- length(results)
  model_names <- names(results[[1]])  ## OBS_UN, CC_UN, MICE_UN, FULL_UN
  
  ## reconstruct SE matrices for betas from b_sd
  store_se <- lapply(model_names, function(nm) init_mat(px * py, N))
  names(store_se) <- model_names
  
  for (k in seq_len(N)) {
    resk <- results[[k]]
    for (nm in model_names) {
      b_sd_k <- resk[[nm]]$b_sd
      if (is.null(b_sd_k) || length(b_sd_k) == 0) {
        b_sd_k <- rep(NA_real_, px * py)
      }
      store_se[[nm]][, k] <- b_sd_k
    }
  }
  
  beta_se_long <- do.call(
    rbind,
    lapply(names(store_se), function(nm) {
      mk_long(store_se[[nm]], par_names_beta, nm, "SE")
    })
  )
  
  ## Join SE with beta_long
  beta_merged <- beta_long %>%
    left_join(beta_se_long, by = c("Simulation", "Parameter", "Model")) %>%
    mutate(
      Bias   = Estimate - True,
      lower  = Estimate - 1.96 * SE,
      upper  = Estimate + 1.96 * SE,
      covered = (True >= lower & True <= upper)
    )
  
  beta_summary <- beta_merged %>%
    group_by(setting, rho12, rho13, rho23, rho_frob, Model, Parameter) %>%
    summarise(
      bias      = mean(Bias,    na.rm = TRUE),
      variance  = var(Estimate, na.rm = TRUE),
      coverage  = mean(covered, na.rm = TRUE),
      .groups   = "drop"
    )
  
  beta_summary_all[[lab]] <- beta_summary
  
  ## -----------------------------
  ## 3b. Bias & variance for random-effect SDs
  ## -----------------------------
  sigmaRE_summary <- sigmaRE_long %>%
    mutate(Bias = Estimate - True) %>%
    group_by(setting, rho12, rho13, rho23, rho_frob, Model, Parameter) %>%
    summarise(
      bias     = mean(Bias,    na.rm = TRUE),
      variance = var(Estimate, na.rm = TRUE),
      .groups  = "drop"
    )
  sigmaRE_summary_all[[lab]] <- sigmaRE_summary
  
  ## -----------------------------
  ## 3c. Bias & variance for residual SD
  ## -----------------------------
  sigmae_summary <- sigmae_long %>%
    mutate(Bias = Estimate - True) %>%
    group_by(setting, rho12, rho13, rho23, rho_frob, Model, Parameter) %>%
    summarise(
      bias     = mean(Bias,    na.rm = TRUE),
      variance = var(Estimate, na.rm = TRUE),
      .groups  = "drop"
    )
  sigmae_summary_all[[lab]] <- sigmae_summary
  
  ## -----------------------------
  ## 3d. Bias & variance for correlation parameters
  ## -----------------------------
  corr_summary <- corr_long %>%
    mutate(Bias = Estimate - True) %>%
    group_by(setting, rho12, rho13, rho23, rho_frob, Model, Parameter) %>%
    summarise(
      bias     = mean(Bias,    na.rm = TRUE),
      variance = var(Estimate, na.rm = TRUE),
      .groups  = "drop"
    )
  corr_summary_all[[lab]] <- corr_summary
}

## Bind all settings together
beta_summary_all_df    <- bind_rows(beta_summary_all)
beta_long_all_df       <- bind_rows(beta_long_all)
sigmaRE_summary_all_df <- bind_rows(sigmaRE_summary_all)
sigmae_summary_all_df  <- bind_rows(sigmae_summary_all)
corr_summary_all_df    <- bind_rows(corr_summary_all)

## ------------------------------------------------
## 4. Visualizations (x-axis: rho_frob)
## ------------------------------------------------

## 4a. Betas: bias, RE-variance, coverage vs rho_frob

beta_summary_all_df$Parameter <- factor(
  beta_summary_all_df$Parameter,
  levels = par_names_beta
)

beta_summary_all_df <- beta_summary_all_df %>%
  mutate(
    Param_index = as.integer(str_remove(Parameter, "Beta")),
    Outcome     = paste0("Outcome ", ceiling(Param_index / px)),
    Covariate   = paste0("X", ((Param_index - 1) %% px) + 1)
  )

## |Bias(β)| vs correlation strength
p_beta_bias <- ggplot(
  beta_summary_all_df,
  aes(x = rho_frob, y = abs(bias), color = Model, group = Model)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Outcome + Covariate, ncol = px) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Frobenius norm of off-diagonal correlations",
    y = expression("|Bias(" * beta * ")|"),
    title = expression("Absolute bias of " * beta *
                       " estimates across methods and correlation settings")
  )

## Relative efficiency RE = Var(model) / Var(FULL_UN) for β
beta_RE <- beta_summary_all_df %>%
  group_by(setting, rho_frob, Outcome, Covariate) %>%
  mutate(
    var_full = variance[Model == "FULL_UN"],
    RE       = variance / var_full
  ) %>%
  ungroup()

p_beta_var <- ggplot(
  beta_RE,
  aes(x = rho_frob, y = log(RE), color = Model, group = Model)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Outcome + Covariate, ncol = px) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Frobenius norm of off-diagonal correlations",
    y = expression(log(RE(beta))),
    title =  expression("Relative efficiency of " * beta *
                 " estimates across methods and correlation settings")
  )

## Coverage for β
p_beta_cov <- ggplot(
  beta_summary_all_df,
  aes(x = rho_frob, y = coverage, color = Model, group = Model)
) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Outcome + Covariate, ncol = px) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Frobenius norm of off-diagonal correlations",
    y = expression("Empirical coverage of 95% CI for " * beta),
    title = expression("Coverage of nominal 95% Wald CIs for " * beta * 
    " by method and correlation setting")
  )

## 4b. Random-effect SDs: bias & RE-variance

sigmaRE_summary_all_df$Parameter <- factor(
  sigmaRE_summary_all_df$Parameter,
  levels = par_names_sigmaRE
)

p_sigmaRE_bias <- ggplot(
  sigmaRE_summary_all_df,
  aes(x = rho_frob, y = abs(bias), color = Model, group = Model)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Parameter, scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Frobenius norm of off-diagonal correlations",
    y = expression("|Bias(" * sigma[RE] * ")|"),
    title = expression("Absolute bias of random-effect " * sigma * " estimates")
  )

sigmaRE_RE <- sigmaRE_summary_all_df %>%
  group_by(setting, rho_frob, Parameter) %>%
  mutate(
    var_full = variance[Model == "FULL_UN"],
    RE       = variance / var_full
  ) %>%
  ungroup()

p_sigmaRE_var <- ggplot(
  sigmaRE_RE,
  aes(x = rho_frob, y = log(RE), color = Model, group = Model)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Parameter, scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Frobenius norm of off-diagonal correlations",
    y = expression(log(RE(sigma[RE]))),
    title = expression("Relative efficiency of random-effect " * sigma * " estimates")
  )

## 4c. Residual SD: bias & RE-variance

p_sigmae_bias <- ggplot(
  sigmae_summary_all_df,
  aes(x = rho_frob, y = abs(bias), color = Model, group = Model)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Frobenius norm of off-diagonal correlations",
    y = expression("|Bias(" * sigma[e] * ")|"),
    title = expression("Absolute bias of residual " * sigma[e] * " estimates")
  )

sigmae_RE <- sigmae_summary_all_df %>%
  group_by(setting, rho_frob) %>%
  mutate(
    var_full = variance[Model == "FULL_UN"],
    RE       = variance / var_full
  ) %>%
  ungroup()

p_sigmae_var <- ggplot(
  sigmae_RE,
  aes(x = rho_frob, y = log(RE), color = Model, group = Model)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Frobenius norm of off-diagonal correlations",
    y = expression(log(RE(sigma[e]))),
    title = expression("Relative efficiency of residual " * sigma[e] * " estimates")
  )

## 4d. Correlations (rhos): bias & RE-variance

corr_summary_all_df <- corr_summary_all_df %>%
  filter(Parameter %in% c("rho_1_2", "rho_1_3", "rho_2_3")) %>%
  mutate(
    i      = as.integer(substr(Parameter, 5, 5)),
    j      = as.integer(substr(Parameter, 7, 7)),
    offdiag = (i != j),
    Parameter_plot = paste0("rho[", i, ",", j, "]")
  )

corr_offdiag <- corr_summary_all_df %>% filter(offdiag)

p_corr_bias <- ggplot(
  corr_offdiag,
  aes(x = rho_frob, y = abs(bias), color = Model, group = Model)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Parameter_plot, labeller = label_parsed) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Frobenius norm of off-diagonal correlations",
    y = expression("|Bias(" * rho * ")|"),
    title = expression("Absolute bias of off-diagonal " * rho * " estimates")
  )

corr_RE <- corr_offdiag %>%
  group_by(setting, rho_frob, Parameter_plot) %>%
  mutate(
    var_full = variance[Model == "FULL_UN"],
    RE       = variance / var_full
  ) %>%
  ungroup()

p_corr_var <- ggplot(
  corr_RE,
  aes(x = rho_frob, y = log(RE), color = Model, group = Model)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Parameter_plot, labeller = label_parsed) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Frobenius norm of off-diagonal correlations",
    y = expression(log(RE(rho))),
    title = expression("Relative efficiency of off-diagonal " * sigma[e] * " estimates")
  )

## ------------------------------------------------
## 5. Inspect
## ------------------------------------------------
print(p_beta_bias)
print(p_beta_var)
print(p_beta_cov)

print(p_sigmaRE_bias)
print(p_sigmaRE_var)

gridExtra::grid.arrange(p_sigmae_bias, p_sigmae_var, ncol = 1)
gridExtra::grid.arrange(p_corr_bias, p_corr_var, ncol = 1)
