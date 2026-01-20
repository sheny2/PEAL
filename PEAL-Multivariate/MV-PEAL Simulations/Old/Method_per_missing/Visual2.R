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
## 0. Global info
## ------------------------------------------------
H  <- 5          # number of site-specific random effects
px <- 4          # number of covariates per outcome
py <- 3          # number of outcomes

par_names_beta    <- paste0("Beta", 1:(px * py))
par_names_sigmaRE <- c("sigma_u", paste0("sigma_v_", 1:H))
par_name_sigma_e  <- "sigma_e"
corr_param_names  <- paste0("rho_", rep(1:py, each = py), "_", rep(1:py, times = py))

## Helper: init matrix and mk_long (matching sim script)
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
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

## ------------------------------------------------
## 1. Recreate miss_grid to attach missingness info
## ------------------------------------------------
miss_grid <- data.frame(
  setting  = paste0("miss_", 1:8),
  p1_base  = seq(0.1, 0.2, length.out = 8),
  p2_base  = seq(0.1, 0.4, length.out = 8),
  p3_base  = seq(0.1, 0.6, length.out = 8)
) %>%
  mutate(
    missing_base_avg = round((p1_base + p2_base + p3_base) / 3, 2)
  )

## ------------------------------------------------
## 2. Locate all simulation summary files
## ------------------------------------------------
files <- list.files(pattern = "^AllOut_miss_.*\\.rds$")
if (length(files) == 0L) stop("No AllOut_miss_*.rds files found in working directory.")

## Containers for all settings
beta_summary_all    <- list()
beta_long_all       <- list()
sigmaRE_summary_all <- list()
sigmae_summary_all  <- list()
corr_summary_all    <- list()

## ------------------------------------------------
## 3. Loop over each missingness setting, compute summaries
## ------------------------------------------------
for (f in files) {
  cat("Processing", f, "...\n")
  obj <- readRDS(f)
  
  # Get setting label from file name, e.g. "AllOut_miss_1.rds" -> "miss_1"
  lab <- gsub("^AllOut_|\\.rds$", "", f)
  
  # Match to miss_grid
  mg <- miss_grid %>% filter(setting == lab)
  if (nrow(mg) != 1) stop("Cannot match setting label in miss_grid for file: ", f)
  
  missing_avg <- mg$missing_base_avg
  
  # Extract long dfs from simulation output
  beta_long    <- obj$beta_long    %>%
    mutate(setting = lab, missing_base_avg = missing_avg)
  sigmaRE_long <- obj$sigmaRE_long %>%
    mutate(setting = lab, missing_base_avg = missing_avg)
  sigmae_long  <- obj$sigmae_long  %>%
    mutate(setting = lab, missing_base_avg = missing_avg)
  corr_long    <- obj$corr_long    %>%
    mutate(setting = lab, missing_base_avg = missing_avg)
  
  beta_long_all[[lab]] <- beta_long
  
  ## -----------------------------
  ## 3a. Coverage for betas using b_sd
  ## -----------------------------
  results     <- obj$results
  N           <- length(results)
  model_names <- names(results[[1]])  # e.g., OBS_UN, CC_UN, MICE_UN, FULL_UN
  
  # Reconstruct SE matrices for betas from results (b_sd)
  store_se <- lapply(model_names, function(nm) init_mat(px * py, N))
  names(store_se) <- model_names
  
  for (k in seq_len(N)) {
    resk <- results[[k]]
    for (nm in model_names) {
      b_sd_k <- resk[[nm]]$b_sd
      # Fallback if b_sd missing
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
  
  # Join SE with beta_long (by Simulation, Parameter, Model)
  beta_merged <- beta_long %>%
    left_join(beta_se_long, by = c("Simulation", "Parameter", "Model")) %>%
    mutate(
      Bias   = Estimate - True,
      lower  = Estimate - 1.96 * SE,
      upper  = Estimate + 1.96 * SE,
      covered = (True >= lower & True <= upper)
    )
  
  beta_summary <- beta_merged %>%
    group_by(setting, missing_base_avg, Model, Parameter) %>%
    summarise(
      bias      = mean(Bias, na.rm = TRUE),
      variance  = var(Estimate, na.rm = TRUE),
      coverage  = mean(covered, na.rm = TRUE),
      .groups   = "drop"
    )
  
  beta_summary_all[[lab]] <- beta_summary
  
  ## -----------------------------
  ## 3b. Bias & variance for random-effects SDs (σ_RE)
  ## -----------------------------
  sigmaRE_summary <- sigmaRE_long %>%
    mutate(Bias = Estimate - True) %>%
    group_by(setting, missing_base_avg, Model, Parameter) %>%
    summarise(
      bias     = mean(Bias, na.rm = TRUE),
      variance = var(Estimate, na.rm = TRUE),
      .groups  = "drop"
    )
  sigmaRE_summary_all[[lab]] <- sigmaRE_summary
  
  ## -----------------------------
  ## 3c. Bias & variance for residual SD (σ_e)
  ## -----------------------------
  sigmae_summary <- sigmae_long %>%
    mutate(Bias = Estimate - True) %>%
    group_by(setting, missing_base_avg, Model, Parameter) %>%
    summarise(
      bias     = mean(Bias, na.rm = TRUE),
      variance = var(Estimate, na.rm = TRUE),
      .groups  = "drop"
    )
  sigmae_summary_all[[lab]] <- sigmae_summary
  
  ## -----------------------------
  ## 3d. Bias & variance for correlation parameters (ρ)
  ## -----------------------------
  corr_summary <- corr_long %>%
    mutate(Bias = Estimate - True) %>%
    group_by(setting, missing_base_avg, Model, Parameter) %>%
    summarise(
      bias     = mean(Bias, na.rm = TRUE),
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
## 4. Prepare factor / labeling structure
## ------------------------------------------------

beta_summary_all_df <- beta_summary_all_df %>%
  mutate(
    Param_index = as.integer(str_remove(Parameter, "Beta")),
    Outcome     = paste0("Outcome ", ceiling(Param_index / px)),
    Covariate   = paste0("X", ((Param_index - 1) %% px) + 1)
  )

beta_summary_all_df$Parameter <- factor(beta_summary_all_df$Parameter,
                                        levels = par_names_beta)

sigmaRE_summary_all_df$Parameter <- factor(sigmaRE_summary_all_df$Parameter,
                                           levels = par_names_sigmaRE)

## For correlations, identify off-diagonal ρ_ij and create a parsed label
corr_summary_all_df <- corr_summary_all_df %>%
  mutate(
    i       = as.integer(substr(Parameter, 5, 5)),
    j       = as.integer(substr(Parameter, 7, 7)),
    offdiag = (i != j),
    Parameter_plot = paste0("rho[", i, ",", j, "]")
  )

corr_offdiag <- corr_summary_all_df %>%
  filter(offdiag, Parameter %in% c("rho_1_2", "rho_1_3", "rho_2_3"))

## ------------------------------------------------
## 5. Plots
## ------------------------------------------------

## 5a. Betas: bias, RE-variance, coverage vs missingness

# |Bias(β)| vs missingness
p_beta_bias <- ggplot(beta_summary_all_df,
                      aes(x = missing_base_avg,
                          y = abs(bias),
                          color = Model,
                          group = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Outcome + Covariate, ncol = px) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Average baseline missingness rate",
    y = expression("|Bias(" * beta * ")|"),
    title = expression("Absolute bias of " * beta * " estimates across methods and missingness")
  )

# Relative efficiency RE = Var(model) / Var(FULL_UN)
beta_RE <- beta_summary_all_df %>%
  group_by(setting, missing_base_avg, Outcome, Covariate) %>%
  mutate(
    var_full = variance[Model == "FULL_UN"],
    RE       = variance / var_full
  ) %>%
  ungroup()

p_beta_var <- ggplot(beta_RE,
                     aes(x = missing_base_avg,
                         y = log(RE),
                         color = Model,
                         group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Outcome + Covariate, ncol = px) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Average baseline missingness rate",
    y = expression(log("RE(" * beta * ")")),
    title = expression("Relative efficiency of " * beta * " estimates vs. FULL[UN]")
  )

# Coverage for β (95% CI)
p_beta_cov <- ggplot(beta_summary_all_df,
                     aes(x = missing_base_avg,
                         y = coverage,
                         color = Model,
                         group = Model)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Outcome + Covariate, ncol = px) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Average baseline missingness rate",
    y = expression("Empirical coverage of 95% CI for " * beta),
    title = expression("Coverage of nominal 95% Wald CIs for " * beta * " by method")
  )

## 5b. Random-effect SDs σ_RE: bias & RE-variance

p_sigmaRE_bias <- ggplot(sigmaRE_summary_all_df,
                         aes(x = missing_base_avg,
                             y = abs(bias),
                             color = Model,
                             group = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Parameter) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Average baseline missingness rate",
    y = expression("|Bias(" * sigma[RE] * ")|"),
    title = expression("Absolute bias of random-effect " * sigma * " estimates")
  )

sigmaRE_RE <- sigmaRE_summary_all_df %>%
  group_by(setting, missing_base_avg, Parameter) %>%
  mutate(
    var_full = variance[Model == "FULL_UN"],
    RE       = variance / var_full
  ) %>%
  ungroup()

p_sigmaRE_var <- ggplot(sigmaRE_RE,
                        aes(x = missing_base_avg,
                            y = log(RE),
                            color = Model,
                            group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Parameter) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Average baseline missingness rate",
    y = expression(log("RE(" * sigma[RE] * ")")),
    title = expression("Relative efficiency of random-effect " * sigma * " estimates vs. FULL[UN]")
  )

## 5c. Residual SD σ_e: bias & RE-variance

p_sigmae_bias <- ggplot(sigmae_summary_all_df,
                        aes(x = missing_base_avg,
                            y = abs(bias),
                            color = Model,
                            group = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Average baseline missingness rate",
    y = expression("|Bias(" * sigma[e] * ")|"),
    title = expression("Absolute bias of residual " * sigma[e] * " estimates")
  )

sigmae_RE <- sigmae_summary_all_df %>%
  group_by(setting, missing_base_avg) %>%
  mutate(
    var_full = variance[Model == "FULL_UN"],
    RE       = variance / var_full
  ) %>%
  ungroup()

p_sigmae_var <- ggplot(sigmae_RE,
                       aes(x = missing_base_avg,
                           y = log(RE),
                           color = Model,
                           group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Average baseline missingness rate",
    y = expression(log("RE(" * sigma[e] * ")")),
    title = expression("Relative efficiency of residual " * sigma[e] *
                         " estimates vs. FULL[UN]")
  )

## 5d. Correlations ρ: bias & RE-variance

p_corr_bias <- ggplot(corr_offdiag,
                      aes(x = missing_base_avg,
                          y = abs(bias),
                          color = Model,
                          group = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Parameter_plot, labeller = label_parsed) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Average baseline missingness rate",
    y = expression("|Bias(" * rho * ")|"),
    title = expression("Absolute bias of off-diagonal " * rho * " estimates")
  )

corr_RE <- corr_offdiag %>%
  group_by(setting, missing_base_avg, Parameter_plot) %>%
  mutate(
    var_full = variance[Model == "FULL_UN"],
    RE       = variance / var_full
  ) %>%
  ungroup()

p_corr_var <- ggplot(corr_RE,
                     aes(x = missing_base_avg,
                         y = log(RE),
                         color = Model,
                         group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ Parameter_plot, labeller = label_parsed) +
  scale_color_brewer(palette = "Dark2") +
  theme_sim() +
  labs(
    x = "Average baseline missingness rate",
    y = expression(log("RE(" * rho * ")")),
    title = expression("Relative efficiency of off-diagonal " * rho *
                         " estimates vs. FULL[UN]")
  )

## ------------------------------------------------
## 6. Print / arrange selected plots
## ------------------------------------------------

print(p_beta_bias)
print(p_beta_var)
print(p_beta_cov)

print(p_sigmaRE_bias)
print(p_sigmaRE_var)

gridExtra::grid.arrange(p_sigmae_bias, p_sigmae_var, ncol = 1)
gridExtra::grid.arrange(p_corr_bias, p_corr_var, ncol = 1)
