rm(list=ls())
# ==============================================================================
# File: analyze_sim_pattern_levels.R
# Usage: Rscript analyze_sim_pattern_levels.R
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(gridExtra)
  library(stringr)
})

# ----------------------------------------------------------------
# 1. SETUP & TRUTH DEFINITIONS
# ----------------------------------------------------------------
H <- 10; px <- 2; py <- 3

# True Beta (Fixed)
beta_true_mat <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)
beta_true_vec <- as.vector(beta_true_mat)
names(beta_true_vec) <- paste0("Beta", 1:(px*py))

# True Variances (Fixed)
sigma_u_true <- 1
sigma_v_true <- 3 
sigma_e_true <- 3

# ----------------------------------------------------------------
# 2. READ DATA & PREP
# ----------------------------------------------------------------
out_dir <- "results"
files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)

if(length(files) == 0) {
  warning("No .rds files found. Ensure 'results/' exists.")
} else {
  cat(sprintf("Reading %d files...\n", length(files)))
  dat <- do.call(rbind, lapply(files, readRDS))
  
  # --- Deduplicate if necessary ---
  dat <- dat %>% distinct(SimID, Scenario, Model, .keep_all = TRUE)
  
  # --- CHANGE 1: RENAME OBSERVED -> PATTERN ---
  # Do this before creating factors
  dat$Model[dat$Model == "Observed"] <- "Pattern"
  
  # --- Parse Scenario to get True Rho Values ---
  rho_vals <- str_match(dat$Scenario, "Rho_([0-9.-]+)_([0-9.-]+)_([0-9.-]+)")
  
  if(all(is.na(rho_vals))) {
    # Fallback defaults
    dat$True_rho12 <- -0.3; dat$True_rho13 <- 0.3; dat$True_rho23 <- 0.8
    dat$Scenario_F <- factor(dat$Scenario, levels=sort(unique(dat$Scenario)))
    level_labels <- setNames(unique(dat$Scenario), unique(dat$Scenario))
  } else {
    dat$True_rho12 <- as.numeric(rho_vals[,2])
    dat$True_rho13 <- as.numeric(rho_vals[,3])
    dat$True_rho23 <- as.numeric(rho_vals[,4])
    
    # --- Map to "Correlation Levels" ---
    dat <- dat %>% mutate(rho_sum = abs(True_rho12) + abs(True_rho13) + abs(True_rho23))
    
    # Identify unique scenarios and sort by strength
    unique_scens <- unique(dat$Scenario[order(dat$rho_sum)])
    
    # Create mapping: "Rho_0.1_..." -> "Level 1"
    level_labels <- seq_along(unique_scens)
    names(level_labels) <- unique_scens
    
    # Set Factor Levels
    dat$Scenario_F <- factor(dat$Scenario, levels = unique_scens)
  }
  
  # --- CHANGE 2: UPDATE FACTOR LEVELS ---
  # Ensure Model Factor Order: Full -> Pattern -> MICE -> CC
  dat$Model <- factor(dat$Model, levels = c("Full", "Pattern", "MICE", "CC"))
}

# ----------------------------------------------------------------
# 3. METRICS COMPUTATION
# ----------------------------------------------------------------
if(exists("dat")) {
  
  # --- Helper for X-Axis Labels ---
  clean_scenario_labels <- function(x) {
    if(exists("level_labels")) {
      ifelse(x %in% names(level_labels), level_labels[x], x)
    } else {
      x
    }
  }
  
  format_greek <- function(x) {
    x <- gsub("^Beta([0-9]+)$", "beta[\\1]", x)
    x <- gsub("^sigma_([a-z]+)$", "sigma[\\1]", x)
    x <- gsub("^rho([0-9]+)$", "rho[\\1]", x)
    return(x)
  }
  
  # --- BETA METRICS ---
  beta_cols <- grep("^Beta[0-9]+$", names(dat), value = TRUE)
  se_cols   <- grep("^SE[0-9]+$",   names(dat), value = TRUE)
  
  long_est <- dat %>%
    dplyr::select(SimID, Scenario_F, Model, all_of(beta_cols)) %>%
    pivot_longer(cols = all_of(beta_cols), names_to = "Parameter", values_to = "Estimate")
  
  long_se <- dat %>%
    dplyr::select(SimID, Scenario_F, Model, all_of(se_cols)) %>%
    pivot_longer(cols = all_of(se_cols), names_to = "Parameter_SE", values_to = "SE") %>%
    mutate(Parameter = sub("SE", "Beta", Parameter_SE)) %>%
    dplyr::select(-Parameter_SE)
  
  beta_long <- left_join(long_est, long_se, by = c("SimID", "Scenario_F", "Model", "Parameter"))
  
  beta_metrics <- beta_long %>%
    group_by(Scenario_F, Model, Parameter) %>%
    summarise(
      Bias     = mean(Estimate) - beta_true_vec[unique(Parameter)],
      Var_Est  = var(Estimate, na.rm=TRUE),
      Coverage = mean((Estimate - 1.96*SE) <= beta_true_vec[unique(Parameter)] &
                        (Estimate + 1.96*SE) >= beta_true_vec[unique(Parameter)], na.rm=TRUE),
      .groups  = "drop"
    )
  
  # --- CHANGE 3: BETA RE (Use Pattern) ---
  beta_re <- beta_metrics %>%
    dplyr::select(Scenario_F, Parameter, Model, Var_Est) %>%
    pivot_wider(names_from = Model, values_from = Var_Est) %>%
    mutate(
      RE_CC       = log(CC / Full),
      RE_MICE     = log(MICE / Full),
      RE_Pattern  = log(Pattern / Full), # Changed from Observed
      RE_Full     = 0
    ) %>%
    pivot_longer(starts_with("RE_"), names_to="Model", values_to="RE") %>%
    mutate(Model = sub("RE_", "", Model)) %>%
    dplyr::select(Scenario_F, Parameter, Model, RE)
  
  beta_final <- left_join(beta_metrics, beta_re, by = c("Scenario_F", "Parameter", "Model")) %>%
    mutate(Facet_Label = format_greek(Parameter),
           Facet_Label = factor(Facet_Label, levels = format_greek(names(beta_true_vec))))
  
  # --- SIGMA & RHO METRICS ---
  sigma_long <- dat %>%
    dplyr::select(Scenario_F, Model, sigma_u, sigma_v, sigma_e, rho12, rho13, rho23, 
                  True_rho12, True_rho13, True_rho23) %>%
    pivot_longer(cols = c(sigma_u, sigma_v, sigma_e, rho12, rho13, rho23), 
                 names_to = "Parameter", values_to = "Estimate") %>%
    mutate(
      True_Val = case_when(
        Parameter == "sigma_u" ~ sigma_u_true,
        Parameter == "sigma_v" ~ sigma_v_true,
        Parameter == "sigma_e" ~ sigma_e_true,
        Parameter == "rho12"   ~ True_rho12,
        Parameter == "rho13"   ~ True_rho13,
        Parameter == "rho23"   ~ True_rho23
      )
    )
  
  sigma_metrics <- sigma_long %>%
    group_by(Scenario_F, Model, Parameter) %>%
    summarise(
      Bias    = mean(Estimate) - mean(True_Val),
      Var_Est = var(Estimate, na.rm=TRUE),
      .groups = "drop"
    )
  
  # --- CHANGE 4: SIGMA RE (Use Pattern) ---
  sigma_re <- sigma_metrics %>%
    dplyr::select(Scenario_F, Parameter, Model, Var_Est) %>%
    pivot_wider(names_from = Model, values_from = Var_Est) %>%
    mutate(
      RE_CC       = log(CC / Full),
      RE_MICE     = log(MICE / Full),
      RE_Pattern  = log(Pattern / Full), # Changed from Observed
      RE_Full     = 0
    ) %>%
    pivot_longer(starts_with("RE_"), names_to="Model", values_to="RE") %>%
    mutate(Model = sub("RE_", "", Model)) %>%
    dplyr::select(Scenario_F, Parameter, Model, RE)
  
  sigma_final <- left_join(sigma_metrics, sigma_re, by = c("Scenario_F", "Parameter", "Model")) %>%
    mutate(Facet_Label = format_greek(Parameter))
  
  # ----------------------------------------------------------------
  # 4. PLOTS
  # ----------------------------------------------------------------
  if(!dir.exists("plots")) dir.create("plots")
  
  # --- Professional Theme ---
  theme_professional <- theme_bw(base_size = 14) + 
    theme(
      text = element_text(family = "serif"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5), 
      strip.background = element_rect(fill = "grey95", color = "grey20"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.box.margin = margin(t = -5),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face="bold")
    )
  
  # --- CHANGE 5: COLORS (Use Pattern) ---
  my_colors <- c("Full"="black", "Pattern"="#1b9e77", "MICE"="#377eb8", "CC"="#e41a1c")
  
  # 1. BETA BIAS
  p1 <- ggplot(beta_final, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") + 
    geom_line(linewidth=0.7) + geom_point(size=2) +
    facet_wrap(~Facet_Label, labeller = label_parsed) + 
    scale_color_manual(values=my_colors) +
    scale_x_discrete(labels = clean_scenario_labels) +
    labs(title=expression(bold("Bias of Regression Coefficients")),
         x="Correlation Strength Level (Low -> High)", y="Bias") + theme_professional
  ggsave("plots/1_Beta_Bias.pdf", p1, width=12, height=8)
  
  # 2. BETA COVERAGE
  p2 <- ggplot(beta_final, aes(x=Scenario_F, y=Coverage, color=Model, group=Model)) +
    geom_hline(yintercept=0.95, linetype="dashed", color="grey50") + 
    geom_line(linewidth=0.7) + geom_point(size=2) +
    facet_wrap(~Facet_Label, labeller = label_parsed) + 
    coord_cartesian(ylim=c(0.7, 1.0)) + 
    scale_color_manual(values=my_colors) +
    scale_x_discrete(labels = clean_scenario_labels) +
    labs(title=expression(bold("Coverage Probability (95% CI)")),
         x="Correlation Strength Level (Low -> High)", y="Coverage Rate") + theme_professional
  ggsave("plots/2_Beta_Coverage.pdf", p2, width=12, height=8)
  
  # 3. BETA LOG RE
  p3 <- ggplot(beta_final, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") + 
    geom_line(linewidth=0.7) + geom_point(size=2) +
    facet_wrap(~Facet_Label, labeller = label_parsed) + 
    scale_color_manual(values=my_colors) +
    scale_x_discrete(labels = clean_scenario_labels) +
    labs(title=expression(bold("Log Relative Variance")), 
         subtitle = "Value > 0 indicates higher variance than Full Data",
         x="Correlation Strength Level (Low -> High)", y="Log Variance Ratio") + theme_professional
  ggsave("plots/3_Beta_LogRE.pdf", p3, width=12, height=8)
  
  # 4 & 5. SIGMA & RHO METRICS
  sigma_subset <- sigma_final %>% filter(grepl("sigma", Parameter))
  rho_subset   <- sigma_final %>% filter(grepl("rho", Parameter))
  
  # Plot Helpers
  plot_bias <- function(data, title) {
    ggplot(data, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
      geom_hline(yintercept=0, linetype="dashed", color="grey50") + 
      geom_line() + geom_point() +
      facet_wrap(~Facet_Label, labeller = label_parsed) + 
      scale_color_manual(values=my_colors) +
      scale_x_discrete(labels=NULL) +
      labs(title=title, x=NULL, y="Bias") + theme_professional
  }
  
  plot_re <- function(data) {
    ggplot(data, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
      geom_hline(yintercept=0, linetype="dashed", color="grey50") + 
      geom_line() + geom_point() +
      facet_wrap(~Facet_Label, labeller = label_parsed) + 
      scale_color_manual(values=my_colors) +
      scale_x_discrete(labels=clean_scenario_labels) +
      labs(title="Log Relative Variance", x="Correlation Strength Level (Low -> High)", y="Log Var Ratio") + 
      theme_professional
  }
  
  ggsave("plots/4_Sigma_Metrics.pdf", 
         grid.arrange(plot_bias(sigma_subset, expression(bold("Variance Components Bias"))), 
                      plot_re(sigma_subset), nrow=2, heights=c(1, 1.2)), width=12, height=10)
  
  ggsave("plots/5_Rho_Metrics.pdf", 
         grid.arrange(plot_bias(rho_subset, expression(bold("Correlation Bias"))), 
                      plot_re(rho_subset), nrow=2, heights=c(1, 1.2)), width=12, height=10)
  
  cat("Analysis complete. Professional plots (Pattern/Levels) saved in 'plots/'.\n")
}