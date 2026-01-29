rm(list=ls())
# ==============================================================================
# File: analyze_sim_specs.R
# Usage: Rscript analyze_sim_specs.R
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

# True Beta
beta_true_mat <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)
beta_true_vec <- as.vector(beta_true_mat)
names(beta_true_vec) <- paste0("Beta", 1:(px*py))

# True Variance & Rho
true_sigma_vals <- c(
  "sigma_u" = 1, "sigma_v" = 3, "sigma_e" = 3,
  "rho12" = -0.3, "rho13" = 0.3, "rho23" = 0.8
)

# ----------------------------------------------------------------
# 2. READ DATA & CLEANING
# ----------------------------------------------------------------
out_dir <- "results"
files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)

if(length(files) == 0) {
  warning("No .rds files found in 'results/'.")
} else {
  cat(sprintf("Reading %d files...\n", length(files)))
  dat <- do.call(rbind, lapply(files, readRDS))
  
  # --- Remove Duplicates (Fixes many-to-many join warnings) ---
  dat <- dat %>% distinct(SimID, Scenario, Model, .keep_all = TRUE)
  
  # --- Extract Missingness Level ---
  dat$Miss_Level <- as.numeric(str_extract(dat$Scenario, "[0-9]+$"))
  dat$Scenario_F <- factor(dat$Scenario, levels = paste0("Miss_Lev_", sort(unique(dat$Miss_Level))))
  
  # --- Handle Specifications (Models) ---
  # Assuming the 'Model' column contains: "independence", "exchangeable", "unstructured"
  # We set "unstructured" as the reference level
  avail_models <- unique(dat$Model)
  cat("Specifications found:", paste(avail_models, collapse=", "), "\n")
  
  # Set factor order: Unstructured first (baseline), then others
  dat$Model <- factor(dat$Model, levels = c("unstructured", "exchangeable", "independence"))
}

# ----------------------------------------------------------------
# 3. METRICS COMPUTATION
# ----------------------------------------------------------------
if(exists("dat")) {
  
  # --- Greek Label Formatter ---
  format_greek <- function(x) {
    x <- gsub("^Beta([0-9]+)$", "beta[\\1]", x)
    x <- gsub("^sigma_([a-z]+)$", "sigma[\\1]", x)
    x <- gsub("^rho([0-9]+)$", "rho[\\1]", x)
    return(x)
  }
  
  clean_scenario <- function(x) gsub("Miss_Lev_", "", x)
  
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
  
  # --- BETA RE (Ratio of Variance vs Unstructured) ---
  # RE = Var(Spec) / Var(Unstructured)
  beta_re <- beta_metrics %>%
    dplyr::select(Scenario_F, Parameter, Model, Var_Est) %>%
    pivot_wider(names_from = Model, values_from = Var_Est)
  
  # Dynamic calculation based on available columns
  # We try to use 'unstructured' as the baseline denominator
  if("unstructured" %in% names(beta_re)) {
    cols_to_calc <- setdiff(names(beta_re), c("Scenario_F", "Parameter", "unstructured"))
    for(m in cols_to_calc) {
      beta_re[[paste0("RE_", m)]] <- beta_re[[m]] / beta_re[["unstructured"]]
    }
    beta_re$RE_unstructured <- 1 # Baseline
    
    beta_re <- beta_re %>%
      dplyr::select(Scenario_F, Parameter, starts_with("RE_")) %>%
      pivot_longer(cols = starts_with("RE_"), names_to = "Model", values_to = "RE") %>%
      mutate(Model = sub("RE_", "", Model))
    
    beta_final <- left_join(beta_metrics, beta_re, by = c("Scenario_F", "Parameter", "Model"))
  } else {
    beta_final <- beta_metrics %>% mutate(RE = NA)
  }
  
  beta_final <- beta_final %>%
    mutate(Facet_Label = factor(format_greek(Parameter), levels = format_greek(names(beta_true_vec))))
  
  # --- SIGMA & RHO METRICS ---
  sigma_long <- dat %>%
    dplyr::select(Scenario_F, Model, all_of(names(true_sigma_vals))) %>%
    pivot_longer(cols = all_of(names(true_sigma_vals)), names_to = "Parameter", values_to = "Estimate") %>%
    mutate(True_Val = true_sigma_vals[Parameter])
  
  sigma_metrics <- sigma_long %>%
    group_by(Scenario_F, Model, Parameter) %>%
    summarise(
      Bias    = mean(Estimate) - mean(True_Val),
      Var_Est = var(Estimate, na.rm=TRUE),
      .groups = "drop"
    )
  
  # --- SIGMA RE ---
  sigma_re <- sigma_metrics %>%
    dplyr::select(Scenario_F, Parameter, Model, Var_Est) %>%
    pivot_wider(names_from = Model, values_from = Var_Est)
  
  if("unstructured" %in% names(sigma_re)) {
    cols_to_calc <- setdiff(names(sigma_re), c("Scenario_F", "Parameter", "unstructured"))
    for(m in cols_to_calc) {
      sigma_re[[paste0("RE_", m)]] <- sigma_re[[m]] / sigma_re[["unstructured"]]
    }
    sigma_re$RE_unstructured <- 1
    
    sigma_re <- sigma_re %>%
      dplyr::select(Scenario_F, Parameter, starts_with("RE_")) %>%
      pivot_longer(cols = starts_with("RE_"), names_to = "Model", values_to = "RE") %>%
      mutate(Model = sub("RE_", "", Model))
    
    sigma_final <- left_join(sigma_metrics, sigma_re, by = c("Scenario_F", "Parameter", "Model"))
  } else {
    sigma_final <- sigma_metrics %>% mutate(RE = NA)
  }
  
  sigma_final <- sigma_final %>% mutate(Facet_Label = format_greek(Parameter))
  
  # ----------------------------------------------------------------
  # 4. PLOTS
  # ----------------------------------------------------------------
  if(!dir.exists("plots")) dir.create("plots")
  
  # Professional Theme
  theme_professional <- theme_bw(base_size = 14) + 
    theme(
      text = element_text(family = "serif"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey20"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face="bold")
    )
  
  # Colors for Specifications
  # Unstructured = Black (Baseline), Exchangeable = Blue, Independence = Red
  my_colors <- c("unstructured"="black", "exchangeable"="#377eb8", "independence"="#e41a1c")
  
  # 1. BETA BIAS
  p1 <- ggplot(beta_final, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") + 
    geom_line(linewidth=0.7) + geom_point(size=2) +
    facet_wrap(~Facet_Label, scales="free_y", labeller = label_parsed) + 
    scale_color_manual(values=my_colors) +
    scale_x_discrete(labels = clean_scenario) +
    labs(title=expression(bold("Bias of Regression Coefficients") ~ (hat(beta) - beta)),
         x="Missingness Scenario", y="Bias") + theme_professional
  ggsave("plots/1_Beta_Bias.pdf", p1, width=12, height=8)
  
  # 2. BETA COVERAGE
  p2 <- ggplot(beta_final, aes(x=Scenario_F, y=Coverage, color=Model, group=Model)) +
    geom_hline(yintercept=0.95, linetype="dashed", color="grey50") + 
    geom_line(linewidth=0.7) + geom_point(size=2) +
    facet_wrap(~Facet_Label, labeller = label_parsed) + 
    coord_cartesian(ylim=c(0.7, 1.0)) + 
    scale_color_manual(values=my_colors) +
    scale_x_discrete(labels = clean_scenario) +
    labs(title=expression(bold("Coverage Probability (95% CI)")),
         x="Missingness Scenario", y="Coverage Rate") + theme_professional
  ggsave("plots/2_Beta_Coverage.pdf", p2, width=12, height=8)
  
  # 3. BETA RE (Variance Ratio)
  p3 <- ggplot(beta_final, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
    geom_hline(yintercept=1, linetype="dashed", color="grey50") + 
    geom_line(linewidth=0.7) + geom_point(size=2) +
    facet_wrap(~Facet_Label, labeller = label_parsed) + 
    scale_color_manual(values=my_colors) +
    scale_x_discrete(labels = clean_scenario) +
    labs(title=expression(bold("Relative Efficiency (Variance Ratio)")), 
         subtitle = expression(Var(hat(beta)[spec]) / Var(hat(beta)[unstructured])),
         x="Missingness Scenario", y="Variance Ratio") + theme_professional
  ggsave("plots/3_Beta_RE.pdf", p3, width=12, height=8)
  
  # 4 & 5. SIGMA & RHO
  sigma_subset <- sigma_final %>% filter(grepl("sigma", Parameter))
  rho_subset   <- sigma_final %>% filter(grepl("rho", Parameter))
  
  # Plot Helpers
  plot_bias <- function(d, t) {
    ggplot(d, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
      geom_hline(yintercept=0, linetype="dashed", color="grey50") + 
      geom_line() + geom_point() +
      facet_wrap(~Facet_Label, scales="free_y", labeller = label_parsed) + 
      scale_color_manual(values=my_colors) +
      scale_x_discrete(labels=NULL) +
      labs(title=t, x=NULL, y="Bias") + theme_professional
  }
  
  plot_re <- function(d) {
    ggplot(d, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
      geom_hline(yintercept=1, linetype="dashed", color="grey50") + 
      geom_line() + geom_point() +
      facet_wrap(~Facet_Label, scales="free_y", labeller = label_parsed) + 
      scale_color_manual(values=my_colors) +
      scale_x_discrete(labels=clean_scenario) +
      labs(title="Relative Efficiency (Variance Ratio)", x="Missingness Scenario", y="Var Ratio") + 
      theme_professional
  }
  
  ggsave("plots/4_Sigma_Metrics.pdf", 
         grid.arrange(plot_bias(sigma_subset, expression(bold("Variance Components Bias"))),
                      plot_re(sigma_subset), nrow=2, heights=c(1,1.2)), width=12, height=10)
  
  ggsave("plots/5_Rho_Metrics.pdf", 
         grid.arrange(plot_bias(rho_subset, expression(bold("Correlation Bias"))),
                      plot_re(rho_subset), nrow=2, heights=c(1,1.2)), width=12, height=10)
  
  cat("Analysis complete. Plots saved in 'plots/'.\n")
}