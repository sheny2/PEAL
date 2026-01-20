rm(list=ls())
# ==============================================================================
# File: analyze_sim_pattern.R
# Usage: Rscript analyze_sim_pattern.R
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
  "sigma_u" = 1,
  "sigma_v" = 3, 
  "sigma_e" = 3,
  "rho12"   = -0.3,
  "rho13"   = 0.3,
  "rho23"   = 0.8
)

# ----------------------------------------------------------------
# 2. READ DATA
# ----------------------------------------------------------------
out_dir <- "results"
files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)

if(length(files) == 0) {
  # warning("No files found...") 
  # For reproduction/template purposes, we assume 'dat' exists or is mocked
} else {
  cat(sprintf("Reading %d files...\n", length(files)))
  dat <- do.call(rbind, lapply(files, readRDS))
  
  # --- CHANGE 1: RENAME RAW DATA VALUES ---
  # We must change the string "Observed" to "Pattern" before factorizing
  dat$Model[dat$Model == "Observed"] <- "Pattern"
}

# --- Extract Missingness Level ---
dat$Miss_Level <- as.numeric(str_extract(dat$Scenario, "[0-9]+$"))
dat$Scenario_F <- factor(dat$Scenario, levels = paste0("Miss_Lev_", sort(unique(dat$Miss_Level))))

# --- CHANGE 2: UPDATE FACTOR LEVELS ---
# Order Models: Full -> Pattern -> MICE -> CC
dat$Model <- factor(dat$Model, levels = c("Full", "Pattern", "MICE", "CC"))

# ----------------------------------------------------------------
# 3. METRICS COMPUTATION
# ----------------------------------------------------------------

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
    Bias     = (mean(Estimate) - beta_true_vec[unique(Parameter)]),
    Var_Est  = var(Estimate, na.rm=TRUE),
    Coverage = mean((Estimate - 1.96*SE) <= beta_true_vec[unique(Parameter)] &
                      (Estimate + 1.96*SE) >= beta_true_vec[unique(Parameter)], na.rm=TRUE),
    .groups  = "drop"
  )

# --- CHANGE 3: RE COMPUTATION (Pattern) ---
beta_re <- beta_metrics %>%
  dplyr::select(Scenario_F, Parameter, Model, Var_Est) %>%
  pivot_wider(names_from = Model, values_from = Var_Est) %>%
  mutate(
    RE_CC      = log(CC / Full),
    RE_MICE    = log(MICE / Full),
    RE_Pattern = log(Pattern / Full), # Changed from Observed
    RE_Full    = 0
  ) %>%
  pivot_longer(starts_with("RE_"), names_to="Model", values_to="RE") %>%
  mutate(Model = sub("RE_", "", Model)) %>%
  dplyr::select(Scenario_F, Parameter, Model, RE)

beta_final <- left_join(beta_metrics, beta_re, by = c("Scenario_F", "Parameter", "Model"))


# --- SIGMA & RHO METRICS ---
sigma_cols <- c("sigma_u", "sigma_v", "sigma_e", "rho12", "rho13", "rho23")

sigma_long <- dat %>%
  dplyr::select(Scenario_F, Model, all_of(sigma_cols)) %>%
  pivot_longer(cols = all_of(sigma_cols), names_to = "Parameter", values_to = "Estimate") %>%
  mutate(True_Val = true_sigma_vals[Parameter])

sigma_metrics <- sigma_long %>%
  group_by(Scenario_F, Model, Parameter) %>%
  summarise(
    Bias     = (mean(Estimate) - mean(True_Val)),
    Var_Est  = var(Estimate, na.rm=TRUE),
    .groups  = "drop"
  )

# --- CHANGE 4: RE COMPUTATION SIGMA (Pattern) ---
sigma_re <- sigma_metrics %>%
  dplyr::select(Scenario_F, Parameter, Model, Var_Est) %>%
  pivot_wider(names_from = Model, values_from = Var_Est) %>%
  mutate(
    RE_CC      = log(CC / Full),
    RE_MICE    = log(MICE / Full),
    RE_Pattern = log(Pattern / Full), 
    RE_Full    = 0
  ) %>%
  pivot_longer(starts_with("RE_"), names_to="Model", values_to="RE") %>%
  mutate(Model = sub("RE_", "", Model)) %>%
  dplyr::select(Scenario_F, Parameter, Model, RE)

sigma_final <- left_join(sigma_metrics, sigma_re, by = c("Scenario_F", "Parameter", "Model"))

# ----------------------------------------------------------------
# 4. DATA PREP FOR PROFESSIONAL PLOTTING
# ----------------------------------------------------------------

format_greek <- function(x) {
  x <- gsub("Beta([0-9]+)", "beta[\\1]", x)
  x <- gsub("sigma_([a-z]+)", "sigma[\\1]", x)
  x <- gsub("rho([0-9]+)", "rho[\\1]", x)
  return(x)
}

# Apply to Beta Data
beta_final <- beta_final %>%
  mutate(
    Facet_Label = format_greek(Parameter),
    Facet_Label = factor(Facet_Label, levels = unique(format_greek(names(beta_true_vec))))
  )

# Apply to Sigma Data
sigma_final <- sigma_final %>%
  mutate(
    Facet_Label = format_greek(Parameter),
    Facet_Label = factor(Facet_Label, levels = unique(format_greek(sigma_cols)))
  )

clean_scenario <- function(x) gsub("Miss_Lev_", "", x)

# ----------------------------------------------------------------
# 5. PLOTS
# ----------------------------------------------------------------
if(!dir.exists("plots")) dir.create("plots")

# Professional Theme
theme_professional <- theme_bw(base_size = 14) + 
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.background = element_rect(fill = "grey95", color = "grey50"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.title = element_text(face="bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


my_colors <- c("Full"="black", "Pattern"="#1b9e77", "MICE"="#377eb8", "CC"="#e41a1c")

# --- 1. BETA BIAS ---
p1 <- ggplot(beta_final, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") +
  geom_line(linewidth=0.8) + 
  geom_point(size=2) +
  facet_wrap(~Facet_Label, labeller = label_parsed) + 
  scale_color_manual(values=my_colors) +
  scale_x_discrete(labels = clean_scenario) +
  labs(
    title = expression(bold("Bias of Regression Coefficients")),
    x = "Missingness Level",
    y = "Bias"
  ) +
  theme_professional
ggsave("plots/1_Beta_Bias.pdf", p1, width=12, height=8)

# --- 2. BETA COVERAGE ---
p2 <- ggplot(beta_final, aes(x=Scenario_F, y=Coverage, color=Model, group=Model)) +
  geom_hline(yintercept=0.95, linetype="dashed", color="grey40") +
  geom_line(linewidth=0.8) + 
  geom_point(size=2) +
  facet_wrap(~Facet_Label, labeller = label_parsed) + 
  coord_cartesian(ylim=c(0, 1.0)) + 
  scale_color_manual(values=my_colors) +
  scale_x_discrete(labels = clean_scenario) +
  labs(
    title = expression(bold("Coverage Probability (95% CI)")),
    x = "Missingness Level",
    y = "Coverage Rate"
  ) +
  theme_professional
ggsave("plots/2_Beta_Coverage.pdf", p2, width=12, height=8)

# --- 3. BETA RE (Log Variance Ratio) ---
p3 <- ggplot(beta_final, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") +
  geom_line(linewidth=0.8) + 
  geom_point(size=2) +
  facet_wrap(~Facet_Label, labeller = label_parsed) + 
  scale_color_manual(values=my_colors) +
  scale_x_discrete(labels = clean_scenario) +
  labs(
    title = expression(bold("Log Relative Variance")),
    x = "Missingness Level",
    y = "Log Ratio",
    caption = "Value > 0 indicates efficiency loss compared to Full Data"
  ) +
  theme_professional
ggsave("plots/3_Beta_LogRE.pdf", p3, width=12, height=8)

# --- 4. SIGMA METRICS ---
sigma_subset <- sigma_final %>% 
  filter(grepl("sigma", Parameter)) 

p4a <- ggplot(sigma_subset, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") + 
  geom_line(linewidth=0.7) + geom_point(size=1.5) +
  facet_wrap(~Facet_Label, labeller = label_parsed) + 
  scale_color_manual(values=my_colors) +
  scale_x_discrete(labels = clean_scenario) +
  theme_professional + 
  labs(title=expression(bold("Variance Comp. Bias") ~ (hat(sigma) - sigma)), x=NULL, y="Bias")

p4b <- ggplot(sigma_subset, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") + 
  geom_line(linewidth=0.7) + geom_point(size=1.5) +
  facet_wrap(~Facet_Label, labeller = label_parsed) + 
  scale_color_manual(values=my_colors) +
  scale_x_discrete(labels = clean_scenario) +
  theme_professional + 
  labs(title="Log Relative Variance", x="Missingness Level", y="Log Ratio")

ggsave("plots/4_Sigma_Metrics.pdf", grid.arrange(p4a, p4b, nrow=2), width=10, height=10)

# --- 5. RHO METRICS ---
rho_subset <- sigma_final %>% 
  filter(grepl("rho", Parameter))

p5a <- ggplot(rho_subset, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") + 
  geom_line(linewidth=0.7) + geom_point(size=1.5) +
  facet_wrap(~Facet_Label, labeller = label_parsed) + 
  scale_color_manual(values=my_colors) +
  scale_x_discrete(labels = clean_scenario) +
  theme_professional + 
  labs(title=expression(bold("Correlation Bias") ~ (hat(rho) - rho)), x=NULL, y="Bias")

p5b <- ggplot(rho_subset, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") + 
  geom_line(linewidth=0.7) + geom_point(size=1.5) +
  facet_wrap(~Facet_Label, labeller = label_parsed) + 
  scale_color_manual(values=my_colors) +
  scale_x_discrete(labels = clean_scenario) +
  theme_professional + 
  labs(title="Log Relative Variance", x="Missingness Level", y="Log Ratio")

ggsave("plots/5_Rho_Metrics.pdf", grid.arrange(p5a, p5b, nrow=2), width=10, height=10)

cat("Analysis complete. Professional plots (Pattern) saved in 'plots/'.\n")