rm(list=ls())
# ==============================================================================
# File: analyze_sim_specs_rho_synthesis.R
# Usage: Rscript analyze_sim_specs_rho_synthesis.R
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(gridExtra)
  library(stringr)
  library(grid)
})

# ----------------------------------------------------------------
# 1. SETUP & TRUTH DEFINITIONS
# ----------------------------------------------------------------
H <- 10; px <- 2; py <- 3

beta_true_mat <- matrix(seq(-3, 3, length.out = px * py), nrow = px, ncol = py)
beta_true_vec <- as.vector(beta_true_mat)
names(beta_true_vec) <- paste0("Beta", 1:(px*py))

# ----------------------------------------------------------------
# 2. READ DATA
# ----------------------------------------------------------------
out_dir <- "results"
files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)

if(length(files) == 0) {
  stop("No result files found in 'results/'.")
}

cat(sprintf("Reading %d files...\n", length(files)))
dat <- do.call(rbind, lapply(files, readRDS))
dat <- dat %>% distinct(SimID, Scenario, Model, .keep_all = TRUE)

# --- SCENARIO PARSING (CORRELATION LEVELS) ---
# Extract numerical values from string like "rho_0.05_0.10_0.15"
# We handle potential case variations (Rho/rho)
rho_vals <- str_match(dat$Scenario, "(?i)rho_([0-9.-]+)_([0-9.-]+)_([0-9.-]+)")

if(all(is.na(rho_vals))) {
  # Fallback if regex fails: treat scenarios as discrete factors sorted alphabetically
  dat$Scenario_F <- factor(dat$Scenario, levels=sort(unique(dat$Scenario)))
  levels(dat$Scenario_F) <- paste0("Level ", 1:length(levels(dat$Scenario_F)))
} else {
  dat$True_Rho12 <- as.numeric(rho_vals[,2])
  dat$True_Rho13 <- as.numeric(rho_vals[,3])
  dat$True_Rho23 <- as.numeric(rho_vals[,4])
  
  # Compute strength proxy (Sum of absolute correlations)
  dat <- dat %>% mutate(rho_sum = abs(True_Rho12) + abs(True_Rho13) + abs(True_Rho23))
  
  # Identify unique scenarios and sort by strength
  unique_scens <- unique(dat$Scenario[order(dat$rho_sum)])
  
  dat$Scenario_F <- factor(dat$Scenario, levels = unique_scens)
  
  # Rename levels to "Level 1" ... "Level N" for clean plotting
  levels(dat$Scenario_F) <- seq_along(levels(dat$Scenario_F))
}

# --- CAPITALIZATION & MODEL ORDERING ---
# Convert "unstructured" -> "Unstructured", etc.
dat$Model <- str_to_title(dat$Model)
dat$Model <- factor(dat$Model, levels = c("Unstructured", "Exchangeable", "Independence"))

# ----------------------------------------------------------------
# 3. METRICS COMPUTATION (BETA ONLY)
# ----------------------------------------------------------------
beta_cols <- grep("^Beta[0-9]+$", names(dat), value = TRUE)
se_cols   <- grep("^SE[0-9]+$",   names(dat), value = TRUE)

long_est <- dat %>% dplyr::select(SimID, Scenario_F, Model, all_of(beta_cols)) %>%
  pivot_longer(cols = all_of(beta_cols), names_to = "Parameter", values_to = "Estimate")

long_se <- dat %>% dplyr::select(SimID, Scenario_F, Model, all_of(se_cols)) %>%
  pivot_longer(cols = all_of(se_cols), names_to = "Parameter_SE", values_to = "SE") %>%
  mutate(Parameter = sub("SE", "Beta", Parameter_SE)) %>% dplyr::select(-Parameter_SE)

beta_long <- left_join(long_est, long_se, by = c("SimID", "Scenario_F", "Model", "Parameter"))

beta_metrics <- beta_long %>%
  group_by(Scenario_F, Model, Parameter) %>%
  summarise(
    Bias    = mean(Estimate) - beta_true_vec[unique(Parameter)],
    Var_Est = var(Estimate, na.rm=TRUE),
    .groups = "drop"
  )

# RE Calculation: Log(Var_Model / Var_Unstructured)
beta_re <- beta_metrics %>%
  dplyr::select(Scenario_F, Parameter, Model, Var_Est) %>%
  pivot_wider(names_from = Model, values_from = Var_Est)

if("Unstructured" %in% names(beta_re)) {
  beta_re <- beta_re %>%
    mutate(
      RE_Exchangeable = log(Exchangeable / Unstructured),
      RE_Independence = log(Independence / Unstructured),
      RE_Unstructured = 0 # Log(1)
    ) %>%
    dplyr::select(Scenario_F, Parameter, starts_with("RE_")) %>%
    pivot_longer(starts_with("RE_"), names_to="Model", values_to="RE") %>%
    mutate(Model = sub("RE_", "", Model))
} else {
  beta_re <- beta_metrics %>% select(Scenario_F, Parameter, Model) %>% mutate(RE = NA)
}

beta_final <- left_join(beta_metrics, beta_re, by = c("Scenario_F", "Parameter", "Model"))

# ----------------------------------------------------------------
# 4. PLOT PREP
# ----------------------------------------------------------------
format_greek <- function(x) {
  gsub("Beta([0-9]+)", "beta[\\1]", x)
}

get_outcome_label <- function(param_name) {
  num <- as.numeric(gsub("Beta", "", param_name))
  if (num %in% 1:2) return("Outcome~1")
  if (num %in% 3:4) return("Outcome~2")
  if (num %in% 5:6) return("Outcome~3")
  return("")
}

beta_final <- beta_final %>%
  mutate(
    Greek_Label = format_greek(Parameter),
    Outcome_Label = sapply(Parameter, get_outcome_label),
    Final_Label_Str = paste(Outcome_Label, Greek_Label, sep="~':'~"),
    Facet_Label = factor(Final_Label_Str, levels = unique(Final_Label_Str[order(as.numeric(gsub("Beta", "", Parameter)))]))
  )

# Color Scheme
my_colors <- c("Unstructured"="#1b9e77",  # Green
               "Exchangeable"="#377eb8",  # Blue
               "Independence"="#e41a1c")  # Red

# Theme
theme_clean <- theme_bw(base_size = 11) + 
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill="grey95"),
    strip.text = element_text(face="bold", size=14),
    legend.title = element_text(face="bold", size=12),
    legend.text = element_text(size=12),
    plot.title = element_text(face="bold", size=12, hjust=0.5),
    panel.grid.minor = element_blank()
  )

theme_bottom <- theme_clean + theme(
  axis.text.x = element_text(angle=0, hjust=0.5), axis.ticks.x = element_line(),
  legend.position = "bottom"
)

# ----------------------------------------------------------------
# 5. GENERATE PLOTS
# ----------------------------------------------------------------

# 1. BIAS
p_bias <- ggplot(beta_final, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line() + geom_point(size=2) + 
  facet_wrap(~Facet_Label, nrow=1, labeller=label_parsed) +
  scale_color_manual(values=my_colors) + 
  ylim(-0.05, 0.05) +
  theme_clean + 
  labs(title="Fixed-Effects Bias", x=NULL, y="Bias")

# 2. RELATIVE EFFICIENCY (Log Variance Ratio)
p_re <- ggplot(beta_final, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line() + geom_point(size=2) + 
  facet_wrap(~Facet_Label, nrow=1, labeller=label_parsed) +
  scale_color_manual(values=my_colors) + 
  theme_bottom + 
  labs(title="Log Relative Variance (vs Unstructured)", 
       subtitle="Value > 0 indicates higher variance (worse efficiency) than Unstructured",
       color = "Correlation Structure",
       x="Correlation Strength Level (Low -> High)", y="Log Variance Ratio")

# ----------------------------------------------------------------
# 6. SYNTHESIS
# ----------------------------------------------------------------
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
mylegend <- g_legend(p_re)
p_re <- p_re + theme(legend.position="none")

lay_mat <- rbind(c(1,1),
                 c(2,2),
                 c(3,3))

if(!dir.exists("plots")) dir.create("plots")
pdf("plots/Specs-Rho-Comparison-Synthesis.pdf", width=14, height=10)
grid.arrange(p_bias, p_re, mylegend, 
             layout_matrix = lay_mat,
             heights = c(1, 1.2, 0.2))
dev.off()

cat("Synthesis plot saved to 'plots/Specs-Rho-Comparison-Synthesis.pdf'\n")