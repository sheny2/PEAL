rm(list=ls())
# ==============================================================================
# File: analyze_sim_pattern_levels_synthesis.R
# Usage: Rscript analyze_sim_pattern_levels_synthesis.R
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

# --- FIX: Define these explicitly for use in case_when later ---
sigma_u_true <- 1
sigma_v_true <- 3
sigma_e_true <- 3

true_sigma_vals <- c(
  "sigma_u" = sigma_u_true, 
  "sigma_v" = sigma_v_true, 
  "sigma_e" = sigma_e_true,
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
  stop("No result files found in 'results/'.")
}

cat(sprintf("Reading %d files...\n", length(files)))
dat <- do.call(rbind, lapply(files, readRDS))
dat <- dat %>% distinct(SimID, Scenario, Model, .keep_all = TRUE)
dat$Model[dat$Model == "Observed"] <- "Pattern"

# --- Extract Correlation Levels ---
rho_vals <- str_match(dat$Scenario, "Rho_([0-9.-]+)_([0-9.-]+)_([0-9.-]+)")

if(all(is.na(rho_vals))) {
  dat$Scenario_F <- factor(dat$Scenario, levels=sort(unique(dat$Scenario)))
  levels(dat$Scenario_F) <- 1:length(levels(dat$Scenario_F))
  
  # Default True Rho values if regex fails
  dat$True_rho12 <- -0.3; dat$True_rho13 <- 0.3; dat$True_rho23 <- 0.8
} else {
  # Parse true values
  dat$True_rho12 <- as.numeric(rho_vals[,2])
  dat$True_rho13 <- as.numeric(rho_vals[,3])
  dat$True_rho23 <- as.numeric(rho_vals[,4])
  
  dat <- dat %>% mutate(rho_sum = abs(True_rho12) + abs(True_rho13) + abs(True_rho23))
  unique_scens <- unique(dat$Scenario[order(dat$rho_sum)])
  dat$Scenario_F <- factor(dat$Scenario, levels = unique_scens)
  
  # Rename levels to "Level 1" ... "Level 8"
  levels(dat$Scenario_F) <- seq_along(levels(dat$Scenario_F))
}

dat$Model <- factor(dat$Model, levels = c("Full", "Pattern", "MICE", "CC"))

# ----------------------------------------------------------------
# 3. METRICS COMPUTATION
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
    Bias     = (mean(Estimate) - beta_true_vec[unique(Parameter)]),
    Var_Est  = var(Estimate, na.rm=TRUE),
    Coverage = mean((Estimate - 1.96*SE) <= beta_true_vec[unique(Parameter)] &
                      (Estimate + 1.96*SE) >= beta_true_vec[unique(Parameter)], na.rm=TRUE),
    .groups  = "drop"
  )

beta_re <- beta_metrics %>%
  dplyr::select(Scenario_F, Parameter, Model, Var_Est) %>%
  pivot_wider(names_from = Model, values_from = Var_Est) %>%
  mutate(RE_CC = log(CC/Full), RE_MICE = log(MICE/Full), RE_Pattern = log(Pattern/Full), RE_Full = 0) %>%
  pivot_longer(starts_with("RE_"), names_to="Model", values_to="RE") %>%
  mutate(Model = sub("RE_", "", Model)) %>% dplyr::select(Scenario_F, Parameter, Model, RE)

beta_final <- left_join(beta_metrics, beta_re, by = c("Scenario_F", "Parameter", "Model"))

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
  summarise(Bias = mean(Estimate) - mean(True_Val), Var_Est = var(Estimate, na.rm=TRUE), .groups = "drop")

sigma_re <- sigma_metrics %>%
  dplyr::select(Scenario_F, Parameter, Model, Var_Est) %>%
  pivot_wider(names_from = Model, values_from = Var_Est) %>%
  mutate(RE_CC = log(CC/Full), RE_MICE = log(MICE/Full), RE_Pattern = log(Pattern/Full), RE_Full = 0) %>%
  pivot_longer(starts_with("RE_"), names_to="Model", values_to="RE") %>%
  mutate(Model = sub("RE_", "", Model)) %>% dplyr::select(Scenario_F, Parameter, Model, RE)

sigma_final <- left_join(sigma_metrics, sigma_re, by = c("Scenario_F", "Parameter", "Model"))

# ----------------------------------------------------------------
# 4. PLOT PREP
# ----------------------------------------------------------------
format_greek <- function(x) {
  x <- gsub("Beta([0-9]+)", "beta[\\1]", x)
  x <- gsub("sigma_([a-z]+)", "sigma[\\1]", x)
  x <- gsub("rho([0-9]+)", "rho[\\1]", x)
  return(x)
}

# --- CUSTOM OUTCOME LABEL LOGIC ---
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

sigma_final$Facet_Label <- factor(format_greek(sigma_final$Parameter), 
                                  levels=unique(format_greek(c("sigma_u","sigma_v","sigma_e","rho12","rho13","rho23"))))

my_colors <- c("Full"="black", "Pattern"="#1b9e77", "MICE"="#377eb8", "CC"="#e41a1c")

theme_clean <- theme_bw(base_size = 11) + 
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill="grey95"),
    strip.text = element_text(face="bold", size=14), # Large Titles
    plot.title = element_text(face="bold", size=12, hjust=0.5),
    panel.grid.minor = element_blank()
  )

theme_bottom <- theme_clean + theme(
  axis.text.x = element_text(angle=0, hjust=0.5), axis.ticks.x = element_line(),
  legend.position = "bottom"
)

# ----------------------------------------------------------------
# 5. SUB-PLOT GENERATION
# ----------------------------------------------------------------

# --- BETA ---
p_beta_bias <- ggplot(beta_final, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line() + geom_point(size=1) + facet_wrap(~Facet_Label, nrow=1, labeller=label_parsed) +
  coord_cartesian(ylim=c(-0.1, 0.02)) + 
  scale_color_manual(values=my_colors) + theme_clean + labs(title="Fixed-effects: Bias", x=NULL, y="Bias")

p_beta_cov <- ggplot(beta_final, aes(x=Scenario_F, y=Coverage, color=Model, group=Model)) +
  geom_hline(yintercept=0.95, linetype="dashed", color="grey50") +
  geom_line() + geom_point(size=1) + facet_wrap(~Facet_Label, nrow=1, labeller=label_parsed) +
  scale_color_manual(values=my_colors) + coord_cartesian(ylim=c(0.5, 1.0)) +
  theme_clean + labs(title="Fixed-effects: Coverage", x=NULL, y="Cov")

p_beta_re <- ggplot(beta_final, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line() + geom_point(size=1) + facet_wrap(~Facet_Label, nrow=1, labeller=label_parsed) +
  scale_color_manual(values=my_colors) + theme_bottom + 
  labs(title="Fixed-effects: Log Relative Variance", 
       color = "Outcome Missingness Handling Method",
       x="Correlation Strength Level (Weak -> Strong)", y="Log RE")

# --- SIGMA ---
sig_sub <- sigma_final %>% filter(grepl("sigma", Parameter))
p_sig_bias <- ggplot(sig_sub, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line() + geom_point(size=1) + facet_wrap(~Facet_Label, nrow=1, labeller=label_parsed) +
  coord_cartesian(ylim=c(-0.5, 0.5)) + 
  scale_color_manual(values=my_colors) + theme_clean + labs(title="Variance: Bias", x=NULL, y="Bias")

p_sig_re <- ggplot(sig_sub, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line() + geom_point(size=1) + facet_wrap(~Facet_Label, nrow=1, labeller=label_parsed) +
  scale_color_manual(values=my_colors) + theme_bottom +
  labs(title="Variance: Log Relative Variance", x="Correlation Strength Level", y="Log RE")

# --- RHO ---
rho_sub <- sigma_final %>% filter(grepl("rho", Parameter))
p_rho_bias <- ggplot(rho_sub, aes(x=Scenario_F, y=Bias, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  coord_cartesian(ylim=c(-0.01, 0.01)) + 
  geom_line() + geom_point(size=1) + facet_wrap(~Facet_Label, nrow=1, labeller=label_parsed) +
  scale_color_manual(values=my_colors) + theme_clean + labs(title="Correlation: Bias", x=NULL, y="Bias")

p_rho_re <- ggplot(rho_sub, aes(x=Scenario_F, y=RE, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line() + geom_point(size=1) + facet_wrap(~Facet_Label, nrow=1, labeller=label_parsed) +
  scale_color_manual(values=my_colors) + theme_bottom +
  labs(title="Correlation: Log Relative Variance", x="Correlation Strength Level", y="Log RE")

# ----------------------------------------------------------------
# 6. SYNTHESIS
# ----------------------------------------------------------------
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
mylegend <- g_legend(p_beta_re)

p_beta_re <- p_beta_re + theme(legend.position="none")
p_sig_re <- p_sig_re + theme(legend.position="none")
p_rho_re <- p_rho_re + theme(legend.position="none")

lay_mat <- rbind(c(1,1), c(2,2), c(3,3),
                 c(4,5),
                 c(6,7),
                 c(8,8))

if(!dir.exists("plots")) dir.create("plots")
pdf("plots/Method-CorrLevel.pdf", width=14, height=16)
grid.arrange(p_beta_bias, p_beta_cov, p_beta_re, 
             p_sig_bias, p_rho_bias, 
             p_sig_re, p_rho_re, 
             mylegend, 
             layout_matrix = lay_mat,
             heights = c(1, 1, 1.2, 1, 1.2, 0.2)
             # top = textGrob("MV-PEAL (Pattern) vs MICE vs CC vs Full Under Varying Correlation Levels", 
             #                gp=gpar(fontsize=18, fontface="bold"))
             )
dev.off()

cat("Synthesis plot saved to 'plots/Method-CorrLevel.pdf'\n")