# ==============================================================================
# mice_diagnostic_viz_sitelocal.R
# Distributional performance of SiteLocal MICE across scenarios
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

summary_tbl <- readRDS("results/mice_diag_summary.rds")

# ----------------------------------------------------------------
# PREP
# ----------------------------------------------------------------
scen_order <- paste0("Miss_Lev_", 1:8)
TRUE_CORR  <- c(rho12 = -0.3, rho13 = 0.3, rho23 = 0.8)

sl <- summary_tbl %>%
  filter(Strategy == "SiteLocal") %>%
  mutate(Scenario = factor(Scenario, levels = scen_order),
         ScenNum  = as.integer(gsub("Miss_Lev_", "", as.character(Scenario))))

# -- Missingness rates
miss_long <- sl %>%
  dplyr::select(ScenNum, matches("^miss_rate_Y[123]__mean$")) %>%
  pivot_longer(-ScenNum, names_to = "Outcome", values_to = "miss_rate",
               names_pattern = "miss_rate_(Y[123])__mean")

# -- Mean bias
bias_long <- sl %>%
  dplyr::select(ScenNum, matches("^mean_bias_Y[123]__(mean|sd)$")) %>%
  pivot_longer(-ScenNum,
               names_to = c("Outcome", ".value"),
               names_pattern = "mean_bias_(Y[123])__(mean|sd)")

# -- Variance ratio
var_long <- sl %>%
  dplyr::select(ScenNum, matches("^var_ratio_Y[123]__(mean|sd)$")) %>%
  pivot_longer(-ScenNum,
               names_to = c("Outcome", ".value"),
               names_pattern = "var_ratio_(Y[123])__(mean|sd)")

# -- Correlation bias
cor_bias_long <- sl %>%
  dplyr::select(ScenNum, matches("^cor_bias_rho[0-9]+__(mean|sd)$")) %>%
  pivot_longer(-ScenNum,
               names_to = c("Pair", ".value"),
               names_pattern = "cor_bias_(rho[0-9]+)__(mean|sd)") %>%
  mutate(PairLabel = recode(Pair,
                            rho12 = "rho(Y1,Y2) = -0.30",
                            rho13 = "rho(Y1,Y3) =  0.30",
                            rho23 = "rho(Y2,Y3) =  0.80"
  ))

# -- Absolute correlation recovered
cor_imp_long <- sl %>%
  dplyr::select(ScenNum, matches("^cor_imp_rho[0-9]+__mean$")) %>%
  pivot_longer(-ScenNum, names_to = "Pair", values_to = "cor_imp",
               names_pattern = "cor_imp_(rho[0-9]+)__mean") %>%
  mutate(TrueCorr  = TRUE_CORR[Pair],
         PairLabel = recode(Pair,
                            rho12 = "rho(Y1,Y2) = -0.30",
                            rho13 = "rho(Y1,Y3) =  0.30",
                            rho23 = "rho(Y2,Y3) =  0.80"
         ))

# ----------------------------------------------------------------
# THEME & PALETTES
# ----------------------------------------------------------------
theme_diag <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f0f4f8", color = "grey70"),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 10, color = "grey45"),
    legend.position  = "bottom"
  )

outcome_colors <- c(Y1 = "#2166ac", Y2 = "#4dac26", Y3 = "#d01c8b")
outcome_shapes <- c(Y1 = 16, Y2 = 17, Y3 = 15)
pair_colors    <- c(
  "rho(Y1,Y2) = -0.30" = "#e08214",
  "rho(Y1,Y3) =  0.30" = "#542788",
  "rho(Y2,Y3) =  0.80" = "#1b7837"
)

x_scale <- scale_x_continuous(breaks = 1:8, labels = paste0("L", 1:8))

# ================================================================
# PLOT 1 — Missingness rates
# ================================================================
p_miss <- ggplot(miss_long, aes(x = ScenNum, y = miss_rate,
                                color = Outcome, shape = Outcome, group = Outcome)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  x_scale +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
  scale_color_manual(values = outcome_colors) +
  scale_shape_manual(values = outcome_shapes) +
  labs(title    = "Missingness Rate per Outcome",
       subtitle = "Sequential logistic mechanism; L1 = mildest, L8 = most severe",
       x = NULL, y = "Mean % Missing") +
  theme_diag

# ================================================================
# PLOT 2 — Mean bias
# ================================================================
p_bias <- ggplot(bias_long, aes(x = ScenNum, y = mean,
                                color = Outcome, group = Outcome)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.6) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Outcome),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Outcome), size = 2.8) +
  x_scale +
  scale_color_manual(values = outcome_colors) +
  scale_fill_manual(values  = outcome_colors) +
  scale_shape_manual(values = outcome_shapes) +
  labs(title    = "Mean Bias  [E(Y_imp) - E(Y_true)]",
       subtitle = "Ribbon = +/-1 SD across reps; dashed = zero bias",
       x = NULL, y = "Mean Bias") +
  theme_diag

# ================================================================
# PLOT 3 — Variance ratio
# ================================================================
p_var <- ggplot(var_long, aes(x = ScenNum, y = mean,
                              color = Outcome, group = Outcome)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.6) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Outcome),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Outcome), size = 2.8) +
  x_scale +
  scale_y_continuous(trans  = "log10",
                     breaks = c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5),
                     labels = c("0.70","0.80","0.90","1.00","1.10","1.20","1.50")) +
  scale_color_manual(values = outcome_colors) +
  scale_fill_manual(values  = outcome_colors) +
  scale_shape_manual(values = outcome_shapes) +
  labs(title    = "Variance Ratio  [Var(Y_imp) / Var(Y_true)]",
       subtitle = "Log scale; dashed = perfect recovery at 1.0",
       x = NULL, y = "Variance Ratio (log)") +
  theme_diag

# ================================================================
# PLOT 4 — Correlation bias
# ================================================================
p_cor_bias <- ggplot(cor_bias_long, aes(x = ScenNum, y = mean,
                                        color = PairLabel, group = PairLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.6) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = PairLabel),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  x_scale +
  scale_color_manual(values = pair_colors) +
  scale_fill_manual(values  = pair_colors) +
  labs(title    = "Correlation Bias  [rho_imp - rho_true]",
       subtitle = "Dashed = zero bias; ribbon = +/-1 SD across reps",
       x = "Scenario", y = "Correlation Bias",
       color = "Pair", fill = "Pair") +
  theme_diag

# ================================================================
# PLOT 5 — Absolute correlation recovered
# ================================================================
p_cor_abs <- ggplot(cor_imp_long, aes(x = ScenNum, y = cor_imp,
                                      color = PairLabel, group = PairLabel)) +
  geom_hline(aes(yintercept = TrueCorr, color = PairLabel),
             linetype = "dotted", linewidth = 1.0) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  x_scale +
  scale_color_manual(values = pair_colors) +
  labs(title    = "Recovered Correlation",
       subtitle = "Dotted line = true rho per pair; solid = imputed estimate",
       x = "Scenario", y = "Imputed rho",
       color = "Pair") +
  theme_diag

# ================================================================
# ASSEMBLE
# ================================================================
full_fig <- p_miss / p_bias / p_var / (p_cor_bias + p_cor_abs) +
  plot_annotation(
    title    = "MICE (Site-Local) — Distributional Recovery of 3-Level Correlated Outcomes",
    subtitle = "8 missingness scenarios | m = 5 imputations | PMM within site",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey45")
    )
  )

ggsave("results/mice_sitelocal_perf.pdf",
       full_fig, width = 12, height = 18, device = cairo_pdf)

ggsave("results/mice_sitelocal_perf.png",
       full_fig, width = 12, height = 18, dpi = 150)

cat("Saved: results/mice_sitelocal_perf.pdf / .png\n")