## ------------------------------------------------
## MV-PEAL Multivariate Simulation Summary Script
## ------------------------------------------------

## 0. Setup ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(readr)
})

## ------------------------------------------------
## 1. Load and combine simulation outputs
## ------------------------------------------------

## Adjust pattern/path if needed
files_allout <- list.files(
  pattern = "^AllOut_rho_.*\\.rds$",
  full.names = TRUE
)

print(files_allout)

# Helper to extract the rho label from file name
get_rho_label <- function(path) {
  # e.g. "AllOut_rho_0.05_0.10_0.15.rds" --> "rho_0.05_0.10_0.15"
  basename(path) |>
    str_remove("^AllOut_") |>
    str_remove("\\.rds$")
}

# Read all objects and attach label
sim_list   <- purrr::map(files_allout, readRDS)
rho_labels <- purrr::map_chr(files_allout, get_rho_label)
names(sim_list) <- rho_labels

# Combine components across rho settings
beta_long_all <- purrr::map2_dfr(
  sim_list, names(sim_list),
  ~ dplyr::mutate(.x$beta_long, rho_label = .y)
)

sigmaRE_long_all <- purrr::map2_dfr(
  sim_list, names(sim_list),
  ~ dplyr::mutate(.x$sigmaRE_long, rho_label = .y)
)

sigmae_long_all <- purrr::map2_dfr(
  sim_list, names(sim_list),
  ~ dplyr::mutate(.x$sigmae_long, rho_label = .y)
)

corr_long_all <- purrr::map2_dfr(
  sim_list, names(sim_list),
  ~ dplyr::mutate(.x$corr_long, rho_label = .y)
)

RE_beta_all <- purrr::map2_dfr(
  sim_list, names(sim_list),
  ~ dplyr::mutate(.x$RE_beta, rho_label = .y)
)

glimpse(beta_long_all)

## ------------------------------------------------
## 2. Helper parsing for parameters & rho grid
## ------------------------------------------------

# From your simulation: px = 4, py = 3
px <- 4
py <- 3

# Recover outcome index and covariate index for nicer faceting.
beta_long_all <- beta_long_all %>%
  mutate(
    BetaIndex = as.integer(str_remove(Parameter, "Beta")),
    Outcome   = paste0("Outcome ", ceiling(BetaIndex / px)),
    Covariate = paste0("X", ((BetaIndex - 1) %% px) + 1)
  )

sigmaRE_long_all <- sigmaRE_long_all %>%
  mutate(
    RE_Component = case_when(
      Parameter == "sigma_u" ~ "Site RE (u)",
      str_detect(Parameter, "^sigma_v_") ~ "Patient RE (v_h)",
      TRUE ~ Parameter
    )
  )

# For corr parameters, extract i,j indices and keep off-diagonals
corr_long_all <- corr_long_all %>%
  mutate(
    i = as.integer(substr(Parameter, 5, 5)),
    j = as.integer(substr(Parameter, 7, 7))
  ) %>%
  filter(i != j) %>%
  mutate(
    Pair = paste0("rho_", i, j)
  )

# Build rho grid from labels: "rho_0.05_0.10_0.15"
rho_grid_df <- tibble(rho_label = rho_labels) %>%
  distinct() %>%
  mutate(parts = str_split(rho_label, "_")) %>%
  mutate(
    rho12   = purrr::map_dbl(parts, ~ as.numeric(.x[2])),
    rho13   = purrr::map_dbl(parts, ~ as.numeric(.x[3])),
    rho23   = purrr::map_dbl(parts, ~ as.numeric(.x[4])),
    rho_avg = (rho12 + rho13 + rho23) / 3
  ) %>%
  select(-parts)

print(rho_grid_df)

## ------------------------------------------------
## 3. Plot: correlation scenarios used
## ------------------------------------------------

p_rho_grid <- ggplot(rho_grid_df, aes(x = rho12, y = rho23, label = rho_label)) +
  geom_point() +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw() +
  labs(
    x = expression(rho[12]),
    y = expression(rho[23]),
    title = "True correlation scenarios used in the simulation"
  )

print(p_rho_grid)
# ggsave("rho_grid_scenarios.png", p_rho_grid, width = 6, height = 4)

## ------------------------------------------------
## 4. Fixed effects: bias & variability
## ------------------------------------------------

# 4.1 Monte Carlo summaries for beta
beta_summ <- beta_long_all %>%
  group_by(rho_label, Model, Parameter, Outcome, Covariate) %>%
  summarise(
    True      = first(True),
    mean_est  = mean(Estimate, na.rm = TRUE),
    bias      = mean_est - True,
    var_est   = var(Estimate, na.rm = TRUE),
    mse       = var_est + bias^2,
    mc_se     = sd(Estimate, na.rm = TRUE) / sqrt(sum(!is.na(Estimate))),
    .groups   = "drop"
  ) %>%
  left_join(rho_grid_df, by = "rho_label")

head(beta_summ)

# 4.1.1 Bias across correlation scenarios
p_beta_bias_rho <- beta_summ %>%
  mutate(
    rho_label = fct_reorder(rho_label, rho_avg)
  ) %>%
  ggplot(aes(x = rho_label, y = bias, color = Model, group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  facet_wrap(~ Outcome + Covariate, 
             # scales = "free_y"
             ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Correlation scenario",
    y = "Monte Carlo bias",
    title = "Fixed-effect bias across correlation scenarios and working correlations"
  )

print(p_beta_bias_rho)
# ggsave("beta_bias_across_rho.png", p_beta_bias_rho, width = 10, height = 8)

# 4.2 Distribution of beta estimates 
p_beta_bias_rho_band <- beta_summ %>%
  mutate(
    rho_label = fct_reorder(rho_label, rho_avg)
  ) %>%
  ggplot(aes(x = rho_label, y = bias, color = Model, group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = bias - 2 * mc_se,
                      ymax = bias + 2 * mc_se),
                  position = position_dodge(width = 0.3),
                  linewidth = 0.4, size = 0.1) +
  facet_wrap(~ Outcome + Covariate) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Correlation scenario",
    y = "Monte Carlo bias (with 95% MC CI)",
    title = "Fixed-effect bias across correlation scenarios and working correlations"
  )

print(p_beta_bias_rho_band)

# ggsave(paste0("p_beta_bias_rho_band.png"), p_beta_bias_rho_band, width = 10, height = 6)

## ------------------------------------------------
## 5. Variance components
## ------------------------------------------------

# 5.1 Random-effect SD bias
sigmaRE_summ <- sigmaRE_long_all %>%
  group_by(rho_label, Model, Parameter, RE_Component) %>%
  summarise(
    True     = first(True),
    mean_est = mean(Estimate, na.rm = TRUE),
    bias     = mean_est - True,
    var_est  = var(Estimate, na.rm = TRUE),
    mse      = var_est + bias^2,
    .groups  = "drop"
  ) %>%
  left_join(rho_grid_df, by = "rho_label")

p_sigmaRE_bias <- sigmaRE_summ %>%
  mutate(rho_label = fct_reorder(rho_label, rho_avg)) %>%
  ggplot(aes(x = rho_label, y = bias, color = Model, group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  facet_wrap(~ RE_Component, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Correlation scenario",
    y = "Monte Carlo bias",
    title = "Bias of random-effect SDs across correlation scenarios"
  )

print(p_sigmaRE_bias)
# ggsave("sigmaRE_bias_across_rho.png", p_sigmaRE_bias, width = 10, height = 6)

# 5.2 Residual SD bias
sigmae_summ <- sigmae_long_all %>%
  group_by(rho_label, Model, Parameter) %>%
  summarise(
    True     = first(True),
    mean_est = mean(Estimate, na.rm = TRUE),
    bias     = mean_est - True,
    var_est  = var(Estimate, na.rm = TRUE),
    mse      = var_est + bias^2,
    .groups  = "drop"
  ) %>%
  left_join(rho_grid_df, by = "rho_label")

p_sigmae_bias <- sigmae_summ %>%
  mutate(rho_label = fct_reorder(rho_label, rho_avg)) %>%
  ggplot(aes(x = rho_label, y = bias, color = Model, group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Correlation scenario",
    y = "Monte Carlo bias",
    title = "Bias of residual SD across correlation scenarios"
  )

print(p_sigmae_bias)
# ggsave("sigmae_bias_across_rho.png", p_sigmae_bias, width = 8, height = 5)

## ------------------------------------------------
## 6. Correlation parameters
## ------------------------------------------------

# 6.1 Monte Carlo summaries for off-diagonal correlations
corr_summ <- corr_long_all %>%
  group_by(rho_label, Model, Pair, Parameter, i, j) %>%
  summarise(
    True     = first(True),
    mean_est = mean(Estimate, na.rm = TRUE),
    bias     = mean_est - True,
    var_est  = var(Estimate, na.rm = TRUE),
    mse      = var_est + bias^2,
    .groups  = "drop"
  ) %>%
  left_join(rho_grid_df, by = "rho_label")

p_corr_bias <- corr_summ %>%
  mutate(rho_label = fct_reorder(rho_label, rho_avg)) %>%
  ggplot(aes(x = rho_label, y = bias, color = Model, group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  facet_wrap(~ Pair, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Correlation scenario",
    y = "Monte Carlo bias",
    title = "Bias of correlation parameters (off-diagonals)"
  )

print(p_corr_bias)
# ggsave("corr_bias_across_rho.png", p_corr_bias, width = 10, height = 6)

# 6.2 True vs mean estimated correlation
p_corr_true_vs_est <- ggplot(corr_summ, aes(x = True, y = mean_est, color = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.8) +
  facet_wrap(~ rho_label) +
  theme_bw() +
  labs(
    x = "True correlation",
    y = "Mean estimated correlation",
    title = "True vs mean estimated correlation by scenario and working correlation"
  )

print(p_corr_true_vs_est)
# ggsave("corr_true_vs_est.png", p_corr_true_vs_est, width = 10, height = 6)

## ------------------------------------------------
## 7. Relative efficiency of working correlation assumptions
## ------------------------------------------------

RE_beta_all2 <- RE_beta_all %>%
  left_join(rho_grid_df, by = "rho_label") %>%
  mutate(
    Parameter_index = as.integer(str_remove(Parameter, "Beta")),
    Outcome   = paste0("Outcome ", ceiling(Parameter_index / px)),
    Covariate = paste0("X", ((Parameter_index - 1) %% px) + 1)
  )

glimpse(RE_beta_all2)

# 7.1 Exchangeable vs Unstructured
p_RE_ex <- RE_beta_all2 %>%
  mutate(rho_label = fct_reorder(rho_label, rho_avg)) %>%
  ggplot(aes(x = rho_label, y = RE_EX_vs_UN, group = Covariate, color = Covariate)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() +
  geom_point() +
  facet_wrap(~ Outcome) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Correlation scenario",
    y = "Var(EX) / Var(UN)",
    title = "Relative efficiency: Exchangeable vs Unstructured (by outcome and covariate)"
  )

print(p_RE_ex)
# ggsave("RE_EX_vs_UN.png", p_RE_ex, width = 10, height = 6)

# 7.2 Independence vs Unstructured
p_RE_in <- RE_beta_all2 %>%
  mutate(rho_label = fct_reorder(rho_label, rho_avg)) %>%
  ggplot(aes(x = rho_label, y = RE_IN_vs_UN, group = Covariate, color = Covariate)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() +
  geom_point() +
  facet_wrap(~ Outcome) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Correlation scenario",
    y = "Var(IN) / Var(UN)",
    title = "Relative efficiency: Independence vs Unstructured (by outcome and covariate)"
  )

print(p_RE_in)


gridExtra::grid.arrange(p_RE_ex, p_RE_in, ncol = 1)

## ------------------------------------------------
## 8. Optional: session info
## ------------------------------------------------

print(sessionInfo())

