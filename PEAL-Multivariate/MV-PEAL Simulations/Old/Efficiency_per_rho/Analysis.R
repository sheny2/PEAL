library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

rho_grid <- list(
  c(-0.05, -0.10, 0.15), 
  c(-0.15, -0.20, 0.25), 
  c(-0.25, -0.30, 0.35), 
  c(-0.35, -0.40, 0.45), 
  c(-0.45, -0.50, 0.55), 
  c(-0.55, -0.60, 0.65),
  c(-0.65, -0.70, 0.75), 
  c(-0.75, -0.80, 0.85)
)

rho_grid <- list(
c(0.05, 0.10, 0.15),
c(0.05, 0.15, 0.25),
c(0.1, 0.20, 0.35),
c(0.1, 0.25, 0.45),
c(0.15, 0.30, 0.55),
c(0.15, 0.35, 0.65),
c(0.2, 0.45, 0.75),
c(0.2, 0.50, 0.85)
)

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

# Use the same labels as in your simulation run
rho_labels <- sapply(rho_grid, function(r) {
  paste0("rho_", paste(sprintf("%0.2f", r), collapse = "_"))
})

rho_meta <- data.frame(
  setting = rho_labels,
  rho12   = sapply(rho_grid, `[`, 1),
  rho13   = sapply(rho_grid, `[`, 2),
  rho23   = sapply(rho_grid, `[`, 3),
  scenario_id = seq_along(rho_labels) 
)
rho_meta

RE_all <- map_dfr(rho_meta$setting, function(lbl) {
  tab <- readRDS(paste0("RE_table_beta_", lbl, ".rds"))
  tab$setting <- lbl
  tab
})

# Attach rho12, rho13, rho23, and scenario_id
RE_all <- RE_all %>%
  left_join(rho_meta, by = "setting")


RE_long <- RE_all %>%
  select(Parameter, setting, scenario_id, rho12, rho13, rho23,
         RE_EX_vs_UN, RE_IN_vs_UN) %>%
  pivot_longer(
    cols      = c(RE_EX_vs_UN, RE_IN_vs_UN),
    names_to  = "Comparison",
    values_to = "RE"
  ) %>%
  mutate(
    Comparison = dplyr::recode(Comparison,
                               RE_EX_vs_UN = "EX vs UN",
                               RE_IN_vs_UN = "IN vs UN"),
    # build labels from rho_meta (all numeric), ordered from weak to strong
    scenario_label = factor(
      setting,
      levels = rho_meta$setting,
      labels = with(rho_meta,
                    paste0("(",
                           sprintf("%.2f", rho12), ", ",
                           sprintf("%.2f", rho13), ", ",
                           sprintf("%.2f", rho23), ")"))
    )
  )


ggplot(RE_long,
       aes(x = scenario_label, y = RE, group = Comparison, color = Comparison)) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.7, alpha = 0.5) +
    geom_line() +
  geom_point() +
  facet_wrap(~ Parameter, ncol = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    x = "Correlation scenario (rho12, rho13, rho23)",
    y = "Relative efficiency vs UN (Variance ratio)",
    title = "Relative Efficiency of EX and IN vs UN across correlation settings"
  )



## ---- 1. Collect mean bias for betas over sims, per setting & model ----
Bias_beta_all <- map_dfr(rho_meta$setting, function(lbl) {
  out <- readRDS(paste0("AllOut_", lbl, ".rds"))  # contains beta_long
  
  out$beta_long %>%
    mutate(Bias = Estimate - True) %>%
    group_by(Model, Parameter) %>%   # average over simulations
    summarise(
      Bias = mean(Bias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(setting = lbl)
})

# Attach rho info & scenario labels
Bias_beta_all <- Bias_beta_all %>%
  left_join(rho_meta, by = "setting") %>%
  mutate(
    scenario_label = factor(
      setting,
      levels = rho_meta$setting,
      labels = with(
        rho_meta,
        paste0("(",
               sprintf("%.2f", rho12), ", ",
               sprintf("%.2f", rho13), ", ",
               sprintf("%.2f", rho23), ")")
      )
    )
  )

## ---- 2. Plot: mean bias of betas vs correlation scenario ----
ggplot(Bias_beta_all,
       aes(x = scenario_label, y = Bias,
           group = Model, color = Model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Parameter, ncol = 4, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    x = "Correlation scenario (rho12, rho13, rho23)",
    y = "Mean Bias of beta",
    title = "Bias of beta estimates across correlation settings"
  )


## ---- 3. Collect mean bias for random-effect sigmas ----
Bias_sigmaRE_all <- map_dfr(rho_meta$setting, function(lbl) {
  out <- readRDS(paste0("AllOut_", lbl, ".rds"))  # contains sigmaRE_long
  
  out$sigmaRE_long %>%
    mutate(Bias = Estimate - True) %>%
    group_by(Model, Parameter) %>%
    summarise(
      Bias = mean(Bias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(setting = lbl)
})

Bias_sigmaRE_all <- Bias_sigmaRE_all %>%
  left_join(rho_meta, by = "setting") %>%
  mutate(
    scenario_label = factor(
      setting,
      levels = rho_meta$setting,
      labels = with(
        rho_meta,
        paste0("(",
               sprintf("%.2f", rho12), ", ",
               sprintf("%.2f", rho13), ", ",
               sprintf("%.2f", rho23), ")")
      )
    )
  )

## ---- 4. Plot: mean bias of random-effect sigmas vs correlation scenario ----
ggplot(Bias_sigmaRE_all,
       aes(x = scenario_label, y = Bias,
           group = Model, color = Model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Parameter, ncol = 4, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    x = "Correlation scenario (rho12, rho13, rho23)",
    y = "Mean Bias of random-effect SDs",
    title = "Bias of random-effect standard deviations across correlation settings"
  )

## ---- 5. Collect mean bias for residual sigma_e ----
Bias_sigmae_all <- map_dfr(rho_meta$setting, function(lbl) {
  out <- readRDS(paste0("AllOut_", lbl, ".rds"))  # contains sigmae_long
  
  out$sigmae_long %>%
    mutate(Bias = Estimate - True) %>%
    group_by(Model, Parameter) %>%
    summarise(
      Bias = mean(Bias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(setting = lbl)
})

Bias_sigmae_all <- Bias_sigmae_all %>%
  left_join(rho_meta, by = "setting") %>%
  mutate(
    scenario_label = factor(
      setting,
      levels = rho_meta$setting,
      labels = with(
        rho_meta,
        paste0("(",
               sprintf("%.2f", rho12), ", ",
               sprintf("%.2f", rho13), ", ",
               sprintf("%.2f", rho23), ")")
      )
    )
  )

## ---- 6. Plot: mean bias of residual sigma vs correlation scenario ----
ggplot(Bias_sigmae_all,
       aes(x = scenario_label, y = Bias,
           group = Model, color = Model)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    x = "Correlation scenario (rho12, rho13, rho23)",
    y = "Mean Bias of residual SD (sigma_e)",
    title = "Bias of residual standard deviation across correlation settings"
  )

