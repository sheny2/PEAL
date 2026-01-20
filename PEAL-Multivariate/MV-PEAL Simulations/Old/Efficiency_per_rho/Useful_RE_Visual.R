## ------------------------------------------------
## MV-PEAL Multivariate Simulation – Relative Efficiency Plots
## ------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
})

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
## 1. Load AllOut_* objects
## ------------------------------------------------

files_allout <- list.files(
  pattern = "^AllOut_rho_.*\\.rds$",
  full.names = TRUE
)

get_rho_label <- function(path) {
  basename(path) |>
    str_remove("^AllOut_") |>
    str_remove("\\.rds$")
}

sim_list   <- purrr::map(files_allout, readRDS)
rho_labels <- purrr::map_chr(files_allout, get_rho_label)
names(sim_list) <- rho_labels

## Extract RE tables from each file
RE_beta_all <- purrr::map2_dfr(
  sim_list, names(sim_list),
  ~ mutate(.x$RE_beta, rho_label = .y)
)

## ------------------------------------------------
## 2. Build rho metadata: rho12, rho13, rho23 + Frobenius norm
## ------------------------------------------------

rho_grid_df <- tibble(rho_label = rho_labels) %>%
  mutate(parts = str_split(rho_label, "_")) %>%
  mutate(
    rho12 = map_dbl(parts, ~ as.numeric(.x[2])),
    rho13 = map_dbl(parts, ~ as.numeric(.x[3])),
    rho23 = map_dbl(parts, ~ as.numeric(.x[4])),
    rho_strength = sqrt(rho12^2 + rho13^2 + rho23^2),
    scenario_label = paste0("(",
                            sprintf("%.2f", rho12), ", ",
                            sprintf("%.2f", rho13), ", ",
                            sprintf("%.2f", rho23), ")")
  ) %>%
  select(-parts)

## ------------------------------------------------
## 3. Attach outcome / covariate labels and rho info
## ------------------------------------------------

px <- 4

RE_all2 <- RE_beta_all %>%
  left_join(rho_grid_df, by = "rho_label") %>%
  mutate(
    Param_index = as.integer(str_remove(Parameter, "Beta")),
    Outcome     = paste0("Outcome ", ceiling(Param_index / px)),
    Covariate   = paste0("X", ((Param_index - 1) %% px) + 1)
  )

## ------------------------------------------------
## 4. Build long RE table for both plots
## ------------------------------------------------

RE_long <- RE_all2 %>%
  select(rho_label, scenario_label, rho_strength,
         Outcome, Covariate,
         RE_EX_vs_UN, RE_IN_vs_UN) %>%
  pivot_longer(
    cols = c(RE_EX_vs_UN, RE_IN_vs_UN),
    names_to = "Comparison",
    values_to = "RE_value"
  ) %>%
  mutate(
    Comparison = recode(Comparison,
                        RE_EX_vs_UN = "EX vs UN",
                        RE_IN_vs_UN = "IN vs UN")
  )



p_RE_categorical <- ggplot(RE_long,
                           aes(x = scenario_label, y = log(RE_value),
                               group = Comparison, color = Comparison)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed",
             linewidth = 0.7, alpha = 0.5) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Outcome + Covariate, ncol = 4) +
  theme_sim() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "Correlation scenario (rho12, rho13, rho23)",
    y = "Relative efficiency vs UN (log variance ratio)",
    title = "Relative Efficiency of EX and IN vs UN across correlation scenarios"
  )

p_RE_categorical



p_RE_strength <- RE_long %>%
  ggplot(aes(x = rho_strength, y = log(RE_value),
             group = Covariate, color = Covariate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  facet_grid(Comparison ~ Outcome) +
  theme_sim() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "Correlation strength (Frobenius norm)",
    y = "Relative efficiency vs UN (log variance ratio)",
    title = "Relative Efficiency vs Correlation Strength\n(Exchangeable & Independence vs Unstructured)"
  )

p_RE_strength
