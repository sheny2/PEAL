## ------------------------------------------------
## MV-PEAL Multivariate Simulation – Relative Efficiency vs Missingness
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
## 1. Load AllOut_* objects (missingness settings)
## ------------------------------------------------

files_allout <- list.files(
  pattern = "^AllOut_miss_.*\\.rds$",  # e.g. AllOut_miss_1.rds, ...
  full.names = TRUE
)

get_setting_label <- function(path) {
  basename(path) |>
    str_remove("^AllOut_") |>
    str_remove("\\.rds$")  # -> "miss_1", "miss_2", ...
}

sim_list   <- purrr::map(files_allout, readRDS)
setting_labels <- purrr::map_chr(files_allout, get_setting_label)
names(sim_list) <- setting_labels

## Extract RE tables from each file
RE_beta_all <- purrr::map2_dfr(
  sim_list, names(sim_list),
  ~ dplyr::mutate(.x$RE_beta, setting = .y)
)

## ------------------------------------------------
## 2. Rebuild missingness grid metadata (must match simulation file)
## ------------------------------------------------

miss_grid_df <- data.frame(
  setting  = paste0("miss_", 1:8),
  p1_base  = seq(0.1, 0.8, length.out = 8),
  p2_base  = seq(0.1, 0.8, length.out = 8),
  p3_base  = seq(0.1, 0.8, length.out = 8)
) %>%
  mutate(
    missing_base_avg = round((p1_base + p2_base + p3_base) / 3, 2)
  )

miss_grid_df

## ------------------------------------------------
## 3. Attach outcome / covariate labels and missingness info
## ------------------------------------------------

px <- 4   # from simulation

RE_all2 <- RE_beta_all %>%
  left_join(miss_grid_df, by = "setting") %>%
  mutate(
    Param_index = as.integer(str_remove(Parameter, "Beta")),
    Outcome     = paste0("Outcome ", ceiling(Param_index / px)),
    Covariate   = paste0("X", ((Param_index - 1) %% px) + 1)
  )

## ------------------------------------------------
## 4. Build long RE table for plotting
## ------------------------------------------------

RE_long <- RE_all2 %>%
  select(setting, missing_base_avg,
         Outcome, Covariate,
         RE_EX_vs_UN, RE_IN_vs_UN) %>%
  pivot_longer(
    cols      = c(RE_EX_vs_UN, RE_IN_vs_UN),
    names_to  = "Comparison",
    values_to = "RE_value"
  ) %>%
  mutate(
    Comparison = dplyr::recode(Comparison,
                               RE_EX_vs_UN = "EX vs UN",
                               RE_IN_vs_UN = "IN vs UN")
  )

## ------------------------------------------------
## 5. Plot 1: RE vs average missingness (faceted by Covariate & Outcome)
## ------------------------------------------------

p_RE_missing <- ggplot(
  RE_long,
  aes(x = missing_base_avg, y = log(RE_value),
      group = Comparison, color = Comparison)
) +
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
    x = "Average missingness probability across outcomes",
    y = "Relative efficiency vs UN (log variance ratio)",
    title = "Relative Efficiency of EX and IN vs UN across missingness levels"
  )

print(p_RE_missing)
# ggsave("RE_vs_missing_faceted.png", p_RE_missing, width = 10, height = 7)

## ------------------------------------------------
## 6. Plot 2: RE vs missingness, faceted by Comparison × Outcome
## ------------------------------------------------

p_RE_missing_strength <- RE_long %>%
  ggplot(aes(x = missing_base_avg, y = log(RE_value),
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
    x = "Average missingness probability across outcomes",
    y = "Relative efficiency vs UN (log variance ratio)",
    title = "Relative Efficiency vs Average Missingness\n(Exchangeable & Independence vs Unstructured)"
  )

print(p_RE_missing_strength)
# ggsave("RE_vs_missing_by_comparison.png", p_RE_missing_strength, width = 10, height = 7)
