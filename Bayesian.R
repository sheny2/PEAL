
# Fit the model
three_level_brms <- brm(
  formula = Y ~ X1 + X2 + X3 + X4 + X5  + (1 | site) + (1 | site:patient), # Random intercepts
  data = rearranged_data,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "Intercept"),
    prior(cauchy(0, 2), class = "sd")
  ),
  iter = 2000,
  warmup = 500,
  chains = 2,
  cores = getOption("mc.cores", 7),
  control = list(adapt_delta = 0.95)
)

# Summarize the model
summary(three_level_brms)

# Plot the model diagnostics
plot(three_level_brms)


