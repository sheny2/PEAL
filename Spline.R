library(splines)
library(ggplot2)

x <- seq(0, 10, length.out = 100)
y <- cos(x) + rnorm(100, sd = 0.3)
data <- data.frame(x, y)

# Linear model without spline
lm_linear <- lm(y ~ x, data = data)

# Model with natural cubic spline (3 degrees of freedom)
lm_spline <- lm(y ~ ns(x, df = 3), data = data)

summary(lm_linear)
summary(lm_spline)

# Compare fits
data$pred_linear <- predict(lm_linear)
data$pred_spline <- predict(lm_spline)

ggplot(data, aes(x, y)) +
  geom_point() +
  geom_line(aes(y = pred_linear), color = "red", linewidth = 1) +
  geom_line(aes(y = pred_spline), color = "blue", linewidth = 1) +
  labs(title = "Linear vs Spline Regression",
       subtitle = "Red: Linear, Blue: Natural Spline (df=3)")


# Using bs() for B-splines
lm_bspline <- lm(y ~ bs(x, degree = 3, knots = c(2, 5, 8)), data = data)

# Plotting
data$pred_bspline <- predict(lm_bspline)
ggplot(data, aes(x, y)) +
  geom_point() +
  geom_line(aes(y = pred_bspline), color = "green", size = 1) +
  labs(title = "B-spline Regression")




#############################

# Create simulated data
n <- 200
x <- seq(0, 10, length.out = n)

# True underlying function with a change point at x=6
true_function <- function(x) {
  ifelse(x < 6,
         1.5 * sin(0.8 * x),  # First regime
         0.5 * x - 1.5)        # Second regime (linear)
}

# Generate response with noise
y <- true_function(x) + rnorm(n, sd = 0.3)

# Create dataframe
sim_data <- data.frame(x = x, y = y, true_y = true_function(x))

# 1. Simple linear model (will miss the nonlinearity)
linear_model <- lm(y ~ x, data = sim_data)

# 2. Polynomial model (degree 3)
poly_model <- lm(y ~ poly(x, 3), data = sim_data)

# 3. Natural cubic spline with 5 degrees of freedom
spline_model <- lm(y ~ ns(x, df = 5), data = sim_data)

# 4. Spline with knot at the change point (x=6)
spline_knot_model <- lm(y ~ ns(x, knots = 6), data = sim_data)

# Generate predictions
sim_data$pred_linear <- predict(linear_model)
sim_data$pred_poly <- predict(poly_model)
sim_data$pred_spline <- predict(spline_model)
sim_data$pred_spline_knot <- predict(spline_knot_model)

# Create plot
ggplot(sim_data, aes(x = x)) +
  geom_point(aes(y = y), alpha = 0.5) +                   # Observed data
  geom_line(aes(y = true_y), color = "black", size = 1.2) + # True function
  geom_line(aes(y = pred_linear), color = "red", size = 1) + # Linear
  geom_line(aes(y = pred_poly), color = "blue", size = 1) +   # Polynomial
  geom_line(aes(y = pred_spline), color = "green", size = 1) + # Spline (auto knots)
  geom_line(aes(y = pred_spline_knot), color = "purple", size = 1, linetype = "dashed") + # Spline (knot at 6)
  labs(title = "Comparison of Regression Approaches",
       subtitle = "Black: True, Red: Linear, Blue: Cubic, Green: Spline (auto knots), Purple: Spline (knot at 6)",
       x = "Predictor (x)",
       y = "Response (y)") +
  theme_minimal()

