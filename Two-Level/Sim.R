library(ggplot2)
library(gridExtra)

set.seed(20250412)
# Set up parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Simulation parameters
n_sim <- 50
K <- 5
px <- 5
pz <- px + 1

s2 <- rep(2, K)
V = diag(runif(n = pz, 0, 3))

# True parameters
true_beta <- c(1:px)/px*10
true_sigma_u <- sqrt(diag(V))
true_sigma_e <- sqrt(s2[1])

# Run simulations with bias calculation
sim_results <- foreach(i = 1:n_sim, .combine = rbind,
                       .packages = c("Matrix", "mvtnorm")) %dopar% {

                         nn <- sample(1000, K, replace = TRUE)
                         X <- matrix(NA, sum(nn), px)
                         Z <- matrix(NA, sum(nn), pz)
                         id.site <- as.character(rep(1:K, nn))
                         Y <- rep(NA, sum(nn))
                         u <- matrix(NA, K, pz)

                         source("DLMM-2.R")

                         for(ii in 1:K){
                           X[id.site==ii,] <- matrix(rbinom(nn[ii]*px, size=1, prob = 0.3), nn[ii], px)
                           Z[id.site==ii,] <- cbind(1, X[id.site==ii,])
                           ui <- c(mvtnorm::rmvnorm(1, rep(0,pz), V))
                           u[ii,] <- ui
                           ei <- rnorm(nn[ii], 0, sqrt(s2[ii]))
                           Y[id.site==ii] <- X[id.site==ii,]%*%true_beta + Z[id.site==ii,] %*%ui + ei
                         }

                         # Fit model
                         ShXYZ <- lmm.get.summary2(Y, X, Z, id.site)
                         fit02.dlmm <- tryCatch({
                           lmm.fit2(ShXYZ = ShXYZ, corstr = 'independence', reml = T)
                         }, error = function(e) NULL)

                         if(is.null(fit02.dlmm)) return(NULL)

                         # source("dlmm-engine.R")
                         # SiXYZ <- lmm.get.summary(Y, X, Z, weights = NULL, id.site)
                         # fit02.dlmm <- tryCatch({
                         #   lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=T, common.s2=T, hessian=T)
                         # }, error = function(e) NULL)
                         # if(is.null(fit02.dlmm)) return(NULL)

                         # Calculate bias
                         est_beta <- fit02.dlmm$b
                         est_sigma_u <- diag(sqrt(fit02.dlmm$V))
                         est_sigma_e <- sqrt(fit02.dlmm$s2)

                         bias_beta <- est_beta - true_beta
                         bias_sigma_u <- est_sigma_u - true_sigma_u
                         bias_sigma_e <- est_sigma_e - true_sigma_e

                         data.frame(
                           sim = i,
                           param = c(paste0("beta", 1:px), paste0("sigma_u", 1:pz), "sigma_e"),
                           true = c(true_beta, true_sigma_u, true_sigma_e),
                           est = c(est_beta, est_sigma_u, est_sigma_e),
                           bias = c(bias_beta, bias_sigma_u, bias_sigma_e)
                         )
                       }

stopCluster(cl)

# Process results
sim_results <- (sim_results)
sim_results$param <- factor(sim_results$param,
                            levels = c(paste0("beta", 1:px), paste0("sigma_u", 1:pz), "sigma_e"))

# Calculate mean bias and confidence intervals
bias_summary <- do.call(rbind, lapply(split(sim_results, sim_results$param), function(x) {
  data.frame(
    param = unique(x$param),
    mean_bias = mean(x$bias),
    lower = quantile(x$bias, 0.025),
    upper = quantile(x$bias, 0.975)
  )
}))

# Create bias visualization
bias_plot <- ggplot(bias_summary, aes(x = param, y = mean_bias)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Bias in Parameter Estimates with 95% CI",
       x = "Parameter",
       y = "Bias (Estimate - True)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Relative bias plot
bias_summary$rel_bias <- bias_summary$mean_bias / rep(c(true_beta, true_sigma_u, true_sigma_e),
                                                      each = nrow(bias_summary)/length(unique(bias_summary$param)))

rel_bias_plot <- ggplot(bias_summary, aes(x = param, y = rel_bias)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Relative Bias (Bias/True Value)",
       x = "Parameter",
       y = "Relative Bias") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots
# grid.arrange(bias_plot, rel_bias_plot, ncol = 1)
bias_plot

# Alternative: Color by parameter type
sim_results$param_type <- ifelse(grepl("beta", sim_results$param), "Fixed Effects",
                                 ifelse(grepl("sigma_u", sim_results$param), "Random Effects SD",
                                        "Error SD"))

ggplot(sim_results, aes(x = param, y = bias, fill = param_type)) +
  geom_jitter(aes(x = param, y = bias), color = "GRAY", size = 0.5, alpha = 0.3) +
  geom_boxplot(alpha = 0.7) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Distribution of Bias",
       x = "Parameter",
       y = "Bias (Estimate - True)",
       fill = "Parameter Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
