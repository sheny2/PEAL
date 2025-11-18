#############################################
## Parallel simulation for peal.fit.RI_mv  ##
## Unstructured correlation, complete data ##
#############################################

suppressPackageStartupMessages({
  library(doParallel)
  library(doRNG)
  library(foreach)
  library(data.table)
  library(dplyr)
  library(MASS)
  library(Matrix)   # nearPD
  library(reshape2)
  library(ggplot2)
})

## --------------------------------------------------
## Engine (peal.fit.RI_mv and generate_Zhv_matrix)
## --------------------------------------------------
source("PEAL_Engine_Multi-RI_Rho.R")

## --------------------------------------------------
## Global settings & true parameters
## --------------------------------------------------
N      <- 200   # number of simulations
H      <- 5     # number of sites
px     <- 6
p_bin  <- 3
p_cont <- px - p_bin
py     <- 3     # number of outcomes

X_cols <- paste0("X", 1:px)
Y_cols <- paste0("Y", 1:py)

set.seed(1)

## Fixed effects: keep fixed across simulations
beta_true <- matrix(runif(px * py, -3, 3), nrow = px, ncol = py)

## Variance components (true)
sigma_u_true       <- 0.5
sigma_v_hosp_true  <- runif(H, min = 0.5, max = 1)  # hospital-specific SD, fixed across sims
sigma_e_true       <- 0.3

## Correlation structure (non-exchangeable)
rho12 <- 0.20
rho13 <- 0.30
rho23 <- 0.60

corr_mat_true <- matrix(c(
  1,     rho12, rho13,
  rho12, 1,     rho23,
  rho13, rho23, 1
), nrow = py, ncol = py, byrow = TRUE)

# Ensure SPD
if (min(eigen(corr_mat_true, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
  corr_mat_true <- as.matrix(nearPD(corr_mat_true, corr = TRUE)$mat)
}

sigma_e_vec <- rep(sigma_e_true, py)
De          <- diag(sigma_e_vec, py)
Sigma_e     <- De %*% corr_mat_true %*% De

## True parameter vectors for bias evaluation
true_beta_vec    <- c(beta_true)
true_sigmaRE_vec <- c(sigma_u_true, sigma_v_hosp_true)  # site + hospital SDs
true_sigmae      <- sigma_e_true
true_corr_vec    <- as.vector(corr_mat_true)            # column-major

par_names_beta    <- paste0("Beta", 1:(px * py))
par_names_sigmaRE <- c("sigma_u", paste0("sigma_v_", 1:H))
par_names_corr    <- paste0("rho[", rep(1:py, each = py), ",", rep(1:py, times = py), "]")

## --------------------------------------------------
## One simulation replicate: DGP + fit
## --------------------------------------------------
run_one <- function(seed) {
  set.seed(seed)
  
  ## ----- DGP -----
  # Number of patients per site
  m_hosp <- sample(50:75, H, replace = TRUE)
  nn     <- rep(m_hosp, times = 1)
  
  id.hosp <- rep(1:H, times = m_hosp)
  id.pat  <- sequence(nn)
  n_visits <- sample(1:20, sum(nn), replace = TRUE)
  
  # Expand IDs to visits
  id.visit         <- sequence(n_visits)
  id.hosp.expanded <- rep(id.hosp, times = n_visits)
  id.pat.expanded  <- rep(id.pat,  times = n_visits)
  
  # Random effects
  u_h  <- rnorm(H, mean = 0, sd = sigma_u_true)
  v_hi <- rnorm(sum(nn), mean = 0, sd = rep(sigma_v_hosp_true, times = m_hosp))
  
  u_h_patient   <- rep(u_h, times = m_hosp)
  u_h_expanded  <- rep(u_h_patient, times = n_visits)
  v_hi_expanded <- rep(v_hi, times = n_visits)
  
  # Covariates
  Nobs   <- sum(n_visits)
  X_bin  <- matrix(rbinom(Nobs * p_bin, size = 1, prob = 0.3), nrow = Nobs, ncol = p_bin)
  X_cont <- matrix(rnorm(Nobs * p_cont, mean = 0, sd = 1),    nrow = Nobs, ncol = p_cont)
  X_hij  <- cbind(X_bin, X_cont)
  
  # Multivariate residuals
  epsilon_hij <- MASS::mvrnorm(n = Nobs, mu = rep(0, py), Sigma = Sigma_e)
  
  # Outcomes
  Y_hij <- matrix(0, nrow = Nobs, ncol = py)
  for (k in 1:py) {
    Y_hij[, k] <- X_hij %*% beta_true[, k] + u_h_expanded + v_hi_expanded + epsilon_hij[, k]
  }
  
  # Build data frame
  three_lvl_dat <- data.table(
    site    = id.hosp.expanded,
    patient = id.pat.expanded,
    visit   = id.visit,
    X_hij, Y_hij
  ) %>% as.data.frame()
  
  data.table::setnames(three_lvl_dat,
                       c("site","patient","visit", paste0("X",1:px), paste0("Y",1:py)))
  
  visit_count <- three_lvl_dat %>%
    group_by(site, patient) %>%
    summarise(total_visits = n(), .groups = "drop")
  
  rearranged_data <- merge(three_lvl_dat, visit_count,
                           by = c("site","patient")) %>%
    arrange(site, total_visits, patient) %>%
    mutate(site = factor(site))
  
  ## ----- Prepare inputs for peal.fit.RI_mv -----
  Y <- as.matrix(rearranged_data[, Y_cols])
  X <- as.matrix(rearranged_data[, X_cols])
  
  Z <- vector("list", H)
  for (i in 1:H) {
    count_mat <- rearranged_data %>%
      dplyr::filter(site == i) %>%
      dplyr::group_by(site, patient) %>%
      dplyr::summarise(n_hi = n(), .groups = "drop")
    Z[[i]] <- generate_Zhv_matrix(count_mat)
  }
  
  id.site <- rearranged_data$site
  
  ## ----- Fit model (unstructured correlation) -----
  fit <- tryCatch(
    peal.fit.RI_mv(
      Y = Y, X = X, Z = Z, id.site = id.site,
      weights = NULL, reml = TRUE,
      corstr = "unstructured",
      estimate_rho = TRUE,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(
      b        = rep(NA_real_, px * py),
      sigmaRE  = rep(NA_real_, H + 1),   # sigma_u + H hospital SDs
      sigma_e  = NA_real_,
      corr_vec = rep(NA_real_, py * py)
    ))
  }
  
  # Extract parameters
  b_hat        <- as.numeric(fit$b)
  sigmaRE_hat  <- sqrt(fit$theta * fit$s2)  # length H+1
  sigmae_hat   <- sqrt(fit$s2)
  Corr_hat     <- fit$Corr
  corr_vec_hat <- as.vector(Corr_hat)
  
  list(
    b        = b_hat,
    sigmaRE  = sigmaRE_hat,
    sigma_e  = sigmae_hat,
    corr_vec = corr_vec_hat
  )
}

## --------------------------------------------------
## Parallel run
## --------------------------------------------------
num_cores <-  max(1, parallelly::availableCores())
cl <- makeCluster(num_cores)
registerDoParallel(cl)
registerDoRNG(12345)

results <- foreach(k = 1:N,
                   .packages = c("data.table","dplyr","MASS","Matrix")) %dopar% {
  source("PEAL_Engine_Multi-RI_Rho.R")
  run_one(sample(1:1e6, 1))
}

stopCluster(cl)

## --------------------------------------------------
## Stack results into matrices
## --------------------------------------------------
b_mat        <- sapply(results, `[[`, "b")          # (px*py) x N
sigmaRE_mat  <- sapply(results, `[[`, "sigmaRE")    # (H+1) x N
sigmae_vec   <- sapply(results, `[[`, "sigma_e")    # length N
corr_mat_hat <- sapply(results, `[[`, "corr_vec")   # (py*py) x N

colnames(b_mat)       <- paste0("sim", 1:N)
colnames(sigmaRE_mat) <- paste0("sim", 1:N)
names(sigmae_vec)     <- paste0("sim", 1:N)
colnames(corr_mat_hat)<- paste0("sim", 1:N)

## --------------------------------------------------
## Make long data frames & bias summaries
## --------------------------------------------------
mk_long <- function(mat, par_names, model_label = "PEAL", value_nm = "Estimate") {
  # melt matrix to get row (ParamIndex) + col (Simulation)
  df <- reshape2::melt(mat,
                       varnames = c("ParamIndex", "Simulation"),
                       value.name = value_nm)
  df$Parameter <- factor(par_names[df$ParamIndex], levels = par_names)
  df$Model     <- model_label
  df
}

## Fixed effects
beta_long <- mk_long(b_mat, par_names_beta, "PEAL", "Estimate")
beta_long$True <- rep(true_beta_vec, times = N)
beta_long <- beta_long %>%
  mutate(Bias = Estimate - True)

beta_summary <- beta_long %>%
  group_by(Parameter) %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias_mean  = mean(Bias, na.rm = TRUE),
    Bias_sd    = sd(Bias, na.rm = TRUE),
    .groups = "drop"
  )

## Random effects SDs
sigmaRE_long <- mk_long(sigmaRE_mat, par_names_sigmaRE, "PEAL", "Estimate")
sigmaRE_long$True <- rep(true_sigmaRE_vec, times = N)
sigmaRE_long <- sigmaRE_long %>%
  mutate(Bias = Estimate - True)

sigmaRE_summary <- sigmaRE_long %>%
  group_by(Parameter) %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias_mean  = mean(Bias, na.rm = TRUE),
    Bias_sd    = sd(Bias, na.rm = TRUE),
    .groups = "drop"
  )

## Residual SD
sigmae_df <- data.frame(
  Simulation = factor(names(sigmae_vec)),
  Estimate   = as.numeric(sigmae_vec)
)
sigmae_df$True <- true_sigmae
sigmae_df <- sigmae_df %>%
  mutate(Bias = Estimate - True)

sigmae_summary <- sigmae_df %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias_mean  = mean(Bias, na.rm = TRUE),
    Bias_sd    = sd(Bias, na.rm = TRUE)
  )

## Correlations
corr_long <- mk_long(corr_mat_hat, par_names_corr, "PEAL", "Estimate")
corr_long$True <- rep(true_corr_vec, times = N)
corr_long <- corr_long %>%
  mutate(Bias = Estimate - True)

corr_summary <- corr_long %>%
  group_by(Parameter) %>%
  summarise(
    True_Value = first(True),
    Mean_Est   = mean(Estimate, na.rm = TRUE),
    Bias_mean  = mean(Bias, na.rm = TRUE),
    Bias_sd    = sd(Bias, na.rm = TRUE),
    .groups = "drop"
  )

## --------------------------------------------------
## Boxplots of parameter estimate bias
## --------------------------------------------------

## 1. Fixed effects β
FE_plot <- ggplot(beta_long, aes(x = Parameter, y = Bias)) +
  geom_boxplot(outlier.alpha = 0.2, fill = "skyblue", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Bias distribution of fixed effects (β) across 100 simulations") +
  ylab("Bias (Estimate - True)") +
  xlab("Parameter")

## 2. Random-effects SDs
RE_plot <- ggplot(sigmaRE_long, aes(x = Parameter, y = Bias)) +
  geom_boxplot(outlier.alpha = 0.2, fill = "lightgreen", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Bias distribution of random-effect SDs across 100 simulations") +
  ylab("Bias (Estimate - True)") +
  xlab("Parameter")

## 3. Residual SD σ_e (single box)
Resid_plot <- ggplot(sigmae_df, aes(x = "sigma_e", y = Bias)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "orange", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  ggtitle("Bias distribution of residual SD (sigma_e) across 100 simulations") +
  ylab("Bias (Estimate - True)") +
  xlab("")

## 4. Correlation entries
Corr_plot <- ggplot(corr_long, aes(x = Parameter, y = Bias)) +
  geom_boxplot(outlier.alpha = 0.2, fill = "plum", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Bias distribution of correlation estimates across 100 simulations") +
  ylab("Bias (Estimate - True)") +
  xlab("Correlation entry")

## Inspect numeric summaries in console
beta_summary
sigmaRE_summary
sigmae_summary
corr_summary





corr_long$Group     <- "Correlation"
sigmae_df$Group     <- "Residual_SD"
sigmaRE_long$Group  <- "RandomEffect_SD"

var_df <- rbind(
  corr_long[, c("Parameter","Simulation","Estimate","Bias","Group")] %>% filter(Parameter %in% c("rho[1,2]","rho[1,3]","rho[2,3]")),
  cbind("Parameter"="sigma_e",
        sigmae_df[, c("Simulation","Estimate","Bias","Group")]),
  sigmaRE_long[, c("Parameter","Simulation","Estimate","Bias","Group")]
)

var_df$Parameter <- factor(var_df$Parameter,
                           levels = unique(var_df$Parameter))

Var_plot <- ggplot(var_df,
                   aes(x = Parameter, y = Bias, fill = Group)) +
  geom_boxplot(outlier.alpha = 0.2, alpha = 0.6, 
               # box.color = NA, 
               # median.colour = NA
  ) +
  scale_fill_manual(values = c(
    "Correlation"       = "#1E88E5",  # blue
    "Residual_SD"       = "#D81B60",  # pink/red
    "RandomEffect_SD"   = "#FFC107"   # amber/yellow
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  ggtitle("Bias of Correlation, Residual SD, and Random-Effect SD") +
  ylab("Bias (Estimate − True)") +
  xlab("Parameter")


gridExtra::grid.arrange(FE_plot, RE_plot, Resid_plot, Corr_plot, nrow = 2)
gridExtra::grid.arrange(FE_plot, Var_plot, nrow = 2)

