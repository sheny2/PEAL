## ===========================================
## Parallel simulation: bias over 100 runs
## ===========================================

suppressPackageStartupMessages({
  library(doParallel)
  library(doRNG)
  library(foreach)
  library(data.table)
  library(dplyr)
  library(MASS)
  library(Matrix)   # for nearPD
  library(reshape2)
  library(ggplot2)
})

## ---- Engine (peal.em.fit.RI_mv lives here)
source("PEAL_Engine_Missing_Multi-RI3-Rho.R")

## -------------------------
## Core settings & true parameters
## -------------------------
N      <- 200  # number of simulations
H      <- 5    # number of sites
px     <- 6
p_bin  <- 3
p_cont <- px - p_bin
py     <- 3   # number of outcomes

X_cols <- paste0("X", 1:px)
Y_cols <- paste0("Y", 1:py)

set.seed(1)
beta_true <- matrix(runif(px * py, -3, 3), nrow = px, ncol = py)

sigma_u_true      <- 0.5
sigma_v_hosp_true <- runif(H, min = 0.5, max = 1)   # fixed across sims
sigma_e_true      <- 0.3

# Correlation (non-exchangeable)
rho12 <- 0.20
rho13 <- 0.30
rho23 <- 0.60
corr_mat_true <- matrix(c(
  1,     rho12, rho13,
  rho12, 1,     rho23,
  rho13, rho23, 1
), nrow = py, byrow = TRUE)

# Ensure SPD
if (min(eigen(corr_mat_true, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
  corr_mat_true <- as.matrix(nearPD(corr_mat_true, corr = TRUE)$mat)
}

sigma_e_vec <- rep(sigma_e_true, py)
De          <- diag(sigma_e_vec, py)
Sigma_e     <- De %*% corr_mat_true %*% De

# True parameter vectors for comparison
true_beta_vec    <- c(beta_true)
true_sigmaRE_vec <- c(sigma_u_true, sigma_v_hosp_true)  # RE SDs
true_sigmae      <- sigma_e_true
true_corr_vec    <- as.vector(corr_mat_true)            # column-major

par_names_beta    <- paste0("Beta", 1:(px * py))
par_names_sigmaRE <- c("sigma_u", paste0("sigma_v_", 1:H))
par_names_corr    <- paste0("rho[", rep(1:py, each = py), ",", rep(1:py, times = py), "]")

## -------------------------
## One replicate: DGP + missingness + fit
## -------------------------
run_one <- function(seed) {
  set.seed(seed)
  
  ## --- DGP ---
  m_hosp <- sample(50:75, H, replace = TRUE)
  nn     <- rep(m_hosp, times = 1)
  id.hosp <- rep(seq_len(H), times = m_hosp)
  id.pat  <- sequence(nn)
  n_visits <- sample(1:15, sum(nn), replace = TRUE)
  
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
  X_bin  <- matrix(rbinom(Nobs * p_bin, 1, 0.3), nrow = Nobs, ncol = p_bin)
  X_cont <- matrix(rnorm(Nobs * p_cont, 0, 1),   nrow = Nobs, ncol = p_cont)
  X_hij  <- cbind(X_bin, X_cont)
  
  # Residuals
  epsilon_hij <- MASS::mvrnorm(n = Nobs, mu = rep(0, py), Sigma = Sigma_e)
  
  # Outcomes
  Y_hij <- matrix(0, nrow = Nobs, ncol = py)
  for (k in 1:py) {
    Y_hij[, k] <- X_hij %*% beta_true[, k] + u_h_expanded + v_hi_expanded + epsilon_hij[, k]
  }
  
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
  
  rearranged_data <- merge(three_lvl_dat, visit_count, by = c("site","patient")) %>%
    arrange(site, total_visits, patient) %>%
    mutate(site = factor(site))
  
  ## --- MAR missingness (your Option 3 logic) ---
  rearranged_data_with_missing <- rearranged_data
  n <- nrow(rearranged_data)
  
  p1_base <- 0.2
  p2_base <- 0.3
  p3_base <- 0.5
  
  # Y3 depends on Y1,Y2 observed
  Y3_miss_prob <- ifelse(!is.na(rearranged_data$Y1) & !is.na(rearranged_data$Y2),
                         p3_base + 0.2,
                         p3_base - 0.2)
  Y3_miss_prob <- pmin(pmax(Y3_miss_prob, 0), 1)
  Y3_missing   <- rbinom(n, 1, Y3_miss_prob) == 1
  rearranged_data_with_missing$Y3[Y3_missing] <- NA
  
  # Y1, Y2 depend on Y3 observed
  Y1_miss_prob <- ifelse(!is.na(rearranged_data_with_missing$Y3),
                         p1_base + 0.2,
                         p1_base - 0.1)
  Y2_miss_prob <- ifelse(!is.na(rearranged_data_with_missing$Y3),
                         p2_base + 0.2,
                         p2_base - 0.1)
  Y1_miss_prob <- pmin(pmax(Y1_miss_prob, 0), 1)
  Y2_miss_prob <- pmin(pmax(Y2_miss_prob, 0), 1)
  
  Y1_missing <- rbinom(n, 1, Y1_miss_prob) == 1
  Y2_missing <- rbinom(n, 1, Y2_miss_prob) == 1
  
  rearranged_data_with_missing$Y1[Y1_missing] <- NA
  rearranged_data_with_missing$Y2[Y2_missing] <- NA
  
  # Ensure no row with all Y's missing
  all_y_missing <- apply(
    rearranged_data_with_missing[, paste0("Y",1:py), drop = FALSE],
    1, function(x) all(is.na(x))
  )
  if (any(all_y_missing)) {
    for (i in which(all_y_missing)) {
      y_to_keep <- sample(paste0("Y",1:py), 1)
      rearranged_data_with_missing[i, y_to_keep] <-
        rearranged_data[i, y_to_keep]
    }
  }
  
   ## --- Fit PEAL EM (unstructured) ---
  fit <- tryCatch(
    peal.em.fit.RI_mv(
      # data       = rearranged_data,
      data       = rearranged_data_with_missing,
      X_cols     = X_cols,
      Y_cols     = Y_cols,
      site_col   = "site",
      patient_col= "patient",
      corstr_init= "unstructured",
      em_iter    = 1,      
      reml       = TRUE,
      verbose    = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(
      b        = rep(NA_real_, px * py),
      sigmaRE  = rep(NA_real_, H + 1),
      sigma_e  = NA_real_,
      corr_vec = rep(NA_real_, py * py)
    ))
  }
  
  # Fixed effects
  b_hat <- as.numeric(fit$history$init$b)
  
  # Random effects SDs: theta * s2 with the initial estimates before EM
  sigmaRE_hat <- sqrt(fit$history$init$theta * fit$history$init$s2)
  
  # Residual SD
  sigmae_hat <- sqrt(fit$s2)
  
  # Correlation matrix
  Corr_hat <- fit$Corr
  corr_vec_hat <- as.vector(Corr_hat)
  
  list(
    b        = b_hat,
    sigmaRE  = sigmaRE_hat,
    sigma_e  = sigmae_hat,
    corr_vec = corr_vec_hat
  )
}

## -------------------------
## Parallel run
## -------------------------
num_cores <-  max(1, parallelly::availableCores())
cl <- makeCluster(num_cores)
registerDoParallel(cl)
registerDoRNG(12345)

results <- foreach(k = 1:N,
                   .packages = c("data.table","dplyr","MASS","Matrix")) %dopar% {
                     source("PEAL_Engine_Missing_Multi-RI3-Rho.R")
                     run_one(sample(1:1e6, 1))
                   }

stopCluster(cl)

## -------------------------
## Stack results into matrices
## -------------------------
b_mat        <- sapply(results, `[[`, "b")          # (px*py) x N
sigmaRE_mat  <- sapply(results, `[[`, "sigmaRE")    # (H+1) x N
sigmae_vec   <- sapply(results, `[[`, "sigma_e")    # length N
corr_mat_hat <- sapply(results, `[[`, "corr_vec")   # (py*py) x N

colnames(b_mat)       <- paste0("sim", 1:N)
colnames(sigmaRE_mat) <- paste0("sim", 1:N)
names(sigmae_vec)     <- paste0("sim", 1:N)
colnames(corr_mat_hat)<- paste0("sim", 1:N)

## -------------------------
## Make long data frames & bias summaries
## -------------------------
mk_long <- function(mat, par_names, model_label = "PEAL", value_nm = "Estimate") {
  df <- reshape2::melt(mat,
                       varnames = c("ParamIndex", "Simulation"),
                       value.name = value_nm)
  # ParamIndex is the row index (1..nrow(mat))
  df$Parameter <- factor(par_names[df$ParamIndex], levels = par_names)
  df$Model <- model_label
  df
}

# Beta
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

# Random effects SDs
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

# Residual SD
sigmae_df <- data.frame(
  Simulation = factor(paste0("sim", 1:N)),
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

# Correlations
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

## -------------------------
## Example plots: bias distributions
## -------------------------

# Fixed effects bias (density)

FE_plot <- ggplot(beta_long, aes(x = Parameter, y = Bias)) +
  geom_boxplot(outlier.alpha = 0.2, fill = "skyblue", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Bias of Fixed Effects") +
  ylab("Bias (Estimate − True)") +
  xlab("Parameter")


# Random effects SD bias
RE_plot <- ggplot(sigmaRE_long, aes(x = Parameter, y = Bias)) +
  geom_boxplot(outlier.alpha = 0.2, fill = "lightgreen", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Bias of Random Effect SDs") +
  ylab("Bias (Estimate − True)") +
  xlab("Parameter")


Resid_plot <- ggplot(sigmae_df, aes(x = "sigma_e", y = Bias)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "orange", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  ggtitle("Bias of Residual SD") +
  ylab("Bias (Estimate − True)") +
  xlab("")


Corr_plot <- ggplot(corr_long, aes(x = Parameter, y = Bias)) +
  geom_boxplot(outlier.alpha = 0.2, fill = "plum", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Bias of Correlation Estimates") +
  ylab("Bias (Estimate − True)") +
  xlab("Correlation Entry")



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
