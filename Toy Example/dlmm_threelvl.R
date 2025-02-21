
#### DLMM for a three-level implementation
# rm(list=ls())
require(tidyverse)
require(lme4)
require(nlme)
require(mvtnorm)

source('dlmm_threelvlengine.R')

set.seed(2025)


# Parameters
H <- 10  # Number of hospitals
m <- 50  # Patients per hospital

beta0 <- 5
px <- 5
beta <- c(5, 3, -2, 4, 1)  # Fixed effect coefficients for 5 covariates
sigma_u <- 2  # Hospital-level variance
sigma_v <- 1.5  # Patient-level variance
sigma_e <- 1  # Error variance

three_lvl_dat <- data.frame()


for (h in 1:H) {
  u_h <- rnorm(1, mean = 0, sd = sigma_u)
  for (i in 1:m) {
    v_hi <- rnorm(1, mean = 0, sd = sigma_v)
    n_i <- sample(1:10, 1)  # Random number of visits for each patient
    for (j in 1:n_i) {
      X_hij <- rbinom(px, size = 1, prob = 0.3)  # binary covariates
      epsilon_hij <- rnorm(1, mean = 0, sd = sigma_e)
      y_hij <- beta0 + sum(X_hij * beta) + u_h + v_hi + epsilon_hij

      three_lvl_dat <- rbind(three_lvl_dat, c(h, i, j, X_hij, y_hij))
    }
  }
}

# Add column names
colnames(three_lvl_dat) <- c("site", "patient", "visit", paste0("X", 1:5), "Y")

head(three_lvl_dat)



# Calculate the number of visits per patient per site
visit_count <- three_lvl_dat %>%
  group_by(site, patient) %>%
  summarise(total_visits = n(), .groups = "drop")

# Reorder data (as a preprocess probably)
rearranged_data <- merge(three_lvl_dat, visit_count, by = c("site", "patient")) %>%
  arrange(site, total_visits, patient)


# Lmer
fit03.lmer <- lmer(Y ~ X1 + X2 + X3 + X4 + X5 + (1|site/patient), data = rearranged_data,
                   # control=lmerControl(optimizer="nloptwrap",
                   #                     optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))
                   # control=lmerControl(optimizer="nloptwrap",
                   #                       optCtrl=list(algorithm="NLOPT_LN_COBYLA",
                   #                                    xtol_rel=1e-6,
                   #                                    xtol_abs=1e-10,
                   #                                    ftol_abs=1e-10)))
                   control = lmerControl(optimizer="bobyqa"))
                   # control = lmerControl(optimizer="Nelder_Mead"))

summary(fit03.lmer)


fit03.lmer <- lmer(Y ~ X1 + X2 + X3 + X4 + X5 + (1|site/patient),
                   data = rearranged_data,
                   control = lmerControl(optimizer="bobyqa"))
summary(fit03.lmer)



# XYZ
Y <- rearranged_data$Y

X <- as.matrix(rearranged_data[, paste0("X", 1:px)])
X <- cbind(1, X)

count_mat = rearranged_data %>%
  group_by(site, patient) %>%
  summarise(n_hi = n(), .groups = 'drop')
Z <- generate_Z_matrix(count_mat, H)

id.site <- rearranged_data$site

ShXYZ <- lmm.get.summary3(Y, X, Z, id.site)


# DLMM
fit03.dlmm = lmm.fit3(Y = NULL, X = NULL, Z = NULL,
                      id.site = NULL, weights = NULL,
           pooled = F, reml = T,
           common.s2 = T,
           ShXYZ = ShXYZ,  # only need summary stats
           corstr = 'independence',
           mypar.init = NULL)

fit03.dlmm.pool = lmm.fit3(Y = Y, X = X, Z = Z,
                           id.site = id.site, weights = NULL,
                 pooled = T, reml = T,
                 common.s2 = T,
                 ShXYZ = NULL,
                 corstr = 'independence',
                 mypar.init = NULL)

# Compare fixed effect estimates
cbind(fit03.dlmm$b, summary(fit03.lmer)$coef[,1])


# sigma_u sigma_v sigma_e est est
# sqrt(diag(fit03.dlmm$V))
# sqrt(fit03.dlmm$s2)

dlmm_sigma_u <- sqrt(diag(fit03.dlmm$V))[1]
dlmm_sigma_v <- sqrt(diag(fit03.dlmm$V))[2]
dlmm_sigma_e <- sqrt(fit03.dlmm$s2)

c(dlmm_sigma_u, dlmm_sigma_v, dlmm_sigma_e)


var_comp <- as.data.frame(VarCorr(fit03.lmer))
lmm_sigma_v <- var_comp$sdcor[1]
lmm_sigma_u <- var_comp$sdcor[2]
lmm_sigma_e <- var_comp$sdcor[3]

c(lmm_sigma_u, lmm_sigma_v, lmm_sigma_e)


# pool no pool
plot(fit03.dlmm$b, fit03.dlmm.pool$b)

# prediction
plot(predict(fit03.lmer), fit03.dlmm.pool$Yihat)
