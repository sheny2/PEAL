library(doParallel)
library(foreach)
library(data.table)
library(dplyr)


# Parameters
H <- 5  # of sites
m_hosp <- sample(500:510, H, replace = T) # of patients (1k to 3k later)


px <- 9  # of covariates
p_bin <- 5  # of binary X
p_cont <- px - p_bin  # of continuous X


beta0 = 5 # intercept
beta <- c(2, 4, 6, 8, 10, 3, 5, 7, 9)  # Fixed effects for covariates

sigma_e <- 1  # error variance
sigma_u <- 3 # site-level variance
sigma_v_hosp <- runif(H, min = 1, max = 5)  # Varying sigma_v by hospital


# Generate data
nn <- rep(m_hosp, times = 1)  # Number of patients per hospital
id.hosp <- rep(1:H, times = m_hosp)   # Hospital ID
id.pat <- sequence(nn)                # Patient ID
n_visits <- sample(1:30, sum(nn), replace = TRUE)  # number of visits of patients

# Expand hospital and patient IDs for visits
id.visit <- sequence(n_visits)
id.hosp.expanded <- rep(id.hosp, times = n_visits)
id.pat.expanded <- rep(id.pat, times = n_visits)

# Random effects
u_h <- rnorm(H, mean = 0, sd = sigma_u)  # Hospital effects
v_hi <- rnorm(sum(nn), mean = 0, sd = rep(sigma_v_hosp, times = m_hosp))  # Patient effects varying by hospital

# Expansion of hospital effects
u_h_patient <- rep(u_h, times = m_hosp)
u_h_expanded <- rep(u_h_patient, times = n_visits)

# Expansion of patient effects
v_hi_expanded <- rep(v_hi, times = n_visits)

# covariates
X_bin <- matrix(rbinom(sum(n_visits) * p_bin, size = 1, prob = 0.3), nrow = sum(n_visits), ncol = p_bin)
X_cont <- matrix(rnorm(sum(n_visits) * p_cont, mean = 0, sd = 1), nrow = sum(n_visits), ncol = p_cont)
X_hij <- cbind(X_bin, X_cont)  # Combine binary & continuous covariates

epsilon_hij <- rnorm(sum(n_visits), mean = 0, sd = sigma_e)

# Compute outcome
y_hij <- beta0 + X_hij %*% beta + u_h_expanded + v_hi_expanded + epsilon_hij


three_lvl_dat <- data.table(
  site = id.hosp.expanded,
  patient = id.pat.expanded,
  visit = id.visit,
  X_hij,
  Y = y_hij
) %>% data.frame()

# Assign column names for covariates
setnames(three_lvl_dat, c("site", "patient", "visit", paste0("X", 1:px), "Y"))

# Display first few rows
print(head(three_lvl_dat))



# Calculate the number of visits per patient per site
visit_count <- three_lvl_dat %>%
  dplyr::group_by(site, patient) %>%
  dplyr::summarise(total_visits = n(), .groups = "drop")

# Reorder data (as a preprocess probably)
rearranged_data <- merge(three_lvl_dat, visit_count,
                         by = c("site", "patient")) %>%
  arrange(site, total_visits, patient) %>%
  mutate(site = factor(site))


# XYZ
Y <- rearranged_data$Y

X <- as.matrix(rearranged_data[, paste0("X", 1:px)])
X <- cbind(1, X)


Z <- list()

for(i in 1:H){
  count_mat = rearranged_data %>%
    filter(site == i) %>%
    group_by(site, patient) %>%
    dplyr::summarise(n_hi = n(), .groups = 'drop')

  Z[[i]] <- (generate_Zhv_matrix(count_mat))
}



id.site <- rearranged_data$site

ShXYZ <- lmm.get.summary3(Y, X, Z, id.site)


# Define number of cores for parallel execution
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)


bench::mark(lmm.fit3(Y = NULL, X = NULL, Z = NULL,
                       id.site = NULL, weights = NULL,
                       pooled = F, reml = T,
                       common.s2 = T,
                       ShXYZ = ShXYZ,  # only need summary stats
                       corstr = 'independence',
                       mypar.init = NULL))

bench::mark(foreach(k = 1, .packages = c("data.table", "dplyr")) %dopar% {

  source("DLMM_Engine_Efficient.R")

  lmm.fit3(Y = NULL, X = NULL, Z = NULL,
           id.site = NULL, weights = NULL,
           pooled = F, reml = T,
           common.s2 = T,
           ShXYZ = ShXYZ,  # only need summary stats
           corstr = 'independence',
           mypar.init = NULL)
})


# Run simulations in parallel using foreach
result <- foreach(k = 1, .packages = c("data.table", "dplyr")) %dopar% {

  source("DLMM_Engine_Efficient.R")

  lmm.fit3(Y = NULL, X = NULL, Z = NULL,
           id.site = NULL, weights = NULL,
           pooled = F, reml = T,
           common.s2 = T,
           ShXYZ = ShXYZ,  # only need summary stats
           corstr = 'independence',
           mypar.init = NULL)
}

stopCluster(cl)





fit03.dlmm <- result[[1]]

fit03.dlmm$b
sqrt(fit03.dlmm$V)
sqrt(fit03.dlmm$s2)



true_beta = c(beta0, beta)
true_beta
true_sigma = c(sigma_u, sigma_v_hosp)
true_sigma
c(sigma_e)




plot(fit03.dlmm$b, true_beta)
points(fit03.dlmm$b[1], true_beta[1], col = "red")
abline(a = 0, b = 1, col = "blue", lwd = 2)

plot(sqrt(fit03.dlmm$V), true_sigma)
points(sqrt(fit03.dlmm$V)[1], true_sigma[1], col = "red")
abline(a = 0, b = 1, col = "blue", lwd = 2)

