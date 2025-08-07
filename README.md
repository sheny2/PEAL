## PEAL: Privacy-preserving Efficient Aggregation of Longitudinal Data
This repository contains an R implementation of the PEAL (Privacy-preserving Efficient Aggregation of Longitudinal data) algorithm. PEAL is a novel, one-shot distributed algorithm designed to fit three-level linear mixed-effects models on longitudinal data without sharing individual patient data (IPD).

The key features of PEAL are:

* Privacy-Preserving: Operates exclusively on site-level summary statistics, ensuring patient-level data remains local. 

* Communication-Efficient: Requires only a single round of communication from participating sites. 

* Lossless: Achieves results that are statistically identical to a traditional pooled analysis where all IPD is centralized. 

This tutorial will guide you through simulating a multi-site dataset and using the PEAL engine to fit a model, then comparing the results to a standard linear mixed model fit with `lme4`.


## Part 1: PEAL with Random Intercepts (RI)
In this part, we'll fit a model with a random intercept for each site and patient-specific random intercepts within each site.
You'll need the PEAL engine script. Make sure the path to the engine file is correct.

```
source("../PEAL Engine/PEAL_Engine_RI.R")
```

We've prepared a simulated three-level dataset ready for analysis
```
three_lvl_dat <- readRDS("three_lvl_dat.rds")
head(three_lvl_dat)
```

### Data Preparation
PEAL operates on hierarchical data. Our dataset, `three_lvl_dat`, contains records of patient visits nested within sites. 
We first perform a simple preprocessing step to reorder the data, which can sometimes help with model fitting, though it's not strictly required by PEAL.

```ruby
# Calculate the number of visits per patient to use for reordering
visit_count <- three_lvl_dat %>%
  dplyr::group_by(site, patient) %>%
  dplyr::summarise(total_visits = n(), .groups = "drop")

# Merge visit counts and reorder the data by site and patient
rearranged_data <- merge(three_lvl_dat, visit_count, by = c("site", "patient")) %>%
  arrange(site, total_visits, patient) %>%
  mutate(site = factor(site))
```


### Step 1.1: Generate Summary Statistics
This is the core of PEAL's privacy-preserving approach. Each site generates summary statistics from its local data. No individual patient data is shared. Here, we simulate this process by looping through our combined dataset. The function 

```ruby
# Define outcome Y, fixed-effects matrix X, and number of sites H
Y <- rearranged_data$Y
px <- 6 # Number of covariates
X <- as.matrix(rearranged_data[, paste0("X", 1:px)])
H <- n_distinct(rearranged_data$site)

# Generate the random-effects design matrix (Z) for each site
Z <- list()
for (i in 1:H) {
  count_mat <- rearranged_data %>%
    filter(site == i) %>%
    group_by(site, patient) %>%
    dplyr::summarise(n_hi = n(), .groups = 'drop')
  
  Z[[i]] <- (generate_Zhv_matrix(count_mat))
}

# Get the list of summary statistics
id.site <- rearranged_data$site
ShXYZ <- peal.get.summary(Y, X, Z, id.site)
```

### Step 1.2: Fit the PEAL Model
The central server uses only the aggregated summary statistics (ShXYZ) to fit the model. Notice that Y, X, and Z are set to NULL.

```ruby
# Fit the model using the RI engine
fit.peal <- peal.fit.RI(
  Y = NULL, X = NULL, Z = NULL, id.site = NULL,
  pooled = F,
  reml = T,
  ShXYZ = ShXYZ # Only summary stats are needed!
)
```


(optional) To prove the lossless property, we run a standard pooled analysis with lme4 and compare the results.

```ruby
# Fit a standard LMM for comparison
fit.lmer <- lmer(
  Y ~ -1 + X1 + X2 + X3 + X4 + X5 + X6 +
    (1 | site) + # Site-level random intercept
    (0 + dummy(site, "1") | site:patient) +
    (0 + dummy(site, "2") | site:patient) +
    (0 + dummy(site, "3") | site:patient),
  data = rearranged_data,
  control = lmerControl(optimizer = "optimx", optCtrl = list(method = "L-BFGS-B"))
)

# Compare Fixed Effects
peal_beta <- c(fit.peal$b)
lmm_beta <- fixef(fit.lmer)
print(cbind(lmm_beta, peal_beta))

# Compare Random Effects (Standard Deviations)
summ.lmer <- summary(fit.lmer)
lmm_sigma <- c(
  site = sqrt(summ.lmer$varcor$site),
  diag(sqrt(summ.lmer$varcor$site.patient.2)),
  diag(sqrt(summ.lmer$varcor$site.patient.1)),
  diag(sqrt(summ.lmer$varcor$site.patient)),
  residual = summ.lmer$sigma
)
peal_sigma <- c(sqrt(fit.peal$V), sqrt(fit.peal$s2))
print(cbind(lmm_sigma, peal_sigma))
```

The output shows that the estimates for both fixed effects and random effect standard deviations are virtually identical, confirming PEAL's lossless nature.
