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

```r
source("../PEAL Engine/PEAL_Engine_RI.R")
```

We've prepared a simulated three-level dataset ready for analysis
```r
three_lvl_dat <- readRDS("three_lvl_dat.rds")
head(three_lvl_dat)
```

| | site | patient | visit | X1 | X2 | X3 | X4 | X5 | X6 | Y |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 1 | 1 | 1 | 0 | 0 | -0.2333860 | 0.5258316 | 0.8517787 | 9.713744 |
| 2 | 1 | 1 | 2 | 0 | 0 | 0 | -1.7007881 | -2.3441436 | 0.2597312 | -16.215685 |
| 3 | 1 | 1 | 3 | 1 | 0 | 0 | -0.6558381 | -0.1574295 | -1.1642077 | -5.654015 |
| 4 | 1 | 1 | 4 | 0 | 0 | 0 | -0.6304705 | -1.4508508 | -1.4120594 | -15.497097 |
| 5 | 1 | 1 | 5 | 0 | 1 | 0 | -0.9323276 | -0.3451917 | -1.2283032 | -7.509270 |
| 6 | 1 | 2 | 1 | 1 | 1 | 1 | 0.3219037 | -1.6069905 | -0.5368642 | -3.596303 |


### Data Preparation
PEAL operates on hierarchical data. Our dataset, `three_lvl_dat`, contains records of patient visits nested within sites. 
We first perform a simple preprocessing step to reorder the data, which can sometimes help with model fitting, though it's not strictly required by PEAL.

```r
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

```r
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

```r
# Fit the model using the RI engine
fit.peal <- peal.fit.RI(
  Y = NULL, X = NULL, Z = NULL, id.site = NULL,
  pooled = F,
  reml = T,
  ShXYZ = ShXYZ # Only summary stats are needed!
)
```


(optional) To prove the lossless property, we run a standard pooled analysis with lme4 and compare the results.

```r
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

| Fixed Effects | lmm_beta | peal_beta |
|---|---|---|
| X1 | 1.013306 | 1.013307 |
| X2 | 2.093979 | 2.093980 |
| X3 | 2.938045 | 2.938046 |
| X4 | 4.008884 | 4.008884 |
| X5 | 5.023196 | 5.023196 |
| X6 | 6.007726 | 6.007727 |

| Random Effect | lmm_sigma | peal_sigma |
|---|---|---|
| site | 0.1094553 | 0.1094358 |
| site1:patient | 0.9705814 | 0.9705763 |
| site2:patient | 0.9562715 | 0.9562770 |
| site3:patient | 0.8411468 | 0.8411524 |
| residual | 0.7070876 | 0.7070875 |


The output shows that the estimates for both fixed effects and random effect standard deviations are virtually identical, confirming PEAL's lossless nature.




## Part 2: PEAL with Random Slopes (RS)
Now, we extend the model to include random slopes, allowing the effect of a covariate (e.g., X6) to vary by patient.


### Step 2.1: Generate Summary Statistics for RS Model
The process is similar, but we first need to define the design matrices for the random slope components.
```r
# Load the Random Slope engine
source("../PEAL Engine/PEAL_Engine_RS.R")

# Define formulas for the fixed effects and random slopes
formula <- ~ X1 + X2 + X3 + X4 + X5 + X6 - 1
X <- model.matrix(formula, data = rearranged_data)

formula_RS <- ~ X6 - 1
X_RS <- model.matrix(formula_RS, data = rearranged_data)

# Construct the Z matrix list and get summaries
Z <- construct_Z_list(rearranged_data, colnames(X_RS))
m_h_all <- (rearranged_data %>% group_by(site) %>% dplyr::summarise(n_distinct(patient)))[, 2]
ShXYZ <- peal.get.summary(Y, X, Z, id.site, m_h_all = m_h_all)
```

### Step 2.2: Fit the PEAL RS Model
We use the peal.fit.RS function, again providing only the summary statistics.

```r
# Fit the RS model
fit.peal.rs <- peal.fit.RS(
  Y = NULL, X = NULL, Z = NULL, id.site = NULL,
  pooled = F,
  reml = T,
  ShXYZ = ShXYZ
)
```

### Step 2.3: Compare with `lme4` for RS
We again verify the results against a corresponding lme4 model.

```r
# Fit the equivalent LMM with random slopes
fit.lmer.rs <- lmer(
  Y ~ -1 + X1 + X2 + X3 + X4 + X5 + X6 +
       (1 | site) +                          
       (0 + dummy(site, "1")|| site:patient) +  (0 +  covar_composite_1 || site:patient) +  
       (0 + dummy(site, "2")|| site:patient)  + (0 +  covar_composite_2 || site:patient) +  
       (0 + dummy(site, "3")|| site:patient)  + (0 +  covar_composite_3 || site:patient),   
  data = rearranged_data %>%
    mutate(
      covar_composite_1 = dummy(site, "1") * X6,
      covar_composite_2 = dummy(site, "2") * X6,
      covar_composite_3 = dummy(site, "3") * X6
    ),
  control = lmerControl(
    optimizer = "optimx",
    optCtrl = list(method = "L-BFGS-B")
  )
)
```

| Fixed Effects | lmm_beta | peal_beta |
|---|---|---|
| X1 | 1.025053 | 1.025055 |
| X2 | 2.104015 | 2.104018 |
| X3 | 2.936898 | 2.936901 |
| X4 | 3.999226 | 3.999225 |
| X5 | 5.029522 | 5.029522 |
| X6 | 6.012435 | 6.012434 |

| Random Effect | lmm_sigma | peal_sigma |
|---|---|---|
| site RI | 0.1112932 | 0.1113393 |
| pt RI by site1 | 0.9590522 | 0.9590176 |
| pt RI by site2 | 0.9623781 | 0.9624495 |
| pt RI by site3 | 0.8379339 | 0.8378841 |
| pt RS by site1 | 0.3219149 | 0.3219251 |
| pt RS by site2 | 0.2795292 | 0.2795442 |
| pt RS by site3 | 0.2733985 | 0.2733912 |
| residual | 0.6379580 | 0.6379587 |

Once again, the outputs for both fixed and random effects match perfectly.



## 📊 Visualizing the Results
A scatter plot is the best way to visualize the one-to-one correspondence of the estimates. The code below generates plots comparing the PEAL and lme4 estimates for both fixed effects and random effect standard deviations.
<img width="924" height="607" alt="Screenshot 2025-08-06 at 7 47 50 PM" src="https://github.com/user-attachments/assets/f464f661-803d-47ab-bcb9-a4bc62abe623" />

