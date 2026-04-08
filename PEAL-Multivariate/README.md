##  MV-PEAL: Multivariate Privacy-preserving Efficient Aggregation of Longitudinal Data
MV-PEAL is an extension of the original PEAL framework, designed to handle multiple correlated outcomes (multivariate) across three-level longitudinal hierarchies (visits nested within patients, nested within sites).

MV-PEAL is specifically built to address real-world data challenges, such as Missing At Random (MAR) data, while maintaining the core PEAL principles:

* Privacy-Preserving: No individual patient data (IPD) leaves the local site.

* Multivariate Outcomes: Simultaneously models multiple Y outcomes and their correlations.

* Computational Efficiency: Utilizes Rcpp integration (`peal_core.cpp`) for high-performance matrix operations.

<img width="1442" height="669" alt="Screenshot 2026-04-08 at 11 40 27 AM" src="https://github.com/user-attachments/assets/9f42fb3b-be02-4bcb-99b8-3589290d9a79" />


##  Setup and Requirements
To use the MV-PEAL engine, ensure you have the following files in your working directory:

```r
library(tidyverse)
library(lme4)
library(Rcpp)
source("/MV-PEAL Engine/MV-PEAL.R")
```

### Part 1: Data Preparation and Missingness

In multivariate clinical studies, it is common for patients to miss certain visits or for specific lab results (outcomes) to be unavailable. MV-PEAL is robust to these missingness patterns.

Our simulated dataset contains 5 sites with multiple patients and visits. It includes two covariates ($X_1, X_2$) and three distinct outcomes ($Y_1, Y_2, Y_3$).

```r
# Load the simulated multivariate data
sim_data <- readRDS("simulated_mv_data.rds")
mv_dat <- sim_data$observed_data

head(mv_dat)
```
data preview is like
| site | patient | visit | X1 | X2 | Y1 | Y2 | Y3 | total_visits |
|---|---|---|---|---|---|---|---|---|
|1	|4	|1	|1	|-1.6815387	|1.127510	|-3.806824	|-9.609663	|2
|1	|4	|2	|0	|-1.6632470	|10.544026	|-1.232500	|NA	|2
|1	|12	|1	|1	|0.6766099	|-12.698192	|NA	|NA	|2
|1	|12	|2	|0	|1.0738382	|-4.701681	|NA	|NA	|2
|1	|25	|1	|1	|0.7611508	|-12.991691	|1.794230	|15.056124	|2
|1	|25	|2	|0	|-0.1535701	|2.291546	|1.723773	|NA	|2



### Part 2: Fitting the MV-PEAL Model
The `peal.fit` function is the main interface. 
It allows you to specify the fixed effects, multiple outcomes, and the correlation structure for the error terms.

```r
# Define column names
X_names <- c("X1", "X2")
Y_names <- c("Y1", "Y2", "Y3")

# Fit the Multivariate Model
fit_mv <- peal.fit(
    data = mv_dat,
    X_cols = X_names,
    Y_cols = Y_names,
    site_col = "site",
    patient_col = "patient",
    corstr = "unstructured", # Captures correlations between Y1, Y2, and Y3
    reml = TRUE,
    use_rcpp = TRUE          # Enables C++ acceleration
)
```

### Part 3: Extracting Parameter Estimates

Once the MV-PEAL model is fitted, you can extract various parameter estimates to understand the relationships between covariates and the multiple outcomes, as well as the underlying variance structures.

3.1 Fixed Effects ($\beta$)

The fixed effects represent the average relationship between your predictors ($X$) and each outcome ($Y$) across the entire population. In MV-PEAL, these are stored in the b component of the fit object.

```r
# Extract fixed effect coefficients
fixed_effects <- fit_mv$b
rownames(fixed_effects) <- X_names
colnames(fixed_effects) <- Y_names

print(fixed_effects)
```

3.2 Variance Components

MV-PEAL decomposes the total variance into site-level, patient-level, and residual components. These are typically extracted as standard deviations ($\sigma$):

* Site-level SD ($\sigma_u$): Captures variability between different clinical sites.

* Patient-level SD ($\sigma_v$): Captures variability between patients within the same site.

* Residual SD ($\sigma_e$): The remaining unexplained error for each outcome.

```r
# Site-level standard deviation (global)
sigma_u_est <- sqrt(fit_mv$theta[1] * fit_mv$s2)

# Patient-level standard deviations (per site)
sigma_v_hosp_est <- sqrt(fit_mv$theta[2:(H+1)] * fit_mv$s2)

# Residual standard deviation
sigma_e_est <- sqrt(fit_mv$s2)
```

3.3 Outcome Correlations ($\rho$)

A primary advantage of MV-PEAL is its ability to estimate how the different outcomes ($Y_1, Y_2, Y_3$) correlate with one another. 
By using `corstr = "unstructured"`, you can extract the full correlation matrix.

```r
# Extract the correlation matrix between outcomes
cor_matrix <- fit_mv$Corr
rownames(cor_matrix) <- colnames(cor_matrix) <- Y_names

# Extract specific pairwise correlations (e.g., Y2 and Y1, Y3 and Y1, Y3 and Y2)
rhos_est <- c(cor_matrix[2, 1], cor_matrix[3, 1], cor_matrix[3, 2]) 
```

