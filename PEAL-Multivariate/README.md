##  MV-PEAL: Multivariate Privacy-preserving Efficient Aggregation of Longitudinal Data
MV-PEAL is an extension of the original PEAL framework, designed to handle multiple correlated outcomes (multivariate) across three-level longitudinal hierarchies (visits nested within patients, nested within sites).

MV-PEAL is specifically built to address real-world data challenges, such as Missing At Random (MAR) data, while maintaining the core PEAL principles:

* Privacy-Preserving: No individual patient data (IPD) leaves the local site.

* Multivariate Outcomes: Simultaneously models multiple Y outcomes and their correlations.

* Computational Efficiency: Utilizes Rcpp integration (`peal_core.cpp`) for high-performance matrix operations.

[MV-PEAL_pipeline.pdf](https://github.com/user-attachments/files/26573464/MV-PEAL_pipeline.pdf)


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



