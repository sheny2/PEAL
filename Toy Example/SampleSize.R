library(Matrix)
# Function to generate the record count matrix
generate_record_count <- function(data) {

  counts <- table(data[,"n_hi"])

  result_matrix <- cbind(as.numeric(names(counts)), as.numeric(counts))

  colnames(result_matrix) <- c("n_hi", "frequency")

  return(result_matrix)
}


# Function to generate the Z_hv matrix
generate_Zhv_matrix <- function(data) {

  record_count_matrix <- generate_record_count(data)

  diagonal_blocks <- list()

  for (i in 1:nrow(record_count_matrix)) {
    n_hi <- record_count_matrix[i, "n_hi"]
    frequency <- record_count_matrix[i, "frequency"]

    identity_block <- diag(frequency)
    ones_vector <- matrix(1, nrow = n_hi, ncol = 1)

    kronecker_product <- kronecker(identity_block, ones_vector)

    diagonal_blocks[[i]] <- kronecker_product
  }

  # Combine all blocks along the diagonal to form the big matrix
  big_matrix <- do.call(Matrix::bdiag, diagonal_blocks)

  return(as.matrix(big_matrix))
}



# Example
data <- data.frame(
  patient = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  n_hi = c(2, 3, 8, 5, 6, 6, 3, 2, 3)
)

data

count_mat <- generate_record_count(data)
count_mat

Zhv <- generate_Zhv_matrix(data)
Zhv

N_h = sum(data[,2])

Z_h <- cbind(rep(1, N_h), Zhv)
Z_h









library(tidyverse)
# Parameters
H <- 3  # Number of hospitals
m <- 10   # Patients per hospital
n <- 5   # Visits per patient

beta0  <- 5
beta1 <- 5  # Fixed effect coefficients
sigma_u <- 3  # Hospital-level variance
sigma_v <- 2  # Patient-level variance
sigma_e <- 1  # Error variance


three_lvl_dat <- data.frame()

# DGP
# set.seed(123)
for (h in 1:H) {
  u_h <- rnorm(1, mean = 0, sd = sigma_u)  # Random hospital effect
  for (i in 1:m) {
    v_hi <- rnorm(1, mean = 0, sd = sigma_v)  # Random patient effect
    for (j in 1:n) {
      X_hij <- rnorm(1, mean = 0, sd = 1)
      epsilon_hij <- rnorm(1, mean = 0, sd = sigma_e)  # Random error
      y_hij <- beta0 + sum(X_hij * beta1) + u_h + v_hi + epsilon_hij  # Outcome

      three_lvl_dat <- rbind(three_lvl_dat, c(h, i, j, X_hij, y_hij))
    }
  }
}

colnames(three_lvl_dat) <- c("site", "patient", "visit", "X", "outcome")


site1_summary = three_lvl_dat %>%
  filter(site == "1") %>%
  group_by(patient) %>%
  summarise(n_hi = n())

Z_hv = generate_Zhv_matrix(site1_summary)

N_h = sum(site1_summary[,2])

Z_h <- cbind(rep(1, N_h), Z_hv)
Z_h








# Function to generate the record count matrix for a single hospital
generate_record_count <- function(data) {
  counts <- table(data[, "n_hi"])
  result_matrix <- cbind(as.numeric(names(counts)), as.numeric(counts))
  colnames(result_matrix) <- c("n_hi", "frequency")
  return(result_matrix)
}

# Function to generate Z_hv matrix for a single hospital
generate_Zhv_matrix <- function(data) {

  record_count_matrix <- generate_record_count(data)
  diagonal_blocks <- list()

  for (i in 1:nrow(record_count_matrix)) {
    n_hi <- record_count_matrix[i, "n_hi"]
    frequency <- record_count_matrix[i, "frequency"]

    identity_block <- diag(frequency)
    ones_vector <- matrix(1, nrow = n_hi, ncol = 1)

    kronecker_product <- kronecker(identity_block, ones_vector)
    diagonal_blocks[[i]] <- kronecker_product
  }

  big_matrix <- do.call(Matrix::bdiag, diagonal_blocks)
  return(as.matrix(big_matrix))
}


# Function to create the giant block diagonal matrix Z for the entire dataset
generate_Z_matrix <- function(data, H) {

  # List to store Z_h for each hospital
  Z_h_list <- list()

  for (h in 1:H) {
    site_data <- subset(data, site == h)

    N_h <- sum(site_data[, "n_hi"])

    Z_hv <- generate_Zhv_matrix(site_data)

    Z_h <- cbind(rep(1, N_h), Z_hv)

    Z_h_list[[h]] <- Z_h
  }

  Z <- do.call(Matrix::bdiag, Z_h_list)
  return(as.matrix(Z))
}

# Example Usage

# Simulated data for multiple hospitals
data <- data.frame(
  site = rep(1:3, each = 3),  # Hospitals 1, 2, 3
  patient = rep(1:3, times = 3),  # Patients
  n_hi = c(2, 3, 8, 5, 6, 6, 3, 2, 3)  # Number of visits per patient
)

# Number of hospitals
H <- length(unique(data$site))

# Generate Z matrix for the entire dataset
Z <- generate_Z_matrix(data, H)
print(Z)







