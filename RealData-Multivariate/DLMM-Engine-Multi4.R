
library(tidyverse)
library(data.table)
library(Matrix)

create_block_matrix <- function(Sigma, M, rho) {
  N <- nrow(Sigma)
  block_matrix <- kronecker(matrix(rho, M, M), Sigma) - kronecker(diag(M), rho*Sigma) + kronecker(diag(M), Sigma)
  return(block_matrix)
}

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
  big_matrix <- cbind(1, big_matrix)
  return(as.matrix(big_matrix))
}

# Function to get summary stats from each site for distributed LMM
lmm.get.summary3 <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(Y))
  X <- as.matrix(X)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  py <- ncol(Y)

  ShXYZ <- list()
  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    # wth <- weights[id.site == sh]
    Xh <- X[id.site == sh, ]
    Yh <- matrix(Y[id.site == sh, ], ncol = 1)
    Zh <- Z[h][[1]]

    # Xh <- kronecker(matrix(1, py, 1), Xh)
    # Zh <- kronecker(matrix(1, py, 1), Zh)
    Xh <- kronecker(diag(py), Xh)
    Zh <- kronecker(diag(py), Zh)

    ShX  <- t(Xh) %*% Xh
    ShXZ <- t(Xh) %*% Zh
    ShXY <- t(Xh) %*% Yh
    ShZ  <- t(Zh) %*% Zh
    ShZY <- t(Zh) %*% Yh
    ShY  <- sum(Yh^2)
    Nh <- sum(id.site == sh)


    ShXYZ[[sh]] <- list(ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ = ShZ, ShZY = ShZY, ShY = ShY, Nh = Nh)
  }

  return(ShXYZ)
}



# Function to fit 3-level DLMM
lmm.profile03_multivariate <- function(par, pooled = FALSE, reml = TRUE,
                                                        Y, X, Z, id.site, weights = NULL,
                                                        ShXYZ, M, rcpp = FALSE) {
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  pz <- ncol(Z) / length(id.site.uniq)

  lpterm1 <- lpterm2 <- remlterm <- 0
  bterm1 <- matrix(0, px, px)
  bterm2 <- rep(0, px)
  Wh <- list()
  N <- 0

  # Extract variance parameters
  sigma_u2 <- par[1]  # Random effect variance
  sigma_vh2 <- par[2:(length(id.site.uniq) + 1)]  # Variance per site
  sigma2 <- tail(par, 1)  # Residual variance
  rho <- par[length(par)]  # Exchangeable correlation

  # Construct the exchangeable residual covariance structure
  Sigma_residual <- sigma2 * ((1 - rho) * diag(M) + rho * matrix(1, M, M))

  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    ShX  <- ShXYZ[[sh]]$ShX
    ShXZ <- ShXYZ[[sh]]$ShXZ
    ShXY <- ShXYZ[[sh]]$ShXY
    ShZ  <- ShXYZ[[sh]]$ShZ
    ShZY <- ShXYZ[[sh]]$ShZY
    ShY  <- ShXYZ[[sh]]$ShY
    Nh   <- ShXYZ[[sh]]$Nh

    N <- N + Nh
    pzh <- ncol(ShZ)

    # Create the variance structure for the group
    V <- diag(c(sigma_u2, rep(sigma_vh2[h], (pzh - 1))), pzh)
    Vinv <- diag(1/diag(V))

    # Modify A matrix to include Sigma_residual
    A <- diag(1, pzh * M) + kronecker(Sigma_residual, ShZ %*% V) / sigma2

    logdet <- as.numeric(determinant(A, logarithm = TRUE)$modulus) + Nh * log(sigma2)
    lpterm1 <- lpterm1 + logdet

    L_Wh <- chol(sigma2 * Vinv + ShZ)
    Wh[[h]] <- chol2inv(L_Wh)

    bterm1 <- bterm1 + (ShX - ShXZ %*% Wh[[h]] %*% t(ShXZ)) / sigma2
    bterm2 <- bterm2 + (ShXY - ShXZ %*% Wh[[h]] %*% ShZY) / sigma2
    lpterm2 <- lpterm2 + (ShY - t(ShZY) %*% Wh[[h]] %*% ShZY) / sigma2
  }

  b = solve(bterm1, bterm2)
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b)

  if (reml) {
    remlterm <- determinant(bterm1, logarithm = TRUE)$modulus
    lp <- -(lpterm1 + qterm + remlterm) / 2
  } else {
    lp <- -(lpterm1 + (1 + log(qterm * 2 * pi / N)) * N) / 2
  }

  s2p = qterm / (N)

  res <- list(lp = lp, b = b, sigma2 = sigma2, s2p = s2p,
              allterms = list(lpterm1 = lpterm1, lpterm2 = lpterm2,
                              qterm = qterm, remlterm = remlterm,
                              bterm1 = bterm1, bterm2 = bterm2))
  return(res)
}



lmm.fit3_multivariate <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                                  pooled = FALSE, reml = TRUE, common.s2 = TRUE,
                                  ShXYZ = list(), corstr = 'exchangeable',
                                  mypar.init = NULL, hessian = FALSE, verbose = TRUE, M) {
  if (pooled) {
    id.site.uniq <- unique(id.site)
    K <- length(id.site.uniq)
    px <- ncol(X)
    pz <- ncol(Z) / K
    ShXYZ <- lmm.get.summary3(Y, X, Z, weights = weights, id.site = id.site)
  } else {
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    K <- length(ShXYZ)
    pz <- K + 1
  }

  if (is.null(mypar.init)) {
    mypar.init <- rep(1, pz + 2)  # Add rho
    if (verbose) cat('Default mypar.init (var comp) = ', mypar.init, '\n')
  }

  fn <- function(parameter) {
    return(-lmm.profile03_multivariate(par = parameter, pooled = FALSE, reml, Y, X, Z, id.site, weights, ShXYZ, M)$lp)
  }

  res <- optim(mypar.init, fn, method = "L-BFGS-B",
               hessian = T,
               lower = rep(1e-6, length(mypar.init)))

  mypar <- res$par
  res.profile <- lmm.profile03_multivariate(par = mypar, pooled = FALSE, reml, Y, X, Z, id.site, weights, ShXYZ, M)

  return(list(b = res.profile$b, V = mypar[1:pz], s2 = mypar[pz + 1],
              rho = mypar[pz + 2], res = res, res.profile = res.profile))
}

