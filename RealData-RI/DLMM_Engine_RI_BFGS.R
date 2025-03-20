library(tidyverse)
library(data.table)
library(lme4)
library(nlme)
library(Matrix)
library(minqa)

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
  if (is.null(weights)) weights <- rep(1, length(Y))
  X <- as.matrix(X)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)

  ShXYZ <- list()
  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    wth <- weights[id.site == sh]
    Xh <- X[id.site == sh, ]
    Yh <- Y[id.site == sh]
    Zh <- Z[h][[1]]

    ShX  <- t(Xh * wth) %*% Xh
    ShXZ <- t(Xh * wth) %*% Zh
    ShXY <- t(Xh * wth) %*% Yh
    ShZ  <- t(Zh * wth) %*% Zh
    ShZY <- t(Zh * wth) %*% Yh
    ShY  <- sum(Yh^2 * wth)
    Nh <- sum(id.site == sh)

    ShXYZ[[sh]] <- list(ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ = ShZ, ShZY = ShZY, ShY = ShY, Nh = Nh)
  }

  return(ShXYZ)
}

# Function to profile out the residual variance s2
lmm.profile03 <- function(par, pooled = FALSE, reml = TRUE,
                          Y, X, Z, id.site, weights = NULL,
                          ShXYZ, rcpp = FALSE) {
  if (pooled) {
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z) / length(id.site.uniq)
  } else {
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    pz <- length(id.site.uniq) + 1
  }

  lpterm1 <- lpterm2 <- remlterm <- 0
  bterm1 <- matrix(0, px, px)
  bterm2 <- rep(0, px)
  Wh <- list()
  N <- 0

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
    sigma_u2 <- par[1]
    sigma_vh2 <- par[1 + h]
    s2 <- tail(par, 1)
    V <- diag(c(sigma_u2, rep(sigma_vh2, (pzh - 1))), pzh)
    # V is diagonal
    Vinv <- diag(1/diag(V))

    # Compute the log-determinant using the Cholesky factorization
    A <- diag(1, pzh) + ShZ %*% V / s2
    logdet <- as.numeric(determinant(A, logarithm = TRUE)$modulus) + Nh * log(s2)
    lpterm1 <- lpterm1 + logdet

    # Compute Wh[[h]] using Cholesky decomposition of (Vinv + ShZ)
    L_Wh <- chol(s2 * Vinv + ShZ)
    Wh[[h]] <- chol2inv(L_Wh)

    bterm1 <- bterm1 + (ShX - ShXZ %*% Wh[[h]] %*% t(ShXZ)) / s2
    bterm2 <- bterm2 + (ShXY - ShXZ %*% Wh[[h]] %*% ShZY) / s2
    lpterm2 <- lpterm2 + (ShY - t(ShZY) %*% Wh[[h]] %*% ShZY) / s2
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

  res <- list(lp = lp, b = b, s2 = s2, s2p = s2p,
              allterms = list(lpterm1 = lpterm1, lpterm2 = lpterm2,
                              qterm = qterm, remlterm = remlterm,
                              bterm1 = bterm1, bterm2 = bterm2))
  return(res)
}


# Function to fit 3-level DLMM
lmm.fit3 <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                     pooled = FALSE, reml = TRUE, common.s2 = TRUE,
                     ShXYZ = list(), corstr = 'independence',
                     mypar.init = NULL, hessian = FALSE, verbose = TRUE) {
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

  if (common.s2) {
    ns <- 1

    if (is.null(mypar.init)) {
      mypar.init <- rep(1, pz + 1)
      if (verbose) cat('Default mypar.init (var comp) = ', mypar.init, '\n')
    }

    fn <- function(parameter) {
      return(-lmm.profile03(par = parameter, pooled = F, reml, Y, X, Z, id.site, weights, ShXYZ)$lp)
    }


    res <- optim(mypar.init, fn, method = "L-BFGS-B", lower = rep(1e-6, length(mypar.init)))
    if (verbose) cat(ifelse(res$convergence == 0, "Convergence Reached", "Non-convergence!"), '\n',
                     "The number of function evaluations used is ", res$counts[1], '\n')


    mypar <- res$par
    res.profile <- lmm.profile03(par = mypar, pooled = FALSE, reml, Y, X, Z, id.site, weights, ShXYZ)
    s2p <- res.profile$s2p
    s2 <- tail(mypar,1)
    V <- mypar[1:pz]
  }

  return(list(b = res.profile$b, V = V, s2 = s2,
              s2p = s2p,
              res = res, res.profile = res.profile))
}
