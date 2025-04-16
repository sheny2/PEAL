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
    Xh <- kronecker(diag(py), Xh) # larger format, dimension times py
    Zh <- kronecker(diag(py), Zh)

    ShX  <- crossprod(Xh, Xh)
    ShXZ <- crossprod(Xh, Zh)
    ShXY <- crossprod(Xh, Yh)
    ShZ  <- crossprod(Zh, Zh)
    ShZY <- crossprod(Zh, Yh)
    ShY  <- sum(Yh^2)
    Nh <- sum(id.site == sh)


    ShXYZ[[sh]] <- list(ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ = ShZ, ShZY = ShZY, ShY = ShY, Nh = Nh)
  }

  return(ShXYZ)
}


lmm.profile03Multi <- function(par, pooled = FALSE, reml = TRUE,
                          Y, X, Z, id.site, weights = NULL,
                          ShXYZ, corstr, rcpp = FALSE, py = NULL) {
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

    if(corstr == 'independence'){
      sigma_u2 = par[1]
      sigma_vh2 = par[1 + h]
      s2 = tail(par, 1)

      a_h = pzh/py
      V <- diag(c(sigma_u2, rep(sigma_vh2, (a_h - 1))), a_h)
      V <- create_block_matrix(Sigma = V, M = py, rho = 0)
    }

    else if(corstr == 'exchangeable'){
      sigma_u2 = par[1]
      sigma_vh2 = par[1 + h]
      s2 = tail(par, 3)[1]
      rho_v = tail(par, 3)[2]
      rho = tail(par, 3)[3]

      a_h = pzh/py
      V <- diag(c(sigma_u2, rep(sigma_vh2, (a_h - 1))), a_h)
      V <- create_block_matrix(Sigma = V, M = py, rho = 0)
      V <- replace(V, V == 0, rho)

      # D <- diag(sqrt(sigma_vh2), pzh)
      # a_h = pzh/py
      # D[1,1] = sqrt(sigma_u2)
      # D[1+a_h,1+a_h] = sqrt(sigma_u2)
      # D[1+a_h*2,1+a_h*2] = sqrt(sigma_u2)
      # # a_h = pzh/py
      # # D <- diag(c(sigma_u2, rep(sigma_vh2, (a_h - 1))), a_h)
      # # D <- create_block_matrix(Sigma = D, M = py, rho = 0)
      # # D <- sqrt(D)
      # V <- D %*% (matrix(rho, pzh, pzh) + diag(1 - rho, pzh)) %*% D
    }

    Vinv <- solve(V)

    # Compute the log-determinant
    A <- diag(1, ncol(V)) + ShZ %*% V / s2
    logdet <- as.numeric(determinant(A, logarithm = TRUE)$modulus) + py * Nh * log(s2)
    lpterm1 <- lpterm1 + logdet

    # Compute Wh[[h]] using Cholesky decomposition of (Vinv + ShZ)
    # L_Wh <- chol(s2 * Vinv + ShZ)
    # Wh[[h]] <- chol2inv(L_Wh)
    Wh[[h]] <- solve(s2 * Vinv + ShZ)

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

  res <- list(lp = lp, b = b, s2 = s2,
              allterms = list(lpterm1 = lpterm1, lpterm2 = lpterm2,
                              qterm = qterm, remlterm = remlterm,
                              bterm1 = bterm1, bterm2 = bterm2))
  return(res)
}


# Function to fit 3-level DLMM
lmm.fit3Multi <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                     pooled = FALSE, reml = TRUE, common.s2 = TRUE,
                     ShXYZ = list(), corstr = 'independence',
                     mypar.init = NULL, hessian = FALSE, verbose = TRUE, py = NULL) {
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

    if (is.null(mypar.init)) {
      if(corstr == 'independence'){
        mypar.init <- c(rep(1, pz+1))
        if (verbose) cat('Default mypar.init (var comp) = ', mypar.init, '\n')
        lower_bounds <- c(rep(1e-6, pz+1))
        upper_bounds <- c(rep(Inf, pz+1))
      }else if(corstr == 'exchangeable'){
        mypar.init <- c(rep(1, pz+1), 0.1, 0.1) # one for RE one for error
        if (verbose) cat('Default mypar.init (var comp + rho) = ', mypar.init, '\n')
        lower_bounds <- c(rep(1e-6, pz+1), rep(-1+1e-6,2))
        upper_bounds <- c(rep(Inf, pz+1), rep(1-1e-6,2))
      }
    }

    fn <- function(parameter) {
      return(-lmm.profile03Multi(par = parameter, pooled = FALSE,
                                 reml, Y, X, Z, id.site, weights, ShXYZ, corstr, py = py)$lp)
    }


    res <- optim(mypar.init, fn, method = "L-BFGS-B",
                 hessian = T,
                 # control = list(maxit = 500),
                 lower = lower_bounds, upper = upper_bounds
                 # lower = rep(1e-6, length(mypar.init))
    )

    if (verbose) cat(ifelse(all(res$convergence == 0, eigen(res$hessian)$value > 0),
                            "Convergence Reached", "Non-convergence!"), 'and',
                     ifelse(all(eigen(res$hessian)$value > 0),
                            "Hessian PD", "Hessian not PD"), '\n',
                     "The number of function evaluations used is ", res$counts[1], '\n')


    mypar <- res$par
    res.profile <- lmm.profile03Multi(par = mypar, pooled = FALSE,
                                      reml, Y, X, Z, id.site, weights, ShXYZ, corstr, py = py)

    s2 <- mypar[pz+1]
    V <- mypar[1:pz]
    rho_v = ifelse(corstr == 'exchangeable', tail(mypar, 2), NA)
    rho = ifelse(corstr == 'exchangeable', tail(mypar, 1), NA)
  }

  return(list(b = res.profile$b, V = V, s2 = s2,
              rho = rho, rho_v = rho_v, # rho_v and rho
              res = res, res.profile = res.profile))
}
