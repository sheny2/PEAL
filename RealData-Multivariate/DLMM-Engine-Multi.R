library(tidyverse)
library(data.table)
library(Matrix)

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

  ShXYZ <- list()
  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    wth <- weights[id.site == sh]
    Xh <- X[id.site == sh, ]
    Yh <- Y[id.site == sh, ]
    Zh <- Z[h][[1]]

    ShX  <- t(Xh) %*% Xh
    ShXZ <- t(Xh * wth) %*% Zh
    ShXY <- t(Xh * wth) %*% Yh
    ShZ  <- t(Zh * wth) %*% Zh
    ShZY <- t(Zh * wth) %*% Yh
    ShY  <- t(Yh) %*% (Yh * wth)  # ShY is now a matrix of size (py by py)
    Nh <- sum(id.site == sh)

    ShXYZ[[sh]] <- list(ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ = ShZ, ShZY = ShZY, ShY = ShY, Nh = Nh)
  }

  return(ShXYZ)
}


# Function to profile out the residual variance s2
lmm.profile03 <- function(par, pooled = FALSE, reml = TRUE,
                          Y, X, Z, id.site, weights = NULL,
                          ShXYZ, corstr, rcpp = FALSE) {
  if (pooled) {
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z) / length(id.site.uniq)
  } else {
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    pz <- length(id.site.uniq) + 1
    py = ncol(ShXYZ[[1]]$ShY)
  }

  lpterm1 <- lpterm2 <- remlterm <- 0
  bterm1 <- matrix(0, px, px)
  # bterm2 <- rep(0, px)
  bterm2 <- matrix(0, px, py)              # bterm2 is now a matrix
  lpterm2 <- matrix(0, py, py)
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
      V <- diag(c(sigma_u2, rep(sigma_vh2, (pzh - 1))), pzh)
    }else if(corstr == 'exchangeable'){
      sigma_u2 = par[1]
      sigma_vh2 = par[1 + h]
      s2 = tail(par, 2)[1]
      rho = tail(par, 2)[2]
      D <- diag(sqrt(sigma_vh2), pzh)
      D[1,1] = sqrt(sigma_u2)
      V <- D %*% (matrix(rho, pzh, pzh) + diag(1 - rho, pzh)) %*% D
    }

    Vinv <- solve(V, diag(nrow(V)))  # Inverse of V

    # Compute the log-determinant with numerical stability
    A <- diag(1, pzh) + ShZ %*% V / s2
    logdet <- as.numeric(determinant(A, logarithm = TRUE)$modulus) + Nh * log(s2)
    lpterm1 <- lpterm1 + logdet


    # Compute Wh[[h]] using Cholesky decomposition of (Vinv + ShZ)
    L_Wh <- chol(s2 * Vinv + ShZ)
    Wh[[h]] <- chol2inv(L_Wh)
    # Wh[[h]] <- solve(Vinv + ShZ, diag(nrow(Vinv)))

    bterm1 <- bterm1 + (ShX - ShXZ %*% Wh[[h]] %*% t(ShXZ)) / s2
    # bterm2 <- bterm2 + (ShXY - ShXZ %*% Wh[[h]] %*% ShZY)

    for (k in 1:py) {
      bterm2[, k] <- bterm2[, k] + (ShXY[, k] - ShXZ %*% Wh[[h]] %*% ShZY[, k]) / s2
    }

    lpterm2 <- lpterm2 + (ShY - t(ShZY) %*% Wh[[h]] %*% ShZY) / s2
  }

  # L <- chol(bterm1)
  # b <- backsolve(L, forwardsolve(t(L), bterm2))
  b <- solve(bterm1, bterm2)

  # qterm <- sum(diag(lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b))

  qterm <- 0
  for (k in 1:py) {
    qterm <- qterm + as.numeric(lpterm2[k, k] - 2 * sum(bterm2[, k] * b[, k]) + t(b[, k]) %*% bterm1 %*% b[, k])
  }


  if (reml) {
    remlterm <- determinant(bterm1, logarithm = TRUE)$modulus
    s2 <- qterm / (N*py - px* py)
    lp <- -(lpterm1 + qterm + remlterm) / 2
  } else {
    s2 <- qterm / N
    lp <- -(lpterm1 + (1 + log(qterm * 2 * pi / N)) * N) / 2
  }

  res <- list(lp = lp, b = b, s2 = s2,
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
    H <- length(ShXYZ)
    pz <- H + 1
  }

  if (common.s2) {

    if (is.null(mypar.init)) {
      if(corstr == 'independence'){
        mypar.init <- c(rep(1, pz+1))
      }else if(corstr == 'exchangeable'){
        mypar.init <- c(rep(1, pz+1), 0.1)
      }
      if (verbose) cat('Default mypar.init (var comp) = ', mypar.init, '\n')
    }

    fn <- function(parameter) {
      return(-lmm.profile03(par = parameter, pooled = F, reml, Y, X, Z, id.site, weights, ShXYZ, corstr)$lp)
    }


    res <- optim(mypar.init, fn, method = "L-BFGS-B", lower = rep(1e-6, length(mypar.init)))
    if (verbose) cat(ifelse(res$convergence == 0, "Convergence Reached", "Non-convergence!"), '\n',
                     "The number of function evaluations used is ", res$counts[1], '\n')


    mypar <- res$par
    res.profile <- lmm.profile03(par = mypar, pooled = FALSE, reml, Y, X, Z, id.site, weights, ShXYZ, corstr)
    s2 <- tail(mypar,1)
    V <- mypar[1:pz]
  }

  return(list(b = res.profile$b, V = V, s2 = s2,
              res = res, res.profile = res.profile))
}

