library(tidyverse)
library(data.table)
library(lme4)
library(nlme)
library(Matrix)
library(minqa)

lmm.profile.multi <- function(V, s2, rho, pooled = FALSE, reml = TRUE,
                              Y, X, Z, id.site, weights = NULL,
                              ShXYZ, m = 3, rcpp = FALSE) {
  if (pooled) {
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z) / length(id.site.uniq)
  } else {
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    pz <- dim(ShXYZ[[1]]$ShXZ)[2]  # Should be px + 1 (for intercept)
  }

  N <- 0
  lpterm1 <- lpterm2 <- remlterm <- 0
  bterm1 <- matrix(0, px*m, px*m)
  bterm2 <- matrix(0, px*m, 1)

  # Create exchangeable correlation matrix
  R <- matrix(rho, m, m)
  diag(R) <- 1
  R_inv <- solve(R)
  V_inv <- solve(V)

  for (h in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[h]
    ShX  <- ShXYZ[[sh]]$ShX      # [px × px]
    ShXZ <- ShXYZ[[sh]]$ShXZ     # [px × pz × m]
    ShXY <- ShXYZ[[sh]]$ShXY     # [px × m]
    ShZ  <- ShXYZ[[sh]]$ShZ      # [pz × pz]
    ShZY <- ShXYZ[[sh]]$ShZY     # [pz × m]
    ShY  <- ShXYZ[[sh]]$ShY      # [m × m]
    Nh   <- ShXYZ[[sh]]$Nh

    N <- N + Nh * m

    A_inv = kronecker(diag(Nh), R_inv)

    # # need to mordify below
    # # Compute the log-determinant using the Cholesky factorization
    # A <- diag(1, pz) + ShZ %*% V / s2
    # logdet <- as.numeric(determinant(A, logarithm = TRUE)$modulus) + Nh * log(s2)
    # lpterm1 <- lpterm1 + logdet
    #
    # # Compute Wh[[h]] using Cholesky decomposition of (Vinv + ShZ)
    # L_Wh <- chol(s2 * Vinv + ShZ)
    # Wh <- chol2inv(L_Wh)
    #
    # bterm1 <- bterm1 + (ShX - ShXZ %*% Wh %*% t(ShXZ)) / s2
    # bterm2 <- bterm2 + (ShXY - ShXZ %*% Wh %*% ShZY) / s2
    # lpterm2 <- lpterm2 + (ShY - t(ShZY) %*% Wh %*% ShZY) / s2
  }

  # Estimate fixed effects
  b <- solve(bterm1, bterm2)
  b_matrix <- matrix(b, nrow = px, ncol = m)

  # Compute profile likelihood
  qterm <- as.numeric(lpterm2 - 2 * t(b) %*% bterm2 + t(b) %*% bterm1 %*% b)

  if (reml) {
    remlterm <- determinant(bterm1, logarithm = TRUE)$modulus
    lp <- -(lpterm1 + qterm + remlterm)/2
  } else {
    lp <- -(lpterm1 + (1 + log(qterm * 2 * pi / N)) * N)/2
  }

  return(list(lp = lp, b = b_matrix, s2 = s2, rho = rho))
}



# Modified fitting function
lmm.fit.multi <- function(Y, X, Z, id.site, weights = NULL,
                          pooled = FALSE, reml = TRUE,
                          ShXYZ = list(), m = 1,
                          corstr = "exchangeable",
                          mypar.init = NULL, hessian = FALSE, verbose = TRUE) {

  if (pooled) {
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z)
    K <- length(id.site.uniq)
    ShXYZ <- lmm.get.summary.multi(Y, X, Z, id.site)
  } else {
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    pz <- dim(ShXYZ[[1]]$ShXZ)[2]
    K <- length(ShXYZ)
  }

  # Objective function
  fn <- function(mypar) {
    V <- matrix(0, pz, pz)
    diag(V) <- mypar[1:pz]

    if (corstr == "exchangeable") {
      rho <- mypar[pz + 1]
      s2 <- mypar[pz + 2]
    } else if (corstr == "independence") {
      rho <- 0
      s2 <- mypar[pz + 1]
    }

    rho <- max(min(rho, 0.99), -0.99) # range

    res <- lmm.profile.multi(V, s2, rho, pooled, reml,
                             Y, X, Z, id.site, weights, ShXYZ, m)
    return(-res$lp)
  }

  # Initial values
  if (is.null(mypar.init)) {
    if (corstr == "exchangeable") {
      mypar.init <- c(rep(1, pz), 0.3, 1)
    } else if (corstr == "independence") {
      mypar.init <- c(rep(1, pz), 1)
    }
  }

  # Optimization
  lower <- rep(1e-6, length(mypar.init))
  if (corstr == "exchangeable") {
    lower[pz + 1] <- -0.99
  }

  res <- optim(mypar.init, fn, method = "L-BFGS-B", hessian = hessian, lower = lower)

  # Extract results
  V <- matrix(0, pz, pz)
  diag(V) <- res$par[1:pz]

  if (corstr == "exchangeable") {
    rho <- res$par[pz + 1]
    s2 <- res$par[pz + 2]
  } else {
    rho <- 0
    s2 <- res$par[pz + 1]
  }

  res.profile <- lmm.profile.multi(V, s2, rho, pooled, reml,
                                   Y, X, Z, id.site, weights, ShXYZ, m)

  return(list(b = res.profile$b, V = V, s2 = s2, rho = rho,
              res = res, res.profile = res.profile))
}


# Summary statistics function
lmm.get.summary.multi <- function(Y, X, Z, id.site) {
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  pz <- ncol(Z)
  m <- ncol(Y)

  ShXYZ <- list()
  for(hh in seq_along(id.site.uniq)) {
    sh <- id.site.uniq[hh]
    Xh <- X[id.site == sh, ]
    Zh <- Z[id.site == sh, ]
    Yh <- Y[id.site == sh, ]  # nh × m matrix

    # Initialize summary statistics
    ShX <- t(Xh) %*% Xh       # px × px
    ShZ <- t(Zh) %*% Zh       # pz × pz
    ShXZ <- array(0, dim = c(px, pz, m))  # px × pz × m
    ShXY <- matrix(0, px, m)  # px × m
    ShZY <- matrix(0, pz, m)  # pz × m
    ShY <- matrix(0, m, m)    # m × m (sum of Y'Y across all patients)

    # Compute cross-products correctly for multivariate case
    for(j in 1:m) {
      ShXZ[,,j] <- t(Xh) %*% Zh          # X'Z for outcome j
      ShXY[,j] <- t(Xh) %*% Yh[,j]       # X'y for outcome j
      ShZY[,j] <- t(Zh) %*% Yh[,j]       # Z'y for outcome j
    }

    # Compute Y'Y correctly (sum of outer products)
    ShY <- t(Yh) %*% Yh                  # m × m

    nh <- sum(id.site == sh)
    ShXYZ[[sh]] <- list(ShX = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ = ShZ, ShZY = ShZY, ShY = ShY, Nh = nh)
  }
  return(ShXYZ)
}


# lmm.get.summary.multi <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL){
#   X <- as.matrix(X)
#   Z <- as.matrix(Z)
#   id.site <- as.character(id.site)
#   id.site.uniq <- unique(id.site)
#   px <- ncol(X)
#   pz <- ncol(Z)
#
#   ShXYZ <- list()
#   for(hh in seq_along(id.site.uniq)){
#     sh = id.site.uniq[hh]
#     Xh <- X[id.site == sh, ]
#     Zh <- Z[id.site == sh, ]
#     Yh <- Y[id.site == sh]
#
#     ShX  = t(Xh) %*% Xh
#     ShXZ = t(Xh) %*% Zh
#     ShXY = t(Xh) %*% Yh
#     ShZ  = t(Zh) %*% Zh
#     ShZY = t(Zh) %*% Yh
#     ShY  = sum(Yh ^ 2)
#     Nh <- sum(id.site == sh)
#     ShXYZ[[sh]] <- list(ShX  = ShX, ShXZ = ShXZ, ShXY = ShXY,
#                         ShZ  = ShZ, ShZY = ShZY, ShY  = ShY, Nh = Nh)
#   }
#
#   return(ShXYZ)
# }
