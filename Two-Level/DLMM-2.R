library(tidyverse)
library(data.table)
library(lme4)
library(nlme)
library(Matrix)
library(minqa)


lmm.get.summary2 <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL){
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  pz <- ncol(Z)

  ShXYZ <- list()
  for(hh in seq_along(id.site.uniq)){
    sh = id.site.uniq[hh]
    Xh <- X[id.site == sh, ]
    Zh <- Z[id.site == sh, ]
    Yh <- Y[id.site == sh]

    ShX  = t(Xh) %*% Xh
    ShXZ = t(Xh) %*% Zh
    ShXY = t(Xh) %*% Yh
    ShZ  = t(Zh) %*% Zh
    ShZY = t(Zh) %*% Yh
    ShY  = sum(Yh ^ 2)
    Nh <- sum(id.site == sh)
    ShXYZ[[sh]] <- list(ShX  = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ  = ShZ, ShZY = ShZY, ShY  = ShY, Nh = Nh)
  }

  return(ShXYZ)
}



lmm.profile02 <- function(V, s2, pooled = FALSE, reml = TRUE,
                          Y, X, Z, id.site, weights = NULL,
                          ShXYZ, rcpp = FALSE) {
  if (pooled) {
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z) / length(id.site.uniq)
  } else {
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    pz <- ncol(ShXYZ[[1]]$ShXZ)
  }

  lpterm1 <- lpterm2 <- remlterm <- 0
  bterm1 <- matrix(0, px, px)
  bterm2 <- rep(0, px)
  Vinv <- diag(1/diag(V))
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

    # Compute the log-determinant using the Cholesky factorization
    A <- diag(1, pz) + ShZ %*% V / s2
    logdet <- as.numeric(determinant(A, logarithm = TRUE)$modulus) + Nh * log(s2)
    lpterm1 <- lpterm1 + logdet

    # Compute Wh[[h]] using Cholesky decomposition of (Vinv + ShZ)
    L_Wh <- chol(s2 * Vinv + ShZ)
    Wh <- chol2inv(L_Wh)

    bterm1 <- bterm1 + (ShX - ShXZ %*% Wh %*% t(ShXZ)) / s2
    bterm2 <- bterm2 + (ShXY - ShXZ %*% Wh %*% ShZY) / s2
    lpterm2 <- lpterm2 + (ShY - t(ShZY) %*% Wh %*% ShZY) / s2
  }

  b = solve(bterm1, bterm2)
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b)

  if (reml) {
    remlterm <- determinant(bterm1, logarithm = TRUE)$modulus
    lp <- -(lpterm1 + qterm + remlterm) / 2
  } else {
    lp <- -(lpterm1 + (1 + log(qterm * 2 * pi / N)) * N) / 2
  }

  res <- list(lp = lp, b = b, s2 = s2)
  return(res)
}


# Function to fit 2-level DLMM
lmm.fit2 <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                     pooled = FALSE, reml = TRUE, common.s2 = TRUE,
                     ShXYZ = list(), corstr = 'independence',
                     mypar.init = NULL, hessian = FALSE, verbose = TRUE) {
  if (pooled) {
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z)
    K <- length(id.site.uniq)
    ShXYZ <- lmm.get.summary2(Y, X, Z, weights, id.site)
  } else {
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    pz <- ncol(ShXYZ[[1]]$ShXZ)
    K <- length(ShXYZ)
  }

  if (common.s2) {
    ns <- 1

    if(is.null(mypar.init)){
      if(corstr == 'independence'){
        mypar.init <- c(rep(0.5, pz), rep(0.5, ns))
      }else if(corstr == 'exchangeable'){
        mypar.init <- c(rep(0.5, pz), 0.1, rep(0.5, ns))
      }else if(corstr == 'unstructured'){
        mypar.init <- c(rep(0.5, pz), rep(0.1, pz * (pz - 1) / 2), rep(0.5, ns))
      }
      cat('default mypar.init (var comp) = ', mypar.init, '\n')
    }

    fn <- function(mypar){
      if(corstr == 'independence'){
        V <- diag(mypar[1 : pz], pz)
        s2 <- (mypar[-c(1 : pz)])
      }else if(corstr == 'exchangeable'){

      }
      return(-lmm.profile02(V, s2, pooled=F, reml, Y, X, Z, id.site, weights, ShXYZ)$lp)
    }


    res <- optim(mypar.init, fn, method = "L-BFGS-B",
                 hessian = T,
                 lower = rep(1e-6, length(mypar.init)))

    if (verbose) cat(ifelse(all(res$convergence == 0, eigen(res$hessian)$value > 0),
                            "Convergence Reached", "Non-convergence!"), 'and',
                     ifelse(all(eigen(res$hessian)$value > 0),
                            "Hessian PD", "Hessian not PD"), '\n',
                     "The number of function evaluations used is ", res$counts[1], '\n')

    mypar <- res$par
    if(corstr == 'independence'){
      V <- diag(mypar[1 : pz], pz)
      s2 <- mypar[- c(1 : pz)]
    }else if(corstr == 'exchangeable'){

    }

    res.profile <- lmm.profile02(V, s2, pooled = FALSE, reml, Y, X, Z, id.site, weights, ShXYZ)
  }

  return(list(b = res.profile$b, V = V, s2 = s2,
              res = res, res.profile = res.profile))
}
