
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
  bterm1 <- matrix(0, px, px)   # bterm1 is still px x px
  bterm2 <- matrix(0, px, py)   # bterm2 is now px x py
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
      V <- diag(c(sigma_u2, rep(sigma_vh2, (pzh - 1))), pzh)
      s2 = tail(par, 1)
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
    # log-determinant
    A <- diag(1, pzh) + ShZ %*% V / s2
    logdet <- as.numeric(determinant(A, logarithm = TRUE)$modulus) + Nh * log(s2)
    lpterm1 <- lpterm1 + logdet

    Wh[[h]] = solve(s2 * Vinv + ShZ, diag(nrow(ShZ)))
    # L_Wh <- chol(s2 * Vinv + ShZ)
    # Wh[[h]] <- chol2inv(L_Wh)

    bterm1 <- bterm1 + (ShX - ShXZ %*% Wh[[h]] %*% t(ShXZ)) / s2
    bterm2 <- bterm2 + (ShXY - ShXZ %*% Wh[[h]] %*% ShZY) / s2     # bterm2 (px x py) now

    # lpterm2: must take trace of (ShY - t(ShZY) %*% Wh[[h]] %*% ShZY)
    # that expression is py by py, so we do sum(diag(...)).
    M <- (ShY - t(ShZY) %*% Wh[[h]] %*% ShZY) / s2
    lpterm2 <- lpterm2 + sum(diag(M))
  }

  b <- solve(bterm1, bterm2)

  # qterm is the final sum-of-squares piece:
  #  lpterm2 - 2 * sum(bterm2 * b) + trace(t(b) %*% bterm1 %*% b).
  # sum(bterm2 * b) does elementwise multiplication => scalar
  # t(b) %*% bterm1 %*% b => (py x py), so we take sum(diag(...)).
  tb_bterm1_b <- t(b) %*% bterm1 %*% b  # py by py
  qterm <- lpterm2 - 2 * sum(bterm2 * b) + sum(diag(tb_bterm1_b))

  if (reml) {
    remlterm <- as.numeric(determinant(bterm1, logarithm = TRUE)$modulus)
    lp <- -(lpterm1 * py + qterm + remlterm * py) / 2
  } else {
    lp <- -(lpterm1 + (1 + log(qterm * 2 * pi / N)) * N) / 2
  }

  res <- list(
    lp = lp,
    b = b,       # (px x py)
    s2 = s2,
    allterms = list(
      lpterm1 = lpterm1, lpterm2 = lpterm2,
      qterm   = qterm, remlterm= remlterm,
      bterm1  = bterm1, bterm2  = bterm2
    )
  )
  return(res)
}


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

    fn <- function(parameter) {
      return(-lmm.profile03(par = parameter, pooled = F, reml, Y, X, Z, id.site, weights, ShXYZ, corstr)$lp)
    }


    if (is.null(mypar.init)) {
      if(corstr == 'independence'){
        mypar.init <- c(rep(1, pz+1))
        if (verbose) cat('Default mypar.init (var comp) = ', mypar.init, '\n')
        lower_bounds <- c(rep(1e-6, pz+1))
        upper_bounds <- c(rep(Inf, pz+1))
      }else if(corstr == 'exchangeable'){
        mypar.init <- c(rep(1, pz+1), 0.1)
        if (verbose) cat('Default mypar.init (var comp + rho) = ', mypar.init, '\n')
        lower_bounds <- c(rep(1e-6, pz+1), -1+1e-6)
        upper_bounds <- c(rep(Inf, pz+1), 1-1e-6)
      }
    }

    res <- optim(mypar.init, fn, method = "L-BFGS-B",
                 lower = lower_bounds, upper = upper_bounds)
    if (verbose) cat(ifelse(res$convergence == 0,
                            "Convergence Reached", "Non-convergence!"), '\n',
                     "The number of function evaluations used is ", res$counts[1], '\n')


    mypar <- res$par
    res.profile <- lmm.profile03(par = mypar, pooled = FALSE, reml, Y, X, Z, id.site, weights, ShXYZ, corstr)
    s2 <- mypar[pz+1]
    V <- mypar[1:pz]
    rho = tail(mypar, 1)
  }

    return(list(b = res.profile$b, V = V, s2 = s2,
                rho = rho,
                res = res, res.profile = res.profile))
}

