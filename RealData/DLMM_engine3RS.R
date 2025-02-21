library(tidyverse)
require(data.table)
require(lme4)
require(nlme)


############ Preprocessing
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






## get summary stats from each site for distributed lmm
lmm.get.summary3 <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL, m_h_all = NULL){
    if(is.null(weights)) weights <- rep(1, length(Y))
    X <- as.matrix(X)
    id.site <- as.character(id.site)
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    # Z <- as.matrix(Z)
    # pz <- ncol(Z)

    ShXYZ <- list()
    for(h in seq_along(id.site.uniq)){
      sh = id.site.uniq[h]
      wth = weights[id.site == sh]
      Xh <- X[id.site == sh, ]
      Yh <- Y[id.site == sh]

      Zh <- Z[h][[1]]
      # Zh <- Z[id.site == sh, ]
      # non_zero_columns <- colSums(Zh != 0) > 0
      # Zh <- Zh[, non_zero_columns, drop = FALSE]

      ShX  = t(Xh*wth) %*% Xh
      ShXZ = t(Xh*wth) %*% Zh
      ShXY = t(Xh*wth) %*% Yh
      ShZ  = t(Zh*wth) %*% Zh
      ShZY = t(Zh*wth) %*% Yh
      ShY  = sum(Yh ^ 2 *wth)
      Nh <- sum(id.site == sh)
      mh = m_h_all[h,]

      ShXYZ[[sh]] <- list(ShX  = ShX, ShXZ = ShXZ, ShXY = ShXY,
                          ShZ  = ShZ, ShZY = ShZY, ShY  = ShY, Nh = Nh, mh = mh)
    }

    return(ShXYZ)
  }





############
## further profile out the residual var s2, used if common.s2=T (deemed as usual situation)
lmm.profile03 <- function(par,                # var components para, = original V / s2    # s2,
                          pooled = F, reml = T,
                          Y, X, Z, id.site, weights=NULL,  # if pooled ipd is available
                          ShXYZ,             # if pooled ipd is not available, use summary stats
                          rcpp = F){
  if(pooled == T){
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z) /length(id.site.uniq)
  }else{
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    K <- length(ShXYZ)
    pz <- K + 1 + K
  }

  if(rcpp == F){
    lpterm1 = lpterm2 = remlterm = 0
    bterm1 <- matrix(0, px, px)  # sum_i Xi' \Sigma_i^-1 Xi
    bterm2 <- rep(0, px)         # sum_i Xi' \Sigma_i^-1 Yi
    # Vinv <- solve(V)
    Wh <- uh <- varuh <- varuh_post <- list()  # save this for each subject for BLUP

    N <- 0
    for(h in seq_along(id.site.uniq)){
      sh <- id.site.uniq[h]

      ShX  <- ShXYZ[[sh]]$ShX
      ShXZ <- ShXYZ[[sh]]$ShXZ
      ShXY <- ShXYZ[[sh]]$ShXY
      ShZ <- ShXYZ[[sh]]$ShZ
      ShZY <- ShXYZ[[sh]]$ShZY
      ShY  <- ShXYZ[[sh]]$ShY
      Nh <- ShXYZ[[sh]]$Nh
      mh <- ShXYZ[[sh]]$mh

      N <- N + Nh

      pzh = ncol(ShZ)
      par_h = c(par[1], par[1+h], par[1+K+h])
      V <- diag(c(par_h[1],
                  rep(par_h[2], mh),
                  rep(par_h[3], ((px-1)*mh))
                  ), pzh)

      # V <- diag(c(par_h[1], rep(par_h[2], (pzh-1))), pzh) # V, Theta differ by h
      Vinv <- solve(V)
      # pzh = ncol(ShZ)
      # V <- diag(c(par[1], rep(par[2], (pzh-1))), pzh) # V
      # Vinv <- solve(V)

      logdet <- log(det(diag(1, pzh) + ShZ %*% V))    # -log(det(diag(wti))) omitted as wti are the fixed weights
      lpterm1 <- lpterm1 + logdet

      Wh[[h]] <- solve(Vinv + ShZ)
      bterm1 <- bterm1 + (ShX - ShXZ %*% Wh[[h]] %*% t(ShXZ)) # / s2i
      bterm2 <- bterm2 + (ShXY - ShXZ %*% Wh[[h]] %*% ShZY) # / s2i
      lpterm2 <- lpterm2 + (ShY - t(ShZY) %*% Wh[[h]] %*% ShZY) #/ s2i
    }

    b <- solve(bterm1, bterm2)
    ## quadratic term: = (Yi-Xi*b)'\Sigma_i^-1 (Yi-Xi*b)
    qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b )
    if(reml == T){
      remlterm = log(det(bterm1))
      s2 <- qterm / (N-px)
      lp <- - (lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b + remlterm) / 2
      # lp <- - (lpterm1 + (1 + log(qterm * 2 * pi / (N - px))) * (N - px)) / 2     # Bates2015JSS eq(41)
    }else{
      s2 <- qterm / N
      lp <- - (lpterm1 + (1 + log(qterm * 2 * pi / N)) * N) / 2               # Bates2015JSS eq(34)
    }

    lk <- -(lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2 # To be used in calculating mAIC
    # lk <- - (lpterm1 + (1+log(qterm*2*pi/N))*N) / 2

    # Loop to estimate uh
    # uh = V Z_i^{T} Sigma_{i}^{-1} (Y_i - X_i \beta)
    # for(h in seq_along(id.site.uniq)){ # re-call this loop
    #   sh <- id.site.uniq[h]
    #
    #   ShX  <- ShXYZ[[sh]]$ShX
    #   ShXZ <- ShXYZ[[sh]]$ShXZ
    #   ShXY <- ShXYZ[[sh]]$ShXY
    #   ShZ <- ShXYZ[[sh]]$ShZ
    #   ShZY <- ShXYZ[[sh]]$ShZY
    #   ShY  <- ShXYZ[[sh]]$ShY
    #   Nh <- ShXYZ[[sh]]$Nh
    #
    #   pzh = ncol(ShZ)
    #   par_h = c(par[1], par[1+h])
    #   V <- diag(c(par_h[1], rep(par_h[2], (pzh-1))), pzh) # V, Theta differ by h
    #   Vinv <- solve(V)
    #
    #   Wh[[h]] <- solve(Vinv + ShZ)
    #
    #   s2h <- s2   # [ii]  # sigma_i^2
    #   uhterm1 <- (ShZY - ShZ %*% Wh[[h]] %*% ShZY) / s2h
    #   uhterm2 <-  ((t(ShXZ) - ShZ %*% Wh[[h]] %*% t(ShXZ)) / s2h) %*% as.matrix(b)
    #   uh[[h]] <- V %*% as.numeric(uhterm1 - uhterm2) * s2h #
    #
    #   vterm1 <- V %*% ((ShZ - ShZ %*% Wh[[h]] %*% ShZ) / s2h) %*% V * (s2h ^ 2)    #
    #   vterm2 <- V %*% ((t(ShXZ) - ShZ %*% Wh[[h]] %*% t(ShXZ)) / s2h)  * s2h   #
    #   varuh[[h]] <- (V * s2h) - vterm1 + (vterm2 %*% solve(bterm1, t(vterm2)))     # 20210111
    #   varuh_post[[h]] <- vterm1 - (vterm2 %*% solve(bterm1, t(vterm2)))            # 20210111
    # }

    res <- list(lp = lp, b = b,
                s2 = s2,
                # uh = uh,
                # lk = lk,
                # varuh = varuh,
                # varuh_post = varuh_post,
                allterms = list(lpterm1 = lpterm1,
                                lpterm2 = lpterm2,
                                qterm = qterm,
                                remlterm = remlterm,
                                # remlterm.i = remlterm.i,
                                bterm1 = bterm1,
                                bterm2 = bterm2)
                )
  }
  return(res)
}




# fit 3 level DLMM
lmm.fit3 <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                     pooled = F, reml = T,
                     common.s2 = T,      # common reShdual var across sites
                     ShXYZ = list(),
                     corstr = 'independence',
                     mypar.init = NULL,
                     hessian = F,
                     verbose){
  if(pooled == T){
    id.site.uniq <- unique(id.site)
    K <- length(id.site.uniq)
    px <- ncol(X)
    pz <- ncol(Z) / K
    ShXYZ <- lmm.get.summary3(Y, X, Z, weights = weights, id.site = id.site)
  }else{
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    K <- length(ShXYZ)
    pz <- 1 + K + K
  }

  ## profile out s2
  if(common.s2 == T){
    ns <- 1 # number of s2 para

    if(is.null(mypar.init)){
      # mypar.init <- c(1, 1)
      mypar.init <- rep(1, pz)

      cat('default mypar.init (var comp) = ', mypar.init, '\n')
    }


    fn <- function(parameter){
      # V <- diag(c(mypar[1], rep(mypar[2], (pz-1))), pz) # Note: V is really Theta here.
      # V <- diag(mypar[1 : pz], pz)

      return(-lmm.profile03(par = parameter, pooled=F, reml, Y, X, Z, id.site, weights, ShXYZ = ShXYZ)$lp)
    }

    # res <- optim(mypar.init, fn, hessian = F, method = c("Nelder-Mead"))
    res <- minqa::bobyqa(mypar.init, fn, lower=rep(1e-6, length(mypar.init)), control=list(maxfun=1e5))

    cat(res$msg, " The number of function evaluations used is ", res$feval)

    # Check for convergence
    if (res$ierr != 0) {
      warning(paste0(res$msg))
    }

    mypar <- res$par
    # V <- diag(c(mypar[1], rep(mypar[2], (pz-1))), pz)
    # V <- diag(mypar[1 : pz], pz)

    res.profile <- lmm.profile03(par = mypar, pooled=F, reml, Y, X, Z, id.site, weights, ShXYZ)
    s2 <- res.profile$s2
    V <- mypar * s2             # scale back

  }


  # ## Inference (Wald test statistic)
  # vd <- diag(solve(res.profile$allterms$bterm1))
  # if(common.s2==T) vd <- diag(solve(res.profile$allterms$bterm1 / s2))  # scale back
  # wald <- res.profile$b / sqrt(vd)
  #
  # ## 95% CI for fixed effects
  # lb <- res.profile$b -  1.96 * sqrt(vd)
  # ub <- res.profile$b +  1.96 * sqrt(vd)
  #
  # ## Marginal AIC: OKAY to use if the main interest is to model fixed population effects with a reasonable correlation structur (Kneib & Greven (2010))
  # mAIC <- 2 * res.profile$lk + 2 * (px + (length(mypar) - ns))
  #
  #
  # ## Prediction
  # uhhat <- as.matrix(do.call(rbind, lapply(seq_len(K), function(kk) {
  #
  #   ll <- length(which(id.site == id.site.uniq[kk]))
  #   dm <- matrix(1, nrow = ll, ncol = pz)
  #
  #   sweep(dm, MARGIN = 2, res.profile$uh[[kk]], `*`)
  #
  # })))
  #
  # uhhat_subj <- as.matrix(do.call(rbind, lapply(seq_len(K), function(kk) {
  #   t(res.profile$uh[[kk]])
  # })))
  #
  # if(pooled == T){
  #   Yhat <- X %*% as.matrix(res.profile$b) # population-level
  #   uhhat_all = uhhat[, rep(1:ncol(uhhat), K)]
  #   Yihat <- Yhat +  rowSums(Z * uhhat_all) # subject-level
  # } else {
  #   Yhat <- Yihat <- NULL
  # }
  # #kk = 6; whch <- which(id.site == id.site.uniq[kk]); length(whch)
  # #sqrt(mean((Yhat[whch] - Y[whch])^2)); sqrt(mean((Yihat[whch] - Y[whch])^2))


  return(list(b = res.profile$b,
              # b.sd = sqrt(vd),     # sd of fixed effect est
              # wald = wald,   # Wald-test statistic
              # lb = lb,       # lower-bound
              # ub = ub,       # uppper-bound
              # XvX = res.profile$allterms$bterm1, # X^{T}V^{-1}X
              # uh = res.profile$uh, # BLUP of random effects
              # uhM = uhhat_subj,
              # varuh = res.profile$varuh,  # Variance (based on prediction error)
              # varuh_post = res.profile$varuh_post, # posterior
              # Yhat = Yhat,         # population-level prediction (WOT random effects) XB
              # Yihat = Yihat,       # subject-specific prediction XB + Zu
              # mAIC = mAIC,
              V = V,
              s2 = s2,
              res = res, res.profile = res.profile))
}


