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



## get summary stats from each site for distributed lmm
lmm.get.summary3 <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL){
  if(is.null(weights)) weights <- rep(1, length(Y))
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  pz <- ncol(Z)

  ShXYZ <- list()
  for(h in seq_along(id.site.uniq)){
    sh = id.site.uniq[h]
    wth = weights[id.site == sh]
    Xh <- X[id.site == sh, ]
    Yh <- Y[id.site == sh]

    Zh <- Z[id.site == sh, ]
    non_zero_columns <- colSums(Zh != 0) > 0
    Zh <- Zh[, non_zero_columns, drop = FALSE]

    # site_summary = three_lvl_dat %>%
    #   filter(site == sh) %>%
    #   group_by(patient) %>%
    #   summarise(n_hi = n())
    #
    # Zh <- generate_Zhv_matrix(site_summary)

    ShX  = t(Xh*wth) %*% Xh
    ShXZ = t(Xh*wth) %*% Zh
    ShXY = t(Xh*wth) %*% Yh
    ShZ  = t(Zh*wth) %*% Zh
    ShZY = t(Zh*wth) %*% Yh
    ShY  = sum(Yh ^ 2 *wth)
    Nh <- sum(id.site == sh)
    ShXYZ[[sh]] <- list(ShX  = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ  = ShZ, ShZY = ShZY, ShY  = ShY, Nh = Nh)
  }

  return(ShXYZ)
}



## further profile out the residual var s2, used if common.s2=T (deemed as usual situation)
lmm.profile03 <- function(V,                # var components para, = original V / s2    # s2,
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
    pz <- ncol(ShXYZ[[1]]$ShXZ)
  }

  if(rcpp == F){
    lpterm1 = lpterm2 = remlterm = 0
    bterm1 <- matrix(0, px, px)  # sum_i Xi' \Sigma_i^-1 Xi
    bterm2 <- rep(0, px)         # sum_i Xi' \Sigma_i^-1 Yi
    Vinv <- solve(V)
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

      N <- N + Nh

      # tmp <- log(det(diag(1, pz) + SiZ %*% V))  # improve
      # if(is.na(tmp)) cat(diag(V), '...', s2i, '...', V[1,], '\n')
      # logdet <- ni * log(s2i) + log(det(diag(1, pz) + SiZ %*% V / s2i))
      logdet <- log(det(diag(1, pz) + ShZ %*% V))    # -log(det(diag(wti))) omitted as wti are the fixed weights
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

    # lk <- -(lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2 # To be used in calculating mAIC
    lk <- - (lpterm1 + (1+log(qterm*2*pi/N))*N) / 2

    # Loop to estimate uh
    # uh = V Z_i^{T} Sigma_{i}^{-1} (Y_i - X_i \beta)
    for(h in seq_along(id.site.uniq)){ # re-call this loop
        sh <- id.site.uniq[h]

        ShX  <- ShXYZ[[sh]]$ShX
        ShXZ <- ShXYZ[[sh]]$ShXZ
        ShXY <- ShXYZ[[sh]]$ShXY
        ShZ <- ShXYZ[[sh]]$ShZ
        ShZY <- ShXYZ[[sh]]$ShZY
        ShY  <- ShXYZ[[sh]]$ShY
        Nh <- ShXYZ[[sh]]$Nh

      s2h <- s2   # [ii]  # sigma_i^2
      uhterm1 <- (ShZY - ShZ %*% Wh[[h]] %*% ShZY) / s2h
      uhterm2 <-  ((t(ShXZ) - ShZ %*% Wh[[h]] %*% t(ShXZ)) / s2h) %*% as.matrix(b)
      uh[[h]] <- V %*% as.numeric(uhterm1 - uhterm2) * s2h #

      vterm1 <- V %*% ((ShZ - ShZ %*% Wh[[h]] %*% ShZ) / s2h) %*% V * (s2h ^ 2)    #
      vterm2 <- V %*% ((t(ShXZ) - ShZ %*% Wh[[h]] %*% t(ShXZ)) / s2h)  * s2h   #
      # varuh[[ii]] <- (V * s2i) - vterm1 + (vterm2 %*% t(vterm2) / lpterm1)        # ? lpterm1 is logdet
      # varuh_post[[ii]] <- vterm1 - (vterm2 %*% t(vterm2) / lpterm1)               # ? lpterm1 is logdet
      varuh[[h]] <- (V * s2h) - vterm1 + (vterm2 %*% solve(bterm1, t(vterm2)))     # 20210111
      varuh_post[[h]] <- vterm1 - (vterm2 %*% solve(bterm1, t(vterm2)))            # 20210111
    }


    res <- list(lp = lp, b = b,
                s2 = s2,
                uh = uh,
                lk = lk,
                varuh = varuh,
                varuh_post = varuh_post,
                allterms = list(lpterm1 = lpterm1,
                                lpterm2 = lpterm2,
                                qterm = qterm,
                                remlterm = remlterm,
                                # remlterm.i = remlterm.i,
                                bterm1 = bterm1,
                                bterm2 = bterm2))
  }
  return(res)
}




# ## lmm profile likelihood w.r.t. variance components
# lmm.profile3 <- function(V, s2,             # var components para
#                         pooled = F, reml = T,
#                         Y, X, Z, id.site, weights=NULL,  # if pooled ipd is available
#                         ShXYZ,             # if pooled ipd is not available, use summary stats
#                         rcpp = F){
#   if(pooled == T){
#     id.site.uniq <- unique(id.site)
#     px <- ncol(X)
#     pz <- ncol(Z)
#   }else{
#     id.site.uniq <- names(ShXYZ)
#     px <- ncol(ShXYZ[[1]]$ShX)
#     pz <- ncol(ShXYZ[[1]]$ShXZ)
#   }
#
#   # allows assuming the same reShdual var across sites
#   if(length(s2) == 1) {
#     s2 <- rep(s2, length(id.site.uniq))
#   }
#
#   if(rcpp == F){
#     lpterm1 = lpterm2 = remlterm = 0
#     remlterm.i = rep(NA, length(id.site.uniq))
#
#     bterm1 <- matrix(0, px, px)
#     bterm2 <- rep(0, px)
#     Vinv <- solve(V)
#     Wh <- uh <- varuh <- varuh_post <- list()  # save this for each subject for BLUP
#
#     for(h in seq_along(id.site.uniq)){
#       Sh <- id.site.uniq[h]
#       if(pooled == T){
#         ShXYZ <- lmm.get.summary3(Y, X, Z, weights = weights, id.site = id.site)
#       } #else{
#       ShX  <- ShXYZ[[Sh]]$ShX
#       ShXZ <- ShXYZ[[Sh]]$ShXZ
#       ShXY <- ShXYZ[[Sh]]$ShXY
#       ShZ <- ShXYZ[[Sh]]$ShZ
#       ShZY <- ShXYZ[[Sh]]$ShZY
#       ShY  <- ShXYZ[[Sh]]$ShY
#       Nh <- ShXYZ[[Sh]]$Nh
#       # }
#
#       s2i <- s2[h]  # Sigma_i^2
#
#       tmp <- log(det(diag(1, pz) + ShZ %*% V / s2i))
#
#       if(is.na(tmp)) cat(diag(V), '...', s2i, '...', V[1,], '\n')
#       # logdet <- Nh * log(s2i) + log(det(diag(1, pz) + ShZ %*% V / s2i))
#       logdet <- Nh * log(s2i) + log(det(diag(1, pz) + ShZ %*% V / s2i))
#
#       # log(max(1e-14, det(diag(1,pz)+ShZ%*%V/s2i)))
#       lpterm1 <- lpterm1 + logdet
#
#       Wi[[h]] <- solve(s2i * Vinv + ShZ)
#       bterm1 <- bterm1 + (ShX - ShXZ %*% Wi[[h]] %*% t(ShXZ)) / s2i
#       bterm2 <- bterm2 + (ShXY - ShXZ %*% Wi[[h]] %*% ShZY) / s2i
#       lpterm2 <- lpterm2 + (ShY - t(ShZY) %*% Wi[[h]] %*% ShZY) / s2i
#
#       if(reml == T){
#         # tmp = log(det((ShX-ShXZ%*%Wi[[h]]%*%t(ShXZ))/s2i))
#         # if(is.na(tmp)) cat(Wi[[h]], '...', s2i, '\n')
#         remlterm.i[h] <- log(abs(det((ShX - ShXZ %*% Wi[[h]] %*% t(ShXZ))/s2i))) # abs() to make sure poShtivity
#         # remlterm <- remlterm + log(det((ShX - ShXZ %*% Wi[[h]] %*% t(ShXZ)) / s2i))
#       }
#     }
#
#     b <- solve(bterm1, bterm2)
#     if(reml == T){
#       remlterm = sum(remlterm.i[is.finite(remlterm.i)])
#       lp <- - (lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b + remlterm) / 2
#     }else{
#       lp <- - (lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2
#     }
#
#     lk <- -(lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2 # To be used in calculating mAIC
#
#     # Loop to estimate uh
#     # uh = V Z_i^{T} Sigma_{i}^{-1} (Y_i - X_i \beta)
#     for(h in seq_along(id.site.uniq)){ # re-call this loop
#       Sh <- id.site.uniq[h]
#       if(pooled == T){
#         Xi <- X[id.site == Sh, ]
#         Zi <- Z[id.site == Sh, ]
#         Yi <- Y[id.site == Sh]
#         ShX  <- t(Xi) %*% Xi
#         ShXZ <- t(Xi) %*% Zi
#         ShXY <- t(Xi) %*% Yi
#         ShZ  <- t(Zi) %*% Zi
#         ShZY <- t(Zi) %*% Yi
#         ShY  <- sum(Yi ^ 2)
#         ni <- sum(id.site == Sh)
#       }else{
#         ShX  <- ShXYZ[[Sh]]$ShX
#         ShXZ <- ShXYZ[[Sh]]$ShXZ
#         ShXY <- ShXYZ[[Sh]]$ShXY
#         ShZ <- ShXYZ[[Sh]]$ShZ
#         ShZY <- ShXYZ[[Sh]]$ShZY
#         ShY <- ShXYZ[[Sh]]$ShY
#         Nh <- ShXYZ[[Sh]]$Nh
#       }
#
#       s2h <- s2[h]  # Sigma_i^2
#       uhterm1 <- (ShZY - ShZ %*% Wh[[h]] %*% ShZY) / s2h
#       uhterm2 <-  ((t(ShXZ) - ShZ %*% Wh[[h]] %*% t(ShXZ)) / s2h) %*% as.matrix(b)
#       uh[[h]] <- V %*% as.numeric(uhterm1 - uhterm2)
#
#       vterm1 <- V %*% ((ShZ - ShZ %*% Wh[[h]] %*% ShZ) / s2h) %*% V
#       vterm2 <- V %*% ((t(ShXZ) - ShZ %*% Wh[[h]] %*% t(ShXZ)) / s2h)
#       varuh[[h]] <- V - vterm1 + (vterm2 %*% t(vterm2) / lpterm1)
#       varuh_post[[h]] <- vterm1 - (vterm2 %*% t(vterm2) / lpterm1)
#     }
#
#
#     res <- list(lp = lp, b = b, lk = lk)
#   }
#   # ShY - t(ShZY)%*%Wi%*%ShZY - 2*(t(ShXY)-t(ShZY)%*%Wi%*%t(ShXZ))%*%b + t(b)%*%(ShX-ShXZ%*%Wi%*%t(ShXZ))%*%b
#   return(res)
# }





# fit 3 level DLMM

lmm.fit3 <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                    pooled = F, reml = T,
                    common.s2 = F,      # common reShdual var across sites
                    ShXYZ = list(),
                    corstr = 'independence', # 'exchangeable', 'ar1', 'unstructured'),
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
    pz <- ncol(ShXYZ[[1]]$ShXZ)
    K <- length(ShXYZ)
  }

  ## 20200629: now further profile out s2 if common.s2=T
  if(common.s2 == T){
    ns <- 1 # number of s2 para
    fn <- function(mypar){
        V <- diag(c(mypar[1], rep(mypar[2], (pz-1))), pz)
        # V <- diag(mypar[1 : pz], pz)

      return(-lmm.profile03(V, pooled=F, reml, Y, X, Z, id.site, weights, ShXYZ = ShXYZ)$lp)
    }

    if(is.null(mypar.init)){
        mypar.init <- c(1, 1)
        # mypar.init <- rep(0.5, pz)

      cat('default mypar.init (var comp) = ', mypar.init, '\n')
    }

    # res <- optim(mypar.init, fn, hessian = F, method = c("Nelder-Mead"))
    res <- minqa::bobyqa(mypar.init, fn, lower=rep(1e-6, length(mypar.init)), control=list(maxfun=1e5))

    mypar <- res$par
    V <- diag(c(mypar[1], rep(mypar[2], (pz-1))), pz) # V is really Theta here!
    # V <- diag(mypar[1 : pz], pz)

    res.profile <- lmm.profile03(V = V, pooled=F, reml, Y, X, Z, id.site, weights, ShXYZ)
    s2 <- res.profile$s2
    V <- V * s2             # scale back

  }else{  # if common.s2=F, can't profile out s2 vector
    ns <- K
    fn <- function(mypar){
      V <- diag(c(mypar[1], rep(mypar[2], (pz-1))), pz)
      s2 <- (mypar[-c(1 : 2)])
        # V <- diag(mypar[1 : pz], pz)
        # s2 <- exp(mypar[-c(1 : pz)])

      return(-lmm.profile3(V, s2, pooled=F, reml, Y, X, Z, id.site, weights, ShXYZ)$lp)
    }

    if(is.null(mypar.init)){
        mypar.init <- c(rep(0.5, 2), rep(0.5, ns))

      cat('default mypar.init (var comp) = ', mypar.init, '\n')
    }

    # res <- optim(mypar.init, fn, hessian = hessian)
    res <- bobyqa(mypar.init, fn, lower=rep(1e-6, length(mypar.init)), control=list(maxfun=1e5))

    mypar <- res$par

    V <- diag(c(mypar[1], rep(mypar[2], (pz-1))), pz)
    s2 <- (mypar[-c(1 : 2)])
    # V <- diag(mypar[1 : pz], pz)
    # s2 <- mypar[- c(1 : pz)]

    res.profile <- lmm.profile3(V = V, s2 = s2, pooled, reml, Y, X, Z, id.site, ShXYZ)
  }



    ## New added
    ## Inference (Wald test statistic)
    vd <- diag(solve(res.profile$allterms$bterm1))
    if(common.s2==T) vd <- diag(solve(res.profile$allterms$bterm1 / s2))  # scale back
    wald <- res.profile$b / sqrt(vd)

    ## 95% CI for fixed effects
    lb <- res.profile$b -  1.96 * sqrt(vd)
    ub <- res.profile$b +  1.96 * sqrt(vd)

    ## Marginal AIC: OKAY to use if the main interest is to model fixed population effects with a reasonable correlation structur (Kneib & Greven (2010))
    mAIC <- 2 * res.profile$lk + 2 * (px + (length(mypar) - ns))

    ## Conditional AIC: Vaida & Blanchard (2005) expresShon details in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2572765/pdf/nihms62635.pdf
    ## However a better approach should be in Kneib & Steven (2010) Appendix:B https://arxiv.org/pdf/1803.05664.pdf
    ## assuming V and Sigma^2 are (UN)KNOWN
    # if(pooled == T & common.s2 == T){
    #   c2 <- diag(s2, ncol(V)) * V
    #
    #   #if(common.s2 == T){
    #   #  c2 <- bdiag(lapply(seq_len(K), function(kk) diag(s2[1], ncol(V)) * V))  #
    #   #} else if(common.s2 == F){
    #   #  c2 <- bdiag(lapply(seq_len(K), function(kk) diag(s2[kk], ncol(V)) * V)) # block diagonal
    #   #}
    #
    #   LLinv <- solve(c2, diag(1, ncol(c2)))
    #   c3 <- eigen(LLinv)
    #   Lambda <- c3$vectors %*% diag(sqrt(c3$values))
    #
    #   c11 <- cbind(X, Z)
    #   c12 <- cbind(matrix(0, ncol = px, nrow = nrow(Lambda)), Lambda)
    #   # dim(c11); dim(c12)
    #   M <- rbind(c11, c12)
    #
    #   MtM <- solve(t(M) %*% M)
    #   H1 <- c11 %*% MtM %*% t(c11)
    #   trace <- sum(diag(H1))
    #   cAIC_vb <- 2 * res.profile$lk + 2 * trace
    # } else{
    #   cAIC_vb = NULL
    # }

    ## Prediction
    uhhat <- as.matrix(do.call(rbind, lapply(seq_len(K), function(kk) {

      ll <- length(which(id.site == id.site.uniq[kk]))
      dm <- matrix(1, nrow = ll, ncol = pz)

      sweep(dm, MARGIN = 2, res.profile$uh[[kk]], `*`)

    })))

    uhhat_subj <- as.matrix(do.call(rbind, lapply(seq_len(K), function(kk) {
      t(res.profile$uh[[kk]])
    })))

    if(pooled == T){
      Yhat <- X %*% as.matrix(res.profile$b) # population-level
      uhhat_all = uhhat[, rep(1:ncol(uhhat), K)]
      Yihat <- Yhat +  rowSums(Z * uhhat_all) # subject-level
    } else {
      Yhat <- Yihat <- NULL
    }
    #kk = 6; whch <- which(id.site == id.site.uniq[kk]); length(whch)
    #sqrt(mean((Yhat[whch] - Y[whch])^2)); sqrt(mean((Yihat[whch] - Y[whch])^2))



  return(list(b = res.profile$b,
              b.sd = sqrt(vd),     # sd of fixed effect est
              wald = wald,   # Wald-test statistic
              lb = lb,       # lower-bound
              ub = ub,       # uppper-bound
              XvX = res.profile$allterms$bterm1, # X^{T}V^{-1}X
              uh = res.profile$uh, # BLUP of random effects
              uhM = uhhat_subj,
              varuh = res.profile$varuh,  # Variance (based on prediction error)
              varuh_post = res.profile$varuh_post, # posterior
              Yhat = Yhat,         # population-level prediction (WOT random effects) XB
              Yihat = Yihat,       # subject-specific prediction XB + Zu
              mAIC = mAIC,
              # cAIC_vb = cAIC_vb,
              V = V,
              s2 = s2,
              res = res, res.profile = res.profile))
}


#
# # test
# id.site.uniq <- names(ShXYZ)
# px <- ncol(ShXYZ[[1]]$ShX)
# pz <- ncol(ShXYZ[[1]]$ShXZ)
#
# mypar.init <- rep(0.5, pz)
# V <- diag(mypar.init[1 : pz], pz)
#
# # -lmm.profile3(V, s2 = 1, pooled=F, reml=F, Y, X, Z, id.site, weights = NULL, ShXYZ)$lp
# -lmm.profile03(V, pooled=F, reml=T, ShXYZ = ShXYZ)$lp
#
#
#
#
# V <- diag(c(9, rep(4, (pz-1))), pz)
# # V <- diag(mypar[1 : pz], pz)
#
# res.profile <- lmm.profile03(V = V, pooled=F, reml=T, ShXYZ = ShXYZ)
# s2 <- res.profile$s2
# V <- V * s2

