library(tidyverse)
library(data.table)
library(lme4)
library(nlme)
library(Matrix)
library(minqa)

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
peal.get.summary <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL, m_h_all = NULL){
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

    ShX  <- crossprod((Xh * wth), Xh)
    ShXZ <- crossprod((Xh * wth), Zh)
    ShXY <- crossprod((Xh * wth), Yh)
    ShZ  <- crossprod((Zh * wth), Zh)
    ShZY <- crossprod((Zh * wth), Yh)
    ShY  = sum(Yh ^ 2 *wth)
    Nh <- sum(id.site == sh)
    mh = m_h_all[h,]

    ShXYZ[[sh]] <- list(ShX  = ShX, ShXZ = ShXZ, ShXY = ShXY,
                        ShZ  = ShZ, ShZY = ShZY, ShY  = ShY, Nh = Nh, mh = mh)
  }

  return(ShXYZ)
}


construct_Z_list <- function(data, X_RS) {
  H <- length(unique(data$site))
  Z <- list()

  for (i in 1:H) {
    count_mat <- data %>%
      filter(site == i) %>%
      group_by(site, patient) %>%
      summarise(n_hi = n(), .groups = 'drop')

    df <- data %>%
      filter(site == i) %>%
      dplyr::select(c(patient, all_of(X_RS)))

    # Extract unique patient IDs
    patients <- unique(df$patient)

    # Extract covariate columns only (excluding patient ID)
    covariate_cols <- colnames(df)[2:ncol(df)]

    # Create patient-specific matrices
    patient_matrices <- lapply(patients, function(pat) {
      patient_data <- df %>%
        filter(patient == pat) %>%
        dplyr::select(all_of(covariate_cols))
      as.matrix(patient_data)
    })

    # Create block diagonal matrix
    block_diag_matrix <- as.matrix(bdiag(patient_matrices))

    # Combine with site-level matrix
    Z[[i]] <- cbind(generate_Zhv_matrix(count_mat), block_diag_matrix)
  }

  return(Z)
}



###### New

## updated lmm.profile3 with correlated patient REs
lmm.profile3 <- function(par, pooled = FALSE, reml = TRUE,
                         Y, X, Z, id.site, weights = NULL,
                         ShXYZ, rcpp = FALSE) {
  if (pooled) {
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
  } else {
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    K  <- length(ShXYZ)
  }

  lpterm1 <- lpterm2 <- remlterm <- 0
  bterm1 <- matrix(0, px, px)
  bterm2 <- rep(0, px)
  N <- 0

  for (h in seq_along(id.site.uniq)) {
    sh  <- id.site.uniq[h]
    ShX <- ShXYZ[[sh]]$ShX
    ShXZ <- ShXYZ[[sh]]$ShXZ
    ShXY <- ShXYZ[[sh]]$ShXY
    ShZ <- ShXYZ[[sh]]$ShZ
    ShZY <- ShXYZ[[sh]]$ShZY
    ShY  <- ShXYZ[[sh]]$ShY
    Nh   <- ShXYZ[[sh]]$Nh
    mh   <- ShXYZ[[sh]]$mh  # number of intercept+slope patient blocks

    N <- N + Nh

    ## extract parameters for site h
    sigma_site2   <- par[1]           # site intercept var
    sigma_int2    <- par[1 + h]       # patient intercept var
    sigma_slope2  <- par[1 + K + h]   # patient slope var
    rho1           <- par[1 + 2*K + h] # correlation
    rho2           <- par[1 + 3*K + h] # correlation
    rho3          <- par[1 + 4*K + h] # correlation

    ## build covariance V
    Sigma_patient <- matrix(c(sigma_int2,
                              rho1 * sqrt(sigma_int2 * sigma_slope2),
                              rho1 * sqrt(sigma_int2 * sigma_slope2),
                              sigma_slope2), 2, 2)

    # repeat mh times for patients
    block_patients <- Matrix::bdiag(replicate(c(mh$`n_distinct(patient)`), Sigma_patient, simplify = FALSE))

    # remaining slopes (if any)
    extra_slopes <- pzh <- ncol(ShZ)
    extra_slopes <- pzh - 1 - 2*mh  # adjust if >1 slope per patient
    if (extra_slopes > 0) {
      block_extra <- diag(sigma_slope2, extra_slopes)
      V <- Matrix::bdiag(sigma_site2, block_patients, block_extra)
    } else {
      V <- Matrix::bdiag(sigma_site2, block_patients)
    }
    V <- as.matrix(V)
    V[,1] <- c(sigma_site2, rep(c(rho2 * sqrt(sigma_site2 * sigma_int2), rho3 * sqrt(sigma_site2 * sigma_slope2)),
                                mh$`n_distinct(patient)`))
    V[1,] <- c(sigma_site2, rep(c(rho2 * sqrt(sigma_site2 * sigma_int2), rho3 * sqrt(sigma_site2 * sigma_slope2)),
                                mh$`n_distinct(patient)`))

    reorder_V <- function(V, has_site = TRUE) {
      p <- ncol(V)
      if (has_site) {
        m <- (p - 1) / 2L                  # number of patients
        stopifnot(p == 1 + 2*m)
        idx_site  <- 1L
        idx_ri    <- idx_site + 1L + 2L*(0:(m-1))   # 2,4,6,...
        idx_rs    <- idx_site + 2L + 2L*(0:(m-1))   # 3,5,7,...
        ord <- c(idx_site, idx_ri, idx_rs)
      } else {
        m <- p / 2L
        stopifnot(p == 2*m)
        idx_ri <- 1L + 2L*(0:(m-1))                 # 1,3,5,...
        idx_rs <- 2L + 2L*(0:(m-1))                 # 2,4,6,...
        ord <- c(idx_ri, idx_rs)
      }
      V[ord, ord, drop = FALSE]
    }
    V = reorder_V(V, has_site = TRUE)


    Vinv <- solve(Matrix::nearPD(V)$mat)


    ## likelihood pieces
    A <- diag(1, ncol(V)) + ShZ %*% V
    logdet <- as.numeric(determinant(A, logarithm = TRUE)$modulus)
    lpterm1 <- lpterm1 + logdet

    A_h <- Vinv + ShZ

    bterm1 <- bterm1 + (ShX - ShXZ %*% solve(A_h, t(ShXZ)))
    bterm2 <- bterm2 + (ShXY - ShXZ %*% solve(A_h, ShZY))
    lpterm2 <- lpterm2 + (ShY  - t(ShZY) %*% solve(A_h, ShZY))
  }

  L <- chol(bterm1)
  b <- backsolve(L, forwardsolve(t(L), bterm2))
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b)

  if (reml) {
    s2 <- qterm / (N - px)
    remlterm <- determinant(bterm1 / s2, logarithm = TRUE)$modulus
    lp <- -(lpterm1 + N * log(s2) + qterm/s2 + remlterm) / 2
  } else {
    s2 <- qterm / N
    lp <- -(lpterm1 + (1 + log(qterm * 2 * pi / N)) * N) / 2
  }

  lk <- - (lpterm1 + (1+log(qterm*2*pi/N))*N) / 2

  res <- list(lp = lp, b = b, s2 = s2, lk = lk,
              allterms = list(lpterm1 = lpterm1, lpterm2 = lpterm2,
                              qterm = qterm, remlterm = remlterm,
                              bterm1 = bterm1, bterm2 = bterm2))
  return(res)
}


## updated peal.fit.RS
peal.fit.CorrRS <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                        pooled = FALSE, reml = TRUE, common.s2 = TRUE,
                        ShXYZ = list(), corstr = 'independence',
                        mypar.init = NULL, hessian = FALSE, verbose = TRUE) {

  if (pooled) {
    id.site.uniq <- unique(id.site)
    K <- length(id.site.uniq)
    px <- ncol(X)
    ShXYZ <- lmm.get.summary3(Y, X, Z, weights = weights, id.site = id.site)
  } else {
    id.site.uniq <- names(ShXYZ)
    px <- ncol(ShXYZ[[1]]$ShX)
    K  <- length(ShXYZ)
  }

  ## new parameter dimension
  pz <- 1 + 3*K   # site var + patient intercepts + patient slopes + correlations

  ## init
  if (is.null(mypar.init)) {
    mypar.init <- c(1, rep(1, 2*K), rep(0, K*3))  # variances=1, rhos=0
    if (verbose) cat('Default mypar.init (var comp & rho) = ', mypar.init, '\n')
  }

  fn <- function(parameter) {
    return(-lmm.profile3(par = parameter, pooled = FALSE, reml,
                         Y, X, Z, id.site, weights, ShXYZ)$lp)
  }

  lower <- c(rep(1e-6, 1 + 2*K), rep(-0.99, K*3))
  upper <- c(rep(Inf,   1 + 2*K), rep( 0.99, K*3))

  res <- optim(mypar.init, fn,
               hessian = TRUE,
               method = "L-BFGS-B",
               lower = lower,
               upper = upper)

  if (verbose) cat(ifelse(all(res$convergence == 0), "Convergence Reached", "Non-convergence!"),
                   "and",
                   ifelse(all(eigen(res$hessian)$value > 0), "Hessian PD", "Hessian not PD"), "\n",
                   "The number of function evaluations used is ", res$counts[1], '\n')

  mypar <- res$par
  res.profile <- lmm.profile3(par = mypar, pooled = FALSE, reml,
                              Y, X, Z, id.site, weights, ShXYZ)
  s2 <- res.profile$s2
  V  <- mypar * s2

  ## inference
  vd <- diag(solve(res.profile$allterms$bterm1))
  if (common.s2) vd <- diag(solve(res.profile$allterms$bterm1 / s2))  # scale back
  wald <- res.profile$b / sqrt(vd)

  lb <- res.profile$b - 1.96 * sqrt(vd)
  ub <- res.profile$b + 1.96 * sqrt(vd)

  cat('Done!', "\n")

  ## split out parameters
  sigma_site2 <- V[1]
  sigma_int2  <- V[2:(K+1)]
  sigma_slope2<- V[(K+2):(2*K+1)]
  rho         <- mypar[(2*K+2):(length(mypar))]

  return(list(b = res.profile$b,
              b.sd = sqrt(vd),
              wald = wald,
              lb = lb, ub = ub,
              sigma_site2 = sigma_site2,
              sigma_int2 = sigma_int2,
              sigma_slope2 = sigma_slope2,
              rho = rho,
              V = V, s2 = s2,
              res = res, res.profile = res.profile))
}

