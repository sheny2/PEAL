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


# -- Stage 0: summaries ----------------------------------------------------

lmm_get_summary_multivar <- function(Y, X, Z_list, id.site, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(Y))
  X <- as.matrix(X); Y <- as.matrix(Y)
  id.site <- as.character(id.site)
  uniq <- unique(id.site)

  out <- list()
  for (h in seq_along(uniq)) {
    sh <- uniq[h]
    idx <- id.site == sh
    w   <- weights[idx]
    Xh <- X[idx,,drop=FALSE]
    Yh <- Y[idx,,drop=FALSE]
    Zh <- Z_list[[h]]

    W  <- Diagonal(x = w)
    # Weighted cross prods
    ShX  <- crossprod(Xh, W %*% Xh)
    ShXZ <- crossprod(Xh, W %*% Zh)
    ShXY <- crossprod(Xh, W %*% Yh)
    ShZ  <- crossprod(Zh, W %*% Zh)
    ShZY <- crossprod(Zh, W %*% Yh)
    ShY  <- crossprod(Yh, W %*% Yh)

    out[[sh]] <- list(ShX=ShX, ShXZ=ShXZ, ShXY=ShXY,
                      ShZ=ShZ, ShZY=ShZY, ShY=ShY, Nh = nrow(Xh))
  }
  out
}


# -- helpers ---------------------------------------------------------------

make_Sigma_e <- function(sigma2, rho, py) {
  # Exchangeable residual covariance
  if (rho <= -1/(py-1) || rho >= 1 || sigma2 <= 0) stop("Invalid (sigma2, rho).")
  Sig <- matrix(rho * sigma2, py, py); diag(Sig) <- sigma2
  Sig
}

# Extract leading-site OLS beta as a quick-and-dirty start
leading_site_beta0 <- function(Sh) {
  # Sh: list with ShX, ShXY for the leading site (site "1")
  solve(Sh$ShX, Sh$ShXY)  # px x py
}

# Build residualized summaries given B0
residualize_summaries <- function(ShXYZ, B0) {
  out <- ShXYZ
  for (nm in names(ShXYZ)) {
    Sh <- ShXYZ[[nm]]
    ShZY_res <- Sh$ShZY - t(Sh$ShXZ) %*% B0
    ShY_res  <- Sh$ShY  - t(B0) %*% Sh$ShXY - t(Sh$ShXY) %*% B0 + t(B0) %*% Sh$ShX %*% B0
    out[[nm]]$ShZY_res <- ShZY_res
    out[[nm]]$ShY_res  <- ShY_res
  }
  out
}

# numerically stable log|A| via Cholesky (with jitter)
chol_logdet <- function(A, max_tries = 6, eps = 1e-8) {
  A <- (A + t(A)) / 2
  for (k in 0:max_tries) {
    adj <- if (k == 0) 0 else eps * (10^k)
    ck  <- try(chol(A + Diagonal(nrow(A)) * adj), silent = TRUE)
    if (!inherits(ck, "try-error")) {
      return(list(L = ck, logdet = 2 * sum(log(diag(ck))), jitter = adj))
    }
  }
  stop("Cholesky failed")
}


# # Negative profile loglik using only summaries and B0
# --- NEGATIVE profile loglik (β known) ------------------------------------
negloglik_var <- function(par, ShXYZ_res, py) {
  # par = c(sigma_u2, sigma_v2_site1, ..., sigma_v2_siteK, sigma2, rho)
  K <- length(ShXYZ_res)
  sigma_u2 <- par[1]
  sigma_v2 <- par[2:(K+1)]
  sigma2   <- par[K+2]
  rho      <- par[K+3]

  if (sigma_u2 <= 0 || any(sigma_v2 <= 0) ||
      sigma2 <= 0 || rho <= -1/(py-1) || rho >= 1) return(Inf)

  # Residual covariance and its inverse / logdet
  Sig  <- make_Sigma_e(sigma2, rho, py)
  chS  <- chol_logdet(Matrix(Sig, sparse = FALSE))
  iSig <- chol2inv(chS$L)
  logdet_Sig <- chS$logdet

  total_det <- 0.0   # sum_h [ log|A_h| + p*log|D_h| ]
  total_q   <- 0.0   # sum_h [ tr(iSig * ShY_res) - v_h' A_h^{-1} v_h ]
  Ndot      <- 0L

  for (h in seq_along(ShXYZ_res)) {
    Sh  <- ShXYZ_res[[h]]
    pzh <- ncol(Sh$ShZ)

    # D_h (RE covariance per site h): diag(sigma_u2, sigma_vh2, ..., sigma_vh2)
    Dh  <- diag(c(sigma_u2, rep(sigma_v2[h], pzh - 1)), pzh)
    iDh <- solve(Dh)  # tiny (pzh x pzh)

    # A_h = (D_h^{-1} ⊗ I_p) + (ShZ_h ⊗ Σ^{-1})
    Ah  <- kronecker(iDh, Diagonal(py)) + kronecker(Sh$ShZ, iSig)

    # log|A_h| (via chol) and p*log|D_h|
    chA <- chol_logdet(Ah)
    logdet_Dh <- as.numeric(determinant(Dh, logarithm = TRUE)$modulus)
    total_det <- total_det + chA$logdet + py * logdet_Dh

    # Quadratic:
    # q0_h = tr(Σ^{-1} * S_YY(res))
    q0_h <- sum(iSig * Sh$ShY_res)  # trace(A^T B) trick

    # v_h = vec(S_YZ_res %*% Σ^{-1})  with S_YZ_res = (py x pzh)
    S_YZ_res <- as.matrix(t(Sh$ShZY_res))  # py × pzh
    v <- as.vector(iSig %*% S_YZ_res)

    # Solve Ah^{-1} v without forming the inverse
    tmp <- backsolve(chA$L, forwardsolve(t(chA$L), v))
    quad_h <- q0_h - drop(crossprod(v, tmp))

    total_q <- total_q + quad_h
    Ndot    <- Ndot + Sh$Nh
  }

  # Add residual determinant across all observations: Ndot * log|Σ_e|
  nll <- 0.5 * (total_det + total_q + Ndot * logdet_Sig)
  return(nll)
}



# # Given Sigma_e and D_h’s, assemble Sxx and sxy from summaries (no raw data), then solve for B
solve_beta_global <- function(ShXYZ, Sigma_e, sigma_u2, sigma_v2_vec) {
  K  <- length(ShXYZ)
  py <- ncol(ShXYZ[[1]]$ShY)
  px <- ncol(ShXYZ[[1]]$ShX)

  iSig <- solve(Sigma_e)

  # use base matrix, not Matrix::dgeMatrix, to keep chol/backsolve simple
  Sxx <- matrix(0, nrow = px*py, ncol = px*py)
  sxy <- matrix(0, nrow = px*py, ncol = 1)

  for (h in seq_along(ShXYZ)) {
    Sh  <- ShXYZ[[h]]
    pzh <- ncol(Sh$ShZ)  # Z'Z is pz x pz
    Dh  <- diag(c(sigma_u2, rep(sigma_v2_vec[h], pzh-1)), pzh)
    iDh <- solve(Dh)

    # K: (py*pz) x (py*pz)
    Kh <- kronecker(diag(py), iDh) + kronecker(iSig, Sh$ShZ)
    Kh_inv <- solve(Kh)

    # Blocks
    A   <- kronecker(iSig, Sh$ShX)                # (py*px) x (py*px)
    Bxz <- kronecker(diag(py), Sh$ShXZ)           # (py*px) x (py*pz)
    Bzx <- t(Bxz)                                 # (py*pz) x (py*px)

    # Accumulate Sxx
    Sxx <- Sxx + (A - Bxz %*% Kh_inv %*% Bzx)

    # Right-hand side: use X'Y %*% iSig and Z'Y %*% iSig
    XY_sig <- Sh$ShXY %*% iSig                    # px x py
    ZY_sig <- Sh$ShZY %*% iSig                    # pz x py

    vXY <- as.vector(XY_sig)                      # vec(X'Y Σ^{-1})
    vZY <- as.vector(ZY_sig)                      # vec(Z'Y Σ^{-1})

    sxy <- sxy + (vXY - Bxz %*% Kh_inv %*% vZY)
  }

  # Symmetrize for numerical stability
  Sxx <- (Sxx + t(Sxx)) / 2

  # Cholesky solve (base)
  L <- chol(Sxx)                                  # upper-tri
  y <- forwardsolve(t(L), sxy)                    # solve L' y = sxy
  vecB <- backsolve(L, y)                         # solve L  x = y

  Bhat  <- matrix(vecB, nrow = px, ncol = py)
  VvecB <- chol2inv(L)                            # = solve(Sxx)

  list(B = Bhat, VvecB = VvecB)
}





# -- Stage 1 (“two-stage”): quick β0, then global variance ----

fit_stage1_variances <- function(ShXYZ, py, leading_site_name = NULL, par_init = NULL, verbose = TRUE) {
  if (is.null(leading_site_name)) leading_site_name <- names(ShXYZ)[1]

  # β^(0) from the leading site by OLS on summaries
  B0 <- leading_site_beta0(ShXYZ[[leading_site_name]])

  # Residualized summaries
  Sh_res <- residualize_summaries(ShXYZ, B0)

  K <- length(ShXYZ)
  if (is.null(par_init)) {
    par_init <- c( # sigma_u2, sigma_v2(h=1..K), sigma2, rho
      1, rep(1, K), 1, 0.1
    )
  }

  low_rho <- if (py > 1) -1/(py - 1) + 1e-6 else 0
  upper   <- c(rep(Inf, K+2), 1 - 1e-6)
  lower   <- c(rep(1e-8, K+2), -1 + 1e-6)

  obj <- function(p) negloglik_var(p, Sh_res, py)

  opt <- optim(par_init, obj, method = "L-BFGS-B", lower = lower, upper = upper, hessian = T)

  if (verbose) cat(ifelse(all(opt$convergence == 0, eigen(opt$hessian)$value > 0),
                          "Convergence Reached for Stage 1 variance", "Non-convergence!"), 'and',
                   ifelse(all(eigen(opt$hessian)$value > 0),
                          "Hessian PD", "Hessian not PD"), '\n',
                   "The number of function evaluations used is", opt$counts[1], '\n')

  sigma_u2 <- opt$par[1]
  sigma_v2 <- opt$par[2:(K+1)]
  sigma2   <- opt$par[K+2]
  rho      <- opt$par[K+3]
  Sigma_e  <- make_Sigma_e(sigma2, rho, py)

  list(B0 = B0, sigma_u2 = sigma_u2, sigma_v2 = sigma_v2,
       sigma2 = sigma2, rho = rho, Sigma_e = Sigma_e, opt = opt)
}

# -- Stage 2: global β given Σ_e, D_h -------------------------------------

fit_stage2_beta <- function(ShXYZ, Sigma_e, sigma_u2, sigma_v2) {
  res <- solve_beta_global(ShXYZ, Sigma_e, sigma_u2, sigma_v2)
  # reshape SEs per (px x py)
  se_vec <- sqrt(diag(res$VvecB))
  list(B = res$B, SE_vec = se_vec, VvecB = res$VvecB)
}


# -- One-shot driver -------------------------------------------------------

federated_lmm_two_stage <- function(Y, X, Z_list, id.site, weights = NULL,
                                    leading_site_name = NULL, par_init = NULL, verbose = TRUE) {
  ShXYZ <- lmm_get_summary_multivar(Y, X, Z_list, id.site, weights)
  py <- ncol(Y)

  # Stage 1: quick β0 from leading site, then global variance components via profile likelihood
  vfit <- fit_stage1_variances(ShXYZ, py, leading_site_name, par_init, verbose)
  # cat("rho:", vfit$rho, "\n")

  # Stage 2: β via closed-form GLS using summaries and \hat{Σ_e}, \hat{D}_h
  # bfit <- fit_stage2_beta(ShXYZ, vfit$Sigma_e, vfit$sigma_u2, vfit$sigma_v2)
  bfit <- lmm.profile03(par = c(sqrt(vfit$sigma_u2),sqrt(vfit$sigma_v2)), pooled = FALSE, reml = T, Y = NULL, X = NULL, Z = NULL, id.site, weights = NULL, ShXYZ, corstr = 'independence')


  list(
    # B = bfit$B,
    # B_SE_vec = bfit$SE_vec,       # vectorized SEs for vec(B)
    B = bfit$b,
    Sigma_e = vfit$Sigma_e,
    sigma_u = sqrt(vfit$sigma_u2),
    sigma_v = sqrt(vfit$sigma_v2),
    sigma   = sqrt(vfit$sigma2),
    rho     = vfit$rho,
    B0      = vfit$B0,
    opt_var = vfit$opt
  )
}







#############
lmm.profile03 <- function(par, pooled = FALSE, reml = TRUE,
                          Y, X, Z, id.site, weights = NULL,
                          ShXYZ, corstr = "independence", rcpp = FALSE) {
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
      R <- matrix(rho, pzh, pzh)  # Correlation matrix
      diag(R) <- 1
      R_r = R
      R_r[,1] = R_r[1,] = 0
      R_r[1,1] = 1
      V <- D %*% R_r %*% D
    }

    Vinv <- solve(V, diag(nrow(V)))  # Inverse of V
    # log-determinant
    A <- diag(col(V)) + ShZ %*% V / s2
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
    cat("Use REML version for better prediction")
  }

  res <- list(
    lp = lp,
    b = b       # (px x py)
  )
  return(res)
}

