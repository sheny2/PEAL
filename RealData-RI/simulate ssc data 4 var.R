##### Packages & Read in Data #####

Packages <- c("dplyr", "plyr", "magrittr", "ggplot2", "tidyr", "matrixcalc", "gridExtra", "MASS")
lapply(Packages, library, character.only = TRUE)

# set working directory to current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##### parameters #####

# variable names
variablenames <- c("pFVCt", "pDLCOt", "EFt", "RVSPt")

# number of variables
nmeas <- length(variablenames)

# select B random mcmc samples
B <- 500

# load deidentified data & mcmc est
load("deidentified_dat.Rdata")
load("combmcmc_fvcdlcoefrvsp_intsl.Rdata")


##### parameters from each iteration of MCMC chain #####

solmat <- combmcmc$Sol
covmat <- combmcmc$VCV

# number of fixed effects variables
nfe <- combmcmc$Fixed$nfl

FEsolmat <- solmat[, 1:nfe]

# reformatting solutions by variables
colorder <- c(colnames(FEsolmat)[startsWith(colnames(FEsolmat), "traitpFVCt")],
              colnames(FEsolmat)[startsWith(colnames(FEsolmat), "traitpDLCOt")],
              colnames(FEsolmat)[startsWith(colnames(FEsolmat), "traitEFt")],
              colnames(FEsolmat)[startsWith(colnames(FEsolmat), "traitRVSPt")])

ordFEsolmat <- FEsolmat[, order(match(colnames(FEsolmat), colorder))]

# posterior mean for fixed effects estimates
sol.postmean <- apply(solmat[, 1:nfe], 2, mean)
beta <- sol.postmean[order(match(names(sol.postmean), colorder))]

# posterior mean for variance estimates
covdat <- apply(combmcmc$VCV, 2, mean)

tmat <- covdat[endsWith(names(covdat), "pid")]
Gnames <- matrix(names(tmat), sqrt(length(tmat))) %>% diag

Gvec <- 1:length(Gnames)
ordervec <- c(Gvec[startsWith(Gnames, "traitpFVCt")],
              Gvec[startsWith(Gnames, "traitpDLCOt")],
              Gvec[startsWith(Gnames, "traitEFt")],
              Gvec[startsWith(Gnames, "traitRVSPt")])

G <- matrix(tmat, sqrt(length(tmat)))
G2 <- G[ordervec, ordervec]

tmatr <- covdat[endsWith(names(covdat), "units")]
R2 <- matrix(tmatr, sqrt(length(tmatr)))


#####  create var(Y) as a function of Z (datlist) and R #####

create_sigma <- function(datlist, R2){

  Sigma2 <- NULL

  for (i in 1:nmeas){

    Cmat <- NULL

    for (j in 1:nmeas){
      time1 <- datlist[[i]][, 2]; time2 <- datlist[[j]][, 2]

      if (i == j){
        # diagonal matrices
        cmat <- R2[i, j] * diag(1, length(time1))

      }else{
        # off diagonal matrices
        covmat <- outer(1:length(time1), 1:length(time2), FUN = "paste", sep = ",")
        comb <- cbind(1:length(time1), match(time1, time2, nomatch = NA))

        if(all(is.na(comb[,2]))){
          cmat <- matrix(0, nrow = length(time1), ncol = length(time2))

        }else{
          covmat[covmat %in% apply(comb, 1, paste, collapse = "," )] <- 1
          covmat[covmat != 1] <- 0
          cmat <- (covmat %>% as.data.frame %>% sapply(as.numeric)) * R2[i, j]}}

      if(is.null(dim(cmat))){
        cmat <- matrix(cmat, nrow = 1)}
      Cmat <- cbind(Cmat, cmat)}
    Sigma2 <- rbind(Sigma2, Cmat)

  }
  return(Sigma2)
}

# simulate data

# 500 patients
N = 500

allid <- deidentified.dat$id %>% unique
sampleid <- sample(allid, N)
sampledat <- deidentified.dat %>% filter(id %in% sampleid)

# generate design matrix
sdat <- NULL

library(splines)

for (i in 1:N){

  nsdat <- sampledat %>% filter(id %in% sampleid[i]); times <- nsdat$time
  t1 <- times[!is.na(nsdat$FVC)]; t2 <- times[!is.na(nsdat$DLCO)]
  t3 <- times[!is.na(nsdat$EF)]; t4 <- times[!is.na(nsdat$RVSP)]

  # random effects design matrix
  Zi_1 <- data.frame(int = 1, time = t1) %>% as.matrix()
  Zi_2 <- data.frame(int = 1, time = t2) %>% as.matrix()
  Zi_3 <- data.frame(int = 1, time = t3) %>% as.matrix()
  Zi_4 <- data.frame(int = 1, time = t4) %>% as.matrix()

  datlistZ <- list(Zi_1, Zi_2, Zi_3, Zi_4)
  Zp <- Reduce(direct.sum, datlistZ)

  Xi_1 <- cbind(1, ns(t1, knots = c(10, 30), Boundary.knots= c(0, 40)))
  Xi_2 <- cbind(1, ns(t2, knots = c(10, 30), Boundary.knots= c(0, 40)))
  Xi_3 <- cbind(1, ns(t3, knots = c(10, 30), Boundary.knots= c(0, 40)))
  Xi_4 <- cbind(1, ns(t4, knots = c(10, 30), Boundary.knots= c(0, 40)))

  Xp <- Reduce(direct.sum, list(Xi_1, Xi_2, Xi_3, Xi_4))

  Sigma2 <- create_sigma(datlist = datlistZ, R2 = R2)

  XbetaZb <- Xp %*% beta

  value <- mvrnorm(mu = XbetaZb, Sigma = Zp %*% G2 %*% t(Zp) + Sigma2)

  varvec <- rep(variablenames, times= c(length(t1), length(t2), length(t3), length(t4)))
  timevec <- c(t1, t2, t3, t4)

  newsdat <- data.frame(id = i, varvec, timevec, value)
  sdat <- rbind(sdat, newsdat)
}

ggplot(sdat, aes(timevec, value, group = id)) + geom_line(alpha = 0.5) + facet_wrap(~varvec)



