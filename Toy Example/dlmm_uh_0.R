
#### DLMM for UnitedHealth covid-19 data, 20200609
# rm(list=ls())
require(data.table)
require(lme4)
require(nlme)
require(mvtnorm)

# setwd('mydir')
source('dlmm-engine.R')

####################  toy example to show 'lossless' property  ##################
# set.seed(123)

K <- 20                       # number of sites
nn <- rep(c(100, 200), K/2)    # size of each site
px <- 5                   # covariates
pz <- px+1                     # assume random effects for all X + intercept
b <- c(1:px)/px*2              # true fixed-effect
# s2 <- c(1:K) / K               # random error variance
s2 <- rep(0.5, K)
# V <- diag(c(0:px)/px) / 10          # covariance of random effects, assume independence here
V <- diag(c(rep(2,2+px/2), rep(5,px/2))) / 10

X <- matrix(NA, sum(nn), px)
Z <- matrix(NA, sum(nn), pz)
id.site <- as.character(rep(1:K, nn))
Y <- rep(NA, sum(nn))
u <- matrix(NA, K, pz)
for(ii in 1:K){
  X[id.site==ii,] <- matrix(rbinom(nn[ii]*px, size=1, prob = 0.3), nn[ii], px)
  # matrix(rnorm(nn[ii]*px), nn[ii], px)
  Z[id.site==ii,] <- cbind(1, X[id.site==ii,])
  ui <- c(rmvnorm(1, rep(0,pz), V))
  u[ii,] <- ui
  ei <- rnorm(nn[ii], 0, sqrt(s2[ii]))
  Y[id.site==ii] <- X[id.site==ii,]%*%b + Z[id.site==ii,] %*%ui + ei
}
X <- cbind(1, X)


# sim_data
sim_data = as.data.frame(cbind(id.site, X[,-1], Y))
names(sim_data) = c("id.site", paste0("X", 1:px), "Y")
sim_data[, 1] <- as.factor(sim_data[, 1])
sim_data[, -1] <- lapply(sim_data[, -1], as.numeric)

# lmer fit
# fit00 <- lmer(Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +
#                 (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10|id.site), data = sim_data)
fit00 <- lmer(Y ~ X1 + X2 + X3 + X4 + X5 +
                (X1 + X2 + X3 + X4 + X5|id.site), data = sim_data)

## assume pooled data are all available, ~ 1 min
fit0 <- lmm.fit(Y, X, Z, id.site, pooled=T, reml=T, common.s2=T, hessian = T)
diag(fit0$V)
cbind(fit0$b, summary(fit00)$coef[,1])

VarCorr(fit00)

tmp = summary(fit00)

# fit0$res.profile$allterms$remlterm.i
# fit0$wald
# fit0$ui
# fit0$Yihat

## get summary stats from each site
SiXYZ <- lmm.get.summary(Y, X, Z=X, weights = NULL, id.site)

## use the summary stats to reconstruct the LMM likelihood, ~ 1 min
fit1 <- lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=T, common.s2=T, hessian=T)


as.data.frame(VarCorr(fit00))

sqrt(diag(fit1$V))
sqrt(fit1$s2)

lmm.sigma = c(as.data.frame(VarCorr(fit00))[1:(px+1),5], tail(as.data.frame(VarCorr(fit00))[,5],1))
dlmm.sigma = c(sqrt(diag(fit1$V)), sqrt(fit1$s2))


lmm.sigma
dlmm.sigma
dlmm.sigma - lmm.sigma


## plot to show estimates from pooled and distributed are identical
plot(fit0$b, fit1$b, main='Linear mixed model: fixed-effects estimation',
     xlab='LMM fixed-effects (pooled)',
     ylab='LMM fixed-effects (distributed)')
abline(0,1)

####################  END: toy example ##############################







# data format:
# Y: vector of outcome, continuous numeric
# X: dataframe of covariates, all numeric
# id.site: vector of site id, please make it as.character()

## all data in data.frame
# X <- data.frame(X[,-1])
dd <- data.frame(y=Y, id.site=id.site, data.frame(X[,-1]))
# px <- ncol(X)
nn <- table(id.site)
xn <- names(data.frame(X[,-1]))

fit.lm <- lm(Y~ as.matrix(X))


## LRT test if any single random effect is significant.
## random effect of covariate m at site i = u_im ~ N(0, vm), i=1...K, m=1...q
## for each m, test H0: vm=0 vs H1: vm>0
## this is LRT on boundary, LR ~ 0.5*chisq(df=0) + 0.5*chisq(df=1),
LR <- matrix(NA, px, 2)
# LMM under H0
fit0 <- lmer(as.formula(paste0('y~', paste0(xn, collapse ='+'), '+(1|id.site)')), data = dd)
for(ii in 1:px){
  cat('covariate #', ii, '=', xn[ii], '...')
  # LMM under H1
  fit1i <- lmer(as.formula(paste0('y ~', paste0(xn, collapse ='+'), '+(1|id.site)+', '(0+', xn[ii], '|id.site)')), data = dd)
  LR[ii,1] <- 2*(summary(fit1i)$logLik -  summary(fit0)$logLik)
  LR[ii,2] <- 0.5 * (1-pchisq(LR[ii,1], df=1))
  if(LR[ii, 2] < 0.05/px) cat('pval =', LR[ii, 2], ', R.E. needed \n')
}

## Z: matrix of covariates of random-effects, intercept + significant Xs
Z <- as.matrix(cbind(1, X[, which(LR[,2] < 0.05/px)]))
X <- as.matrix(cbind(1, X))


## LMM by lme4::lmer
# fit1 <- lmer(as.formula(paste0('y~', paste0(xn, collapse ='+'), '+(1|id.site)', paste0(paste0('+(0+', xn[which(LR[,2] < 0.05/px)], '|id.site)'),
#                                                                                        collapse = ''))), data=dd)
fit1 <- lmer(as.formula(paste0('y~', paste0(xn, collapse ='+'), '+(1|id.site)', paste0(paste0('+(0+', xn, '|id.site)'),
                                                                                       collapse = ''))), data=dd)

# fixed-effect
b0 <- summary(fit1)$coef[,1]
# var components
# plot(diag(V), data.frame(summary(fit1)$varcor)$vcov[1:pz]); abline(0,1)


## LMM using my alg
## assume pooled data are all available
fit.lmm00 <- lmm.fit(Y, X, Z, id.site, pooled=T, common.s2 = F, hessian=T)  # assume different residual var
fit.lmm01 <- lmm.fit(Y, X, Z, id.site, pooled=T, common.s2 = T, hessian=T)  # assume common residual var

## get summary stats from each site, DLMM using summary stats to reconstruct the likelihood
SiXYZ <- lmm.get.summary(Y, X, Z, weights = NULL, id.site)
fit.lmm10 <- lmm.fit(SiXYZ = SiXYZ, pooled=F, common.s2 = F, hessian=T)

fit.lmm11 <- lmm.fit(SiXYZ = SiXYZ, pooled=F, common.s2 = T, hessian=T)
fit.lmm11
## fixed-effect very close
cbind(c(0,b), b0, fit.lmm01$b, fit.lmm11$b)

## still testing the variance est of fixed effects
# sd.lmm11 <- sqrt(diag(solve(fit.lmm10$res.profile$allterms$bterm1)))
# sd.1 <- sqrt(diag(summary(fit1)$vcov))
# round(cbind(sd.lmm11, sd.1), 3)

tmp=summary(fit1)

plot(b0, fit.lmm11$b, main='Linear mixed model: fixed-effects estimation',
     xlab='LMM fixed-effects (lme4)',
     ylab='LMM fixed-effects (distlmm)')
abline(0,1)


res <- list(SiXYZ=SiXYZ,
            fit.lmm00 = fit.lmm00,
            fit.lmm01 = fit.lmm01,
            fit.lmm10 = fit.lmm10,
            fit.lmm11 = fit.lmm11)

res$px = px
res$nn = nn
res$xn = xn
res$fit.lm = fit.lm
res$LR = LR


## can you also make some summary stats e.g. a "Table 1" of mean/sd or proportion of covariates and outcome

# res$table1 <- xxx



saveRDS(res, file='dlmm_UH_20200609.rds')



