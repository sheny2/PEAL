## load packages
require(pda)
require(lme4)

## sample data
?LOS
data(LOS)

## split the data to 3 separate sets (patient-level data)
LOS_split <- split(LOS, LOS$site)




X
Y
Z
id.site
function(Y = NULL, X = NULL, Z = NULL, weights = NULL, id.site = NULL){
  if(is.null(weights)) weights <- rep(1, length(Y))
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  pz <- ncol(Z)

  SiXYZ <- list()
  for(ii in seq_along(id.site.uniq)){
    si = id.site.uniq[ii]
    wti = weights[id.site == si]
    Xi <- X[id.site == si, ]
    Zi <- Z[id.site == si, ]
    Yi <- Y[id.site == si]
    # if(any(apply(Xi[,-1], 2, function(a)length(unique(a)))==1))
    #   warning(paste0('singular X in site #', ii, ' detected!'))
    # if(any(apply(Zi[,-1], 2, function(a)length(unique(a)))==1))
    #   warning(paste0('singular Z in site #', ii, ' detected!'))

    SiX  = t(Xi*wti) %*% Xi
    SiXZ = t(Xi*wti) %*% Zi
    SiXY = t(Xi*wti) %*% Yi
    SiZ  = t(Zi*wti) %*% Zi
    SiZY = t(Zi*wti) %*% Yi
    SiY  = sum(Yi ^ 2 *wti)
    ni <- sum(id.site == si)
    SiXYZ[[si]] <- list(SiX  = SiX, SiXZ = SiXZ, SiXY = SiXY,
                        SiZ  = SiZ, SiZY = SiZY, SiY  = SiY, ni = ni)
  }

  return(SiXYZ)
}



