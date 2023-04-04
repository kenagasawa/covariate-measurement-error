### main file
rm(list=ls())
library("dplyr") 
library("splines")

source("functions.R")

seed <- 12
set.seed(seed)

df <- as.data.frame(readRDS("Earnings_NLS72.rds"))
N     <- nrow(df)

# number of folds
L=5
index <- seq(L)
#Split up into I_1,...I_L
folds <- split(sample(N,N,replace=FALSE), as.factor(1:L))

## dim of V and # of regressors for the main regression
tmp   <- CFMat_setup(df)
lenv  <- ncol(tmp$X)
delt  <- matrix(0.3,nrow = ncol(tmp$Rv),ncol=lenv)
tmp   <- CFMatY(df,delt,lenv)
lenp  <- ncol(tmp$Ry)

# vector to store estimates
psi <- numeric(0)

## lambda constants were chosen by cross validation
lambdaV <- 0.5*qnorm(1-0.1/(2*lenp))/sqrt(N)
lambdaY <- 0.1*qnorm(1-0.1/(2*lenp))/N^(1/3)
lambdaA <- 0.1*qnorm(1-0.1/(2*lenp))/N^(1/3)


## cross-fitting
for (l in index){
  id.l  <- folds[[l]]
  
  nvecl <- matrix(0,nrow=N,ncol=lenp)
  mvecl <- matrix(0,nrow=N,ncol=lenp)
  Regma <- matrix(0,nrow=N,ncol=lenp)
  for (ll in index[-l]){
    id.ll  <- folds[[ll]]
    idn   <- unlist( lapply(index[-c(l,ll)],function(x) folds[[x]] ) )
    df.nl <- df[idn,]
    df.ll <- df[id.ll,]
    
    ## estimate V
    dm.v  <- CFMat_setup(df.nl)
    Rv.nl <- dm.v$Rv
    
    deltahat <- matrix(0,nrow = ncol(Rv.nl),ncol=lenv)
    for (vi in seq(lenv)){
      nvec          <- dm.v$X[,vi]*Rv.nl
      deltahat[,vi] <- RMD_stable(nvec,nvec,Rv.nl,p0=5,0,lambda=lambdaV)
    }
    
    ## setup design matrix
    dm.y   <- CFMatY(df.ll,deltahat,lenv)
    Y.nl   <- df.ll[,"log_earn"]
    Ry.nl  <- dm.y$Ry
    Ryd.nl <- dm.y$Ryd
    
    mvecl[id.ll,]  <- Ryd.nl-Ry.nl
    nvecl[id.ll,]  <- Y.nl*Ry.nl
    Regma[id.ll,]   <- Ry.nl
  }
  ## estimation of E[Y|D,W,V]
  mvecl   <- mvecl[-id.l,]
  nvecl   <- nvecl[-id.l,]
  Regma   <- Regma[-id.l,]
  betahat <- RMD_stable(mvecl,nvecl,Regma,p0=5,0,lambda=lambdaY)
  
  ### estimate Riez representor
  rhohat <- RMD_stable(mvecl,nvecl,Regma,p0=5,1,lambda=lambdaA)
  
  ### Compute the influence function
  
  ## estimate Vl
  idnl    <- -folds[[l]]
  df.l    <- df[id.l,]
  df.nl   <- df[idnl,]
  
  dm.v  <- CFMat_setup(df.nl)
  Rv.nl <- dm.v$Rv
  deltahat <- matrix(0,nrow = ncol(Rv.nl),ncol=lenv)
  for (vi in seq(lenv)){
    nvec          <- dm.v$X[,vi]*Rv.nl
    deltahat[,vi] <- RMD_stable(nvec,nvec,Rv.nl,p0=5,0,lambda=lambdaV)
  }
  
  # influence function
  dm.l    <- CFMatY(df.l,deltahat,lenv)
  Y.l     <- df.l[,"log_earn"]
  Ry.l    <- dm.l$Ry
  Ryd.l   <- dm.l$Ryd
  DRy.l   <- dm.l$DRy
  DRyd.l  <- dm.l$DRyd
  X.l     <- dm.l$X
  V.l     <- dm.l$V
  N.l     <- nrow(Ry.l)
  
  mu.l    <- Ry.l%*% betahat
  mu_d.l  <- Ryd.l%*% betahat
  alpha.l <- as.numeric(Ry.l%*%rhohat)
  
  
  ## second correction term
  Dmu   <- matrix(0,nrow=N.l,ncol=lenv)
  Dmud <- matrix(0,nrow=N.l,ncol=lenv)
  term1 <- matrix(0,nrow=N.l,ncol=lenv)
  term2 <- matrix(0,nrow=N.l,ncol=lenv)
  for (vi in seq(lenv)){
    # derivative estimate
    Dmu[,vi]  <- as.numeric(DRy.l[,,vi]%*%betahat)
    Dmud[,vi]<- as.numeric(DRyd.l[,,vi]%*%betahat)
    
    term1[,vi] <- (Dmud[,vi]-Dmu[,vi])*(X.l[,vi]-V.l[,vi])
    term2[,vi] <- -alpha.l*Dmu[,vi]*(X.l[,vi]-V.l[,vi])
  }
  
  
  corr.term <- rowSums(term1 + term2)
  ###
  psi.l <- (mu_d.l - mu.l) + alpha.l*(Y.l-mu.l) + corr.term
  psi   <- c(psi,psi.l)
  
}

## point estimate and standard error
est <- mean( psi )
varest <- sum( (psi -est)^2 )/N
std.dev <- sqrt( varest/N )