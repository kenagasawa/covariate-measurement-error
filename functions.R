## setup design matrix to estimate control function V
CFMat_setup <- function(df){
  D   <- unlist(df[,"SATmed72"])/100
  Dp  <- bs((D-mean(D))/sd(D),df=5)
  
  ## proxy variables
  X   <- cbind(poly(df[,"Math_scale"],degree=3,raw=TRUE),
               poly(df[,"Reading_scale"],degree=3,raw=TRUE),
               poly(scale(df[,"SATmeanApplied"]),degree=2,raw=TRUE),
               poly(scale(df[,"Work_score"]),degree=2,raw=TRUE) )
  ## continuous excluded variables
  Z.cont <- c("SAT","Family_score")
  tmp    <- scale(do.call(cbind, lapply(Z.cont,function(x) bs(df[,x],df=5) ) )  )
  
  ## all the excluded variables
  Z   <- cbind(df[,"FatherCollege"],df[,"MotherCollege"],tmp,
               df[,"FatherCollege"]*tmp)
  
  ## other controls
  W   <- model.matrix(~( df[,"Female"] + df[,"minority"]+
                           scale(df[,"logIncome"]) )^2)
  W   <- W[,-1]
  
  N   <- length(D)
  
  lenw <- ncol(W)
  lenz <- ncol(Z)
  WZ  <- matrix(0,nrow=N,ncol=lenw*lenz)
  Z[,apply(Z,2,function(x) sum(is.na(x)))!=0] <- 0
  
  for (i in seq(lenz)){
    WZ[,((i-1)*lenw+1):(i*lenw)] <- Z[,i]*W
  }
  WZ[,apply(WZ,2,function(x) sum(is.na(x)))!=0] <- 0
  
  WD    <- do.call(cbind,lapply(seq(ncol(Dp)),function(x) Dp[,x]*W ) )
  ZD    <- do.call(cbind,lapply(seq(ncol(Dp)),function(x) Dp[,x]*Z ) )
  Rv    <- cbind(1,Dp,W,Z,WZ,WD,ZD)
  
  out  <- list(Y=df[,"log_earn"],D=D,W=W,Rv=Rv,X=X)
  return(out)
}

## setup design matrix for the main regression
d.icr  <- 1

CFMatY <- function(df,delta,lenv){
  dfv  <- 5
  dm   <- CFMat_setup(df)
  Y    <- dm$Y
  W    <- dm$W
  D    <- dm$D
  D0   <- (D-mean(D))/sd(D)
  D0d  <- (D+d.icr-mean(D))/sd(D)
  N    <- length(Y)
  ## compute V
  Vhat <- matrix(0,nrow = N,ncol = lenv)
  for (vi in seq(lenv)){
    Vhat[,vi] <- dm$Rv%*%delta[,vi]
  }
  V   <- scale(Vhat)
  
  ## setup regressors
  check.na <- apply(V,2,function(x) sum(is.na(x)) )
  if ( sum(check.na>0) >0 ){
    V   <- V[,-which(check.na>0)]
    lenv   <- lenv - sum(check.na>0)
  }
  Vp   <- do.call(cbind, lapply(seq(lenv),function(x) bs(V[,x],df=dfv) ) )
  
  tmp  <- do.call(cbind, lapply( split(t(Vp),seq(ncol(Vp))), function(x) x*W ) )
  WV   <- cbind(W,Vp,tmp)
  
  WVD   <- D0*WV
  WVDd  <- D0d*WV
  Ry    <- cbind(1,D0,WV,WVD) 
  Ryd   <- cbind(1,D0d,WV,WVDd) 
  
  
  #setup array of regressors to estimate derivatives wrt V
  DRy  <- array(0,dim = c(N,ncol(Ry),lenv))
  DRyd <- array(0,dim = c(N,ncol(Ry),lenv))
  
  for (vid in seq(lenv)){
    tmp        <- bs(V[,vid],df=dfv)
    bsknots    <- sort( c(rep(range(V[,vid]),4),attr(tmp,"knots") ) )
    derV       <- matrix(0,nrow=N,ncol=ncol(Vp))
    derV[,((vid-1)*dfv+1):(vid*dfv)] <- splineDesign(bsknots,V[,vid],derivs=1)[,-1,drop=FALSE]
    tmpV       <- do.call(cbind, lapply( split(t(derV),seq(ncol(Vp))), function(x) x*W ) )
    derWV      <- cbind(matrix(0,nrow=N,ncol=ncol(W)), derV,tmpV)
    derWVD     <- D0*derWV
    derWVDd    <- D0d*derWV
    DRy[,,vid] <- cbind(matrix(0,nrow=N,ncol=2),derWV,derWVD)
    DRyd[,,vid]<- cbind(matrix(0,nrow=N,ncol=2),derWV,derWVDd)
  }
  
  out <- list(Ry=Ry,Ryd=Ryd,DRy=DRy,DRyd=DRyd,V=Vhat,X=dm$X)
  return(out)
}

########################
#### ML functions  #####
########################
## The codes below are directly borrowed from R codes of Chernozhukov, Newey, and Singh (2022, ECMA)
tol=1e-6
D_LB=0 #each diagonal entry of \hat{D} lower bounded by D_LB
D_add=.2 #each diagonal entry of \hat{D} increased by D_add. 0.1 for 0, 0,.2 otw
max_iter=20 #10

#Performs theoretical iteration to help tune regularization parameter
get_D <- function(B,mvec,rho_hat){
  n=nrow(B)
  p=ncol(B)
  
  df=matrix(0,p,n)
  for (i in 1:n){
    df[,i]=B[i,]*as.vector(rho_hat %*% B[i,]) - mvec[i,]
  }
  df=df^2
  D2=rowMeans(df)
  
  D=sqrt(D2)
  return(D) #pass around D as vector
}

#Euclidean Norm
two.norm <- function(x){
  return(sqrt(x %*% x))
}

#Calculates ML function for alpha and gamma
RMD_stable<-function(mv,nv,B,p0,is_alpha = TRUE,lambda){
  ###Inputs
  ###M,N,B: matrices
  ###p0: initial parameter for p
  ###is_alpha: True = return estimate of Riez representor, False = return estimate of regression function
  ###Output
  ###rho_hat or beta_hat, a vector to estimate alpha or gamma respectively by taking the dot product with b
  k=1
  
  p=ncol(B)
  n=nrow(B)
  M=colMeans(mv)
  N=colMeans(nv)
  G=t(B)%*%B/n
  # low-dimensional moments
  M_hat0=M[1:p0]
  N_hat0=N[1:p0]
  G_hat0=G[1:p0,1:p0]
  
  # initial estimate
  rho_hat0=solve(G_hat0,M_hat0)
  rho_hat=c(rho_hat0,rep(0,p-ncol(G_hat0)))
  beta_hat0=solve(G_hat0,N_hat0)
  beta_hat=c(beta_hat0,rep(0,p-ncol(G_hat0)))
  
  # moments
  M_hat=M
  N_hat=N
  G_hat=G
  
  # # penalty
  # lambda=c*qnorm(1-alpha_lasso/(2*p))/sqrt(n)
  
  if(is_alpha){ 
    ###########
    # alpha_hat
    ###########
    indx <- (M == 0)
    diff_rho=1
    #Loop through max_iter times or until the change in rho between iterations is less than tol
    while(diff_rho>tol & k<=max_iter){
      
      # previous values
      rho_hat_old=rho_hat+0
      
      # normalization
      D_hat_rho=get_D(B,mv,rho_hat_old)
      D_hat_rho=pmax(D_LB,D_hat_rho)
      D_hat_rho=D_hat_rho+D_add
      
      rho_hat=RMD_lasso(M_hat, G_hat, D_hat_rho,lambda,indx)$coefficients
      
      
      # difference
      diff_rho=two.norm(rho_hat-rho_hat_old)
      k=k+1
    }
    return(rho_hat)
    
  } else { 
    ###########
    # gamma_hat
    ###########
    diff_beta=1
    #Loop through max_iter times or until the change in rho between iterations is less than tol
    while(diff_beta>tol & k<=max_iter){
      
      # previous values
      beta_hat_old=beta_hat+0
      
      # normalization
      D_hat_beta=get_D(B,nv,beta_hat_old)
      D_hat_beta=pmax(D_LB,D_hat_beta)
      D_hat_beta=D_hat_beta+D_add
      
      beta_hat=RMD_lasso(N_hat, G_hat, D_hat_beta, lambda)$coefficients
      
      # difference
      diff_beta=two.norm(beta_hat-beta_hat_old)
      k=k+1
      
    }
    return(beta_hat)
    
  }
}


l.dtzg=0.1


#Lasso learning algorithm to tune the regularization paremeter of our model
RMD_lasso <- function(M, G, D, lambda=0,index=NULL,control = list(maxIter = 1000, optTol = 10^(-5), 
                                                                  zeroThreshold = 10^(-6)), beta.start = NULL) {
  
  p <- ncol(G)
  if (is.null(index)) index <- rep(0,p)
  Gt<-G
  Mt<-M
  l <- 0.1
  L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,...)
  lambda_vec=lambda*L*D #v3: insert D here
  
  if (is.null(beta.start)) {
    beta <- rep(0,p) #vs low-dimensional initialization
  }
  else {
    beta <- beta.start
  }
  wp <- beta
  mm <- 1
  if (is.null(index)){
    indexset <- 1:p
  } else{
    indexset    <- (1:p)[!index]
    beta[index] <- 0
  }
  while (mm < control$maxIter) {
    beta_old <- beta
    for (j in indexset) {
      rho=Mt[j]-Gt[j,]%*%beta+Gt[j,j]*beta[j]
      z=Gt[j,j]
      
      if (sum(is.na(rho)) >= 1) {
        beta[j] <- 0
        next
      }
      if (rho < -1 * lambda_vec[j]) 
        beta[j] <- (rho+lambda_vec[j])/z
      if (abs(rho) <= lambda_vec[j]) 
        beta[j] <- 0
      if (rho > lambda_vec[j]) 
        beta[j] <- (rho-lambda_vec[j])/z
    }
    wp <- cbind(wp, beta)
    if (sum(abs(beta - beta_old), na.rm = TRUE) < control$optTol) {
      break
    }
    mm <- mm + 1
  }
  w <- beta
  w[abs(w) < control$zeroThreshold] <- 0
  return(list(coefficients = w, coef.list = wp, num.it = mm))
}
