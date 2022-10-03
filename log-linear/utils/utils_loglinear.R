# ###################################################################################
# Simulate from log-linear model
# ####################################################################################

sim.loglinear<-function(d, A, B, lam.init, T){
	BurnIn=100
	n=nrow(A)
	nu=matrix(0,n,(T+BurnIn))
	lam=matrix(0,n,(T+BurnIn))
	X=matrix(0,n,(T+BurnIn))
	X[,1]=matrix(rpois(n, lam.init),n,1)
    nu[,1]=log(lam.init)
     
	for(t in 2:(T+BurnIn)){
	nu[,t]=d+A%*%nu[,t-1]+B%*%log(X[,t-1]+1)	
	lam[,t]=exp(nu[,t])	
	X[,t]=matrix(rpois(n, lam.init),n,1)
	}
	return(list(X=X[,(BurnIn+1):(T+BurnIn)],nu=nu[,(BurnIn+1):(T+BurnIn)],lam=lam[,(BurnIn+1):(T+BurnIn)]))
}


loglikpois.ts=function(theta,X,mod,lambda.init){
  T=ncol(X)
  n=nrow(X)
  
  loglik=rep(0,T)
  
  if (mod$corr) {
    m=n*(n-1)/2
    ntheta=length(theta)
    lambda=get.intensity(X,mod,theta[1:(ntheta-m)],lambda.init)
    value=1.0
    for (t in (mod$p+1):T) {
      result=.C("gpoisson",
                as.integer(X[,t]),
                as.integer(rep(0,m)),
                as.integer(n),
                as.integer(m),
                as.double(c(lambda[,t],theta[(ntheta-m+1):ntheta])),
                value=as.double(value))
      loglik[t]=result$value
    }
  } else {
    lambda=get.intensity(X,mod,theta,lambda.init)
    for (t in (mod$p+1):T) {
      loglik[t]=sum(dpois(X[,t],lambda[,t],log=TRUE))
    }
  }
  return(-sum(loglik[(mod$p+1):T]))
}

# ###################################################################################
# Forecast log-linear model
# ####################################################################################

forecast.loglinear<-function(d, A, B, lam.init, T.forecast){
  n=nrow(A)
  nu=matrix(0,n,(T.forecast))
  lam=matrix(0,n,(T.forecast))
  X=matrix(0,n,(T.forecast))
  X[,1]=x.forecast[,1]
  nu[,1]=log(lam.init)
  for(t in 2:(T+BurnIn)){
    nu[,t]=d+A%*%nu[,t-1]+B%*%log(X[,t-1]+1)	
    lam[,t]=exp(nu[,t])	
    X[,t]=matrix(rpois(n, lam.init),n,1)
  }
  return(list(X=X[,(BurnIn+1):(T+BurnIn)],nu=nu[,(BurnIn+1):(T+BurnIn)],lam=lam[,(BurnIn+1):(T+BurnIn)]))
}

# ###################################################################################
# Get intensity
# ####################################################################################

get.intensity = function(X,mod,theta,lambda.init) {
  
  T = ncol(X)
  n = nrow(X)
  
  w=rep(0,n)
  A=matrix(0,n,mod$p)
  B=matrix(0,n,n*mod$q)
  struct=list(w=w,A=A,B=B)
  struct=relist(theta,struct)
  w=struct$w
  A=matrix(apply(struct$A,2,diag),n,n*mod$p)
  B=struct$B
  
  lambda=matrix(NA,n,T)
  if (mod$type=="loglinear") {
    lambda[,1:mod$p] = log(lambda.init)
  } else {
    lambda[,1:mod$p] = lambda.init
  }
  for (t in (mod$p+1):T) {
    r=min(t-1,mod$q)
    if (mod$type=="loglinear") {
      Y=log(X[,(t-1):(t-r)]+1)
    } else {
      Y=X[,(t-1):(t-r)]
    }
    lambda[,t]=w+A%*%lambda[,(t-1):(t-mod$p)]+B[,1:(r*n)]%*%Y
  }
  if (mod$type=="loglinear") lambda=exp(lambda)
  return(lambda)
}
