#setwd("/Users/Yuliya/Dropbox/SV-Counts/Poisson")

###########################################################################
###########################################################################
mpois.ts=function(T,w,A,B,rho,lambda.init,type=c("loglinear","linear")) {
        
  type=match.arg(type)
  
  n = length(w)
  if (class(A)=="matrix") {
    p=1
    if (nrow(A)!=ncol(A)) stop("A must be a square matrix")
    if (nrow(A)!=n) stop("A and w are not compatible")
  } else if (class(A)=="array") {
    p=dim(A)[3]
    if (any(dim(A)[1:2]!=n)) stop("A and w are not compatible")
    dim(A)=c(n,n*p)
  } else {
    stop("A must be a matrix or array")
  }
  if (class(B)=="matrix") {
    q=1
    if (nrow(B)!=ncol(B)) stop("B must be a square matrix")
    if (nrow(B)!=n) stop("B and w are not compatible")
  } else if (class(B)=="array") {
    q=dim(B)[3]
    if (any(dim(B)[1:2]!=n)) stop("B and w are not compatible")
    dim(B)=c(n,n*q)
  } else {
    stop("B must be a matrix or array")
  }
  if (!missing(rho)) {
    if (any(dim(rho)!=n))
      stop("dimensions of rho are not compatible with w, A, B")
  }
  X=matrix(NA,n,T)
  lambda=matrix(NA,n,T)
  e=matrix(0,n,n)

  # initialise observations 
  
  for (t in 1:p) {
    X[,t] = rpois(n,lambda.init)
  }
  
  if (type=="loglinear") {
    lambda.init=log(lambda.init)
  }
  lambda[,1:p]=lambda.init

  for (t in (p+1):T) {
    r=min(t-1,q)
    if (type=="loglinear") {
      Y=log(X[,(t-1):(t-r)]+1)
    } else {
      Y=X[,(t-1):(t-r)]
    }
    lambda[,t]=w+A%*%c(lambda[,(t-1):(t-p)])+B[,1:(n*r)]%*%c(Y)
    if (!missing(rho)) {
      e[,]=0
      e[lower.tri(e,FALSE)]=rpois(n*(n-1)/2,rho[lower.tri(rho,FALSE)])
      Z=apply(e+t(e),2,sum)
    }
    if (type=="loglinear") {
      X[,t]=rpois(n,exp(lambda[,t]))+Z
    } else {
      X[,t]=rpois(n,lambda[,t])+Z
    }
    #e=rnorm(2,0,1)
    #e[2]=rho[1,2]*e[1]+sqrt(1-rho[1,2]^2)*e[2]
    #e=e-1
    #if (type=="loglinear") {
    #  X[,t]=rpois(n,exp(lambda[,t]+e))
    #} else {
    #  X[,t]=rpois(n,lambda[,t]*exp(e))
    #}
    #print(lambda[,t])
  }
  if (type=="loglinear") {
    lambda=exp(lambda)
  }

  result=list(X=X,lambda=lambda)
  return(result)
}
###########################################################################
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
###########################################################################
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
###########################################################################
# recompile gpoisson.c and load the resulting shared library
dyn.load("gpoisson.so")
###########################################################################
# 
# w=c(0.2,0.2)
# A=diag(c(0.5,0.5))
# B=matrix(c(0.35,0,+0.2,+0.4),2,2)
# rho=matrix(c(0,2,2,0),2,2)
# T=600
# lambda.init=c(5,5)
# X=mpois.ts(T,w,A,B,rho,lambda.init,"loglinear")
# w=c(1,1)
# # Y=mpois.ts(T,w,A,B,rho,lambda.init,"linear")
# # par(mfrow=c(2,2))
# # plot(X$lambda[1,],type="l",ylim=c(0,max(X$lambda[1,])+5))
# # points(X$X[1,],col=2)
# # plot(X$lambda[2,],type="l",ylim=c(0,max(X$lambda[2,])+5))
# # points(X$X[2,],col=2)
# # plot(Y$lambda[1,],type="l",ylim=c(0,max(Y$lambda[1,])+5))
# # points(Y$X[1,],col=2)
# # plot(Y$lambda[2,],type="l",ylim=c(0,max(Y$lambda[2,])+5))
# # points(Y$X[2,],col=2)
# # cor(t(X$X))
# # cor(t(Y$X))
# 
# I=101:T
# w=c(0.2,0.2)
# mod=list(p=1,q=1,type="loglinear",corr=TRUE)
# rho=1#rho[lower.tri(rho,diag=FALSE)]
# theta=list(w=w,A=diag(A),B=B,rho=rho)
# theta=unlist(theta)
# lambda=get.intensity(X$X[,I],mod,theta[1:(length(theta)-length(rho))],
#   X$lambda[,101])
# par(mfrow=c(2,1))
# plot(X$lambda[1,I],type="l")
# lines(lambda[1,],col=2)
# plot(X$lambda[2,I],type="l")
# lines(lambda[2,],col=2)
# 
# loglikpois.ts(theta,X$X[,101:T],mod,lambda.init)
# # 
# # theta=list(w=c(0.1,0.1),A=c(0,0),B=matrix(c(0.5,0,0,0.5),2,2),rho=0.2)
# # theta=unlist(theta)
# # res=optim(par=theta,fn=loglikpois.ts,gr=NULL,method="Nelder-Mead",
# #   X=X$X[,101:T],mod=mod,lambda.init=lambda.init,hessian=TRUE,
# #   control=list(parscale=rep(0.1,length(theta))))
# # CI=matrix(2*length(theta),length(theta),2)
# # std.err=sqrt(diag(solve(res$hessian)))
# # CI[,1]=res$par-2*std.err
# # CI[,2]=res$par+2*std.err
# # CI
# # 
# # theta=list(x=c(0.1,0.1),A=c(0,0),B=matrix(c(0.5,0,0,0.5),2,2),rho=0.5)
# # theta=unlist(theta)
# # res=optim(par=theta,fn=loglikpois.ts,gr=NULL,method="Nelder-Mead",
# #   X=X$X[,101:T],mod=mod,lambda.init=lambda.init,hessian=TRUE,
# #   control=list(parscale=c(rep(0.1,length(theta))),maxit=5000))
# # std.err=sqrt(diag(solve(res$hessian)))
# # CI=matrix(2*length(theta),length(theta),2)
# # CI[,1]=res$par-2*std.err
# # CI[,2]=res$par+2*std.err
# # CI
# # 
# # res
