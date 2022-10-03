setwd("/Users/yuliya/Dropbox/CountSSM/git/log-linear")
#dyn.load("gpoisson.so")
source("utils/utils_loglinear.R")

set.seed(10)

# Set parameters for data generation
T=600 # sample size
w=c(1,1) 
A=diag(c(0.5,0.5)) 
B=matrix(c(0.35,0,+0.2,+0.4),2,2) 
rho=matrix(c(0,0.2,0.2,0),2,2) 
lambda.init=c(5,5) 

# Generate data
Data = sim.loglinear(w, A, B, lambda.init, T) # X, nu, lam
X = Data$X # take the counts
# T = ncol(X)
n=nrow(X)

# Put parameters in a list
theta=list(w=w,A=diag(A),B=B,rho=rho)
theta=unlist(theta)

# Define the type of the model
mod=list(p=1,q=1,type="loglinear",corr=TRUE)
I=101:T

# Starting values for optimization
theta=list(w=c(0.1,0.1),A=c(0.5,0.5),B=matrix( c(0.8,0.2,0.2,0.8),2,2),rho=0.2)
theta=unlist(theta)

# Optimize 
res=optim(par=theta,fn=loglikpois.ts,gr=NULL,method="Nelder-Mead",
          X=X,mod=mod,lambda.init=lambda.init,hessian=TRUE,
          control=list(parscale=rep(0.1,length(theta))))

# Confidence Intervals
CI=matrix(2*length(theta),length(theta),2)
std.err=sqrt(diag(solve(res$hessian)))
CI[,1]=res$par-2*std.err
CI[,2]=res$par+2*std.err
CI

