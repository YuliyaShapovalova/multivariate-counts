#################################################################################
# Bivariate analysis of count data with state-space model
#################################################################################

# Load libraries
setwd("/Users/Yuliya/Dropbox/CountSSM/git/ssm")
library("xtable")
library("coda")
library("stats")

# Simulate data 
seed = 42
T = 1000
beta = c(1,1)
phi = matrix(c(0.5,0.3,0.0,0.5),2,2)
seta = matrix(c(0.5,0.3,0.3,0.5),2,2)
Data = simmulticount(T,beta,phi,seta,seed) # x-- counts, lambda -- intensity, h -- latent process

X = Data$x
T=ncol(X)
n=nrow(X)


# Run PMCMC
N=2000 # Number of particles in SMC
M=1000   # Number of MCMC steps; M>20000 recommnded for inference; take smaller number to test if the code is running (e.g., M=10)
initialpars=c(10,10,0.7,0.0,0.0,0.7,5,5,0.8) # initial values for the parameters: beta, phi, seta
const=1 # Calibration step of the proposal variance-covariance matrix
Sig_proposal=const*c(0.95,0.95,0.00912,0.0223,0.7*0.00223,0.7*0.000912,0.1*0.01,0.1*0.01,0.7*0.01) # The steps need to be calibrated to achieve optimal performance

results_mcmc=PMH_corr(T,n,N,M,initialpars,Sig_proposal,X) # thetamh -- parameter samples
#results_mcmc=PMH_corr_inR(T,n,N,M,initialpars,Sig_proposal,X) # pure R code can be used for testing -- very slow

# Acceptance rate, optimal should be between 25 and 40%; for optimal calibration check also ACF, trace plots etc.
for(i in 1:9){      
        print(accrate(results_mcmc$thetamh[,i]))
}

# Save results
save(results_mcmc, file="results_mcmc.RData")


# Compute summary statistics based on the results
BurnIn=M*0.2 # 20% of the samples are discarded for warm-up
thetamh=results_mcmc$thetamh
HPD<-HPDinterval(mcmc(thetamh[BurnIn:M,]))

tb<-cbind(apply(thetamh[BurnIn:M,],2,mean),apply(thetamh[BurnIn:M,],2,median),
          apply(thetamh[BurnIn:M,],2,Mode),HPD[,1],HPD[,2])
xtable(tb,digits=4)

