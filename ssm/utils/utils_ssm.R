######################################################################################
# Functions for (bivariate) analysis of count data
# lambda(t)=exp(beta*h(t))
# phi*h(t-1)+eta(t)
# likfort2count - Fortran function for likelihood evaluation
# PMH_corr - particle Metropolis-Hastings with correlation parameter in the latent process
# pf - likelihood evaluation using R
# PMH_nocorr - particle Metropolis-Hastings without correlation in the latent process
######################################################################################
library("psych")

SigmaEtaAlpha=2.5  # prior hyperparameters
SigmaEtaBeta=0.025 # prior hyperparameters

setpar<-function(...) par(mar=c(2.5,2.5,0.5,0.5),mgp=c(1.4,0.2,0),
                          cex=0.7,bty="l",las=1,lwd=0.5,tck=0.01,...)

likfort2count<-function(n,t,l,x,beta,phi,seta){ 
        dyn.load('utils/likfort2count.so')
        #print(a)
        #likfin=0.0
        n=nrow(x)
        t=ncol(x)
        hpsmT=matrix(c(0),n,l)
        hpT=matrix(c(0),n,l)
        hptraj=matrix(c(0),n,T)
        hpsmtraj=matrix(c(0),n,T)
        hpmean=matrix(c(0),n,T)
        #hpall=matrix(c(0),n,l)
        results <- .Fortran("likfort2count",
                            as.integer(n),
                            as.integer(t),
                            as.integer(l),
                            as.integer(x),
                            as.double(beta),
                            as.double(phi),
                            as.double(seta),
                            likfin=as.double(1),
                            hpsmT=as.double(hpsmT),
                            hpT=as.double(hpT),
                            hptraj=as.double(hptraj),
                            hpsmtraj=as.double(hpsmtraj),
                            hpmean=as.double(hpmean))
        return(list(lik=results$likfin,hpsmT=results$hpsmT, hpT=results$hpT,
                    hptraj=results$hptraj, hpsmtraj=results$hpsmtraj, hpmean=hpmean))
        
}

dpoisson<-function(x, par){
        a=-log(factorial(x))
        b=x*log(par)        
        c=-par
        return(a+b+c)
}

accrate<-function(Likelihoodvector){
        N=length(Likelihoodvector)
        x=sum(ifelse(Likelihoodvector[-1]!=Likelihoodvector[-N],1,0))
        return(x/N)
}

resampleSystematic<-function(w,N){
        # [ indx ] = resampleSystematic[ w, N]
        # Systematic resampling method for particle filtering. 
        # Author: Tiancheng Li,Ref:
        #T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
        # submit to IEEE Signal Processing Magazine, August 2013
        
        # Input:
        #w    the input weight sequence 
        #       N    the desired length of the output sequence[i.e. the desired number of resampled particles]
        # Output:
        #indx the resampled index according to the weight sequence
        
        #if (nargin ==1){
        #N = length(w)
        #}
        M=N
        #M = length(w)
        w = w / sum(w)
        Q = cumsum(w)
        indx = matrix(c(0),1, N)
        T = seq(0,1-1/N,length=N) + runif(1)/N#
        
        i = 1
        j = 1
        while(i<=N && j<=M){
                while (Q[j] < T[i]){
                        j = j + 1
                }
                indx[i] = j
                i = i + 1
        }
        return(indx)
}

rmultnorm<-function(N,mu,S,Cholesky=FALSE){
        n=length(mu)
        Z=matrix(rnorm(n*N),n,N)
        if (Cholesky) {
                M=c(mu)+S%*%Z
        } else {
                M=c(mu)+t(chol(S))%*%Z
        }
        return(M)
}


pf<-function(N,n,theta,x){
        #marginal likelihood estimation with particle filter
        Ibeta=c(1:n)
        Iphi=(n+1):(n+n*n)
        Iseta=(n+n*n+1):(2*n+n*n)
        Irho=2*n+n*n+1
        T=ncol(x)
        hmean=matrix(c(0),n,T+1)
        
        beta=theta[Ibeta]
        phi=matrix(theta[Iphi],n,n)
        seta=diag(theta[Iseta])
        seta[1,2]=theta[Irho]
        seta[2,1]=theta[Irho]
        
        #hp=matrix(c(0),N,T+1)
        hp=array(c(0),dim=c(n,N,T+1))
        w=matrix(c(1),N,T+1)/N
        Lik=rep(0,T+1)
        for(t in 2:(T+1)){
                k=resampleSystematic(w[,t-1],N)
                hmean[,t]=apply(hp[,,t-1],1,mean)
                hp[,,t]=phi%*%hp[,k,t-1]+rmultnorm(N,rep(0,n),seta) 
                #w[,t]=dpoisson(sum(x[,t-1]),apply(theta[Ibeta]*exp(hp[,,t]),2,sum))
                w[,t]=apply(dpoisson(x[,t-1],theta[Ibeta]*exp(hp[,,t]), log=TRUE),2,sum)
                #rescale
                wmax=max(w[,t])
                w[,t]=exp(w[,t]-wmax)
                Lik[t]=Lik[t-1]+wmax+log(sum(w[,t]))-log(N)
                w[,t]=w[,t]/sum(w[,t])
                #hmean[,t]=sum(w[,t]*hp[,,t])
        }
        hpT=hp[,,T]
        return(list(lik=Lik[T+1],hmean=hmean,hpT=hpT))
}

# Partilce Metropolis-Hastings for parameters sampling with correlation
PMH_corr<-function(T,n,N,M,initialpars,Sig,x){
        thetamh=matrix(c(0),M,length(initialpars))
        logLik=rep(0,M)
        a=rep(0,M)
        hp_T=matrix(c(0),n,N)
        hp_T_allM=array(0,c(n,N,M))
        hp_trajM=array(0,c(n,T,M))
        thetamh[1,]=initialpars
        #logLik[1]=pf(N,n,thetamh[1,],x)
        
        Sig_beta=diag(Sig[1:n])
        Sig_phi1=diag(Sig[(n+1):(n+2)])
        Sig_phi2=diag(Sig[(n+3):(n+n*n)])
        Sig_seta=diag(Sig[(n+n*n+1):(2*n+n*n)])
        Sig_rho=Sig[2*n+n*n+1]
        
        Ibeta=c(1:n)
        Iphi=(n+1):(n+n*n)
        Iphi1=(n+1):(n+2)
        Iphi2=(n+3):(n+n*n)
        Iseta=(n+n*n+1):(2*n+n*n)
        Irho=2*n+n*n+1
        
        t=T
        l=N
        
        beta=initialpars[Ibeta]
        phi=matrix(initialpars[Iphi],n,n)
        seta=diag(initialpars[Iseta])
        rho=initialpars[Irho]
        seta[1,2]=rho
        seta[2,1]=rho
        estim=likfort2count(n,t,l,x,beta,phi,seta)
        logLik[1]=estim$lik
        hp_T=estim$hpT
        hp_T_allM[,,1]=matrix(hp_T,n,N)
        hp_trajM[,,1]=matrix(estim$hptraj,n,T)
        #logLik[1]=pf(N,n,initialpars,x)
        
        for(i in 2:M){
                thetac=thetamh[i-1,]
                thetac[Ibeta]=c(rmultnorm(1,c(thetamh[i-1,Ibeta]),Sig_beta))
                #################################
                #Proposal for the scale par beta
                #################################
                if(all(thetac[1:n]>0)){
                        #logLikC=pf(N,n,thetac,x)
                        beta=thetac[Ibeta]
                        phi=matrix(thetac[Iphi],n,n)
                        seta=diag(thetac[Iseta])
                        seta[1,2]=thetac[Irho]
                        seta[2,1]=thetac[Irho]
                        estim=likfort2count(n,t,l,x,beta,phi,seta)
                        logLikC=estim$lik
                        hp_T=estim$hpT
                        hp_T_allM[,,i]=matrix(hp_T,n,N)
                        hp_trajM[,,i]=matrix(estim$hptraj,n,T)
                        a[i]=exp(logLikC-logLik[i-1])*(PriorValue(thetac)/PriorValue(thetamh[i-1,])) 
                        #a[i]=exp(logLikC-logLik[i-1])
                }else{
                        logLikC=-Inf
                        a[i]=0
                }     
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                #a[i]=exp(logLikC-logLik[i-1])
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }else{
                        thetamh[i,]=thetamh[i-1,]
                        logLik[i]=logLik[i-1]
                }
                #################################################
                #Proposal for the persistence par phi, first row
                #################################################
                thetac=thetamh[i,]
                thetac[Iphi1]=c(rmultnorm(1,c(thetamh[i,Iphi1]),Sig_phi1))
                beta=thetac[Ibeta]
                phi=matrix(thetac[Iphi],n,n)
                seta=diag(thetac[Iseta])
                seta[1,2]=thetac[Irho]
                seta[2,1]=thetac[Irho]
                Phiev=eigen(phi)$values
                if(all(abs(Phiev)<0.999)){
                        #logLikC=pf(N,n,thetac,x)
                        #a[i]=exp(logLikC-logLik[i])
                        beta=thetac[Ibeta]
                        phi=matrix(thetac[Iphi],n,n)
                        seta=diag(thetac[Iseta])
                        seta[1,2]=thetac[Irho]
                        seta[2,1]=thetac[Irho]
                        estim=likfort2count(n,t,l,x,beta,phi,seta)
                        logLikC=estim$lik
                        hp_T=estim$hpT
                        hp_T_allM[,,i]=matrix(hp_T,n,N)
                        hp_trajM[,,i]=matrix(estim$hptraj,n,T)
                        a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
                }else{
                        logLikC=-Inf
                        a[i]=0
                }    
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                #a[i]=exp(logLikC-logLik[i])
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }
                
                #################################################
                #Proposal for the persistence par phi, second row
                #################################################
                thetac=thetamh[i,]
                thetac[Iphi2]=c(rmultnorm(1,c(thetamh[i,Iphi2]),Sig_phi2))
                beta=thetac[Ibeta]
                phi=matrix(thetac[Iphi],n,n)
                seta=diag(thetac[Iseta])
                seta[1,2]=thetac[Irho]
                seta[2,1]=thetac[Irho]
                Phiev=eigen(phi)$values
                if(all(abs(Phiev)<0.999)){
                        #logLikC=pf(N,n,thetac,x)
                        #a[i]=exp(logLikC-logLik[i])
                        beta=thetac[Ibeta]
                        phi=matrix(thetac[Iphi],n,n)
                        seta=diag(thetac[Iseta])
                        seta[1,2]=thetac[Irho]
                        seta[2,1]=thetac[Irho]
                        estim=likfort2count(n,t,l,x,beta,phi,seta)
                        logLikC=estim$lik
                        hp_T=estim$hpT
                        hp_T_allM[,,i]=matrix(hp_T,n,N)
                        hp_trajM[,,i]=matrix(estim$hptraj,n,T)
                        a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,]))
                }else{
                        a[i]=0
                        logLikC=-Inf
                }     
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                #a[i]=exp(logLikC-logLik[i])
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }
                
                
                ####################################################
                #Proposal for the variance par seta, variance block
                ####################################################
                thetac=thetamh[i,]
                thetac[Iseta]=c(rmultnorm(1,c(thetamh[i,Iseta]),Sig_seta))
                beta=thetac[Ibeta]
                phi=matrix(thetac[Iphi],n,n)
                seta=diag(thetac[Iseta])
                seta[1,2]=thetac[Irho]
                seta[2,1]=thetac[Irho]
                Setaev=eigen(seta)$values
                if((seta[1,1]*seta[2,2]>seta[1,2]^2)&seta[1,1]>0.0001&seta[2,2]>0.0001){
                        #logLikC=pf(N,n,thetac,x)
                        estim=likfort2count(n,t,l,x,beta,phi,seta)
                        logLikC=estim$lik
                        hp_T=estim$hpT
                        hp_T_allM[,,i]=matrix(hp_T,n,N)
                        hp_trajM[,,i]=matrix(estim$hptraj,n,T)
                        a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
                        #a[i]=exp(logLikC-logLik[i])
                }else{
                        logLikC=-Inf
                        a[i]=0
                }       
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                #a[i]=exp(logLikC-logLik[i])*()
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }
                
                ################################################
                #Proposal for the correlation par in seta
                ################################################
                thetac=thetamh[i,]
                thetac[Irho]=c(rmultnorm(1,c(thetamh[i,Irho]),Sig_rho))
                seta=diag(thetac[Iseta])
                seta[1,2]=thetac[Irho]
                seta[2,1]=thetac[Irho]
                Setaev=eigen(seta)$values
                if((seta[1,1]*seta[2,2]>seta[1,2]^2)&seta[1,1]>0.0001&seta[2,2]>0.0001){
                        #logLikC=pf(N,n,thetac,x)
                        beta=thetac[Ibeta]
                        phi=matrix(thetac[Iphi],n,n)
                        #seta=matrix(thetac[Iseta],n,n)
                        estim=likfort2count(n,t,l,x,beta,phi,seta)
                        logLikC=estim$lik
                        hp_T=estim$hpT
                        hp_T_allM[,,i]=matrix(hp_T,n,N)
                        hp_trajM[,,i]=matrix(estim$hptraj,n,T)
                        a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
                        #a[i]=exp(logLikC-logLik[i])
                }else{
                        logLikC=-Inf
                        a[i]=0
                }       
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                #a[i]=exp(logLikC-logLik[i])*()
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }
                
                
                
                
        }
        return(list(thetamh=thetamh,logLik=logLik, h_TM=hp_T_allM,hp_trajM=hp_trajM))        
}


PMH_corr_log<-function(T,n,N,M,initialpars,Sig,x){
  thetamh=matrix(c(0),M,length(initialpars))
  logLik=rep(0,M)
  a=rep(0,M)
  hp_T=matrix(c(0),n,N)
  hp_T_allM=array(0,c(n,N,M))
  hp_trajM=array(0,c(n,T,M))
  thetamh[1,]=initialpars
  #logLik[1]=pf(N,n,thetamh[1,],x)
  
  Sig_beta=diag(Sig[1:n])
  Sig_phi1=diag(Sig[(n+1):(n+2)])
  Sig_phi2=diag(Sig[(n+3):(n+n*n)])
  Sig_seta=diag(Sig[(n+n*n+1):(2*n+n*n)])
  Sig_rho=Sig[2*n+n*n+1]
  
  Ibeta=c(1:n)
  Iphi=(n+1):(n+n*n)
  Iphi1=(n+1):(n+2)
  Iphi2=(n+3):(n+n*n)
  Iseta=(n+n*n+1):(2*n+n*n)
  Irho=2*n+n*n+1
  
  t=T
  l=N
  
  beta=initialpars[Ibeta]
  phi=matrix(initialpars[Iphi],n,n)
  seta=diag(initialpars[Iseta])
  rho=initialpars[Irho]
  seta[1,2]=rho
  seta[2,1]=rho
  estim=likfort2count(n,t,l,x,beta,phi,seta)
  logLik[1]=estim$lik
  hp_T=estim$hpT
  hp_T_allM[,,1]=matrix(hp_T,n,N)
  hp_trajM[,,1]=matrix(estim$hptraj,n,T)
  #logLik[1]=pf(N,n,initialpars,x)
  
  for(i in 2:M){
    thetac=thetamh[i-1,]
    thetac[Ibeta]=c(rmultnorm(1,c(thetamh[i-1,Ibeta]),Sig_beta))
    #################################
    #Proposal for the scale par beta
    #################################
    if(all(thetac[1:n]>0)){
      #logLikC=pf(N,n,thetac,x)
      beta=thetac[Ibeta]
      phi=matrix(thetac[Iphi],n,n)
      seta=diag(thetac[Iseta])
      seta[1,2]=thetac[Irho]
      seta[2,1]=thetac[Irho]
      estim=likfort2count(n,t,l,x,beta,phi,seta)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      hp_trajM[,,i]=matrix(estim$hptraj,n,T)
      a[i]=exp(logLikC-logLik[i-1])*(PriorValue(thetac)/PriorValue(thetamh[i-1,])) 
      #a[i]=exp(logLikC-logLik[i-1])
    }else{
      logLikC=-Inf
      a[i]=0
    }     
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i-1])
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }else{
      thetamh[i,]=thetamh[i-1,]
      logLik[i]=logLik[i-1]
    }
    #################################################
    #Proposal for the persistence par phi, first row
    #################################################
    thetac=thetamh[i,]
    thetac[Iphi1]=c(rmultnorm(1,c(thetamh[i,Iphi1]),Sig_phi1))
    beta=thetac[Ibeta]
    phi=matrix(thetac[Iphi],n,n)
    seta=diag(thetac[Iseta])
    seta[1,2]=thetac[Irho]
    seta[2,1]=thetac[Irho]
    Phiev=eigen(phi)$values
    if(all(abs(Phiev)<0.999)){
      #logLikC=pf(N,n,thetac,x)
      #a[i]=exp(logLikC-logLik[i])
      beta=thetac[Ibeta]
      phi=matrix(thetac[Iphi],n,n)
      seta=diag(thetac[Iseta])
      seta[1,2]=thetac[Irho]
      seta[2,1]=thetac[Irho]
      estim=likfort2count(n,t,l,x,beta,phi,seta)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      hp_trajM[,,i]=matrix(estim$hptraj,n,T)
      a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
    }else{
      logLikC=-Inf
      a[i]=0
    }    
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i])
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }
    
    #################################################
    #Proposal for the persistence par phi, second row
    #################################################
    thetac=thetamh[i,]
    thetac[Iphi2]=c(rmultnorm(1,c(thetamh[i,Iphi2]),Sig_phi2))
    beta=thetac[Ibeta]
    phi=matrix(thetac[Iphi],n,n)
    seta=diag(thetac[Iseta])
    seta[1,2]=thetac[Irho]
    seta[2,1]=thetac[Irho]
    Phiev=eigen(phi)$values
    if(all(abs(Phiev)<0.999)){
      #logLikC=pf(N,n,thetac,x)
      #a[i]=exp(logLikC-logLik[i])
      beta=thetac[Ibeta]
      phi=matrix(thetac[Iphi],n,n)
      seta=diag(thetac[Iseta])
      seta[1,2]=thetac[Irho]
      seta[2,1]=thetac[Irho]
      estim=likfort2count(n,t,l,x,beta,phi,seta)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      hp_trajM[,,i]=matrix(estim$hptraj,n,T)
      a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,]))
    }else{
      a[i]=0
      logLikC=-Inf
    }     
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i])
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }
    
    
    ####################################################
    #Proposal for the variance par seta, variance block
    ####################################################
    thetac=thetamh[i,]
    thetac[Iseta]=c(rmultnorm(1,c(thetamh[i,Iseta]),Sig_seta))
    beta=thetac[Ibeta]
    phi=matrix(thetac[Iphi],n,n)
    seta=diag(thetac[Iseta])
    seta[1,2]=thetac[Irho]
    seta[2,1]=thetac[Irho]
    Setaev=eigen(seta)$values
    if((seta[1,1]*seta[2,2]>seta[1,2]^2)&seta[1,1]>0.0001&seta[2,2]>0.0001){
      #logLikC=pf(N,n,thetac,x)
      estim=likfort2count(n,t,l,x,beta,phi,seta)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      hp_trajM[,,i]=matrix(estim$hptraj,n,T)
      a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
      #a[i]=exp(logLikC-logLik[i])
    }else{
      logLikC=-Inf
      a[i]=0
    }       
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i])*()
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }
    
    ################################################
    #Proposal for the correlation par in seta
    ################################################
    thetac=thetamh[i,]
    thetac[Irho]=c(rmultnorm(1,c(thetamh[i,Irho]),Sig_rho))
    seta=diag(thetac[Iseta])
    seta[1,2]=thetac[Irho]
    seta[2,1]=thetac[Irho]
    Setaev=eigen(seta)$values
    if((seta[1,1]*seta[2,2]>seta[1,2]^2)&seta[1,1]>0.0001&seta[2,2]>0.0001){
      #logLikC=pf(N,n,thetac,x)
      beta=thetac[Ibeta]
      phi=matrix(thetac[Iphi],n,n)
      #seta=matrix(thetac[Iseta],n,n)
      estim=likfort2count(n,t,l,x,beta,phi,seta)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      hp_trajM[,,i]=matrix(estim$hptraj,n,T)
      a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
      #a[i]=exp(logLikC-logLik[i])
    }else{
      logLikC=-Inf
      a[i]=0
    }       
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i])*()
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }
    
    
    
    
  }
  return(list(thetamh=thetamh,logLik=logLik, h_TM=hp_T_allM,hp_trajM=hp_trajM))        
}

# Partilce Metropolis-Hastings for parameters sampling with correlation,
# implemented completely in R --> ok for testing, but slow for actually being used
PMH_corr_inR<-function(T,n,N,M,initialpars,Sig,x){
  thetamh=matrix(c(0),M,length(initialpars))
  logLik=rep(0,M)
  a=rep(0,M)
  hp_T=matrix(c(0),n,N)
  hp_T_allM=array(0,c(n,N,M))
  hp_trajM=array(0,c(n,T,M))
  thetamh[1,]=initialpars
  #logLik[1]=pf(N,n,thetamh[1,],x)
  
  Sig_beta=diag(Sig[1:n])
  Sig_phi1=diag(Sig[(n+1):(n+2)])
  Sig_phi2=diag(Sig[(n+3):(n+n*n)])
  Sig_seta=diag(Sig[(n+n*n+1):(2*n+n*n)])
  Sig_rho=Sig[2*n+n*n+1]
  
  Ibeta=c(1:n)
  Iphi=(n+1):(n+n*n)
  Iphi1=(n+1):(n+2)
  Iphi2=(n+3):(n+n*n)
  Iseta=(n+n*n+1):(2*n+n*n)
  Irho=2*n+n*n+1
  
  t=T
  l=N
  
  beta=initialpars[Ibeta]
  phi=matrix(initialpars[Iphi],n,n)
  seta=diag(initialpars[Iseta])
  rho=initialpars[Irho]
  seta[1,2]=rho
  seta[2,1]=rho
  estim=pf(N,n,initialpars,x)
  logLik[1]=estim$lik
  hp_T=estim$hpT
  hp_T_allM[,,1]=matrix(hp_T,n,N)
  
  
  for(i in 2:M){
    thetac=thetamh[i-1,]
    thetac[Ibeta]=c(rmultnorm(1,c(thetamh[i-1,Ibeta]),Sig_beta))
    #################################
    #Proposal for the scale par beta
    #################################
    if(all(thetac[1:n]>0)){
      #logLikC=pf(N,n,thetac,x)
      beta=thetac[Ibeta]
      phi=matrix(thetac[Iphi],n,n)
      seta=diag(thetac[Iseta])
      seta[1,2]=thetac[Irho]
      seta[2,1]=thetac[Irho]
      #estim=likfort2count(n,t,l,x,beta,phi,seta)
      estim=pf(N,n,thetac,x)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      a[i]=exp(logLikC-logLik[i-1])*(PriorValue(thetac)/PriorValue(thetamh[i-1,])) 
      #a[i]=exp(logLikC-logLik[i-1])
    }else{
      logLikC=-Inf
      a[i]=0
    }     
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i-1])
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }else{
      thetamh[i,]=thetamh[i-1,]
      logLik[i]=logLik[i-1]
    }
    #################################################
    #Proposal for the persistence par phi, first row
    #################################################
    thetac=thetamh[i,]
    thetac[Iphi1]=c(rmultnorm(1,c(thetamh[i,Iphi1]),Sig_phi1))
    beta=thetac[Ibeta]
    phi=matrix(thetac[Iphi],n,n)
    seta=diag(thetac[Iseta])
    seta[1,2]=thetac[Irho]
    seta[2,1]=thetac[Irho]
    Phiev=eigen(phi)$values
    if(all(abs(Phiev)<0.999)){
      #logLikC=pf(N,n,thetac,x)
      #a[i]=exp(logLikC-logLik[i])
      beta=thetac[Ibeta]
      phi=matrix(thetac[Iphi],n,n)
      seta=diag(thetac[Iseta])
      seta[1,2]=thetac[Irho]
      seta[2,1]=thetac[Irho]
      #estim=likfort2count(n,t,l,x,beta,phi,seta)
      estim=pf(N,n,thetac,x)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
    }else{
      logLikC=-Inf
      a[i]=0
    }    
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i])
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }
    
    #################################################
    #Proposal for the persistence par phi, second row
    #################################################
    thetac=thetamh[i,]
    thetac[Iphi2]=c(rmultnorm(1,c(thetamh[i,Iphi2]),Sig_phi2))
    beta=thetac[Ibeta]
    phi=matrix(thetac[Iphi],n,n)
    seta=diag(thetac[Iseta])
    seta[1,2]=thetac[Irho]
    seta[2,1]=thetac[Irho]
    Phiev=eigen(phi)$values
    if(all(abs(Phiev)<0.999)){
      #logLikC=pf(N,n,thetac,x)
      #a[i]=exp(logLikC-logLik[i])
      beta=thetac[Ibeta]
      phi=matrix(thetac[Iphi],n,n)
      seta=diag(thetac[Iseta])
      seta[1,2]=thetac[Irho]
      seta[2,1]=thetac[Irho]
      #estim=likfort2count(n,t,l,x,beta,phi,seta)
      estim=pf(N,n,thetac,x)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,]))
    }else{
      a[i]=0
      logLikC=-Inf
    }     
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i])
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }
    
    
    ####################################################
    #Proposal for the variance par seta, variance block
    ####################################################
    thetac=thetamh[i,]
    thetac[Iseta]=c(rmultnorm(1,c(thetamh[i,Iseta]),Sig_seta))
    beta=thetac[Ibeta]
    phi=matrix(thetac[Iphi],n,n)
    seta=diag(thetac[Iseta])
    seta[1,2]=thetac[Irho]
    seta[2,1]=thetac[Irho]
    Setaev=eigen(seta)$values
    if((seta[1,1]*seta[2,2]>seta[1,2]^2)&seta[1,1]>0.0001&seta[2,2]>0.0001){
      #logLikC=pf(N,n,thetac,x)
      #estim=likfort2count(n,t,l,x,beta,phi,seta)
      estim=pf(N,n,thetac,x)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
      #a[i]=exp(logLikC-logLik[i])
    }else{
      logLikC=-Inf
      a[i]=0
    }       
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i])*()
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }
    
    ################################################
    #Proposal for the correlation par in seta
    ################################################
    thetac=thetamh[i,]
    thetac[Irho]=c(rmultnorm(1,c(thetamh[i,Irho]),Sig_rho))
    seta=diag(thetac[Iseta])
    seta[1,2]=thetac[Irho]
    seta[2,1]=thetac[Irho]
    Setaev=eigen(seta)$values
    if((seta[1,1]*seta[2,2]>seta[1,2]^2)&seta[1,1]>0.0001&seta[2,2]>0.0001){
      #logLikC=pf(N,n,thetac,x)
      beta=thetac[Ibeta]
      phi=matrix(thetac[Iphi],n,n)
      #seta=matrix(thetac[Iseta],n,n)
      #estim=likfort2count(n,t,l,x,beta,phi,seta)
      estim=pf(N,n,thetac,x)
      logLikC=estim$lik
      hp_T=estim$hpT
      hp_T_allM[,,i]=matrix(hp_T,n,N)
      a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
      #a[i]=exp(logLikC-logLik[i])
    }else{
      logLikC=-Inf
      a[i]=0
    }       
    #a[i]=min(1,exp(logLikC-logLik[i-1]))   
    #a[i]=exp(logLikC-logLik[i])*()
    u=runif(1,0,1)
    if(u<a[i]){
      thetamh[i,]=thetac
      logLik[i]=logLikC
    }
    
    
    
    
  }
  return(list(thetamh=thetamh,logLik=logLik,h_TM=hp_T_allM,hp_trajM=hp_trajM))        
}


simmulticount<-function(T,beta,phi,seta,Seed) {
        if (!missing(Seed))set.seed(Seed)
        BurnIn=500
        n=length(beta)
        #dim(beta)<-c(n,1)
        etatrue=rmultnorm((T+BurnIn),rep(0,n),seta)
        h=matrix(c(0),n,(T+BurnIn))
        lambda=matrix(c(0),n,(T+BurnIn))
        x=matrix(c(0),n,(T+BurnIn))
        h[,1]=beta+etatrue[,1]
        lambda[,1]=exp(h[,1])
        for(i in 1:n){
                x[i,1]=rpois(1,lambda[i,1])
        }
        
        for (t in 2:(T+BurnIn)) { 
                h[,t]=phi%*%h[,t-1]+etatrue[,t]
                lambda[,t]=exp(h[,t])
                for(i in 1:n){
                        x[i,t]=rpois(1,beta*lambda[i,t])
                }
        }
        return(list(x=x[,BurnIn:(T+BurnIn)],lambda=lambda[,(BurnIn:(T+BurnIn))],h=h[,(BurnIn:(T+BurnIn))]))
}


PriorValue=function(theta){
        diag_seq=seq(from=1,to=n*n,by=n+1)
        nondiag_seq=seq(from=1,to=n*n,by=1)
        nondiag_seq=nondiag_seq[-diag_seq]
        Iphi=(n+1):(n+n*n)
        prior=prod(dnorm(theta[1:n],c(0,0),c(25,25)))*prod(dnorm(theta[Iphi[nondiag_seq]],c(0,0),c(0.5,0.5)))*
                prod(dbeta((1+theta[Iphi[diag_seq]])/2, 20, 1.5))*     
                prod(digamma(theta[(n+n*n+1):(2*n+n*n)],SigmaEtaAlpha,SigmaEtaBeta))*
                dunif(theta[(2*n+n*n+1)],min=-1,max=1)
        return(prior)
}

digamma<-function(x,a,b) {
        return(b^a/gamma(a)*x^(-a-1)*exp(-b/x))
}

PMH_nocorr<-function(T,n,N,M,initialpars,Sig,x){
        thetamh=matrix(c(0),M,length(initialpars))
        logLik=rep(0,M)
        a=rep(0,M)
        particles=array((0),c(n*T,N,M))
        partconstr=array(0, c(n,N,T))
        thetamh[1,]=initialpars
        #logLik[1]=pf(N,n,thetamh[1,],x)
        
        Sig_beta=diag(Sig[1:n])
        Sig_phi=diag(Sig[(n+1):(n+n*n)])
        Sig_seta=diag(Sig[(n+n*n+1):(2*n+n*n)])
        
        Ibeta=c(1:n)
        Iphi=(n+1):(n+n*n)
        Iseta=(n+n*n+1):(2*n+n*n)
        
        t=T
        l=N
        
        beta=initialpars[Ibeta]
        phi=matrix(initialpars[Iphi],n,n)
        seta=diag(initialpars[Iseta])
        logLik[1]=likfort2count(n,t,l,x,beta,phi,seta)$lik
        particles
        #logLik[1]=pf(N,n,initialpars,x)
        
        for(i in 2:M){
                thetac=thetamh[i-1,]
                thetac[Ibeta]=c(rmultnorm(1,c(thetamh[i-1,Ibeta]),Sig_beta))
                #################################
                #Proposal for the scale par beta
                #################################
                if(all(thetac[1:n]>0)){
                        logLikC=pf(N,n,thetac,x)
                        #beta=thetac[Ibeta]
                        #phi=matrix(thetac[Iphi],n,n)
                        #seta=matrix(thetac[Iseta],n,n)
                        #logLikC=likfort2count(n,t,l,x,beta,phi,seta)
                        a[i]=exp(logLikC-logLik[i-1])*(PriorValue(thetac)/PriorValue(thetamh[i-1,])) 
                }else{
                        logLikC=-Inf
                        a[i]=0
                }       
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }else{
                        thetamh[i,]=thetamh[i-1,]
                        logLik[i]=logLik[i-1]
                }
                #################################################
                #Proposal for the persistence par phi, first row
                #################################################
                thetac=thetamh[i,]
                thetac[Iphi]=c(rmultnorm(1,c(thetamh[i,Iphi]),Sig_phi))
                Phimh=matrix(thetac[Iphi],n,n)
                Phiev=eigen(Phimh)$values
                if(all(abs(Phiev)<0.999)){
                        logLikC=pf(N,n,thetac,x)
                        #beta=thetac[Ibeta]
                        #phi=matrix(thetac[Iphi],n,n)
                        #seta=matrix(thetac[Iseta],n,n)
                        #logLikC=likfort2count(n,t,l,x,beta,phi,seta)
                        a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
                }else{
                        logLikC=-Inf
                        a[i]=0
                }     
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }
                #####################################
                #Proposal for the variance par seta
                #####################################
                thetac=thetamh[i,]
                thetac[Iseta]=c(rmultnorm(1,c(thetamh[i,Iseta]),Sig_seta))
                if(all(thetac[Iseta]>0.00001)){
                        logLikC=pf(N,n,thetac,x)
                        #beta=thetac[Ibeta]
                        #phi=matrix(thetac[Iphi],n,n)
                        #seta=matrix(thetac[Iseta],n,n)
                        #logLikC=likfort2count(n,t,l,x,beta,phi,seta)
                        a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/PriorValue(thetamh[i,])) 
                }else{
                        logLikC=-Inf
                        a[i]=0
                }     
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }
        }
        return(list(thetamh=thetamh,logLik=logLik))        
}

PMH_IW<-function(T,n,N,M,tau,initialpars,Sig,x){
        thetamh=matrix(c(0),M,length(initialpars))
        logLik=rep(0,M)
        a=rep(0,M)
        
        thetamh[1,]=initialpars
        #logLik[1]=pf(N,n,thetamh[1,],x)
        
        Sig_beta=diag(Sig[1:n])
        Sig_phi1=diag(Sig[(n+1):(n+2)])
        Sig_phi2=diag(Sig[(n+3):(n+n*n)])
        Sig_seta=diag(Sig[(n+n*n+1):(2*n+n*n)])
        Sig_rho=Sig[2*n+n*n+1]
        
        Ibeta=c(1:n)
        Iphi=(n+1):(n+n*n)
        Iphi1=(n+1):(n+2)
        Iphi2=(n+3):(n+n*n)
        Iseta=(n+n*n+1):(2*n+n*n)
        Irho=2*n+n*n+1
        
        t=T
        l=N
        
        beta=initialpars[Ibeta]
        phi=matrix(initialpars[Iphi],n,n)
        seta=diag(initialpars[Iseta])
        rho=initialpars[Irho]
        seta[1,2]=rho
        seta[2,1]=rho
        logLik[1]=likfort2count(n,t,l,x,beta,phi,seta)
        #logLik[1]=pf(N,n,initialpars,x)
        for(i in 2:M){
                thetac=thetamh[i-1,]
                thetac[Ibeta]=c(rmultnorm(1,c(thetamh[i-1,Ibeta]),Sig_beta))
                #################################
                #Proposal for the scale par beta
                #################################
                if(all(thetac[1:n]>0)){
                        #logLikC=pf(N,n,thetac,x)
                        beta=thetac[Ibeta]
                        phi=matrix(thetac[Iphi],n,n)
                        seta=diag(thetac[Iseta])
                        seta[1,2]=thetac[Irho]
                        seta[2,1]=thetac[Irho]
                        logLikC=likfort2count(n,t,l,x,beta,phi,seta)
                        a[i]=exp(logLikC-logLik[i-1])*(PriorValue(thetac)/
                                                               PriorValue(thetamh[i-1,])) 
                        #a[i]=exp(logLikC-logLik[i-1])
                        if(is.nan(a[i])){
                                a[i]=0
                        }
                }else{
                        logLikC=-Inf
                        a[i]=0
                }     
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                #a[i]=exp(logLikC-logLik[i-1])
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }else{
                        thetamh[i,]=thetamh[i-1,]
                        logLik[i]=logLik[i-1]
                }
                #################################################
                #Proposal for the persistence par phi, first row
                #################################################
                thetac=thetamh[i,]
                thetac[Iphi1]=c(rmultnorm(1,c(thetamh[i,Iphi1]),Sig_phi1))
                Phimh=matrix(thetac[Iphi],n,n)
                Phiev=eigen(Phimh)$values
                if(all(abs(Phiev)<0.999)){
                        #logLikC=pf(N,n,thetac,x)
                        #a[i]=exp(logLikC-logLik[i])
                        beta=thetac[Ibeta]
                        phi=matrix(thetac[Iphi],n,n)
                        seta=diag(thetac[Iseta])
                        seta[1,2]=thetac[Irho]
                        seta[2,1]=thetac[Irho]
                        logLikC=likfort2count(n,t,l,x,beta,phi,seta)
                        a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/
                                                             PriorValue(thetamh[i,])) 
                        if(is.nan(a[i])){
                                a[i]=0
                        }
                }else{
                        logLikC=-Inf
                        a[i]=0
                }    
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                #a[i]=exp(logLikC-logLik[i])
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }
                
                #################################################
                #Proposal for the persistence par phi, second row
                #################################################
                thetac=thetamh[i,]
                thetac[Iphi2]=c(rmultnorm(1,c(thetamh[i,Iphi2]),Sig_phi2))
                Phimh=matrix(thetac[Iphi],n,n)
                Phiev=eigen(Phimh)$values
                if(all(abs(Phiev)<0.999)){
                        #logLikC=pf(N,n,thetac,x)
                        #a[i]=exp(logLikC-logLik[i])
                        beta=thetac[Ibeta]
                        phi=matrix(thetac[Iphi],n,n)
                        seta=diag(thetac[Iseta])
                        seta[1,2]=thetac[Irho]
                        seta[2,1]=thetac[Irho]
                        logLikC=likfort2count(n,t,l,x,beta,phi,seta)
                        a[i]=exp(logLikC-logLik[i])*(PriorValue(thetac)/
                                                             PriorValue(thetamh[i,]))
                        if(is.nan(a[i])){
                                a[i]=0
                        }
                }else{
                        a[i]=0
                        logLikC=-Inf
                }     
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                #a[i]=exp(logLikC-logLik[i])
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }
                
                ####################################################
                #Proposal for the variance-covariance matrix
                ####################################################
                thetac=thetamh[i,]
                #thetac[Iseta]=c(rmultnorm(1,c(thetamh[i,Iseta]),Sig_seta))
                beta=thetac[Ibeta]
                phi=matrix(thetac[Iphi],n,n)
                seta=diag(thetac[Iseta])
                seta[1,2]=thetac[Irho]
                seta[2,1]=thetac[Irho]
                seta_prev=seta
                #seta=riwish(tau+3,tau*seta)
                seta=rinvwishart(tau+3,tau*seta_prev)
                thetac[Iseta]=diag(seta)
                thetac[Irho]=seta[1,2]
                Setaev=eigen(seta)$values
                if(all(abs(Setaev)<0.999)){
                        #logLikC=pf(N,n,thetac,x)
                        logLikC=likfort2count(n,t,l,x,beta,phi,seta)
                        a[i]=exp(logLikC-logLik[i])*((PriorValue(thetac)/
                                                             PriorValue(thetamh[i,])))*
                                (det(seta_prev)^((2*tau+9)/2)/(det(seta))^((2*tau+9)/2))*
                                exp((tau/2)*(tr(seta%*%t(seta)))-(tau/2)*tr(seta%*%t(seta)))
                        #a[i]=exp(logLikC-logLik[i])
                        if(is.nan(a[i])){
                                a[i]=0
                        }
                }else{
                        logLikC=-Inf
                        a[i]=0
                }       
                #a[i]=min(1,exp(logLikC-logLik[i-1]))   
                #a[i]=exp(logLikC-logLik[i])*()
                u=runif(1,0,1)
                if(u<a[i]){
                        thetamh[i,]=thetac
                        logLik[i]=logLikC
                }
        }
        return(list(thetamh=thetamh,logLik=logLik))        
}

Mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

