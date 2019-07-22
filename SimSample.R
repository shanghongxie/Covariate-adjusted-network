rm(list=ls(all=TRUE))

library(Rcpp)
library(RcppEigen)


##############################
###  Algorithm parameters  ###
##############################

n=500; 

p=5; q=5; po=(p-1)*p/2;delta=1


### precision matrix
Omega=matrix(0,q,po)
Omega[1:3,1]=-c(0.5,1,1.5)
Omega[1:3,5]=-c(1,0.5,1.5)
Omega[1:3,8]=-c(-1.5,0.5,1)
Omega[1:3,10]=-c(0.5,1.5,-1)

obeta=as.vector(Omega)

## kappa: mean 

Kappa=matrix(0,q,p)
Kappa[1:3,1]=c(0.5,1,1.5)
Kappa[1:3,2]=c(1,0.5,1.5)
Kappa[1:3,3]=c(-1.5,1,0.5)

kbeta=as.vector(Kappa)

### 1/sigma^2
obeta0=rep(5,p)   ### sigma0^2=0.1

### 2nd stage parameter
beta0=0   ### intercept
xbeta=c(1,2,rep(0,q-2))  ## coefficients for covariate X
mbeta=c(1,3,rep(0,p-2))  ## coefficients for biomaker nodes M
oeta=c(1,0,0,0,2,0,0,0,0,0,rep(0,po-10))  ## coefficients for connections



#########################
###  Simulation - CV  ###
#########################

source("SimGenerate.R")

sourceCpp("cNetworkC.cpp")

source("cNetworkR.R") 

tau=1; lambda=NULL; rlambda=NULL; nlambda=10; nfolds=10
t=10; maxit=1e+5;
maxitOut=1e+5; maxitIn=1e+5
thresh=1e-6; threshIn=thresh/10
tau=1


library(MASS)

set.seed(123)

nSim=100

############################
#### Generate data  ########
############################

## m: biomarkers (nodes), x: covariates; y: outcome; mi:mutual information
m=list(); x=list(); y=list(); mi=list()
sdm=list(); sdx=list(); sdmi=list()

rejN=NULL

for (iSim in 1:nSim) {
  cat(iSim)
  tmp=GenDataD(n,p,q,kbeta,obeta, obeta0, delta,beta0,xbeta,mbeta,oeta)
  M=tmp$M; X=tmp$X; Y=tmp$Y; MI=tmp$MI
  sdM=tmp$sdM; sdX=tmp$sdX; sdMI=tmp$sdMI
  rejN=c(rejN, tmp$idv)

  m[[iSim]]=M; x[[iSim]]=X; y[[iSim]]=Y; mi[[iSim]]=MI; sdm[[iSim]]=sdM
  sdx[[iSim]]=sdX; sdmi[[iSim]]=sdMI
  
}
mean(rejN)
rej=mean(rejN) 


################################
### 1st stage estimation  ######
################################

resultO=NULL; resultK=NULL; resultO0=NULL
for (iSim in 1:nSim) {
  cat(iSim)
  
  X=x[[iSim]]; M=m[[iSim]]
  out=LmNetwork_1st(X, M, tau, lambda=NULL, nlambda=10, rlambda=NULL, nfolds=10, foldid=NULL, t=10, maxitOut=maxitOut, maxitIn=maxitIn, thresh=thresh, threshIn=threshIn, itrunc=TRUE, standardize=FALSE, keep.coef=TRUE)
  
  ### precision matrix (omega)
  resultO=cbind(resultO, out$coef.thrI$omega)
  ### mean (kappa)
  resultK=cbind(resultK, out$coef.thrI$kappa)
  resultO0=cbind(resultO0, out$coef.thrI$omega0)
  
}



### resultO: precision matrix
### MSE
mseO=mean(apply((resultO-obeta)^2,2,sum))  

tpO=mean(apply(resultO!=0 & obeta!=0,2,sum))
fpO=mean(apply(resultO!=0 & obeta==0,2,sum)) 
tnO=mean(apply(resultO==0 & obeta==0,2,sum));
fnO=mean(apply(resultO==0 & obeta!=0,2,sum)) 

tpOM=apply(resultO!=0 & obeta!=0,2,sum)
tnOM=apply(resultO==0 & obeta==0,2,sum)
fpOM=apply(resultO!=0 & obeta==0,2,sum)
fnOM=apply(resultO==0 & obeta!=0,2,sum)

mccO=mean((tpOM*tnOM-fpOM*fnOM)/sqrt(tpOM+fpOM)/sqrt(tpOM+fnOM)/sqrt(tnOM+fpOM)/sqrt(tnOM+fnOM))

######### resultK, mean

### MSE
mseK=sum(apply((resultK-kbeta)^2,1,mean))  

tpK=mean(apply(resultK!=0 & kbeta!=0,2,sum))
fpK=mean(apply(resultK!=0 & kbeta==0,2,sum)) 

tnK=mean(apply(resultK==0 & kbeta==0,2,sum));
fnK=mean(apply(resultK==0 & kbeta!=0,2,sum)) 

tpKM=apply(resultK!=0 & kbeta!=0,2,sum)
tnKM=apply(resultK==0 & kbeta==0,2,sum)
fpKM=apply(resultK!=0 & kbeta==0,2,sum)
fnKM=apply(resultK==0 & kbeta!=0,2,sum)

mccK=mean((tpKM*tnKM-fpKM*fnKM)/sqrt(tpKM+fpKM)/sqrt(tpKM+fnKM)/sqrt(tnKM+fpKM)/sqrt(tnKM+fnKM))



### MSE
mse0=sum(apply((resultO0-obeta0)^2,1,mean))  


result1st=rbind(c(mseK,tpK,fpK,tnK,fnK,mccK),c(mseO,tpO,fpO,tnO,fnO,mccO))
rownames(result1st)=c("kappa","omega")
colnames(result1st)=c("SSE","TP","FP","TN","FN","MCC")


##########################################
##### Calculate Mutual Information  ######
##########################################
iSim=1;nSim=100

eta0s=NULL; xbetas=NULL; mbetas=NULL; beta0s=NULL

mists=list(); 
misde=list()

for (iSim in 1:nSim){

  cat(iSim)
 
  ######################################
  ################ use log #############
  Omega=matrix(resultO[,iSim],nrow=q)  
  Omega0=resultO0[,iSim]
   
  mist=matrix(rep(0,n*po),n,po)
  # 
  for (i in 1:n){
    ox=c(x[[iSim]][i,]%*%Omega)
    sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        pho=ox[sr]^2/Omega0[s]/Omega0[r]
        ## if the absolute value of partial correlation is too close to 1, let abs(pho)=0.999
        pho=ifelse(abs(pho)>0.999, sign(pho)*0.999, pho)
        mist[i,sr]=-log(1-pho)/2
        sr=sr+1
      }
    }
  }
  
  mists[[iSim]]=mist
  
  
  ### standardize mutual information
  MIm=apply(mist,2,mean)   
  temM=t(mist)-MIm
  MIsd=sqrt(apply(temM^2,1,mean))
  
  
  for (j in 1:length(MIsd)){
    if (MIsd[j]>0.005){  
      temM[j,]=temM[j,]/MIsd[j]  ### if MIsd too small, then exclude
    }
  }
  temM=t(temM);
  sdmist=temM

  
  index=which(MIsd>0.005)
  
  sdmist=sdmist[,index]
  sdmist=as.matrix(sdmist)
  
  misde[[iSim]]=sdmist
  
}


######################################
###########    2nd stage    ##########
######################################

iSim=1;nSim=100

eta0s=NULL; xbetas=NULL; mbetas=NULL;

for (iSim in 1:nSim){
 
  cat(iSim)
  
  sdmist=misde[[iSim]]
  
  mist=mists[[iSim]]
  
  MIm=apply(mist,2,mean)   ### standardized
  temM=t(mist)-MIm
  MIsd=sqrt(apply(temM^2,1,mean))
  index=which(MIsd>0.005)
  

  ALL=cbind(sdx[[iSim]],sdm[[iSim]],sdmist)
  
  ### use adaptive lasso
  
  out=EnetLm(ALL, y[[iSim]], alpha=1.0, lambda=NULL, nlambda=100, rlambda=NULL, nfolds=10, foldid=NULL, inzero=TRUE, 
             adaptive=TRUE, aini=NULL, isd=TRUE, keep.beta=FALSE, thresh=1e-7, maxit=1e+5)
  
  
  beta=out$Beta0
  
  ## xbeta covariates
  xbetas=cbind(xbetas,beta[1:(ncol(sdx[[iSim]]))])
  
  ## mbeta mediators
  mbetas=cbind(mbetas,beta[(ncol(sdx[[iSim]])+1):(ncol(sdx[[iSim]])+ncol(sdm[[iSim]]))])
  
  ## gamma
  eta0=rep(0,po)
  eta0[index]=beta[(ncol(sdx[[iSim]])+ncol(sdx[[iSim]])+1):(ncol(sdx[[iSim]])+ncol(sdm[[iSim]])+ncol(sdmist))]
  
  eta0s=cbind(eta0s,eta0)
 
}



##############################
###     MI (connections)   ###
##############################

### MSE
etamse=sum(apply((eta0s-oeta)^2,1,mean))  


etatp=mean(apply(eta0s!=0 & oeta!=0,2,sum)); 
etafp=mean(apply(eta0s!=0 & oeta==0,2,sum)) 


etatn=mean(apply(eta0s==0 & oeta==0,2,sum));
etafn=mean(apply(eta0s==0 & oeta!=0,2,sum)) 

etatpM=apply(eta0s!=0 & oeta!=0,2,sum)
etatnM=apply(eta0s==0 & oeta==0,2,sum)
etafpM=apply(eta0s!=0 & oeta==0,2,sum)
etafnM=apply(eta0s==0 & oeta!=0,2,sum)

etamcc=mean((etatpM*etatnM-etafpM*etafnM)/sqrt(etatpM+etafpM)/sqrt(etatpM+etafnM)/sqrt(etatnM+etafpM)/sqrt(etatnM+etafnM))


##############################
###    xbeta (covariates)  ###
##############################

## MSE
xbetamse=sum(apply((xbetas-xbeta)^2,1,mean))  

### true positive xbeta
xbetatp=mean(apply(xbetas!=0 & xbeta!=0,2,sum)); 
### false positive xbeta
xbetafp=mean(apply(xbetas!=0 & xbeta==0,2,sum))

xbetatn=mean(apply(xbetas==0 & xbeta==0,2,sum));
xbetafn=mean(apply(xbetas==0 & xbeta!=0,2,sum)) 

xbetatpM=apply(xbetas!=0 & xbeta!=0,2,sum)
xbetatnM=apply(xbetas==0 & xbeta==0,2,sum)
xbetafpM=apply(xbetas!=0 & xbeta==0,2,sum)
xbetafnM=apply(xbetas==0 & xbeta!=0,2,sum)



xbetamccV=(xbetatpM*xbetatnM-xbetafpM*xbetafnM)/sqrt(xbetatpM+xbetafpM)/sqrt(xbetatpM+xbetafnM)/sqrt(xbetatnM+xbetafpM)/sqrt(xbetatnM+xbetafnM)
xbetamcc=mean(xbetamccV)


### xbeta


## MSE
mbetamse=sum(apply((mbetas-mbeta)^2,1,mean))  #


mbetatp=mean(apply(mbetas!=0 & mbeta!=0,2,sum));

mbetafp=mean(apply(mbetas!=0 & mbeta==0,2,sum)) 

mbetatn=mean(apply(mbetas==0 & mbeta==0,2,sum));
mbetafn=mean(apply(mbetas==0 & mbeta!=0,2,sum)) 


mbetatpM=apply(mbetas!=0 & mbeta!=0,2,sum)
mbetatnM=apply(mbetas==0 & mbeta==0,2,sum)
mbetafpM=apply(mbetas!=0 & mbeta==0,2,sum)
mbetafnM=apply(mbetas==0 & mbeta!=0,2,sum)

mbetamccV=(mbetatpM*mbetatnM-mbetafpM*mbetafnM)/sqrt(mbetatpM+mbetafpM)/sqrt(mbetatpM+mbetafnM)/sqrt(mbetatnM+mbetafpM)/sqrt(mbetatnM+mbetafnM)

mbetamcc=mean(mbetamccV[!is.na(mbetamccV)])

result2nd=rbind(c(etamse,etatp,etafp,etatn,etafn,etamcc),c(xbetamse,xbetatp,xbetafp,xbetatn,xbetafn,xbetamcc),
                c(mbetamse,mbetatp,mbetafp,mbetatn,mbetafn,mbetamcc))

rownames(result2nd)=c("eta","xbeta","mbeta")
colnames(result2nd)=c("SSE","TP","FP","TN","FN","MCC")


