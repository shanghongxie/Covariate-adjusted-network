


#############################
#####  Data Generation  #####
#############################

library(MASS)

##############################################
###  Gen 2:Different varying coefficients  ###
##############################################
## delta: variance of covariate xi
GenDataD <- function (n, p, q, kbeta, obeta, obeta0, delta,beta0,xbeta,mbeta,oeta) {
  
  ## number of undirected edges
  po=(p-1)*p/2
  
  Omega=matrix(obeta,nrow=q)
  Kappa=matrix(kbeta,nrow=q)
  
  ###  Generate X, M, MI  ###
  i=1; X=NULL; M=NULL; MI=NULL;idv=0;dv=0;imax=0;iboth=0

  repeat {
    dv=dv+1
   
    xi=mvrnorm(1,mu=rep(0,q),diag(delta,q))
    
    ox=c(xi%*%Omega); kx=c(xi%*%Kappa)
    
    Omegai=diag(obeta0); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        Omegai[s,r]=-ox[sr]; Omegai[r,s]=-ox[sr]
        sr=sr+1
      }
    }
    
    
    temmi=min(eigen(Omegai)$values)
    if (max(Omegai)>obeta0[1]) {imax=imax+1}
    if (max(Omegai)>obeta0[1] & (temmi<0.05)) {iboth=iboth+1}
    if (temmi<0.05) {idv=idv+1; next}
    
    Sigmai=solve(Omegai)
    
    ### mean matrix
    mui=Sigmai%*%kx
  
    
    ### generate biomarker nodes M
    mi=mvrnorm(1, mui, Sigmai)
    
    ### generate 2nd stage data ###
    ## Mutual Information
    count=0
    mii=numeric((p-1)*p/2); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        pho=ox[sr]^2/obeta0[s]/obeta0[r]
        if (abs(pho)>0.9999){
          count=count+1
        }
        
        pho=ifelse(abs(pho)>0.9999, sign(pho)*0.9999, pho)
        mii[sr]=-log(1-pho)/2
        sr=sr+1
      }
    }
    
    X=rbind(X,xi)
    M=rbind(M,mi)
    MI=rbind(MI,mii)
    
    if (i==n) break
    i=i+1
  }
  
  
  # ###  Scaled for 2nd stage ###
  index=which(apply(MI,2,sum)!=0)
  inM=MI[,index]
  
  sdMI=matrix(0,n,po)
  temM=t(inM)-apply(inM,2,mean)
  temM=temM/sqrt(apply(temM^2,1,mean))
  temM=t(temM)
  sdMI[,index]=temM

  tem=t(X)-apply(X,2,mean)
  tem=tem/sqrt(apply(tem^2,1,mean))
  tem=t(tem)
  sdX=tem
     
  tem=t(M)-apply(M,2,mean)
  tem=tem/sqrt(apply(tem^2,1,mean))
  tem=t(tem)
  sdM=tem
   
  ####  Generate Y  ###
  Y=NULL
  for (i in 1:n) {
    xb=beta0+sum(sdX[i,]*xbeta)+sum(sdM[i,]*mbeta)+sum(sdMI[i,]*oeta)
    Y[i]=rnorm(1,mean=xb,sd=1)
  }
  
  
  return(list(X=X, M=M, Y=Y, MI=MI, sdX=sdX, sdM=sdM, sdMI=sdMI, idv=idv,dv=dv,count=count,imax=imax,iboth=iboth))   
}

