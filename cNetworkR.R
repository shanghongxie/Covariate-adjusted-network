

#######################################
############# 1st stage  ###############
#######################################

LmNetwork_1st=function(X, M, tau=1, lambda=NULL, nlambda=10, rlambda=NULL, nfolds=1, foldid=NULL, t=10,  maxitOut=1e+5, maxitIn=1e+5, thresh=1e-6, threshIn=1e-7,itrunc=TRUE, standardize=FALSE, keep.coef=FALSE){

  ### X: covariates;   M: biomarker nodes;   tau=1: lasso, tau=a where a<1 then sparse group lasso; t: stepsize; itrunc: hard-thresholding or not; standardize: standardize covariate or not;
  ### maxitOut: maximum iterations among columns; maxitIn: maximum iterations inside a column; thresh: convergence threshold among columns: threshIn: convergence threshold inside a column;
  
  n=nrow(X); p=ncol(M); q=ncol(X)
  po=(p-1)*p/2; 
  
  ### scaleC and standardized covariates
  if (!standardize){
    xscale=rep(1, q)
  } else{
    tem=scaleC(X)
    xscale=tem$sd;X=tem$x
    rm(tem)
  }
  
  ### Lambda path 1st stage
  if (is.null(lambda)) {
    maxLambdas=maxLambdaC1(tau,M,X)
    
    lambda.max=max(unlist(maxLambdas))
    lambda.min=ifelse(is.null(rlambda), lambda.max*0.0001, lambda.max*rlambda)  
    lambda=lambda.max*(lambda.min/lambda.max)^(c(0:(nlambda-1))/(nlambda-1))
    
  } else {
    nlambda=length(lambda)
  }
  
  
  ######## Run #######
  
  ####### Full #######
  out=NetLmC(lambda, tau, M, X, t, maxitOut, maxitIn, thresh, threshIn)
  
  ###  Record  ###
 
  out$Kappa=out$Kappa/xscale; out$Omega=out$Omega/xscale
  
  nzeroGA=numeric(nlambda); nzeroGO=numeric(nlambda)
  for (il in 1:nlambda) {
    nzeroGA[il]=sum(apply(matrix(out$Kappa[,il], nrow=q)!=0,2,any))
    nzeroGO[il]=sum(apply(matrix(out$Omega[,il], nrow=q)!=0,2,any))
  }
  fit=data.frame(lambda=lambda, ll=out$ll, obj=out$obj, nzeroGA=nzeroGA, nzeroGO=nzeroGO)
  
  
  if (nfolds==1 & is.null(foldid)) {
    output=list(Kappa=out$Kappa, Omega=-out$Omega, fit=fit, flag=out$flag)  ## , t=out$t
    #class(output)="cNetwork"
    return(output)
    
  } else {
    #####  CV cross-validation  #####
    
    ###  Split data  ###
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=nrow(M)))
    } else {
      nfolds=max(foldid)
    }
    weighti=c(table(foldid))
    
    
    ###  CV estimates  ###
    
    outi=list(); cvRSS=matrix(NA, nrow=nfolds, ncol=nlambda)
    for (i in 1:nfolds) {
      temid=(foldid==i)

      outi[[i]]=NetLmCvC(lambda, tau, M[!temid,], X[!temid,], M[temid,], X[temid,], t, maxitOut, maxitIn, thresh, threshIn) ## revise 11/23/16 replace maxitIn by maxitOut
      cvRSS[i, ]=outi[[i]]$llk
    }
    
    
    outtest=NetLmCvC(lambda, tau, M[!temid,], X[!temid,], M[temid,], X[temid,], t, maxitOut, maxitIn, thresh, threshIn)
    outtest$llk
    
    nfoldi=apply(!is.na(cvRSS), 2, sum)
    cvm=apply(cvRSS, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvRSS, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))
    indexCV=which.min(cvm)
    
    # cvm 
    outCV=data.frame(lambda=lambda, cvm=cvm, cvse=cvse, nzeroGA=nzeroGA, nzeroGO=nzeroGO, stringsAsFactors=FALSE)
    
    if (!itrunc) {
      
      if (keep.coef) {
        outputCV=list(Kappa=out$Kappa, Omega=-out$Omega, Omega0=out$Omega0, fit=fit, CV=outCV, lambda.opt=lambda[indexCV], flag=out$flag)   ## , t=out$t
      } else {
        outputCV=list(Kappa=out$Kappa[, indexCV], Omega=-out$Omega[, indexCV], Omega0=out$Omega0[, indexCV], fit=fit, CV=outCV, lambda.opt=lambda[indexCV], flag=out$flag)  ## , t=out$t
      }
 
      return(outputCV)
      
    } else {
      
      #####  Hard-thresholding group  #####
      
      ordConvert=NULL; sr=1
      for (s in 1:(p-1)) {
        for (r in (s+1):p) {
          ordConvert=rbind(ordConvert, c(s,r,sr))
          sr=sr+1
        }
      }
      ordConvert=rbind(0, ordConvert)
      colnames(ordConvert)=c("s","r","sr")
      ordConvert=as.data.frame(ordConvert)
      
      ###  CV non-zero group  ###
      cvllGn0=list(); il0=indexCV; k=1; cvG.min=rep(NA, nlambda)
      
      repeat {
        Kappai=sapply(outi,function(x,ili){x$Kappa[,ili]}, ili=il0)  ## q*p rows and column is folds
        Omegai=sapply(outi,function(x,ili){x$Omega[,ili]}, ili=il0)
        Omega0i=sapply(outi,function(x,ili){x$Omega0[,ili]}, ili=il0)
        
        normAi=matrix(0, nrow=p, ncol=nfolds); normOi=matrix(0, nrow=po, ncol=nfolds)
        
        for (s in 1:p) {
          normAi[s,]=sqrt(apply(Kappai[1:q+q*(s-1),]^2,2,sum))
        }
        for (s in 1:po) {
          normOi[s,]=sqrt(apply(Omegai[1:q+q*(s-1),]^2,2,sum))
        }
        
        
        na=max(apply(normAi!=0,2,sum)); no=max(apply(normOi!=0,2,sum))
        
        orderAi=apply(normAi,2,order,decreasing=T)
        orderOi=apply(normOi,2,order,decreasing=T)
        
        
        cvllj=matrix(0, nrow=na+1, ncol=no+1); k=1
        for (k in 1:nfolds) {
          temid=(foldid==k)
          
          ordai=0
          if (na != 0) {ordai=c(ordai, orderAi[1:na,k])}
          ordoi=0
          if (no != 0) {ordoi=c(ordoi, orderOi[1:no,k])}
          ordoi=data.frame(sr=ordoi)
          ordoi=merge(ordoi, ordConvert, sort=FALSE, all.x=TRUE)
          ordoi=as.matrix(ordoi)
          
          cvllj=cvllj+weighti[k]*cvThrNetGC(M[temid,],X[temid,], cbind(0.0, matrix(Kappai[,k], nrow=q)), cbind(0.0, matrix(Omegai[,k], nrow=q)), Omega0i[,k], ordai, ordoi)
        }
        cvllGn0[[il0]]=cvllj/sum(weighti); #rm(cvllj)
        cvG.min[il0]=min(cvllGn0[[il0]])
        
        
        ## search the neighbours
        
        il1=c(il0-1, il0+1)
        if (il0 == nlambda) {
          il1=c(il0-2, il0-1)
        } else if (il0 == 1) {
          il1=c(il0+1, il0+2)
        }
        
        for (j in 1:2) {
          if (il1[j]>=1 & il1[j]<=nlambda) {
            if (is.na(cvG.min[il1[j]])) {
              Kappai=sapply(outi,function(x,ili){x$Kappa[,ili]}, ili=il1[j])
              Omegai=sapply(outi,function(x,ili){x$Omega[,ili]}, ili=il1[j])
              Omega0i=sapply(outi,function(x,ili){x$Omega0[,ili]}, ili=il1[j])
              
              normAi=matrix(0, nrow=p, ncol=nfolds); normOi=matrix(0, nrow=po, ncol=nfolds)
              for (s in 1:p) {
                normAi[s,]=sqrt(apply(Kappai[1:q+q*(s-1),]^2,2,sum))
              }
              for (s in 1:po) {
                normOi[s,]=sqrt(apply(Omegai[1:q+q*(s-1),]^2,2,sum))
              }
              
              
              na=max(apply(normAi!=0,2,sum)); no=max(apply(normOi!=0,2,sum))
              
              
              orderAi=apply(normAi,2,order,decreasing=T)
              orderOi=apply(normOi,2,order,decreasing=T)
              
              
              cvllj=matrix(0, nrow=na+1, ncol=no+1); k=1
              for (k in 1:nfolds) {
                temid=(foldid==k)
                
                ordai=0
                if (na != 0) {ordai=c(ordai, orderAi[1:na,k])}
                ordoi=0
                if (no != 0) {ordoi=c(ordoi, orderOi[1:no,k])}
                ordoi=data.frame(sr=ordoi)
                ordoi=merge(ordoi, ordConvert, sort=FALSE, all.x=TRUE)
                ordoi=as.matrix(ordoi)
                
                cvllj=cvllj+weighti[k]*cvThrNetGC(M[temid,],X[temid,], cbind(0.0, matrix(Kappai[,k], nrow=q)), cbind(0.0, matrix(Omegai[,k], nrow=q)), Omega0i[,k], ordai, ordoi)
              }
              cvllGn0[[il1[j]]]=cvllj/sum(weighti); #rm(cvllj)
              cvG.min[il1[j]]=min(cvllGn0[[il1[j]]])
              
            }
          } else {
            break
          }
        }
        # if (il1[j]==1 | il1[j]==nlambda)
        #   break
        if (il0==which.min(cvG.min)) {
          break
        } else {
          il0=which.min(cvG.min)
        }
      } ## end of repeat
      
      
      indexG=which.min(cvG.min)
      
      
      tmp=cvllGn0[[indexG]]
      
      if (min(tmp)<0){
        tem=which((1-abs(tmp)/abs(min(tmp)))<=0 & tmp*min(tmp) >0, arr.ind=TRUE)
      }else{
        tem=which((tmp/min(tmp)-1)<=0 & tmp*min(tmp) >0, arr.ind=TRUE)
      }
     
      indexGij=tem[which.min(tem[,1]+tem[,2]),]
      
      outG=cutBetaG(indexGij[1]-1, indexGij[2]-1, out$Kappa[, indexG], out$Omega[, indexG], out$Omega0[, indexG],p, q)
      
      
      outCV.thr=data.frame(lambda=lambda[indexG], nzeroGA=indexGij[1]-1, nzeroGO=indexGij[2]-1, cvm=cvG.min[indexG], stringsAsFactors=FALSE)
      
      
      
      
      #####  Hard-thresholding individual  #####
      il0=indexG
      indexAg=which(apply(matrix(outG$kappa, nrow=q)!=0,2,any))
      indexOg=which(apply(matrix(outG$omega, nrow=q)!=0,2,any))
      
      Kappai=sapply(outi,function(x,ili){x$Kappa[,ili]}, ili=il0)
      Omegai=sapply(outi,function(x,ili){x$Omega[,ili]}, ili=il0)
      Omega0i=sapply(outi,function(x,ili){x$Omega0[,ili]}, ili=il0)
      
      
      temA=matrix(0, nrow=p*q, ncol=nfolds); temO=matrix(0, nrow=po*q, ncol=nfolds)
      if (length(indexAg)>0) {
        for (s in indexAg) {temA[1:q+(s-1)*q, ]=Kappai[1:q+(s-1)*q, ]}
      }
      if (length(indexOg)>0) {
        for (s in indexOg) {temO[1:q+(s-1)*q, ]=Omegai[1:q+(s-1)*q, ]}
      }
      Kappai=temA; Omegai=temO
      na=max(apply(Kappai!=0,2,sum)); no=max(apply(Omegai!=0,2,sum))
      
      absKappai=abs(Kappai); absOmegai=abs(Omegai)
      orderAi=apply(absKappai,2,order,decreasing=T)
      orderOi=apply(absOmegai,2,order,decreasing=T)
      
      cvllIn0=matrix(0, nrow=na+1, ncol=no+1); nz1=na; nz2=no; k=2
      for (nz1 in 1:(na+1)) {
        for (nz2 in 1:(no+1)) {
          cvllk=matrix(NA, nrow=nfolds, ncol=1)
          for (k in 1:nfolds) {
            temid=(foldid==k)
            
            if (nz1<=(p*q)) {
              indexka=which(absKappai[,k]>absKappai[orderAi[nz1,k],k])
            } else {
              indexka=seq(p*q)
            }
            if (nz2<=(po*q)) {
              indexko=which(absOmegai[,k]>absOmegai[orderOi[nz2,k],k])
            } else {
              indexko=seq(po*q)
            }
            
            kappak=rep(0,p*q); omegak=rep(0,po*q)
            if (length(indexka) > 0) {
              for (k1 in 1:length(indexka)) {
                kappak[indexka]=Kappai[indexka,k]
              }
            }
            if (length(indexko) > 0) {
              for (k2 in 1:length(indexko)) {
                omegak[indexko]=Omegai[indexko,k]
              }
            }
            cvllk[k]=PLC(M[temid,],X[temid,], kappak, omegak, Omega0i[,k])
          }
          cvrawi=cvllk; nfoldi=apply(!is.na(cvrawi), 2, sum)
          cvmi=apply(cvrawi, 2, weighted.mean, w=weighti, na.rm=TRUE)
          cvllIn0[nz1,nz2]=cvmi
        }
      }
      
      
      tmp=cvllIn0
      
      if (min(tmp)<0){
        tem=which((1-abs(tmp)/abs(min(tmp)))<=0 & tmp*min(tmp) >0, arr.ind=TRUE)
      }else{
        tem=which((tmp/min(tmp)-1)<=0 & tmp*min(tmp) >0, arr.ind=TRUE)
      } 
      
      indexIij=tem[which.min(tem[,1]+tem[,2]),]

      kappag=outG$kappa
      kappag[abs(kappag)<=sort(abs(kappag),TRUE)[indexIij[1]]]=0
      
      omegag=outG$omega
      omegag[abs(omegag)<=sort(abs(omegag),TRUE)[indexIij[2]]]=0
      omega0g=outG$omega0
      
      ### make omega consistent with paper
      outI=list(kappa=kappag, omega=-omegag, omega0=omega0g)
      
      outCV.thr2=data.frame(lambda=lambda[indexG], nzeroIA=indexIij[1]-1, nzeroIO=indexIij[2]-1, cvm=min(cvllIn0), stringsAsFactors=FALSE)
      
      
      outG$omega=-outG$omega
      
      if (keep.coef) {
        outputCVG=list(Kappa=out$Kappa, Omega=-out$Omega, coef.thrG=outG, coef.thrI=outI, fit=fit, CV=outCV, CVG=cvllGn0[[indexG]], CVI=cvllIn0, CV.thrG=outCV.thr, CV.thrI=outCV.thr2, lambda.opt=lambda[indexCV], lambda.thr=lambda[indexG], flag=out$flag)  ## , t=out$t
      } else {
        outputCVG=list(Kappa=out$Kappa[, indexCV], Omega=-out$Omega[, indexCV], coef.thrG=outG, coef.thrI=outI, fit=fit, CV=outCV, CVG=cvllGn0[[indexG]], CVI=cvllIn0, CV.thrG=outCV.thr, CV.thrI=outCV.thr2, lambda.opt=lambda[indexCV], lambda.thr=lambda[indexG], flag=out$flag)  ## , t=out$t
      }
      

      return(outputCVG)
      
    }
    
  } # CV
  
 
}





############################
#####  Hard threshold  #####
############################

cutBetaG <- function(na, no, kappa, omega, omega0, p, q) {
  po=(p-1)*p/2; na=na+1; no=no+1 # non-zero+1
  Kappai=matrix(kappa, nrow=q); Omegai=matrix(omega, nrow=q)
  normai=sqrt(apply(Kappai^2,2,sum))
  normoi=sqrt(apply(Omegai^2,2,sum))
  
  orderai=order(normai,decreasing=T)
  orderoi=order(normoi,decreasing=T)
  
  if (na<=p) {
    indexka=which(normai>normai[orderai[na]])
  } else {
    indexka=seq(p)
  }
  if (no<=po) {
    indexko=which(normoi>normoi[orderoi[no]])
  } else {
    indexko=seq(po)
  }
  
  kappak=rep(0,p*q); omegak=rep(0,po*q)
  if (length(indexka) > 0) {
    for (k1 in 1:length(indexka)) {
      kappak[1:q+q*(indexka[k1]-1)]=kappa[1:q+q*(indexka[k1]-1)]
    }
  }
  if (length(indexko) > 0) {
    for (k2 in 1:length(indexko)) {
      omegak[1:q+q*(indexko[k2]-1)]=omega[1:q+q*(indexko[k2]-1)]
    }
  }
  
  return(list(kappa=kappak, omega=omegak, omega0=omega0))
}



###############################
########### 2nd stage #########
###############################
library(Matrix)

IniLm=function(x, y){
  
  N0=nrow(x); p=ncol(x)
  nalambda=10
  
  ### Calculation always based on standardized X and centered y
  tem=scaleC(x)
  xscale=tem$sd; x=tem$x; mx=tem$m
  rm(tem)
  
  my=mean(y); y=y-my
  
  wbeta=rep(1, p)
  ### Lambda path
  lambda_max=maxLambdaLmC(x,y,1.0,wbeta,N0)
  lambda_min=ifelse(N0>=p, lambda_max*0.0001, lambda_max*0.01)
  alambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nalambda-1))/(nalambda-1))
  
  repeat {
    outi=EnetLm(x, y, alpha=0.0, lambda=alambda, keep.beta=T)
    if(!is.null(outi))break
    alambda=alambda*2.0
  }
  
  indexi=ncol(outi$Beta)
  beta0=outi$Beta[, indexi]
  wbeta=1/abs(beta0); sgn=sign(beta0[1:p])
  return(list(wbeta=wbeta, sgn=sgn, lambda=alambda[indexi]))
}




######################
###  2nd Stage    ###
######################


EnetLm=function(x, y, alpha=1.0, lambda=NULL, nlambda=100, rlambda=NULL, nfolds=1, foldid=NULL, inzero=TRUE, weight1=NULL, weight2=NULL, adaptive=FALSE, aini=NULL, isd=FALSE, keep.beta=FALSE, thresh=1e-7, maxit=1e+5) {
  # isd=TRUE; thresh=1e-7; maxit=1e+5
  penalty=ifelse(alpha==1, "Lasso", "Enet")
  
  N0=nrow(x); p=ncol(x)
  
  
  ### Adaptive weights based on Ridge (L2)
  if (adaptive) {
    if (is.null(aini)) 
      aini=IniLm(x, y)
    wbeta=aini$wbeta
    rm(aini)
  } else {
    wbeta=rep(1, p)
  }
  
  if (is.null(weight1)) {
    weight1=rep(1,p)
  }
  
  if (is.null(weight2)) {
    weight2=rep(1,p)
  }
  
  
  ### Lambda path
  if (is.null(lambda)) {
    ilambda=1
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p, 0.0001, 0.01)
    }
    lambda=(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    ilambda=0
    nlambda=length(lambda)
  }
  
  #####  Run  #####
  ## wbeta:  weight for L1 penalty, wbeta2: weight for L2 penalty
  
  out=EnetLmC(x, y, alpha, lambda, nlambda, ilambda, wbeta, weight1, weight2, as.integer(isd), p, N0, thresh, maxit)
  nlambdai=out$nlambda ## number of lambdas
  if (nlambdai==0)
    return(NULL)
  lambdai=out$lambda[1:nlambdai]
  
  out$Beta=Matrix(out$Beta[, 1:nlambdai], sparse=TRUE) 
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]
  out$rsq=out$rsq[1:nlambdai]
  
  if (nfolds==1 & is.null(foldid)) {
    fit=data.frame(lambda=lambdai, rsq=out$rsq, nzero=out$nzero)
    return(list(Beta=out$Beta, fit=fit, penalty=penalty, flag=out$flag))
  } else {
    ### Calculation always based on standardized X and centered y
    tem=scaleC(x)
    xscale=tem$sd; x=tem$x; mx=tem$m
    rm(tem)
    
    my=mean(y); y=y-my
    
    ###  Split data for cross-validation
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=N0))
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid)
    N0i=numeric(nfolds); Nf=numeric(nfolds)
    for (i in 1:nfolds) {
      N0i[i]=sum(tb[-i]); Nf[i]=tb[i]
    }
    weighti=as.vector(tapply(rep(1,N0), foldid, sum))
    
    #####  Cross-validation estimates  #####
    outi=list(); cvRSS=matrix(NA, nrow=nfolds, ncol=nlambdai)
    for (i in 1:nfolds) {
      temid=(foldid==i)
      outi[[i]]=cvEnetLmC(x[!temid, ], y[!temid], alpha, lambdai, nlambdai, wbeta, weight1, weight2, as.integer(isd), N0i[i],p, thresh, maxit, x[temid, ], y[temid], Nf[i])     
      cvRSS[i, 1:outi[[i]]$nlambda]=outi[[i]]$RSSp[1:outi[[i]]$nlambda] ## for ith fold
    }
    
    cvRSS=matrix(cvRSS[, 1:nlambdai], ncol=nlambdai)
    cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); rm(cvRSS) #
    cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))
    
    indexi=which.min(cvm)
    indexij=which(cvm<=(cvm[indexi]+cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*";#temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai, rsq=out$rsq, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi, stringsAsFactors=FALSE)
    
    if (!inzero) {
      rm(outi)
      if (!keep.beta) {
        # lambda.1se=lambdai[indexij]
        return(list(Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, flag=out$flag))
      } else {
        return(list(Beta=out$Beta, fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, flag=out$flag))
      }
    }
    
    #####  Cross-validation hard threshold  #####
    il0=indexi; cvm=list(); cv.min=rep(NA, nlambdai)
    repeat {
      numi=out$nzero[il0]
      Betai=sapply(outi, function(x){x$Beta[, il0]})
      Betao=apply(Betai!=0, 2, sum)
      numi2=min(max(Betao), numi)
      
      if (numi2>0) {
        cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
        for (i in 1:nfolds) {
          Betaj=Betai[, i]; temid=foldid==i
          numj=min(Betao[i], numi)
          if (numj==0) {
            cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i])
          } else {
            temo=rank(-abs(Betaj), ties.method="min")
            temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
            temo=temo[order(temo[, 1]), ]
            cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i])
          }
        }
      } else {
        cvRSS=matrix(NA, nrow=nfolds, ncol=1)
        for (i in 1:nfolds) {
          temid=foldid==i
          cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i])
        }
      }
      
      cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum);rm(cvRSS) #
      cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
      cv.min[il0]=min(cvm[[il0]])
      
      il1=c(il0-1, il0+1)
      for (j in 1:2) {
        if (il1[j]>=1 & il1[j]<=nlambdai) {
          if (is.na(cv.min[il1[j]])) {
            numi=out$nzero[il1[j]]
            Betai=sapply(outi, function(x){x$Beta[, il1[j]]})
            Betao=apply(Betai!=0, 2, sum)
            numi2=min(max(Betao), numi)
            
            if (numi2>0) {
              cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
              for (i in 1:nfolds ){
                Betaj=Betai[, i]; temid=foldid==i
                numj=min(Betao[i], numi)
                if (numj==0) {
                  cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i])
                } else {
                  temo=rank(-abs(Betaj), ties.method="min")
                  temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
                  temo=temo[order(temo[, 1]), ]
                  cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i])
                }
              }
            } else {
              cvRSS=matrix(NA, nrow=nfolds, ncol=1)
              for(i in 1:nfolds) {
                temid=foldid==i
                cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i])
              }
            }
            cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum)
            rm(cvRSS)
            cvm[[il1[j]]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
            cv.min[il1[j]]=min(cvm[[il1[j]]])
          }
        } else {
          break
        }
      }
      if (il1[j]==1 | il1[j]==nlambdai)
        break
      if (il0==which.min(cv.min)) {
        break
      } else {
        il0=which.min(cv.min)
      }
    }
    index0=which.min(cv.min)
    
    Beta0=out$Beta[,index0]
    cuti=which.min(cvm[[index0]])
    Beta0[abs(Beta0)<=sort(abs(Beta0),TRUE)[cuti+1]]=0
    
    temCV0=data.frame(lambda=lambdai[index0],cvm=cv.min[index0],nzero=cuti)
    
    if (!keep.beta) {
      return(list(Beta=out$Beta[, indexi], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, flag=out$flag))
    } else {
      return(list(Beta=out$Beta, Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, flag=out$flag))
    }
  }
}
