// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

/*****  Center and standardize  *****/
// [[Rcpp::export]]
List scaleC(Eigen::MatrixXd X){
  int i, p=X.cols(), N=X.rows();
  Eigen::VectorXd mX(p), sdX(p);
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
    X.col(i)/=sdX(i);
  }
  return List::create(Named("x")=X, Named("sd")=sdX, Named("m")=mX);
}


/*****  Lambda path: max  1st stage *****/
// [[Rcpp::export]]
List maxLambdaC1(double tau, Eigen::MatrixXd M, Eigen::MatrixXd X) {
  // M: biomarker nodes and X: covariates
  int p=M.cols(), n=M.rows(), q=X.cols(), po=p*(p-1)/2;
  int s, r, sr, j, k, ij;
  Eigen::VectorXd dli(q), dliAbs(q);
  Eigen::VectorXd SS(q);
  
  Eigen::VectorXd maxLambdasK(p), maxLambdasO(po);
  Eigen::VectorXd lambda2(2), lambdaF(q+1);
  double temA, temB, temC;
  
  // maxLambda - Kappa(z)
  for (s=0; s<p; ++s) {
    dli=X.transpose()*M.col(s)/n; // add minus
    dliAbs=dli.cwiseAbs();    // the coefficient-wise absolute value of dli
    // sort dliAbs in decreasing order; .data(): A pointer to the first element in the array used internally by the vector.
    // std::sort(s.begin(), s.end(), std::greater<int>());
    
    std::sort(dliAbs.data(), dliAbs.data()+dliAbs.size(), std::greater<double>());
    
    for (j=0; j<q; ++j) {
      for (k=0; k<q; ++k) {
        if (dli[k] > dliAbs[j]) {
          SS[k]=dli[k]-dliAbs[j];
        } else if (dli[k] < -dliAbs[j]) {
          SS[k]=dli[k]+dliAbs[j];
        } else {
          SS[k]=0.0;
        }
      }
      lambdaF[j]=SS.squaredNorm()-pow(dliAbs[j]/tau-dliAbs[j], 2)*q;
    }
    lambdaF[q]=dli.squaredNorm();
    
    ij=0;
    while (1) {
      ++ij;
      if (lambdaF[ij] > 0) break; // for tau=1, lambdaF always > 0.
    }
    
    // tau=1, then ij=0
    temA=pow(tau, 2)*ij-pow(1-tau, 2)*q; // = 0 when tau=1
    temB=0.0; temC=0.0;
    for (j=0; j<ij; ++j) {
      temB-=dliAbs[j]; // 2nd order
      temC+=pow(dliAbs[j], 2); // 3rd order
    } // temB=0 and temC= 0 when tau=1 since ij=0. can't go inside the loop
    temB=temB*tau*2.0; // = 0 when tau=1
    if (temA == 0) goto linearA;  // goto linearA when tau=1
    
    lambda2[0]=(-temB+sqrt(temB*temB-temA*temC*4.0))/(temA*2.0);
    lambda2[1]=(-temB-sqrt(temB*temB-temA*temC*4.0))/(temA*2.0);
    if (ij < q) {
      if (lambda2[0] > (dliAbs[ij]/tau) && lambda2[0] < (dliAbs[ij-1]/tau)) {
        maxLambdasK[s]=lambda2[0];
      } else {
        maxLambdasK[s]=lambda2[1];
      }
    } else {
      if (lambda2[0] > (dliAbs[q-1]/tau)) {
        maxLambdasK[s]=lambda2[0];
      } else {
        maxLambdasK[s]=lambda2[1];
      }
    }
    continue;
    
    linearA:
      maxLambdasK[s]=-temC/temB;
  }
  
  // maxLambda - Omega(z)
  sr=0;
  for (s=0; s<(p-1); ++s) {
    for (r=s+1; r<p; ++r) {
      dli=X.transpose()*(M.col(s).array()*M.col(r).array()+M.col(r).array()*M.col(s).array()).matrix();
      dli/=-n; // add minus
      dliAbs=dli.cwiseAbs(); 
      
      std::sort(dliAbs.data(), dliAbs.data()+dliAbs.size(), std::greater<double>());
      
      for (j=0; j<q; ++j) {
        for (k=0; k<q; ++k) {
          if (dli[k] > dliAbs[j]) {
            SS[k]=dli[k]-dliAbs[j];
          } else if (dli[k] < -dliAbs[j]) {
            SS[k]=dli[k]+dliAbs[j];
          } else {
            SS[k]=0.0;
          }
        }
        lambdaF[j]=SS.squaredNorm()-pow(dliAbs[j]/tau-dliAbs[j], 2)*q;
      }
      lambdaF[q]=dli.squaredNorm();
      
      ij=0;
      while (1) {
        ++ij;
        if (lambdaF[ij] > 0) break;
      }
      
      temA=pow(tau, 2)*ij-pow(1-tau, 2)*q;
      temB=0.0; temC=0.0;
      for (j=0; j<ij; ++j) {
        temB-=dliAbs[j];
        temC+=pow(dliAbs[j], 2);
      }
      temB=temB*tau*2.0;
      if (temA == 0) goto linearO;
      
      lambda2[0]=(-temB+sqrt(temB*temB-temA*temC*4.0))/(temA*2.0);
      lambda2[1]=(-temB-sqrt(temB*temB-temA*temC*4.0))/(temA*2.0);
      if (ij < q) {
        if (lambda2[0] > (dliAbs[ij]/tau) && lambda2[0] < (dliAbs[ij-1]/tau)) {
          maxLambdasO[sr]=lambda2[0];
        } else {
          maxLambdasO[sr]=lambda2[1];
        }
      } else {
        if (lambda2[0] > (dliAbs[q-1]/tau)) {
          maxLambdasO[sr]=lambda2[0];
        } else {
          maxLambdasO[sr]=lambda2[1];
        }
      }
      
      ++sr;
      continue;
      
      linearO:
        maxLambdasO[sr]=-temC/temB;
      ++sr;
      
    }
    
    
    
  }
  return List::create(Named("maxLambdasK")=maxLambdasK, Named("maxLambdasO")=maxLambdasO);
}



///////////////////////////////////////////
/////   Network Model for Diagonals   /////
///////////////////////////////////////////

/*****  Network Model  *****/
// [[Rcpp::export]]
List NetLmC(Eigen::VectorXd lambdas, double tau,
            Eigen::MatrixXd M, Eigen::MatrixXd X, double tt,
            int maxitOut, int maxitIn, double thresh, double threshIn){ 
  
  int p=M.cols(), n=M.rows(), q=X.cols(), nlam=lambdas.size(), po=(p-1)*p/2;
  int j, s, r, sr, il;
  int isAdd=0; // countOut=0, countIn=0, 
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlam);
  double lambda1=0.0, lambda2=0.0, tlambda1=0.0, hinge=0.0;
  
  
  Eigen::VectorXd ll=Eigen::VectorXd::Zero(nlam), obj=Eigen::VectorXd::Zero(nlam);
  Eigen::VectorXd ttt=Eigen::VectorXd::Zero(nlam);

  
  Eigen::MatrixXd OutKappa=Eigen::MatrixXd::Zero(p*q, nlam), OutOmega=Eigen::MatrixXd::Zero(po*q, nlam), OutOmega0=Eigen::MatrixXd::Zero(p, nlam);
 
  for (il=0; il<nlam; ++il) {
    double t=tt;
    Eigen::VectorXi isActiveK=Eigen::VectorXi::Zero(p), isActiveO=Eigen::VectorXi::Zero(po);
    
    Eigen::VectorXd theta0=Eigen::VectorXd::Zero(q), theta1=Eigen::VectorXd::Zero(q), dtheta=Eigen::VectorXd::Zero(q), thetaNew=Eigen::VectorXd::Zero(q), SS=Eigen::VectorXd::Zero(q), dli=Eigen::VectorXd::Zero(q);
 
    double zeta1=0.0;
    
    Eigen::MatrixXd Mr=Eigen::MatrixXd::Zero(n, p);
    Eigen::VectorXd Mi=Eigen::VectorXd::Zero(n), Mi2=Eigen::VectorXd::Zero(n), Mi3=Eigen::VectorXd::Zero(n);// prob=Eigen::VectorXd::Zero(n), xb=Eigen::VectorXd::Zero(n);
    
    double llCurrent=0.0, diffllCurrent=0.0;
    double objCurrent=0.0, objNew=0.0, objOld=0.0, lp=0.0;
    double q1=0.0, q2=0.0, q3=0.0;
    
    Eigen::MatrixXd Kappa=Eigen::MatrixXd::Zero(q, p), Omega=Eigen::MatrixXd::Zero(q, po);
    Eigen::VectorXd Omega0=Eigen::VectorXd::Ones(p);
    
    ///  Initialization  ////
    int countOut=0;
    int countIn=0;
    
    Mr=M;
    //prob=prob.array()+0.5;
    llCurrent=M.squaredNorm()/(2.0*n);
    
    objCurrent=llCurrent; objNew=objCurrent;
    
    
    lambda1=lambdas[il]*tau; lambda2=lambdas[il]*(1-tau)*sqrt(q);
    
    ///  Active set  ///
    
    // Active set - Kappa(z)
    for (s=0; s<p; ++s) {
      if (isActiveK[s] == 0) {
        dli=X.transpose()*Mr.col(s)/n; // add minus
        for (j=0; j<q; ++j) {
          if (dli(j) > lambda1) {
            SS[j]=dli(j)-lambda1;
          } else if (dli(j) < -lambda1) {
            SS[j]=dli(j)+lambda1;
          } else {
            SS[j]=0.0;
          }          
        }
        if (SS.norm() > lambda2) {
          isActiveK[s]=1;
        }
      }
    }
 
    
    // Active set - Omega(z)
    sr=0;
    for (s=0; s<(p-1); ++s) {
      for (r=s+1; r<p; ++r) {
        if (isActiveO[sr] == 0) {
          dli=X.transpose()*(Mr.col(s).array()*M.col(r).array()+Mr.col(r).array()*M.col(s).array()).matrix();
          dli/=-n; // add minus
          for (j=0; j<q; ++j) {
            if (dli(j) > lambda1) {
              SS[j]=dli(j)-lambda1;
            } else if (dli(j) < -lambda1) {
              SS[j]=dli(j)+lambda1;
            } else {
              SS[j]=0.0;
            }          
          }
          if (SS.norm() > lambda2) {
            isActiveO[sr]=1;
          }
        }
        ++sr;
      }
    }
    
    // Outer loop
    Outer:
      while (1) {
        ++countOut;
        //llOld=llCurrent;
        objOld=objCurrent;
        
        ///  Network model  ///
        
        // Update Kappa(z)
        //  t=tt; //
        for (s=0; s<p; ++s) {
          
          if (isActiveK[s]==0) goto nextA;
          
          theta0=Kappa.col(s);
          dli=-X.transpose()*Mr.col(s)/n;
          
          countIn=0;
          while (1) {
            ++countIn;
            
            // Obtain t
            if (t < threshIn) t=threshIn;
            while (1) {
              
              // Update
              theta1=Kappa.col(s)-t*dli;
              tlambda1=lambda1*t;
              for (j=0; j<q; ++j) {
                if (theta1[j] > tlambda1) {
                  SS[j]=theta1[j]-tlambda1;
                } else if (theta1[j] < -tlambda1) {
                  SS[j]=theta1[j]+tlambda1;
                } else {
                  SS[j]=0.0;
                }
              }
              hinge=1.0-(t*lambda2/SS.norm());
              if (hinge > 0.0) {
                theta1=hinge*SS;
              } else {
                theta1.setZero(q);
                if ((Kappa.col(s)).squaredNorm() == 0.0) goto nextA;
              }
              
              dtheta=theta1-Kappa.col(s);
              Mi=X*dtheta;
              
              diffllCurrent=Mi.dot(2.0*Mr.col(s)-Mi/Omega0[s])/(2.0*n)+dtheta.dot(dli)+dtheta.dot(dtheta)/(t*2.0);
              if (diffllCurrent > 0.0) break;
              t=t*0.8;
              if (t < threshIn) goto nextA;
            }
            
            // Backtracking
            thetaNew=theta1+(theta1-theta0)*countIn/(countIn+3.0);
            dtheta=thetaNew-Kappa.col(s);
            
            Mi=X*dtheta;
            llCurrent=llCurrent-Mi.dot(2.0*Mr.col(s)-Mi/Omega0[s])/(2.0*n);
            Mr.col(s)=Mr.col(s)-Mi/Omega0[s];
            lp+=lambda1*(thetaNew.cwiseAbs().sum()-Kappa.col(s).cwiseAbs().sum())+lambda2*(thetaNew.norm()-Kappa.col(s).norm());
            objNew=llCurrent+lp;
            
            Kappa.col(s)=thetaNew;
            theta0=theta1;
            dli=-X.transpose()*Mr.col(s)/n;
            if (countIn >= maxitIn) {
              objCurrent=objNew;
              break;
            } else if (fabs(objNew-objCurrent) < fabs(threshIn*objCurrent)) {
              //fabs(llNew-llCurrent) < fabs(threshIn*llCurrent)
              objCurrent=objNew; //llCurrent=llNew;
              break;
            } else {
              objCurrent=objNew;
            }
            
          }
          
          nextA:;
        } // update Kappa(z)  end of inner loop of Kappa(z)
        
        
        // Update Omega(z) - upper triangular
        //  t=tt; //
        sr=0;
        for (s=0; s<(p-1); ++s) {
          for (r=s+1; r<p; ++r) {
            if (isActiveO[sr]==0) goto nextO;
            
            theta0=Omega.col(sr);
            dli=X.transpose()*(Mr.col(s).array()*M.col(r).array()+Mr.col(r).array()*M.col(s).array()).matrix();
            dli/=n;
            
            countIn=0;
            while (1) {
              ++countIn;
              
              // Obtain t
              if (t < threshIn) t=threshIn;
              while (1) {
                
                // Update
                theta1=Omega.col(sr)-t*dli;
                tlambda1=lambda1*t;
                for (j=0; j<q; ++j) {
                  if (theta1[j] > tlambda1) {
                    SS[j]=theta1[j]-tlambda1;
                  } else if (theta1[j] < -tlambda1) {
                    SS[j]=theta1[j]+tlambda1;
                  } else {
                    SS[j]=0.0;
                  }
                }
                hinge=1.0-(t*lambda2/SS.norm());
                if (hinge > 0.0) {
                  theta1=hinge*SS;
                } else {
                  theta1.setZero(q);
                  if ((Omega.col(sr)).squaredNorm() == 0.0) {
                    goto nextO;
                  }
                }
                
                dtheta=theta1-Omega.col(sr);
                Mi3=X*dtheta;
                Mi=Mi3.array()*M.col(r).array();
                Mi2=Mi3.array()*M.col(s).array();
                
                diffllCurrent=-(Mi.dot(2.0*Mr.col(s)+Mi/Omega0[s])+Mi2.dot(2.0*Mr.col(r)+Mi2/Omega0[r]))/(2.0*n)+dtheta.dot(dli)+dtheta.dot(dtheta)/(t*2.0);
                if (diffllCurrent > 0.0) break;
                t=t*0.8;
                if (t < threshIn) goto nextO;
              }
              
              // Backtracking
              thetaNew=theta1+(theta1-theta0)*countIn/(countIn+3.0);
              
              dtheta=thetaNew-Omega.col(sr);
              
              Mi3=X*dtheta;
              Mi=Mi3.array()*M.col(r).array();
              Mi2=Mi3.array()*M.col(s).array();
              
              llCurrent=llCurrent+(Mi.dot(2.0*Mr.col(s)+Mi/Omega0[s])+Mi2.dot(2.0*Mr.col(r)+Mi2/Omega0[r]))/(2.0*n);
              Mr.col(s)=Mr.col(s)+Mi/Omega0[s];
              Mr.col(r)=Mr.col(r)+Mi2/Omega0[r];
              lp+=lambda1*(thetaNew.cwiseAbs().sum()-Omega.col(sr).cwiseAbs().sum())+lambda2*(thetaNew.norm()-Omega.col(sr).norm());
              objNew=llCurrent+lp;
              
              Omega.col(sr)=thetaNew;
              theta0=theta1;
              dli=X.transpose()*(Mr.col(s).array()*M.col(r).array()+Mr.col(r).array()*M.col(s).array()).matrix();
              dli/=n;
              if (countIn >= maxitIn) {
                objCurrent=objNew;
                break;
              } else if (fabs(objNew-objCurrent) < fabs(threshIn*objCurrent)) {
                //fabs(llNew-llCurrent) < fabs(threshIn*llCurrent)
                objCurrent=objNew;
                break;
              } else {
                objCurrent=objNew;
              }
            }
            
            nextO:;
   
            ++sr;
        } // update Omega(z)
      } // update Omega(z)
        
        
        // Update Omega0 - diagonals
        for (s=0; s<p; ++s) {
          Mi=Omega0[s]*(M.col(s)-Mr.col(s));
          q1=M.col(s).squaredNorm();
          q2=-n;
          q3=-Mi.squaredNorm();
          
          zeta1=(-q2+sqrt(q2*q2-4.0*q1*q3))/(2*q1);
          Mi=M.col(s)-Mi/zeta1;
          llCurrent=llCurrent+zeta1/(2.0*n)*Mi.squaredNorm()-Omega0[s]/(2.0*n)*Mr.col(s).squaredNorm()+0.5*(log(Omega0[s])-log(zeta1));
          
          Mr.col(s)=Mi;
          Omega0[s]=zeta1;
        } // update Omega0
        //objCurrent=llCurrent+lp;
        
        objCurrent=llCurrent+lp;
        
  
        if (fabs(objOld-objCurrent) < fabs(thresh*objOld)){flag[il]=1; break;}
        if (countOut >= maxitOut){flag[il]=0; goto nextL;}
  } // end of outer loop
      
      
      ///  Check active sets  ///
      isAdd=0;
    
    // Active set - Kappa(z)
    for (s=0; s<p; ++s) {
      if (isActiveK[s] == 0) {
        dli=X.transpose()*Mr.col(s)/n; // add minus
        for (j=0; j<q; ++j) {
          if (dli(j) > lambda1) {
            SS[j]=dli(j)-lambda1;
          } else if (dli(j) < -lambda1) {
            SS[j]=dli(j)+lambda1;
          } else {
            SS[j]=0.0;
          }          
        }
        if (SS.norm() > lambda2) {
          isActiveK[s]=1; isAdd=1;
        }
      }
    }
    
    
    
    // Active set - Omega(z)
    sr=0;
    for (s=0; s<(p-1); ++s) {
      for (r=s+1; r<p; ++r) {
        if (isActiveO[sr] == 0) {
          dli=X.transpose()*(Mr.col(s).array()*M.col(r).array()+Mr.col(r).array()*M.col(s).array()).matrix();
          dli/=-n; // add minus
          for (j=0; j<q; ++j) {
            if (dli(j) > lambda1) {
              SS[j]=dli(j)-lambda1;
            } else if (dli(j) < -lambda1) {
              SS[j]=dli(j)+lambda1;
            } else {
              SS[j]=0.0;
            }          
          }
          if (SS.norm() > lambda2) {
            isActiveO[sr]=1; isAdd=1;
          }
        }
        ++sr;
      }
    }
    
    
    if (isAdd == 1) {
      goto Outer;
    }
    
    nextL:
    ll[il]=llCurrent;
    obj[il]=objCurrent;
    
    OutKappa.col(il)=Kappa;
    OutOmega.col(il)=Omega;
    OutOmega0.col(il)=Omega0;
    
    ttt[il]=t;
  } // loop:lambda
  
  return List::create(Named("Kappa")=OutKappa, Named("Omega")=OutOmega, Named("Omega0")=OutOmega0,
                            Named("flag")=flag, Named("ttt")=ttt,
                            Named("lambda")=lambdas, Named("ll")=ll, Named("obj")=obj); // Named("countOut")=countOut, 11/23/16
}




// *********  Network Model cross-validation 1st stage ********* //

// [[Rcpp::export]]
List NetLmCvC(Eigen::VectorXd lambdas, double tau,
              Eigen::MatrixXd M, Eigen::MatrixXd X, 
              Eigen::MatrixXd Mk, Eigen::MatrixXd Xk, double tt,
              int maxitOut, int maxitIn, double thresh, double threshIn){ 
  
  int p=M.cols(), n=M.rows(), nk=Mk.rows(), q=X.cols(), nlam=lambdas.size(), po=(p-1)*p/2;
  int j, s, r, sr, il;
  int isAdd=0; //countOut=0, countIn=0, 11/23/16
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlam);
  
  
  
  Eigen::VectorXd ll=Eigen::VectorXd::Zero(nlam), obj=Eigen::VectorXd::Zero(nlam);
  Eigen::VectorXd llk=Eigen::VectorXd::Zero(nlam);
  
  Eigen::VectorXd ttt=Eigen::VectorXd::Zero(nlam);
  
  Eigen::MatrixXd OutKappa=Eigen::MatrixXd::Zero(p*q, nlam), OutOmega=Eigen::MatrixXd::Zero(po*q, nlam), OutOmega0=Eigen::MatrixXd::Zero(p, nlam);
  
  
  for (il=0; il<nlam; ++il) {
    double t=tt;
    Eigen::VectorXi isActiveK=Eigen::VectorXi::Zero(p), isActiveO=Eigen::VectorXi::Zero(po);

    Eigen::VectorXd theta0=Eigen::VectorXd::Zero(q), theta1=Eigen::VectorXd::Zero(q), dtheta=Eigen::VectorXd::Zero(q), thetaNew=Eigen::VectorXd::Zero(q), SS=Eigen::VectorXd::Zero(q), dli=Eigen::VectorXd::Zero(q);

    double zeta1=0.0;
    
    Eigen::MatrixXd Mr=Eigen::MatrixXd::Zero(n, p), Mkr=Eigen::MatrixXd::Zero(nk,p);
    Eigen::VectorXd Mi=Eigen::VectorXd::Zero(n), Mi2=Eigen::VectorXd::Zero(n), Mi3=Eigen::VectorXd::Zero(n);
    
    double llCurrent=0.0, diffllCurrent=0.0;
    double objCurrent=0.0, objNew=0.0, objOld=0.0, lp=0.0;
    double q1=0.0, q2=0.0, q3=0.0;
    double lambda1=0.0, lambda2=0.0, tlambda1=0.0, hinge=0.0;
    
    
    Eigen::MatrixXd Kappa=Eigen::MatrixXd::Zero(q, p), Omega=Eigen::MatrixXd::Zero(q, po);
    Eigen::VectorXd Omega0=Eigen::VectorXd::Ones(p);
    
    ///  Initialization  ////
    int countOut=0;
    int countIn=0;
    
    Mr=M;
    //prob=prob.array()+0.5;
    llCurrent=M.squaredNorm()/(2.0*n);
    
    objCurrent=llCurrent; objNew=objCurrent;
    
    
    lambda1=lambdas[il]*tau; lambda2=lambdas[il]*(1-tau)*sqrt(q);
   
    
    ///  Active set  ///
    
    // Active set - Kappa(z)
    for (s=0; s<p; ++s) {
      if (isActiveK[s] == 0) {
        dli=X.transpose()*Mr.col(s)/n; // add minus
        for (j=0; j<q; ++j) {
          if (dli(j) > lambda1) {
            SS[j]=dli(j)-lambda1;
          } else if (dli(j) < -lambda1) {
            SS[j]=dli(j)+lambda1;
          } else {
            SS[j]=0.0;
          }          
        }
        if (SS.norm() > lambda2) {
          isActiveK[s]=1;
        }
      }
    }
    //isActiveA.setOnes();
    
    // Active set - Omega(z)
    sr=0;
    for (s=0; s<(p-1); ++s) {
      for (r=s+1; r<p; ++r) {
        if (isActiveO[sr] == 0) {
          dli=X.transpose()*(Mr.col(s).array()*M.col(r).array()+Mr.col(r).array()*M.col(s).array()).matrix();
          dli/=-n; // add minus
          for (j=0; j<q; ++j) {
            if (dli(j) > lambda1) {
              SS[j]=dli(j)-lambda1;
            } else if (dli(j) < -lambda1) {
              SS[j]=dli(j)+lambda1;
            } else {
              SS[j]=0.0;
            }          
          }
          if (SS.norm() > lambda2) {
            isActiveO[sr]=1;
          }
        }
        ++sr;
      }
    }
    
    
    // Outer loop
    Outer:
      while (1) {
        ++countOut;
        //llOld=llCurrent;
        objOld=objCurrent;
        
        ///  Network model  ///
        
      
        
        // Update Kappa(z)
        //   t=tt; //
        for (s=0; s<p; ++s) {
          if (isActiveK[s]==0) goto nextA;
          
          theta0=Kappa.col(s);
          dli=-X.transpose()*Mr.col(s)/n;
          
          countIn=0;
          while (1) {
            ++countIn;
            
            // Obtain t
            if (t < threshIn) t=threshIn;
            while (1) {
              
              // Update
              theta1=Kappa.col(s)-t*dli;
              tlambda1=lambda1*t;
              for (j=0; j<q; ++j) {
                if (theta1[j] > tlambda1) {
                  SS[j]=theta1[j]-tlambda1;
                } else if (theta1[j] < -tlambda1) {
                  SS[j]=theta1[j]+tlambda1;
                } else {
                  SS[j]=0.0;
                }
              }
              hinge=1.0-(t*lambda2/SS.norm());
              if (hinge > 0.0) {
                theta1=hinge*SS;
              } else {
                theta1.setZero(q);
                if ((Kappa.col(s)).squaredNorm() == 0.0) goto nextA;
              }
              
              dtheta=theta1-Kappa.col(s);
              Mi=X*dtheta;
              
              diffllCurrent=Mi.dot(2.0*Mr.col(s)-Mi/Omega0[s])/(2.0*n)+dtheta.dot(dli)+dtheta.dot(dtheta)/(t*2.0);
              if (diffllCurrent > 0.0) break;
              t=t*0.8;
              if (t < threshIn) goto nextA;
            }
            
            // Backtracking
            thetaNew=theta1+(theta1-theta0)*countIn/(countIn+3.0);
            dtheta=thetaNew-Kappa.col(s);
            
            Mi=X*dtheta;
            llCurrent=llCurrent-Mi.dot(2.0*Mr.col(s)-Mi/Omega0[s])/(2.0*n);
            Mr.col(s)=Mr.col(s)-Mi/Omega0[s];
            lp+=lambda1*(thetaNew.cwiseAbs().sum()-Kappa.col(s).cwiseAbs().sum())+lambda2*(thetaNew.norm()-Kappa.col(s).norm());
            objNew=llCurrent+lp;
            
            Kappa.col(s)=thetaNew;
            theta0=theta1;
            dli=-X.transpose()*Mr.col(s)/n;
            if (countIn >= maxitIn) {
              objCurrent=objNew;
              break;
            } else if (fabs(objNew-objCurrent) < fabs(threshIn*objCurrent)) {
              //fabs(llNew-llCurrent) < fabs(threshIn*llCurrent)
              objCurrent=objNew; //llCurrent=llNew;
              break;
            } else {
              objCurrent=objNew;
            }
            
          }
          
          nextA:;
        } // update Kappa(z)  end of inner loop of Kappa(z)
        
        
        // Update Omega(z) - upper triangular
        //  t=tt; //
        sr=0;
        for (s=0; s<(p-1); ++s) {
          for (r=s+1; r<p; ++r) {
            if (isActiveO[sr]==0) goto nextO;
            
            theta0=Omega.col(sr);
            dli=X.transpose()*(Mr.col(s).array()*M.col(r).array()+Mr.col(r).array()*M.col(s).array()).matrix();
            dli/=n;
            
            countIn=0;
            while (1) {
              ++countIn;
              
              // Obtain t
              if (t < threshIn) t=threshIn;
              while (1) {
                
                // Update
                theta1=Omega.col(sr)-t*dli;
                tlambda1=lambda1*t;
                for (j=0; j<q; ++j) {
                  if (theta1[j] > tlambda1) {
                    SS[j]=theta1[j]-tlambda1;
                  } else if (theta1[j] < -tlambda1) {
                    SS[j]=theta1[j]+tlambda1;
                  } else {
                    SS[j]=0.0;
                  }
                }
                hinge=1.0-(t*lambda2/SS.norm());
                if (hinge > 0.0) {
                  theta1=hinge*SS;
                } else {
                  theta1.setZero(q);
                  if ((Omega.col(sr)).squaredNorm() == 0.0) {
                    goto nextO;
                  }
                }
                
                dtheta=theta1-Omega.col(sr);
                Mi3=X*dtheta;
                Mi=Mi3.array()*M.col(r).array();
                Mi2=Mi3.array()*M.col(s).array();
                
                diffllCurrent=-(Mi.dot(2.0*Mr.col(s)+Mi/Omega0[s])+Mi2.dot(2.0*Mr.col(r)+Mi2/Omega0[r]))/(2.0*n)+dtheta.dot(dli)+dtheta.dot(dtheta)/(t*2.0);
                if (diffllCurrent > 0.0) break;
                t=t*0.8;
                if (t < threshIn) goto nextO;
              }
              
              // Backtracking
              thetaNew=theta1+(theta1-theta0)*countIn/(countIn+3.0);
              
              dtheta=thetaNew-Omega.col(sr);
              
              Mi3=X*dtheta;
              Mi=Mi3.array()*M.col(r).array();
              Mi2=Mi3.array()*M.col(s).array();
              
              llCurrent=llCurrent+(Mi.dot(2.0*Mr.col(s)+Mi/Omega0[s])+Mi2.dot(2.0*Mr.col(r)+Mi2/Omega0[r]))/(2.0*n);
              Mr.col(s)=Mr.col(s)+Mi/Omega0[s];
              Mr.col(r)=Mr.col(r)+Mi2/Omega0[r];
              lp+=lambda1*(thetaNew.cwiseAbs().sum()-Omega.col(sr).cwiseAbs().sum())+lambda2*(thetaNew.norm()-Omega.col(sr).norm());
              objNew=llCurrent+lp;
              
              Omega.col(sr)=thetaNew;
              theta0=theta1;
              dli=X.transpose()*(Mr.col(s).array()*M.col(r).array()+Mr.col(r).array()*M.col(s).array()).matrix();
              dli/=n;
              if (countIn >= maxitIn) {
                objCurrent=objNew;
                break;
              } else if (fabs(objNew-objCurrent) < fabs(threshIn*objCurrent)) {
                //fabs(llNew-llCurrent) < fabs(threshIn*llCurrent)
                objCurrent=objNew;
                break;
              } else {
                objCurrent=objNew;
              }
            }
            
            nextO:;
          
            ++sr;
        } // update Omega(z)
      } // update Omega(z)
        
        
        // Update Omega0 - diagonals
        for (s=0; s<p; ++s) {
          Mi=Omega0[s]*(M.col(s)-Mr.col(s));
          q1=M.col(s).squaredNorm();
          q2=-n;
          q3=-Mi.squaredNorm();
          
          zeta1=(-q2+sqrt(q2*q2-4.0*q1*q3))/(2*q1);
          Mi=M.col(s)-Mi/zeta1;
          llCurrent=llCurrent+zeta1/(2.0*n)*Mi.squaredNorm()-Omega0[s]/(2.0*n)*Mr.col(s).squaredNorm()+0.5*(log(Omega0[s])-log(zeta1));
          
          Mr.col(s)=Mi;
          Omega0[s]=zeta1;
        } // update Omega0
        //objCurrent=llCurrent+lp;
        
        objCurrent=llCurrent+lp;
        
        //if (fabs(llOld-llCurrent) < fabs(thresh*llOld)){flag[il]=1; break;}
        if (fabs(objOld-objCurrent) < fabs(thresh*objOld)){flag[il]=1; break;}
        if (countOut >= maxitOut){flag[il]=0; goto nextL;}
  } // end of outer loop
      
      
      ///  Check active sets  ///
      isAdd=0;
    
    // Active set - Kappa(z)
    for (s=0; s<p; ++s) {
      if (isActiveK[s] == 0) {
        dli=X.transpose()*Mr.col(s)/n; // add minus
        for (j=0; j<q; ++j) {
          if (dli(j) > lambda1) {
            SS[j]=dli(j)-lambda1;
          } else if (dli(j) < -lambda1) {
            SS[j]=dli(j)+lambda1;
          } else {
            SS[j]=0.0;
          }          
        }
        if (SS.norm() > lambda2) {
          isActiveK[s]=1; isAdd=1;
        }
      }
    }
    
    
    
    // Active set - Omega(z)
    sr=0;
    for (s=0; s<(p-1); ++s) {
      for (r=s+1; r<p; ++r) {
        if (isActiveO[sr] == 0) {
          dli=X.transpose()*(Mr.col(s).array()*M.col(r).array()+Mr.col(r).array()*M.col(s).array()).matrix();
          dli/=-n; // add minus
          for (j=0; j<q; ++j) {
            if (dli(j) > lambda1) {
              SS[j]=dli(j)-lambda1;
            } else if (dli(j) < -lambda1) {
              SS[j]=dli(j)+lambda1;
            } else {
              SS[j]=0.0;
            }          
          }
          if (SS.norm() > lambda2) {
            isActiveO[sr]=1; isAdd=1;
          }
        }
        ++sr;
      }
    }
    
    
    if (isAdd == 1) {
      goto Outer;
    }
    
    nextL:
    ll[il]=llCurrent;
    obj[il]=objCurrent;
    
    OutKappa.col(il)=Kappa;
    OutOmega.col(il)=Omega;
    OutOmega0.col(il)=Omega0;
    
    
    Mkr=Mk;
    for (s=0; s<p; ++s) {
      //  if (isActiveA[s] == 1) {
      Mkr.col(s)-=Xk*Kappa.col(s)/Omega0[s];
      // }
    }
    
    sr=0;
    for (s=0; s<(p-1); ++s) {
      for (r=s+1; r<p; ++r) {
        Mkr.col(s)+=((Xk*Omega.col(sr)).array()*Mk.col(r).array()).matrix()/Omega0[s];
        Mkr.col(r)+=((Xk*Omega.col(sr)).array()*Mk.col(s).array()).matrix()/Omega0[r];
        ++sr; ///// add this Xiang
      }
    }
    
    for (s=0; s<p; ++s) {
      llk[il]+=Mkr.col(s).squaredNorm()*Omega0[s]/(2.0*nk)-0.5*log(Omega0[s]);
    }
    
   
    ttt[il]=t;
  } // loop:lambda
  
  return List::create(Named("Kappa")=OutKappa, Named("Omega")=OutOmega, Named("Omega0")=OutOmega0,
                            Named("flag")=flag, Named("ttt")=ttt,
                            Named("lambda")=lambdas, Named("ll")=ll, Named("llk")=llk,Named("obj")=obj); //Named("countOut")=countOut, 
}



//////////////////////////////////
/////   Number of non-zero   /////
//////////////////////////////////

/*****  Hard-thresholding - Network Group  *****/
// [[Rcpp::export]]
Eigen::MatrixXd cvThrNetGC(Eigen::MatrixXd M, Eigen::MatrixXd X,
                           Eigen::MatrixXd Kappa, Eigen::MatrixXd Omega, Eigen::VectorXd Omega0,
                           Eigen::VectorXi ordA, Eigen::MatrixXi ordO){ 
  
  int na=ordA.size(), no=ordO.rows();
  int p=M.cols(), n=M.rows();
  int  i, j, s, r, sr;
  Eigen::MatrixXd ll=Eigen::MatrixXd::Zero(na, no);
  Eigen::MatrixXd Mr=Eigen::MatrixXd::Zero(n, p);
  Eigen::VectorXd Mi=Eigen::VectorXd::Zero(n), Mi2=Eigen::VectorXd::Zero(n), Mi3=Eigen::VectorXd::Zero(n);
  
  for (i=0; i<na; ++i) {
    
    if (i==0) {
      Mr=M;
    } else {
      Mr=M;
      for (j=i; j>0; --j) {
        Mr.col(ordA[j]-1)-=X*Kappa.col(ordA[j])/Omega0[ordA[j]-1];
      }
    }
    
    for (j=0; j<no; ++j) {
      if (j==0) {
        for (s=0; s<p; ++s) {
          ll(i, 0)+=Mr.col(s).squaredNorm()*Omega0[s]*0.5/n-0.5*log(Omega0[s]);
        }
      } else {
        s=ordO(j,1); r=ordO(j,2); sr=ordO(j,0);
        
        Mi=X*Omega.col(sr);
        Mi2=Mi.array()*M.col(r-1).array();
        Mi3=Mi.array()*M.col(s-1).array();
        
        ll(i, j)=ll(i, j-1)+(Mi2.dot(2.0*Mr.col(s-1)+Mi2/Omega0[s-1])+Mi3.dot(2.0*Mr.col(r-1)+Mi3/Omega0[r-1]))/(2.0*n);
        Mr.col(s-1)+=Mi2/Omega0[s-1];
        Mr.col(r-1)+=Mi3/Omega0[r-1];
      }
    }
  }
  
  return(ll);
}



/*****  pseudo-likelihood (least squared)  *****/
// [[Rcpp::export]]
double PLC(Eigen::MatrixXd M, Eigen::MatrixXd X,
           Eigen::VectorXd kappa, Eigen::VectorXd omega, Eigen::VectorXd omega0){ 
  
  int p=M.cols(), n=M.rows(), q=X.cols();
  int  j, s, r, sr, po=(p-1)*p/2;
  Eigen::MatrixXd Kappa=Eigen::MatrixXd::Zero(q, p), Omega=Eigen::MatrixXd::Zero(q, po);
  Eigen::MatrixXd Mr=Eigen::MatrixXd::Zero(n, p);
  double ll=0.0;
  
  for (s=0; s<p; ++s) {
    for (j=0; j<q; ++j) {
      Kappa(j,s)=kappa(j+s*q);
    }
  }
  for (s=0; s<po; ++s) {
    for (j=0; j<q; ++j) {
      Omega(j,s)=omega(j+s*q);
    }
  }
  
  Mr=M;
  
  for (s=0; s<p; ++s) {
    Mr.col(s)-=X*Kappa.col(s)/omega0[s];
  }
  sr=0;
  for (s=0; s<(p-1); ++s) {
    for (r=s+1; r<p; ++r) {
      Mr.col(s)+=((X*Omega.col(sr)).array()*M.col(r).array()).matrix()/omega0[s];
      Mr.col(r)+=((X*Omega.col(sr)).array()*M.col(s).array()).matrix()/omega0[r];
      ++sr;
    }
  }
  
  for (s=0; s<p; ++s) {
    ll+=omega0[s]*Mr.col(s).squaredNorm();
  }
  ll/=(2.0*n);
  
  
  // add omega0 part 
  for (s=0; s<p; ++s){
    ll-=0.5*log(omega0[s]);
  }
  
  return(ll);  // it is ll instead of RSS 
}


/*****  pseudo-likelihood in 1st stage, cut  *****/
// [[Rcpp::export]]
double cvTrimPLC(Eigen::MatrixXd M, Eigen::MatrixXd X,
                 Eigen::VectorXd kappa, Eigen::VectorXd omega, Eigen::VectorXd omega0){ 
  
  int p=M.cols(), n=M.rows(), q=X.cols();
  int  j, s, r, sr, po=(p-1)*p/2;
  Eigen::MatrixXd Kappa=Eigen::MatrixXd::Zero(q, p), Omega=Eigen::MatrixXd::Zero(q, po);
  Eigen::MatrixXd Mr=Eigen::MatrixXd::Zero(n, p);
  double ll=0.0;
  
  for (s=0; s<p; ++s) {
    for (j=0; j<q; ++j) {
      Kappa(j,s)=kappa(j+s*q);
    }
  }
  for (s=0; s<po; ++s) {
    for (j=0; j<q; ++j) {
      Omega(j,s)=omega(j+s*q);
    }
  }
  
  Mr=M;
  
  for (s=0; s<p; ++s) {
    Mr.col(s)-=X*Kappa.col(s)/omega0[s];
  }
  sr=0;
  for (s=0; s<(p-1); ++s) {
    for (r=s+1; r<p; ++r) {
      Mr.col(s)+=((X*Omega.col(sr)).array()*M.col(r).array()).matrix()/omega0[s];
      Mr.col(r)+=((X*Omega.col(sr)).array()*M.col(s).array()).matrix()/omega0[r];
      ++sr;
    }
  }
  
  for (s=0; s<p; ++s) {
    ll+=omega0[s]*Mr.col(s).squaredNorm();
  }
  ll/=(2.0*n);
  
  
  // add omega0 part 
  for (s=0; s<p; ++s){
    ll-=0.5*log(omega0[s]);
  }
  
  return(ll);  // it is ll instead of RSS 
}



/*****  pseudo-likelihood (objective function)  *****/
// [[Rcpp::export]]
double ObC(Eigen::MatrixXd M, Eigen::MatrixXd X,
           Eigen::VectorXd kappa, Eigen::VectorXd omega, Eigen::VectorXd omega0, double lambda){ 
  
  int p=M.cols(), n=M.rows(), q=X.cols();
  int  j, s, r, sr, po=(p-1)*p/2;
  Eigen::MatrixXd Kappa=Eigen::MatrixXd::Zero(q, p), Omega=Eigen::MatrixXd::Zero(q, po);
  Eigen::MatrixXd Mr=Eigen::MatrixXd::Zero(n, p);
  double ll=0.0;
  double out=0.0;
  
  for (s=0; s<p; ++s) {
    for (j=0; j<q; ++j) {
      Kappa(j,s)=kappa(j+s*q);
    }
  }
  for (s=0; s<po; ++s) {
    for (j=0; j<q; ++j) {
      Omega(j,s)=omega(j+s*q);
    }
  }
  
  Mr=M;
  
  for (s=0; s<p; ++s) {
    Mr.col(s)-=X*Kappa.col(s)/omega0[s];
  }
  sr=0;
  for (s=0; s<(p-1); ++s) {
    for (r=s+1; r<p; ++r) {
      Mr.col(s)+=((X*Omega.col(sr)).array()*M.col(r).array()).matrix()/omega0[s];
      Mr.col(r)+=((X*Omega.col(sr)).array()*M.col(s).array()).matrix()/omega0[r];
      ++sr;
    }
  }
  
  for (s=0; s<p; ++s) {
    ll+=omega0[s]*Mr.col(s).squaredNorm();
  }
  ll/=(2.0*n);
  
  
  // add omega0 part 
  for (s=0; s<p; ++s){
    ll-=0.5*log(omega0[s]);
  }
  
  
  for (j=0; j<((p+(p-1)*p/2)); ++j) {
    out=out+lambda*(kappa.cwiseAbs().sum()+omega.cwiseAbs().sum());
  }
  out=out+ll;
  
  return(out);  // it is ll instead of RSS 
}





//////////////////////////////////////
/////     2nd stage model       /////
/////////////////////////////////////


/*****  LM: Lambda path (max) inner product <xj,y> *****/
// [[Rcpp::export]]
double maxLambdaLmC(Eigen::MatrixXd X, Eigen::VectorXd y, double alpha, Eigen::VectorXd wbeta, int N0){
  Eigen::VectorXd Li(N0);
  
  Li=(y.transpose()*X).cwiseAbs()/N0;   // <xj,y>/N0
  Li=Li.array()/wbeta.array()/alpha ;
  return Li.maxCoeff();
}



/*****  Used for CV trimming  *****/   //  
// [[Rcpp::export]]
Eigen::VectorXd cvTrimLmC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco, 
                          Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF) {
  int i, j;
  Eigen::VectorXd RSS, xbF=Eigen::VectorXd::Zero(NF); // xb=Eigen::VectorXd::Zero(N),
  
  if(nn2>0){
    RSS.setZero(nn2); //nn2= # of part of data
    
    if(nn==0){
      RSS(0)=yF.squaredNorm();
    }else{
      for(i=0;i<nn;i++){
        j=loco(i); //   index of nonzero beta
        xbF+=XF.col(j)*beta(i); // 
        RSS(i)=(yF-xbF).squaredNorm();
      }
    }
    if(nn2>nn){for(i=nn;i<nn2;i++){RSS(i)=RSS(nn-1);}}
  }else{
    RSS.setZero(1);
    RSS(0)=yF.squaredNorm();
  }
  
  return(RSS);
}



/*****  LM: Enet (L1+L2)  *****/
// [[Rcpp::export]]
List EnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y,
             double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta, Eigen::VectorXd weight1,
             Eigen::VectorXd weight2,int isd, int p, int N0, double thresh, int maxit){ 
  
  int  i, j, it=0, il, iadd, ia=0; 
  double zi, obj0, obj1, b0, db0, objQi, objQj, rss0, rss1, RSS0; // lambda2,
  double lambdaMax;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd RSS(nlambda), RSQ(nlambda);
  double xr, dbMax, thresh2=1e-5;
  Eigen::VectorXd mX(p), sdX(p), di(p);
  
  if (isd == 0) {
    for (i=0;i<p;++i) {
      mX(i)=X.col(i).mean();
      X.col(i)=X.col(i).array()-mX(i);
      sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
      X.col(i)/=sdX(i);
    }
  }
  
  y=y.array()-y.mean();
  
  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLmC(X, y, alpha, wbeta, N0);
    } else {
      lambdaMax=maxLambdaLmC(X, y, 0.001, wbeta, N0);
    }
    lambda=lambda.array()*lambdaMax;
  }
  
  RSS0=y.squaredNorm();
  obj0=RSS0/N0/2.0; rss0=RSS0;
  
  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N0);
  }
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta.array()*weight1.array();
    lambda2=lambda(il)*(1-alpha)*weight2; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){  
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    local:
      while(1){
        ++it;
        
        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N0+beta0(j);
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+1); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j); 
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=lambda2(j)*pow(b0, 2);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+1);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=lambda2(j)*pow(b0, 2);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }
          
          rss0+=db0*(db0*N0+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update
        
        obj0=rss0/N0/2.0+objQj+objQi/2.0;
        
        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}
      }//while
      
      iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(y.dot(X.col(i))/N0);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    // isd=1: already standardized 
    if (isd == 1) {
      Beta.col(il)=beta0; 
    } else {
      Beta.col(il)=beta0.array()/sdX.array();
    }
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;
    
    if(RSQ(il) > 0.999) goto exit;
  }//for lambda
  
  exit:
    return(List::create(Named("Beta")=Beta, Named("flag")=flag, Named("rsq")=RSQ,
                              Named("RSS")=RSS, Named("lambda")=lambda, Named("nlambda")=il));
}


/*****  LM: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvEnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y,
               double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, Eigen::VectorXd weight1, Eigen::VectorXd weight2, int isd, 
               int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF){
  
  int i, j, it=0, il, iadd, ia=0; 
  double zi, obj0, obj1, rss0, rss1, b0, db0, objQi, objQj, RSS0; // lambda2, 
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda); // beta matrix for different lambdas
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXd RSSp(nlambda), RSS(nlambda), RSQ(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xbF(NF);
  double xr, dbMax, thresh2=1e-5;
  Eigen::VectorXd mX=Eigen::VectorXd::Zero(p), sdX=Eigen::VectorXd::Ones(p), di=Eigen::VectorXd::Zero(p);
  List tem;
  
  //   if (isd == 0){
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
    X.col(i)/=sdX(i);
  }
  //     
  //      tem=scaleC(X);
  //      sdX=tem$sd;
  //      X=tem$x;
  
  //       
  //   }
  //   
  //   y=y.array()-y.mean();
  
  RSS0=y.squaredNorm();
  obj0=RSS0/N/2.0; rss0=RSS0;
  
  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N);
  }
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta.array()*weight1.array();
    lambda2=lambda(il)*(1-alpha)*weight2; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha) 
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    
    //it=0;
    local:
      while(1){
        ++it;
        
        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N+beta0(j);
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+1); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j); 
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2)*lambda2(j);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+1);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2)*lambda2(j);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }
          
          rss0+=db0*(db0*N+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update
        
        obj0=rss0/N/2.0+objQj+objQi/2.0;
        
        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}
      }//while
      
      iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(y.dot(X.col(i))/N);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    // if (isd == 1) {
    // Beta.col(il)=beta0;
    // } else {
    Beta.col(il)=beta0.array()/sdX.array();
    // }
    
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;
    
    xbF.setZero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*(beta0(j)/sdX(j));}
    RSSp(il)=(yF-xbF).squaredNorm();
    
    //if(RSQ(il) > 0.999) goto exit;
  }//for lambda
  
  exit:
    return(List::create(Named("Beta")=Beta, Named("flag")=flag,
                        Named("RSS")=RSS, Named("rsq")=RSQ, Named("RSSp")=RSSp, Named("nlambda")=il,
                              Named("sdX")=sdX));
}


