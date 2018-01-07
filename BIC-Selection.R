# This function performs BIC-selection of the tuning parameters lmbda1, lmbda2 in the integrative graphical modeling approach of McGillivray and Michailidis (2018). 

# X: List of M data sets
# M: Number of data sources
# p: Vector of length M, indicating the number of variables for each of the M data sources.
# pcore: Number of core variables
# maxiter: Maximum number of iterations for the ADMM algorithm.
# tol: Threshold tolerance for convergence. 



VerticalIntegration_BIC <- function(X,M,p,pcore,rho,lmbda1,lmbda2,maxiter,tol){
  
  n=dim(X[[1]])[1]
  N1=length(lmbda1)
  N2=length(lmbda2)
  AIC=BIC=counts=matrix(0,N1,N2)
  Scov=list()
  Rcor=list()
  thetas0=list()
  for (m in 1:M){
    Scov[[m]]=((n-1)/n)*cov(X[[m]])	
    thetas0[[m]]=matrix(0,p[m],p[m])
  }
  thetas.est=vector("list",N1*N2)
  
  for (i in 1:N1){
    for (j in 1:N2){
      lmbda_grid=c(lmbda1[i],lmbda2[j])
      if (i>1){
        thetas0=list()
        for (m in 1:M){
          thetas0[[m]]=thetas.est[[j + (i-2)*N2]][[m]]
        }  
        res.thetas=VerticalIntegration_ADMM(X, M,pcore, rho,lmbda=lmbda_grid,maxiter, tol, init=TRUE, thetas0)$thetas.hat
        Rt=list()
        for (m in 1:M){
          Rt[[m]]=res.thetas[[m]]
        }
        thetas.est[[j + (i-1)*N2]]=Rt
      } else {
        res.thetas=VerticalIntegration_ADMM(X, M,pcore, rho,lmbda=lmbda_grid,maxiter, tol, init=TRUE, thetas0)$thetas.hat
        Rt=list()
        for (m in 1:M){
          Rt[[m]]=res.thetas[[m]]
        }
        thetas.est[[j + (i-1)*N2]]=Rt
      }
      AIC[i,j]=0
      BIC[i,j]=0
      counts[i,j]=0
      for (m in 1:M){
        The=res.thetas[[m]]	
        The[abs(The)<1e-3]=0
        upper.entries=The[upper.tri(The)==TRUE]
        ct=length(upper.entries[abs(upper.entries)>0])
        counts[i,j]=counts[i,j]+ct
        if (ct==0){
          AIC[i,j]=Inf
          BIC[i,j]=Inf
        } else {
          AIC[i,j] = AIC[i,j] - n*determinant(The,logarithm=TRUE)$modulus[[1]] + n*tr(Scov[[m]]%*%The) + 2*ct
          AIC[i,j]=AIC[i,j]/M
          BIC[i,j] = BIC[i,j] - n*determinant(The,logarithm=TRUE)$modulus[[1]] + n*tr(Scov[[m]]%*%The) + log(n)*ct
          BIC[i,j]=BIC[i,j]/M
        }
      } # End of m loop
    } # End of j loop
  } # End of i loop
  
  id.AIC=which(AIC==min(AIC), arr.ind=TRUE)
  if (dim(id.AIC)[1]>1){
    id.AIC=id.AIC[1,]
  }
  opt.tune.AIC <- c(lmbda1[id.AIC[1]],lmbda2[id.AIC[2]])
  id1=id.AIC[1]
  id2=id.AIC[2]  
  thetAIC=thetas.est[[id2 + N2*(id1-1)]]
  
  id.BIC=which(BIC==min(BIC), arr.ind=TRUE)
  if (dim(id.BIC)[1]>1){
    id.BIC=id.BIC[1,]
  }
  opt.tune.BIC <- c(lmbda1[id.BIC[1]],lmbda2[id.BIC[2]])
  id1=id.BIC[1]
  id2=id.BIC[2]
  thetBIC=thetas.est[[id2 + N2*(id1-1)]]
  
  return(list(AIC=AIC,tune.AIC=opt.tune.AIC,thetas.hat.AIC=thetAIC,
              BIC=BIC,tune=opt.tune.BIC,thetas.hat.BIC=thetBIC,counts=counts))
  
} 

