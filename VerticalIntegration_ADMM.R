library(expm)
library(psych)
library(spam)


pos.part <- function(x){
  return(max(x,0))
}


soft.thres <- function(x,lmb){
  st=pmax(0,x-lmb)-pmax(0,-x-lmb)
  return(st)
}


soft.thres.mat <-
function(a,lam,penalize.diagonal){ 
	# If penalize.diagonal is FALSE, the diagonal is not penalized. 
 output <- sign(a)*pmax(0, abs(a)-lam)
  if(penalize.diagonal==FALSE){
  	diag(output) <- diag(a)	
  } 
  return(output)
}



sgl.pen <-
function(Q,V,rho,lmbda1,lmbda2,penalize.diagonal){

	lam1 = lmbda1/rho
	lam2 = lmbda2/rho

	soft.thres.Q = Q
	for(v in 1:V){
		soft.thres.Q[[v]] = soft.thres.mat(Q[[v]],lam1,penalize.diagonal=penalize.diagonal)
	}   
	
    pd=dim(Q[[1]])
	norm.soft.thres.Q=matrix(0,pd,pd)
	for(v in 1:V){
		norm.soft.thres.Q = norm.soft.thres.Q + (soft.thres.Q[[v]])^2
	}
	norm.soft.thres.Q = sqrt(norm.soft.thres.Q)

	notshrunk = 1*(norm.soft.thres.Q>lam2)
    norm.soft.thres.Q[notshrunk==0]=1

	output = Q
	for(v in 1:V){
		output[[v]] = soft.thres.Q[[v]]*(1-lam2/norm.soft.thres.Q)
		output[[v]] = output[[v]]*notshrunk
		if (penalize.diagonal==FALSE){
			diag(output[[v]])=diag(Q[[v]])
		}
	}
	return(output)
}



VerticalIntegration_ADMM <- function(X,M,pcore,rho,lmbda,maxiter,tol,init,thetas.init){
  s=1
  conv=FALSE
  n=dim(X[[1]])[1]
  eps=1e-2
  S=vector("list",M)
  R=vector("list",M)
  theta.final=vector("list",M)
  A=B=U=vector("list",M)
  p=rep(0,M)
  nsize=rep(n,M)
  
  for (m in 1:M){
  	ns=nsize[m]
    S[[m]]=((ns-1)/ns)*cov(X[[m]])
    R[[m]]=cov2cor(S[[m]])
    resid=rep(0,M)
    Xm=X[[m]]
    p[m]=dim(Xm)[2]
    A[[m]]=diag(p[m])
    B[[m]]=matrix(0,p[m],p[m])
    U[[m]]=matrix(0,p[m],p[m])
}
  
  if (init==TRUE){
    B=thetas.init
  }
  
	for(m in 1:M){
	for(j in 1:p[m]){
		X[[m]][,j] = X[[m]][,j]-mean(X[[m]][,j])
	}}

  while (s < maxiter & !conv){

    # A-minimization (likelihood part)
    for (m in 1:M){
      Ei=rho*(B[[m]] - U[[m]]) - S[[m]]
      r=eigen(Ei)
      Q=r$vectors
      lmb=r$values
      Lmb=(1/(2*rho))*(lmb + sqrt(lmb^2 + 4*rho))
      A[[m]]=Q%*%diag(Lmb)%*%t(Q)
    }
  
  
  	# B-minimization
    B.old=B   
		Q=list()
		for (m in 1:M){
			Q[[m]]=A[[m]][1:pcore,1:pcore]+U[[m]][1:pcore,1:pcore]
		}
  
  Qlist=sgl.pen(Q,M,rho,lmbda[1],lmbda[2],penalize.diagonal=FALSE)
  for (m in 1:M){
  	B[[m]][1:pcore,1:pcore]=Qlist[[m]]
  }
    
    for (m in 1:M){
      pt=p[[m]]
      At=A[[m]]+U[[m]]
      if (pcore<pt){
          id1=1:pcore
          id2=(pcore+1):pt
          p1=length(id1)
          p2=length(id2)
          B[[m]][id1,id2]=matrix(soft.thres(At[id1,id2],lmbda[1]/rho),p1,p2)
          B[[m]][id2,id1]=matrix(soft.thres(At[id2,id1],lmbda[1]/rho),p2,p1)
          B[[m]][id2,id2]=matrix(soft.thres(At[id2,id2],lmbda[1]/rho),p2,p2)
      }
      diag(B[[m]])=diag(A[[m]]+U[[m]])
    }
    

    # Dual variable updates
    for (m in 1:M){
      U[[m]]=A[[m]]+U[[m]]-B[[m]]
    } # End of m loop
    
    for (m in 1:M){
      resid[m]=norm(A[[m]]-B[[m]],"F")
    }
    
    diffB=0
    for (m in 1:M){
      diffB=diffB + sum((B[[m]]-B.old[[m]])^2)
    }
    Mx=max(c(sum(resid),sqrt(diffB)))/M
    
    if (Mx<tol){
      #print("Converged")
      for (m in 1:M){
        theta.final[[m]]=B[[m]]
      }
      conv=TRUE
      return(list(thetas.hat=theta.final,niter=s))
    } else {
      s=s+1
      conv=FALSE
    }
    
    if (s==maxiter){
      #print("Did not converge")
      for (m in 1:M){
      theta.final[[m]]=B[[m]]
      }
      conv=TRUE
      return(list(thetas.hat=theta.final,niter=s))
    }
  } # End of while loop
  
} # End of function










