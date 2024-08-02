library(RMTL)
MTL_model = function(X,Z,lambda1,lambda2,M,p,type,k=2,opts=list(init = 0, tol = 10^-3, maxIter = 1000)){
  cvfit<-cvMTL(X, Z, type="Classification", Regularization=type,
               Lam2=lambda2, 
               Lam1_seq=lambda1)
  lambda_min = cvfit$Lam1.min
  fit<-MTL(X, Z, type="Classification", Regularization=type,
           Lam2=lambda2, 
           Lam1_seq=lambda_min,k=2,opts=list(init = 0, tol = 10^-3, maxIter = 1000))
  beta = matrix(0,M,p)
  beta[,1]=fit$C
  beta[,-1]=t(fit$W)
  return(beta)
}