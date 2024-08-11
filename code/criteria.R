library(pROC)
# Calculate AUC
cal_auc <- function(y, y_hat)
{
  as.double(pROC::auc(roc(y, y_hat, quiet = T)))
}

# Calculate criteria for predictive performance
predictive_criteria<-function(beta,X_test,Y_test,Z_test,pi,M){
  f1=0
  accuracy=0
  auc=0
  auc_adjust=0
  y_hat=list()
  y_prob=list()
  for(m in 1:M){
    y_hat[[m]]=ifelse((X_test[[m]]%*%beta[m,-1]+beta[m,1])>0,1,0)
    y_prob[[m]]=exp(X_test[[m]]%*%beta[m,-1]+beta[m,1])/(1+exp(X_test[[m]]%*%beta[m,-1]+beta[m,1]))
    TP=sum(y_hat[[m]]!=0 & Y_test[[m]]!=0)
    TPFP=sum(y_hat[[m]]!=0)
    TPFN=sum(Y_test[[m]]!=0)
    precision=TP/TPFP
    recall=TP/TPFN
    if(TPFP!=0 & TP!=0){
      f1=f1+2*precision*recall/(precision+recall)/M
    }
    
    accuracy=accuracy+sum(y_hat[[m]]==Y_test[[m]])/length(Y_test[[m]])/M
    auc=auc+cal_auc(as.vector(Y_test[[m]]),as.vector(y_prob[[m]]))/M
    auc_adjust=auc_adjust+(cal_auc(as.vector(Z_test[[m]]),as.vector(y_prob[[m]]))-pi[[m]]/2)/(1-pi[[m]])/M
  }
  return(c(f1=f1,accuracy=accuracy,auc=auc,auc_adjust=auc_adjust))
}



# Calculate all criteria
cal_criteria<-function(beta_hat,beta_true,M,p0,X_test,Y_test,Z_test,pi){
  predictive_res=predictive_criteria(beta_hat,X_test,Y_test,Z_test,pi,M)
  if(ncol(beta_hat)==(p0+1)){
    beta_hat=beta_hat[,2:(p0+1)]
  }
  CP = 0
  TPR = 0
  FPR = 0
  FNR = 0
  for(m in 1:M){
    TP = sum(I(beta_true[m,]!=0 & beta_hat[m,]!=0))
    FP = sum(I(beta_true[m,]==0 & beta_hat[m,]!=0))
    FN = sum(I(beta_true[m,]!=0 & beta_hat[m,]==0))
    TPFN = sum(I(beta_true[m,]!=0))
    TNFP = sum(I(beta_true[m,]==0))
    TPR = TP/TPFN+TPR
    FPR = FP/TNFP+FPR
  }
  RMSE=sqrt(norm(beta_hat-beta_true,'f')^2/(M*p0))
  return(c(TPR=TPR/M,FPR=FPR/M,RMSE=RMSE,predictive_res))
}
