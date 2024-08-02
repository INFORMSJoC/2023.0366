# I-LR: this method uses the traditional logistic regression with the penalty described in Equation (6).
source('I-PU.R')
R_step<-function(X,Y,beta,M){
  mu_star=list()
  mu_m=list()
  for(m in 1:M){
    temp=X[[m]]%*%beta[m,-1]+beta[m,1]
    mu_star[[m]]=exp(temp)/(1+exp(temp))
    mu_m[[m]]=4*(Y[[m]]-mu_star[[m]])+temp
  }
  return(mu_m)
}

I_LR<-function(beta_hat, X_m, X, Y, I_bd, G, A, 
                      M, N, p, group, lambda1, lambda2, loop, a = 3, rho1 = 1, rho2 = 1){
  p0=p-1
  beta_pre=matrix(0, M, p)
  iter = 1
  mu_m = R_step(X_m, Y, beta_hat, M)
  repeat{
    
    mu_m=R_step(X_m,Y,beta_hat,M)
    for(m in 1:M){
      beta_hat[m,1] = mean(mu_m[[m]]-X_m[[m]]%*%beta_hat[m,-1])
      mu_m[[m]] = mu_m[[m]]-beta_hat[m,1]
    }
    beta_hat[,-1]<-M_step_group(beta_hat[,-1],X,mu_m,I_bd,G,A,sample_size,N,M,p0,group,lambda1,lambda2,a,rho1,rho2)
    
    iter = iter+1
    res = norm(beta_hat-beta_pre,'F')
    #if(iter%%50==0){
    #print(res)
    #}
    
    if(res<10^-4 | iter>loop){
      return(beta_hat)
      break
    }
    beta_pre=beta_hat
  }
}


I_LR_train<- function(beta_hat,X_m,X,Y,Z,X_valid,Y_valid,Z_valid,sample_size,M,pi,N,p,group,lambda1_seq,lambda2_seq,loop,a = 3,rho1=1,rho2=1,beta_true = NULL){
  tune_seq = expand.grid(lambda1_seq, lambda2_seq)
  neg_log_mat=matrix(0, nrow(tune_seq), M)
  beta_ILR_list=list()
  X_bd = as.matrix(bdiag(X_m))# calculate the bdiag of X_m
  X = do.call('rbind', X_m)# Patch X_m by rows
  I_bd = as.matrix(bdiag(genI(sample_size)))
  G = (diag(1, M)+rho1*matrix(1, M, M)/rho2)/(M*rho1 + rho2)# G is define in Remark2
  G_extend = gen_G_extend(G, sample_size, M, N)
  A = solve(diag(N)+1/N*hadamard.prod((X%*%t(X)), G_extend))#### A is the part of the beta update that involves solving for the inverse, and is presented separately to avoid duplicating calculations
  # train model
  for(i in nrow(tune_seq) : 1){
    
    beta_ILR = I_LR(beta_hat,X_m,X,Z, I_bd,G,A,M,N,p,group[-1],
                    tune_seq[i, 1], tune_seq[i, 2], rho1 = 1, rho2 = 1,loop = loop)
    beta_ILR_list[[i]] = beta_ILR
    if (! is.null(beta_true)){  
    cat(i, 'lambda1:', tune_seq[i, 1], 'lambda2:', tune_seq[i, 2], '\n')
    print(cal_criteria(beta_ILR, beta_true, M, p0, X_valid, Y_valid, Z_valid, pi))}
  }
  # select model
  for(i in nrow(tune_seq):1){
    beta_ILR=beta_ILR_list[[i]]
    for(m in 1:M){
      neg_log_mat[i, m] = neg_log_likelihood(X_valid[[m]], Y_valid[[m]], beta_ILR[m,])
    }
  }
  ind_min=which.min(rowMeans(neg_log_mat))
  return(list(beta_ILR_list = beta_ILR_list,
              ind_min = ind_min,
              best_tuning = tune_seq[ind_min,])) 
}

neg_log_likelihood <- function(x,y,beta) {  
  eta <- x%*%beta[-1]+beta[1]
  neg_log_lik <- -sum(-log(1+exp(eta))+y*eta)  
  return(neg_log_lik)  
}  

