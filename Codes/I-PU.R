# calculate soft 
soft = function(x, tau){
  sgn = sign(x)
  x = abs(x)-tau
  x[x < 0] = 0
  return(sgn*x)	
}
entrywise_prod<-function(A,b){
  for(i in 1:length(b)){
    A[i,]=A[i,]*b[i]
  }
  return(A)
}

# Calculate MCP penalty
MCP_penalty<-function(beta,lambda,a){
  return(ifelse(beta <= a*lambda, lambda*beta-beta^2/(2*a), 0.5*a*lambda^2))
}
# Calculate the first derivative of the MCP penalty.
MCP_D<-function(beta,lambda,a){
  return(ifelse(beta<a*lambda,lambda-beta/a,0))
}
# Calculate the CMCP.
CMCP<-function(beta, lambda, a){
  pj= length(beta)
  b = 0.5*a*pj*lambda
  mu=MCP_D(sum(MCP_penalty(abs(beta), lambda,  a)), lambda,b) *MCP_D(abs(beta), lambda, a)
  return(mu)
}
# Calculate the CMCP for each group.
CMCP_G <- function(x,group,lambda, a){
  res=tapply(x, group, CMCP, lambda=lambda, a = a)
  if(is.list(res)){
    return(do.call('c',res))
  }
  else{
    return(res)
  }
}
# Update phi for each group.
up_phi<-function(xi,lambda,rho1,a){
  xi_norm=sqrt(sum(xi^2))
  if(xi_norm<=a*lambda){
    xi=soft(1,lambda/(rho1*xi_norm))/(1-1/(a*rho1))*xi
  }
  return(xi)
}
# This function is used to calculate phi, psi, and Delta.
cal_diff<-function(M,p,beta,type=1){
  # type=1 is used for calculate phi and psi, type=2 is used to calculate Delta
  index=1
  if(type==1){
    result=matrix(0,(M*(M-1)/2),p)
    for(m in 1:(M-1)){
      for(l in (m+1):M){
        result[index,]=beta[m,]-beta[l,]
        index=index+1
      }
    }
  }
  else{
    result=matrix(0,(M*(M-1)/2),M)
    for(m in 1:(M-1)){
      for(l in (m+1):M){
        result[index,m]=1
        result[index,l]=-1
        index=index+1
      }
    }
  }
  return(result)
}
# This function is used to update phi. 
up_phi_group <- function(xi,group,lambda,rho1, a){
  res = tapply(xi, group, up_phi, lambda = lambda, rho1 = rho1, a = a)
  if(is.list(res)){
    return(do.call('c', res))
  }
  else{
    return(res)
  }
}

# Calculate the expectation of y.
E_step<-function(X,Z,beta,M,pi,nl,nu){
  omega = Z 
  mu = Z
  T_m = list()
  for(m in 1:M){
    unlabel = which(Z[[m]] == 0)
    temp = X[[m]]%*%beta[m,-1]+beta[m,1]
    d_m = log((nl[m]+pi[m]*nu[m])/(pi[m]*nu[m]))
    omega[[m]][unlabel] = (exp(temp)/(1+exp(temp)))[unlabel] # calculate the omega defined in (7)
    mu[[m]] = exp(temp+d_m)/(1+exp(temp+d_m))
    T_m[[m]] = 4*(omega[[m]] - mu[[m]]) + temp
  }
  return(T_m)
}
# M-step
M_step_group<-function(beta,X,T_m,I_bd,G,A,sample_size,N,M,p,group,lambda1,lambda2, a ,rho1=1,rho2=1){
  T_m_vec = do.call('c', T_m)
  IXY = 1/N*I_bd%*%(entrywise_prod(X,T_m_vec))
  phi = cal_diff(M,p,beta,1)
  Delta = cal_diff(M,p,beta,2)
  psi = beta
  nu= 0
  upsilon = 0
  r = 1
  k = 1
  while(r>10^-3 & k<3){
    B = IXY+rho1*t(Delta)%*%(phi-nu/rho1)+rho2*psi+upsilon
    C = 1/N*entrywise_prod(X,vec(t(A%*%(rowSums(hadamard.prod(X,t(I_bd)%*%G%*%B))))))
    beta = G%*%(B-I_bd%*%C)
    varpi = beta-upsilon/rho2
    xi = cal_diff(M,p,beta,1)+nu/rho1
    mu = t(apply(beta, MARGIN=1, CMCP_G, group = group, lambda = lambda2, a = a))
    psi = soft(varpi, mu/rho2)
    phi = t(apply(xi, 1, up_phi_group, group = group, lambda = lambda1, rho1 = rho1, a = a))
    nu = nu + rho1*(Delta%*%beta-phi)
    upsilon = upsilon + rho2*(psi-beta)
    r = norm(Delta%*%beta-phi,'F') + norm(psi-beta,'F')
    k = k+1
  }
  colnames(psi)<-NULL
  return(psi)
}

# generate I
genI<-function(sample_size){
  Ilist = list()
  for(i in 1 : length(sample_size)){
    Ilist[[i]]=t(rep(1, sample_size[i]))
  }
  return(Ilist)
}
# generate G_extend
gen_G_extend<-function(G,sample_size,M,N){
  G_extend<-matrix(0,N,N)
  for(m in 1 : M){
    for(l in 1 : M){
      if(m == 1){
        row_index = 1:sample_size[1]
      }
      else{
        row_index = (sum(sample_size[1 : (m-1)])+1):sum(sample_size[1:m])
      }
      if(l == 1){
        col_index = 1 : sample_size[1]
      }
      else{
        col_index = (sum(sample_size[1 : (l-1)])+1):sum(sample_size[1 : l])
      }
      G_extend[row_index, col_index] = G[m, l]
    }
  }
  
  return(G_extend)
}


I_PU <- function(beta_hat,X_m,X,Z,I_bd,G,A,M,pi,N,p,group,lambda1,lambda2,loop,a = 3,rho1=1,rho2=1){
  p0 = p-1
  nl = rep(0,M)
  nu = rep(0,M)
  sample_size = rep(0,M)
  for(m in 1 : M){
    nl[m] = sum(Z[[m]]==1)
    nu[m] = sum(Z[[m]]==0)
    
  }
  beta_pre = matrix(0, M, p)
  iter=1
  repeat{
    T_m = E_step(X_m, Z, beta_hat, M, pi, nl, nu) # use E-step to get mu
    for(m in 1:M){
      beta_hat[m,1] = mean(T_m[[m]]-X_m[[m]]%*%beta_hat[m,-1])
      T_m[[m]] = T_m[[m]]-beta_hat[m,1]
    }
    
    beta_hat[,-1] <- M_step_group(beta_hat[,-1], X, T_m, I_bd, G, A, sample_size, N,
                                   M, p0, group, lambda1, lambda2, a = a, rho1, rho2)
    
    iter = iter+1
    res = norm(beta_hat-beta_pre,'F')
    if(res < 10^-4 | iter > loop){
      return(beta_hat)
      break
    }
    beta_pre = beta_hat
  }
}

I_PU_train<- function(beta_hat,X_m,X,Z,X_valid,Y_valid,Z_valid,sample_size,M,pi,N,p,group,lambda1_seq,lambda2_seq,loop,a = 3,rho1=1,rho2=1,beta_true = NULL){
  tune_seq = expand.grid(lambda1_seq, lambda2_seq)
  dev_mat=matrix(0, nrow(tune_seq), M)
  beta_IPU_list=list()
  X_bd = as.matrix(bdiag(X_m))# Calculate the block diagonal of X_m.
  X = do.call('rbind', X_m)# Concatenate X_m by rows.
  I_bd = as.matrix(bdiag(genI(sample_size)))
  G = (diag(1, M)+rho1*matrix(1, M, M)/rho2)/(M*rho1 + rho2)#G is defined in Remark 2.
  G_extend = gen_G_extend(G, sample_size, M, N)
  A = solve(diag(N)+1/N*hadamard.prod((X%*%t(X)), G_extend))# A is the part of the beta update that involves matrix inversion.
  # train model
  for(i in nrow(tune_seq) : 1){
    beta_IPU = I_PU(beta_hat,X_m,X,Z, I_bd,G,A,M,pi,N,p,group[-1],
                          tune_seq[i, 1], tune_seq[i, 2], rho1 = 1, rho2 = 1,loop = loop)
    beta_IPU_list[[i]] = beta_IPU
    if (! is.null(beta_true)){  
      cat(i, 'lambda1:', tune_seq[i, 1], 'lambda2:', tune_seq[i, 2], '\n')
      print(cal_criteria(beta_IPU, beta_true, M, p0, X_valid, Y_valid, Z_valid, pi))}
  }
  # select model
  for(i in nrow(tune_seq):1){
    beta_IPU=beta_IPU_list[[i]]
    for(m in 1:M){
      dev_mat[i, m] = deviances(X_valid[[m]], Z_valid[[m]], pi[m], beta_IPU[m,])
    }
  }
  ind_min=which.min(rowMeans(dev_mat))
  return(list(beta_IPU_list = beta_IPU_list,
              ind_min = ind_min,
              best_tuning = tune_seq[ind_min,])) 
  
  
}
