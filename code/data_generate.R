library(MASS)
library(Matrix)
library(matrixcalc)
data_generate = function(seed, example, p0){
  # p0: the dimensions without intercept
  # M: the number of datasets
  # sample_size: sample size of each dataset
  # N : total sample size
  # g_size : the size of each group
  # group : group vector without intercept
  # group_all: group vector
  # beta_true: the true coefficients
  # Sigma : covariance of x

  p = p0 + 1
  if (example == 1){
    important_pos = 1:16
    M  = 10
    n  = 200
    nu = rep(n/2, M)
    nl = rep(n/2, M)
    pi = rep(0, M)
    sample_size = rep(n, M)
    N = sum(sample_size)
    g_size = 5
    group = rep(2 : (p0/g_size+1), each = g_size)
    beta_true = matrix(0, M, p0)
    for(m in 1:(M/2)){
      beta_true[m, ] = rep(c(1, -1, 0), c(5, 3, p0-8))
    }
    for(m in (M/2+1) : M){
      beta_true[m, ] = rep(c(-1, 1, 0), c(5, 3, p0-8))
    }
  }
  
  else if (example == 2){
    M=10
    n=200
    important_pos = c(1,2,3,4,5,6,8,9,11,12,14,15,18,20)
    nu=rep(n/2,M)
    nl=rep(n/2,M)
    pi=rep(0,M)
    sample_size=rep(n,M)
    N=sum(sample_size)
    g_size=5
    group           = rep(4:(p0/g_size+2), each=g_size)
    group=c(2,2,3,3,3,group)
    group_all=c(1, group)
    
    beta_true = matrix(0,M,p0)
    beta_true[1, ]=c(rep(c(1, 0.5, -0.5, 0, -0.5),c(2, 1, 1, 1, 3)), rep(0, p0-8))
    beta_true[2, ]=c(rep(c(1, 0.5, -0.5, 0, -0.5),c(2, 1, 1, 1, 3)), rep(0, p0-8))
    beta_true[3, ]=c(rep(c(1, 0.5, -0.5, 0, -0.5),c(2, 1, 1, 1, 3)), rep(0, p0-8))
    beta_true[4, ]=c(rep(c(1, 0.5, -0.5, 0, 0.5),c(2, 1, 1, 1, 3)), rep(0, p0-8))
    beta_true[5, ]=c(rep(c(1, 0.5, -0.5, 0, 0.5),c(2, 1, 1, 1, 3)), rep(0, p0-8))
    beta_true[6, ]=c(rep(c(-1, 0.5, -0.5, 0, 0.5),c(2, 1, 1, 1, 3)), rep(0, p0-8))
    beta_true[7, ]=c(rep(c(-1, 0.5, -0.5, 0, 0.75),c(2, 1, 1, 4, 2)), rep(0, p0-10))
    beta_true[8, ]=c(rep(c(-1, 0.5, -0.5, 0, 0.75),c(2, 1, 1, 4, 2)), rep(0, p0-10))
    beta_true[9, ]=c(rep(c(-1, 0.5, -0.5, 0, 0.75),c(2, 1, 1, 4, 2)), rep(0, p0-10))
    beta_true[10, ]=c(rep(c(-1, 0.5, -0.5, 0, 0.75),c(2, 1, 1, 4, 2)), rep(0, p0-10))
  }
  else if (example == 3){
    important_pos = 1:16
    M = 20
    n = 100
    nu = rep(n/2, M)
    nl = rep(n/2, M)
    pi = rep(0, M)
    sample_size = rep(n, M)
    N = sum(sample_size)
    g_size = 5
    group = rep(2 : (p0/g_size+1), each = g_size)
    beta_true = matrix(0, M, p0)
    
    for(m in 1:(M/2)){
      beta_true[m, ] = rep(c(1, -1, 0), c(5, 3, p0-8))
    }
    for(m in (M/2+1) : M){
      beta_true[m, ] = rep(c(-1, 1, 0), c(5, 3, p0-8))
    }
  }
  
  
  group_all = c(1, group)
  rho         = 0.3
  Sigma       = diag(1, p0)
  for(i in 1 : p0){
    for(j in i : p0){
      Sigma[i, j]  = rho^abs(i-j)*I(group[i] == group[j])
      Sigma[j, i]  = Sigma[i, j]
    }
  }
  rho1 = 1
  rho2 = 1
  X_m  = list()
  Y = list()
  X_test = list()
  Y_test = list()
  Z = list()
  Z_test = list()
  X_valid = list()
  Y_valid = list()
  Z_valid = list()
  set.seed(1000 * seed)
  for(m in 1:M){
    X_all = mvrnorm(n = 50000, rep(0, p0), Sigma)
    beta_all = X_all[, 1:20] %*% beta_true[m,1:20] 
    Y_all = rbinom(50000, 1, exp(beta_all)/(1+exp(beta_all)))
    pi[m] = sum(Y_all)/50000
    total = c(1 : 50000)
    ind_l = sample(which(Y_all == 1), nl[m]+200) 
    ind_u = sample(total[-ind_l], nu[m]+2000)
    ind_train = sort(c(ind_l[1 : nl[m]], ind_u[1 : nu[m]]))
    ind_valid = sort(c(ind_l[(nl[m] + 1) : (nl[m]+100)], ind_u[(nu[m]+1) : (nu[m]+1000)]))
    ind_test = sort(c(ind_l[(nl[m] + 101) : (nl[m]+200)], ind_u[(nu[m]+1001) : (nu[m]+2000)]))
    Y_test[[m]] = Y_all[ind_test]
    X_test[[m]] = X_all[ind_test, ]
    Y_valid[[m]] = Y_all[ind_valid]
    X_valid[[m]] = X_all[ind_valid, ]
    Y[[m]] = Y_all[ind_train]
    X_m[[m]] = X_all[ind_train, ] 
    Y_all[ind_u] = 0
    Z[[m]] = Y_all[ind_train]
    Z_test[[m]] = Y_all[ind_test]
    Z_valid[[m]] = Y_all[ind_valid]
  }

  return (list(M  = M,
               n  = n,
               nu = nu,
               nl = nl,
               pi = pi,
               sample_size = sample_size,
               N = N,
               p = p,
               beta_true = beta_true,
               X_m  = X_m,
               Y = Y,
               X_test = X_test,
               Y_test = Y_test,
               X_valid = X_valid,
               Y_valid = Y_valid,
               Z = Z,
               Z_test = Z_test,
               Z_valid = Z_valid,
               group_all = group_all,
               important_pos = important_pos
  ))
}



