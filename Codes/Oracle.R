### Oracle: In the oracle method, it is assumed that true responses and cluster structures are known, and then a logistic model is constructed

# W is define in 5.1
cal_W = function(p0, beta_true){
  s=rep(0, p0)
  for(i in 1 : p0){
    s[i] = length(unique(beta_true[, i]))
  }
  S = sum(s)
  cums = c(0, cumsum(s))
  W = list()
  for(m in 1 : M){
    temp = matrix(0, p0, S)
    for(i in 1 : p0){
      unique_group = unique(beta_true[, i])
      subgroup = beta_true[, i]
      index = 1
      pos_list = list()
      for (j in 1 : length(unique_group)){
        pos_list[[j]] = which(subgroup == unique_group[j])
      }
      for(j in 1 : length(unique_group)){
        subgroup[pos_list[[j]]] = j
      }
      temp[i,cums[i]+subgroup[m]] = 1
    }
    W[[m]] = temp
  }
  W = (do.call('rbind', W))
  return(W)
}


oracle_model = function(important_coef, Y, X_m, Z, beta_true, p0, M){
  X_bd = as.matrix(bdiag(X_m))
  W = cal_W(p0, beta_true)
  X_or = X_bd %*% W
  X_use = X_or[, 1 : (2*important_coef)]
  Y_vec = do.call('c', Y)
  oracle = glmnet(X_use, Y_vec, family = 'binomial', lambda=0, intercept = F)
  oracle = as.vector(coef(oracle)[2 : (2*important_coef+1)])
  oracle = c(oracle, rep(0, p0-important_coef))
  oracle = W%*%oracle
  oracle = matrix(oracle, M, p0, byrow = T)
  oracle = cbind(rep(0, M), oracle)
  return(oracle)
}
