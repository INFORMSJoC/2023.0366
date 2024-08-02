options(warn=-1)
library(glmnet)
library(ncvreg)
library(PUlasso)
library(grpreg)
# This code includes the following single models, and their estimation results can be selected as initial values for warm start.
# lr_lasso: Logistic regression model with Lasso regularization. We use the "glmnet" package for training.
# lr_group: Logistic regression model with group penalty. We use the "grpreg" package for training.
# pu_lasso: PU model with Lasso regularization. We use the "PUlasso" package for training.
# pu_grplasso: PU model with group Lasso. We use the "PUlasso" package for training.
# pu_mcp: PU model with mcp penalty
# pu_group: PU model with group penalty.


pu_mcp<-function(X, Z, pi, init_value, lambda){
  # init_value: A vector representing an initial point where we start PU_mcp algorithm from.
  # lambda: Tuning parameter.
  nl = sum(Z == 1)
  nu = sum(Z == 0)
  y_hat = Z
  beta_pre = init_value
  beta_hat = init_value
  iter = 1
  d = log((nl+pi*nu)/(pi*nu))
  unlabel = which(Z == 0)
  repeat{
    iter = iter + 1
    temp = X %*% beta_hat[-1]+beta_hat[1]
    y_hat[unlabel] = (exp(temp)/(1+exp(temp)))[unlabel]
    mu_star = exp(temp+d)/(1+exp(temp+d))
    mu = 4*(y_hat-mu_star)+temp
    fit = ncvreg(X, mu, family = 'gaussian', penalty = 'MCP',lambda = lambda)
    beta_hat = coef(fit)
    if (sqrt(sum((beta_hat-beta_pre)^2)) < 10^-3 |iter > 3){
      return (beta_hat)
      break
    }
    beta_pre = beta_hat
  }
  
}

pu_group<-function(X, Z, pi, init_value, lambda, group){
  nl = sum(Z == 1)
  nu = sum(Z == 0)
  y_hat = Z
  beta_pre = init_value
  beta_hat = init_value
  iter = 1
  d = log((nl+pi*nu)/(pi*nu))
  unlabel = which(Z == 0)
  repeat{
    iter = iter + 1
    temp = X %*% beta_hat[-1]+beta_hat[1]
    y_hat[unlabel] = (exp(temp)/(1+exp(temp)))[unlabel]
    mu_star = exp(temp+d)/(1+exp(temp+d))
    mu = 4*(y_hat-mu_star)+temp
    fit = grpreg(X=X, y = mu, family = 'gaussian', penalty = 'grMCP',lambda = lambda, group = group)
    beta_hat = coef(fit)
    if (sqrt(sum((beta_hat-beta_pre)^2))<10^-3 |iter > 10){
      return (beta_hat)
      break
    }
    beta_pre = beta_hat
  }
  
}

# single_model: This function applies each of the six specified models to individual datasets separately
single_model = function(M, Y, X_m, Z, p, pi, group, Y_test, Z_test, X_test){
  # lr_lasso
  beta_lr_lasso = matrix(0, M, p) 
  for(m in 1 : M){
    CV_lr_lasso = cv.glmnet(X_m[[m]], Z[[m]], family='binomial', 
                         nfolds = 5, lambda = seq(0, 0.01, 0.001))
    beta_lr_lasso[m, ] = as.vector(coef(CV_lr_lasso))
  }
  # lr_group
  beta_lr_group = matrix(0, M, p)
  for(m in 1 : M){
    CV_lr_group = cv.grpreg(X = X_m[[m]], y = Z[[m]], family='binomial', 
                          nfolds = 5, penalty = 'grMCP', group = group[-1])
    beta_lr_group[m, ] = as.vector(coef(CV_lr_group))
  }
  # pu_lasso
  beta_pu_lasso = matrix(0, M, p)
  for(m in 1:M){
    CV_pu_lasso = cv.grpPUlasso(X = X_m[[m]], z=Z[[m]], py1=pi[[m]], 
                          lambda = seq(0.01, 0.025, 0.001), initial_coef = beta_lr_lasso[m, ])
    lambda_index = which(CV_pu_lasso$lambda == CV_pu_lasso$lambda.min)
    beta_pu_lasso[m, ] = coef(CV_pu_lasso)[ ,lambda_index]
  }
  
  # pu_mcp
  beta_pu_mcp = matrix(0, M, p)
  lambda_mcp = seq(0.01,0.12,0.002)
  dev_list_mcp = matrix(0, length(lambda_mcp), M)
  for(m in 1 : M){
    beta_mcp_list = list()
    for (i in 1 : length(lambda_mcp)){
      fit_pu_mcp = pu_mcp(X = X_m[[m]], Z = Z[[m]], pi = pi[[m]], init_value = beta_pu_lasso[m,], lambda = lambda_mcp[i])
      beta_mcp_list[[i]] = fit_pu_mcp
      dev_list_mcp[i, m]=deviances(X_test[[m]], Z_test[[m]], pi[m], fit_pu_mcp)
    }
    ind_mcp = which.min(dev_list_mcp[ , m])
    beta_pu_mcp[m, ] = beta_mcp_list[[ind_mcp]]
  }
  
  # pu_grplasso
  beta_pu_grplasso = matrix(0, M, p)
  for(m in 1 : M){
    CV_pu_grplasso = cv.grpPUlasso(X = X_m[[m]], z = Z[[m]], py1 = pi[[m]], 
                                      lambda = seq(0, 0.025, 0.001),
                                      initial_coef = beta_pu_mcp[m, ], group = group[-1])
    lambda_index = which(CV_pu_grplasso$lambda == CV_pu_grplasso$lambda.min)
    beta_pu_grplasso[m, ] = coef(CV_pu_grplasso)[ ,lambda_index]
  }
  # pu_group
  lambda_pu_group = seq(0.01, 0.15, 0.002)
  dev_list_pu_group = matrix(0, length(lambda_pu_group), M)
  beta_pu_group = matrix(0, M, p)
  for(m in 1 : M){
    beta_pu_group_list = list()
    for (i in 1 : length(lambda_pu_group)){
      tryCatch({fit_pu_group <- pu_group(X = X_m[[m]], Z = Z[[m]], pi = pi[[m]],
                                           init_value = beta_pu_grplasso[m, ], lambda = lambda_pu_group[i], group = group[-1])}
               , error = function(e){fit_pu_group = rep(0, p)})
      beta_pu_group_list[[i]] = fit_pu_group
      dev_list_pu_group[i, m] = deviances(X_test[[m]], Z_test[[m]], pi[m], fit_pu_group)
    }
    ind_pu_group = which.min(dev_list_pu_group[, m])
    beta_pu_group[m, ] = beta_pu_group_list[[ind_pu_group]]
  }
  return(list(
    beta_lr_lasso = beta_lr_lasso, 
    beta_lr_group = beta_lr_group,
    beta_pu_lasso = beta_pu_lasso,
    beta_pu_mcp = beta_pu_mcp,
    beta_pu_grplasso = beta_pu_grplasso,
    beta_pu_group = beta_pu_group))
}


