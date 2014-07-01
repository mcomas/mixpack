require(mvtnorm)
# This function generates a mixture with specific
mclust_rmixnorm = function(n, mclust_solution){
  func_pi = mclust_solution$parameters$pro
  func_mean = mclust_solution$parameters$mean
  func_sigma = mclust_solution$parameters$variance$sigma
  rmixnorm(n, 
           func_pi,
           func_mean,
           func_sigma)
}
rmixmod_rmixnorm = function(n, mixmod_solution){
  func_pi = mixmod_solution@bestResult@parameters@proportions
  func_mean = t(mixmod_solution@bestResult@parameters@mean)
  func_sigma = do.call('abind', list(mixmod_solution@bestResult@parameters@variance, along = 3))
  rmixnorm(n, 
           func_pi,
           func_mean,
           func_sigma)
}

rmixnorm = function(n, Pi, Mu, S, labels=F){
  z = sample(x=1:length(Pi), size=n, prob=Pi, replace=T)
  rmn = matrix(0, nrow=n, ncol=nrow(Mu))
  for(i in 1:length(Pi)){
    n_z = sum(z==i)
    if(n_z!=0){
      if(ncol(Mu)==1)
        rmn[z==i,] = rnorm(n_z, mean=Mu[,i], sd=sqrt(S[,,i]))
      else
        rmn[z==i,] = rmvnorm(n_z, mean=Mu[,i], sigma=S[,,i])
    }
  }
  if(labels)
    rmn = cbind(rmn, z)
  rmn
}

dmixnorm = function(x, Pi, Mu, S, part = 1:length(Pi)){
  #z = sample(x=1:length(Pi), size=n, prob=Pi, replace=T)
  #rmn = matrix(0, nrow=n, ncol=nrow(Mu))
  if(is.vector(x)){
    dmn = 0
  }else{
    dmn = rep(0, times=nrow(x))
  }
  for(i in part){
    if(ncol(Mu)==1)
      dmn = dmn + Pi[i] * dnorm(x, mean=Mu[,i], sd=sqrt(S[,,i]))
    else
      dmn = dmn + Pi[i] * dmvnorm(x, mean=Mu[,i], sigma=S[,,i])
  }
  dmn / sum(Pi[part])
}

get_order = function(fittedMean, fittedVariance, mainMean, mainVariance){
  K = ncol(mainMean)
  res = matrix(0, nrow=K, ncol=K)
  for(i in 1:K){
    for(j in 1:K){
      res[i,j] = norm(mainVariance[,,j]-fittedVariance[,,i], type="2") + norm(mainMean[,j]-fittedMean[,i], type = "2")
    }
  }
  return( apply(res, 1, which.min))
}

clr_mixnorm = function(X, Pi, Mu, Sigma){
  log.dnorm = llply( 1:length(Pi), function(i) log(Pi[i]) + dmvnorm(X, mean=Mu[,i], sigma=Sigma[,,i], log=T))
  log.dnorm.mean = Reduce('+', log.dnorm) / length(log.dnorm)
  data.frame(do.call(cbind, llply( log.dnorm, function(comp) comp - log.dnorm.mean)))
}