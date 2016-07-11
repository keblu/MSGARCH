f.theta.mixture = function(K, NbParams, theta){
  
  if (is.vector(theta))
    theta = matrix(theta, nrow = 1)
  nTotalParams = ncol(theta)
  nTheta = nrow(theta)
  Pvector = theta[(NbParams+1):nTotalParams]
  tmp = theta[1,]
  idx = 1
  for(i in 1:(K-1)){
    tmp[(NbParams+idx):(NbParams+idx+K-1)] = Pvector[i]
    idx = idx + K
  }
  
  newTheta = matrix(nrow = nTheta,ncol = length(tmp))
  
  for(i in 1:nTheta){
    Pvector = theta[i,(NbParams+1):nTotalParams]
    idx = 1
    for(j in 1:(K-1)){
      newTheta[i,(NbParams+idx):(NbParams+idx+K-1)] = Pvector[j]
      idx = idx + K
    }
  }
  newTheta[,1:NbParams] = theta[,1:NbParams]
  return(newTheta)
}

f.theta.mixture.reverse = function(K, NbParams, theta){
  
  if (is.vector(theta))
    theta = matrix(theta, nrow = 1)
  nTotalParams = ncol(theta)
  nTheta = nrow(theta)
  newtheta = matrix(nrow = nTheta, ncol =  NbParams + K - 1)
  
  for(i in 1:nTheta){
    Pvector = theta[i,(NbParams+1):nTotalParams]
    newtheta[i,1:NbParams] = theta[i,1:NbParams]
    idx = 1
    for(j in 1:(K-1)){
      newtheta[i,(NbParams+j)] = Pvector[idx]
      idx = idx + K
    }
  }
  return(newtheta)
}