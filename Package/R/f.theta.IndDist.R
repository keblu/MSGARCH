f.theta.RegIndDist = function(K, NbTotalParams, NbParamsModel, theta) {
  
  if (is.vector(theta)) 
    theta = matrix(theta, nrow = 1)
  nTheta = nrow(theta)
  NbDistParams = NbTotalParams - NbParamsModel
  if (NbDistParams[1] == 0) {
    return(theta)
  }
  NbParamsTransition = ncol(theta) - sum(NbParamsModel) - NbDistParams[1]
  newLength = sum(NbParamsModel) + sum(NbDistParams) + NbParamsTransition
  newTheta = matrix(data = NA, nrow = nTheta, ncol = newLength)
  for (i in 1:nTheta) {
    DistParams = theta[i, (sum(NbParamsModel) + 1):(sum(NbParamsModel) + NbDistParams[1])]
    PParams = theta[i, (sum(NbParamsModel) + NbDistParams[1] + 1):length(theta[i, 
      ])]
    
    ind1 = 1
    ind2 = 1
    
    for (j in 1:K) {
      newTheta[i, ind2:(ind2 + NbParamsModel[j] - 1)] = theta[i, ind1:(ind1 + 
        NbParamsModel[j] - 1)]
      newTheta[i, (ind2 + NbParamsModel[j]):(ind2 + NbParamsModel[j] + NbDistParams[1] - 
        1)] = DistParams
      ind1 = ind1 + NbParamsModel[j]
      ind2 = ind2 + NbParamsModel[j] + length(DistParams)
    }
    newTheta[i, ind2:newLength] = PParams
  }
  return(newTheta)
}

f.theta.RegIndDist.reverse = function(K, NbTotalParams, NbParamsModel, theta) {
  
  if (is.vector(theta)) 
    theta = matrix(theta, nrow = 1)
  nTheta = nrow(theta)
  NbParams = length(theta[1, ])
  NbDistParams = NbTotalParams - NbParamsModel
  if (NbDistParams[1] == 0) {
    return(theta)
  }
  newLength = sum(NbParamsModel) + NbDistParams[1] + NbParams - sum(NbTotalParams)
  
  newTheta = matrix(data = NA, nrow = nTheta, ncol = newLength)
  for (i in 1:nTheta) {
    
    PParams = theta[i, (sum(NbTotalParams) + 1):NbParams]
    DistParams = theta[i, (NbParamsModel[1] + 1):(NbParamsModel[1] + NbDistParams[1])]
    
    
    ind1 = 1
    ind2 = 1
    for (j in 1:K) {
      newTheta[i, ind1:(ind1 + NbParamsModel[j] - 1)] = theta[i, ind2:(ind2 + 
        NbParamsModel[j] - 1)]
      ind1 = ind1 + NbParamsModel[j]
      ind2 = ind2 + NbParamsModel[j] + length(DistParams)
    }
    newTheta[i, ind1:(ind1 + length(DistParams) - 1)] = DistParams
    ind1 = ind1 + length(DistParams)
    newTheta[i, ind1:newLength] = PParams
  }
  return(newTheta)
}
