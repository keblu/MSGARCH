f.error = function(message) {
  cat(paste0("try/catch ", message, "\n"))
  return(FALSE)
}


f.process.ctr = function(ctr = list()) {
  con = list(theta0 = NULL, do.init = TRUE, N.mcmc = 10000, N.burn = 5000,
             N.thin = 10, NP = 500, itermax = 500, enhance.theta0 = TRUE)
  con[names(ctr)] = ctr
  return(con)
}

f.check.y = function(y){
  if(is.null(y)){
    stop("y is NULL")
  }
  if(!is.numeric(y)){
    stop("y must be numeric")
  }
  if(all(is.nan(y))){
    stop("nan dectected in y")
  }
  if(!is.null(dim(y))){
    if(any(dim(y) == 1)){
      y = as.vector(as.matrix(y))
    } else {
      stop("y is not a vector")
    }
  }
  y = as.vector(y)
  y = as.matrix(y)
  return(y)
}

f.check.theta = function(spec, theta){
  
  
  if(is.null(theta)){
    stop("theta is NULL")
  }
  
  if(!is.numeric(theta)){
    stop("theta must be a numeric")
  }
  
  if(all(is.nan(theta))){
    stop("nan dectected in theta")
  }
  
  len.theta = length(spec$theta0)
  if (is.vector(theta)) {
    theta = matrix(theta, nrow = 1)
  }
  if(is.data.frame(theta)){
   theta =  data.matrix(theta)
  }
  
  if(dim(theta)[2] != len.theta){
    stop(paste0("Each parameter estimate in theta must be of length ",len.theta))
  }
  
  if (isTRUE(spec$is.shape.ind)) {
    theta = spec$func$f.do.shape.ind(theta)
  }
  
  if (isTRUE(spec$is.mix)) {
    theta = spec$func$f.do.mix(theta)
  }
  return(theta)
}

f.sort.theta = function(spec, theta){
  thetaUncVol = theta
  if(isTRUE(spec$is.shape.ind)){
    theta = spec$func$f.do.shape.ind(theta = theta)
  }
  
  Nbparams = spec$n.params
  Nmodel = length(Nbparams)
  if(Nmodel == 1) {
    return (theta)
  }
  name = spec$name
  unique.spec = unique(name,FALSE)
  paramsLoc = c(0,cumsum(Nbparams))
  tmp = theta
  pos = 1:Nmodel
  for (i in 1:length(unique.spec)){
    postmp = pos
    idx = name == unique.spec[i]
    posidx = pos[idx]
    Nmodelidx = length(posidx)
    idxLoc = paramsLoc[c(idx)]
    idxParams = spec$n.params[c(idx)][1]
    
    unc.vol = MSGARCH::unc.vol(object = spec, thetaUncVol)
    unc.vol.idx = unc.vol[idx]
    unc.vol.sort = sort(unc.vol.idx, index.return	= TRUE)
    for (i in 1:Nmodelidx){
      new.pos = which(unc.vol == unc.vol.sort$x[i])
      postmp[posidx[i]] = new.pos
    }
    
    for (i in 1:Nmodelidx){
      ind = unc.vol.sort$ix[i]
      tmp[(idxLoc[i]+1):(idxLoc[i]+idxParams)] = theta[(idxLoc[ind]+1):(idxLoc[ind]+idxParams)]
    }
    
  }
  
  if(!isTRUE(spec$is.mix)){
    
    p = matrix(nrow = Nmodel, ncol = Nmodel)
    for (i in 1:(Nmodel-1)){
      p[i,1:Nmodel] = tmp[(paramsLoc[Nmodel+1]+1):(paramsLoc[Nmodel+1] + Nmodel)]
    }
    p[Nmodel,] = 1-colSums(matrix(p[1:(Nmodel-1),], ncol = Nmodel ))
    tmpp = p
    for (i in 1:(Nmodel)){
      for (j in 1:(Nmodel)){
        tmpp[i,j] = p[postmp[i],postmp[j]]
      }
    }
    new.p = as.vector(tmpp[1:(Nmodel-1),])
    tmp[(paramsLoc[Nmodel+1]+1):length(tmp)] = new.p
  } else {
    p = rep(0,Nmodel)
    for (i in 1:(Nmodel-1)){
      p[i] = tmp[(paramsLoc[Nmodel+1]+1)]
    }
    p[Nmodel] = 1-sum(p)
    tmpp = p
    
    for (j in 1:(Nmodel)){
      tmpp[j] = p[postmp[j]]
    }
    new.p = tmpp[1:(Nmodel-1)]
    tmp[(paramsLoc[Nmodel+1]+1):length(tmp)] = new.p
  }
  if(isTRUE(spec$is.shape.ind)){
    tmp = spec$func$f.do.shape.ind.reverse(tmp)
  }
  return(tmp)
}

f.enhance.theta = function(spec, theta, y){
  K = spec$K
  l.y = length(y)
  sep = seq(from = 1 , to =  l.y, length.out = 11)
  vol = NULL
  if(spec$K== 1){
    vol = sqrt(var(y))
  }
  for(i in 1:(length(sep)-1)){
    vol[i] = sqrt(var(y[sep[i]:sep[i+1]]))
  }
  vol.goal = quantile(vol, prob = seq(0.1,0.9, length.out = K))
  pos = c(1,cumsum(spec$n.params)+1)
  
  if (isTRUE(spec$is.shape.ind)) {
    theta = spec$func$f.do.shape.ind(theta)
  }
  
  if (isTRUE(spec$is.mix)) {
    theta = spec$func$f.do.mix(theta)
  }
  
  for(i in 1:K){
    f.fun = function(x){
      theta.try = theta
      theta.try[,pos[i]] = x
      if (isTRUE(spec$is.shape.ind)) {
        theta.try = spec$func$f.do.shape.ind.reverse(theta.try)
      }
      
      if (isTRUE(spec$is.mix)) {
        theta.try = spec$func$f.do.mix.reverse(theta.try)
      }
      
      unc.vol = MSGARCH::unc.vol(spec, theta =theta.try)
      
      if (isTRUE(spec$is.shape.ind)) {
        theta.try = spec$func$f.do.shape.ind(theta.try)
      }
      
      if (isTRUE(spec$is.mix)) {
        theta.try = spec$func$f.do.mix(theta.try)
      }
      
      return(unc.vol[i] - vol.goal[i])
    }
    theta[,pos[i]] = uniroot(f.fun,lower = spec$lower[pos[i]], upper = spec$upper[pos[i]])$root
  }
  if(spec$K > 1){
    if(!spec$is.mix){
      pos = pos[length(pos)]
      theta[pos:length(theta)] = (1-0.8)/(K-1)
      for(i in 1:(K-1)){
        theta[pos] = 0.8
        pos = pos + K + 1
      }
    }
  }
  if (isTRUE(spec$is.shape.ind)) {
    theta = spec$func$f.do.shape.ind.reverse(theta)
  }
  
  if (isTRUE(spec$is.mix)) {
    theta = spec$func$f.do.mix.reverse(theta)
  }
  
  return(theta)
}