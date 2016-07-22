
f.sort.theta = function(theta, spec) {
  thetaUncVol = theta
  if (isTRUE(spec$DistRegInd)) {
    theta = MSGARCH::f.theta.RegIndDist(spec$K, spec$NbParams, spec$NbParamsModel, 
      theta)
  }
  
  Nbparams = spec$NbParams
  Nmodel = length(Nbparams)
  if (Nmodel == 1) {
    return(theta)
  }
  name = spec$name
  unique.spec = unique(name, FALSE)
  Nbparams = spec$NbParams
  paramsLoc = c(0, cumsum(Nbparams))
  tmp = theta
  pos = 1:Nmodel
  for (i in 1:length(unique.spec)) {
    postmp = pos
    idx = name == unique.spec[i]
    posidx = pos[idx]
    Nmodelidx = length(posidx)
    idxLoc = paramsLoc[c(idx)]
    idxParams = spec$NbParams[c(idx)][1]
    
    unc.vol = spec$f.unc.vol(thetaUncVol, y = 0)
    unc.vol.idx = unc.vol[idx]
    unc.vol.sort = sort(unc.vol.idx, index.return = TRUE)
    for (i in 1:Nmodelidx) {
      new.pos = which(unc.vol == unc.vol.sort$x[i])
      postmp[posidx[i]] = new.pos
    }
    
    for (i in 1:Nmodelidx) {
      ind = unc.vol.sort$ix[i]
      tmp[(idxLoc[i] + 1):(idxLoc[i] + idxParams)] = theta[(idxLoc[ind] + 1):(idxLoc[ind] + 
        idxParams)]
    }
    
  }
  
  if (!isTRUE(spec$mixture)) {
    
    p = matrix(nrow = Nmodel, ncol = Nmodel)
    for (i in 1:(Nmodel - 1)) {
      p[i, 1:Nmodel] = tmp[(paramsLoc[Nmodel + 1] + 1):(paramsLoc[Nmodel + 
        1] + Nmodel)]
    }
    p[Nmodel, ] = 1 - colSums(matrix(p[1:(Nmodel - 1), ], ncol = Nmodel))
    tmpp = p
    for (i in 1:(Nmodel)) {
      for (j in 1:(Nmodel)) {
        tmpp[i, j] = p[postmp[i], postmp[j]]
      }
    }
    new.p = as.vector(tmpp[1:(Nmodel - 1), ])
    tmp[(paramsLoc[Nmodel + 1] + 1):length(tmp)] = new.p
  } else {
    p = rep(0, Nmodel)
    for (i in 1:(Nmodel - 1)) {
      p[i] = tmp[(paramsLoc[Nmodel + 1] + 1)]
    }
    p[Nmodel] = 1 - sum(p)
    tmpp = p
    
    for (j in 1:(Nmodel)) {
      tmpp[j] = p[postmp[j]]
    }
    new.p = tmpp[1:(Nmodel - 1)]
    tmp[(paramsLoc[Nmodel + 1] + 1):length(tmp)] = new.p
  }
  if (isTRUE(spec$DistRegInd)) {
    tmp = MSGARCH::f.theta.RegIndDist.reverse(spec$K, spec$NbParams, spec$NbParamsModel, 
      tmp)
  }
  return(tmp)
}

f.error = function(message) {
  cat(paste0("try/catch ", message, "\n"))
  return(FALSE)
}


f.process.ctr = function(ctr = list()) {
  con = list(theta0 = NULL, do.init = FALSE, N.mcmc = 5000, N.burn = 1000, N.thin = 10)
  con[names(ctr)] = ctr
  return(con)
}
