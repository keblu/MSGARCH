f.spec = function(models, do.mix = FALSE, do.shape.ind = FALSE) {
  
  ################################################### create the relevant model C++ objects
  K = length(models)  # number of models
  
  if (K == 1) {
    do.mix = FALSE
    do.shape.ind = FALSE
  }
  options(warn = -1)
  if (K > 1) {
    tmp = list()
    for (i in 1:K) {
      tmp[i] = new(models[[i]])
    }
  } else {
    tmp = new(models[[1]])
  }
  options(warn = 0)
  if (K > 1) {
    mod = new(MSgarch, tmp)
  } else {
    mod = tmp
  }
  ################################################### detach for speed
  dist = NULL
  name = mod$name
  for (i in 1:length(name)) {
    dist[i] = stringr::str_sub(name[i], start = stringr::str_locate(name, "_")[i, 
      1] + 1, nchar(name[i]))
  }
  uniqueDist = unique(dist)
  if (isTRUE(do.shape.ind) && length(uniqueDist) > 1) {
    stop("The distribution of each regime must be the same if the distribution are not regime dependent")
  }
  calc_ht = mod$calc_ht
  eval_model = mod$eval_model
  ineq_func.base = mod$ineq_func
  f_sim = mod$f_sim
  n.params = mod$NbParams
  n.params.vol = mod$NbParamsModel
  f_pdf_Rcpp = mod$f_pdf
  f_cdf_Rcpp = mod$f_cdf
  f_rnd_Rcpp = mod$f_rnd
  f_unc_vol_Rcpp = mod$f_unc_vol
  if (K > 1) {
    f_get_Pstate_Rcpp = mod$f_get_Pstate
  } else {
    f_get_Pstate_Rcpp = function(theta, y, PLast) {
      if (!isTRUE(PLast)) {
        out = matrix(1, nrow = length(y), ncol = 1)
      } else {
        out = matrix(1, nrow = 1, ncol = 1)
      }
      return(out)
    }
  }
  NbtotalParams = sum(n.params)
  
  f.do.mix = function(theta) {
    return(f.theta.mixture(K, NbtotalParams, theta))
  }
  
  f.do.mix.reverse = function(theta) {
    return(f.theta.mixture.reverse(K, NbtotalParams, theta))
  }
  
  f.do.shape.ind = function(theta) {
    return(f.theta.RegIndDist(K, n.params, n.params.vol, theta))
  }
  
  f.do.shape.ind.reverse = function(theta) {
    return(f.theta.RegIndDist.reverse(K, n.params, n.params.vol, theta))
  }
  
  # ======================# model simulation #======================#
  
  f.sim.base = function(n, theta, K_ = K, burnin = 500, do.state = FALSE, do.mix = FALSE, 
    do.shape.ind = FALSE) {
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
      draws = matrix(data = NA, nrow = nrow(theta), ncol = n)
      state = matrix(data = NA, nrow = nrow(theta), ncol = n)
    for(i in 1:nrow(theta)){
      tmp = f_sim(n, theta[i,], burnin)
      if(K_ == 1 ){
        draws[i,] = tmp
        state[i,] = rep(1,n)
      } else{
        draws[i,] = tmp$draws
        state[i,] = tmp$state
      }
    }
    out = list()
    out$draws = draws
    if (isTRUE(do.state)) {
      out$state = state
    }
   return(out)
  }
  
  # =========================================================#
  f_Pstate.base = function(theta, y, K_ = K, do.mix = FALSE, do.shape.ind = FALSE) {
    
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    
    out = array(dim = c(nrow(y), nrow(theta), K_ ))
    for(i in 1:nrow(theta)){
      tmp = f_get_Pstate_Rcpp(theta[i,], y, FALSE)
      for(j in 1:K_){
        out[,i,j] = tmp[,j]
      }
    }
    
    
   return(out)
  }
  
  
  f_Plast.base = function(theta, y, K_ = K, do.mix = FALSE, do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    
    out = matrix(data = NA,nrow = nrow(theta), ncol = K_)
    for(i in 1:nrow(theta)){
      out[i,] = f_get_Pstate_Rcpp(theta[i,], y, TRUE)
    }
    
    return(out)
  }
  
  # ======================# model simulation #======================#
  f_unc_vol.base = function(theta, y = 0, do.mix = FALSE, do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }

    for(i in 1:nrow(theta)){
      out = f_unc_vol_Rcpp(theta, y)
    }
  
    
    out = sqrt(out)
    return(out)
  }
  # =========================================================#
  
  # ======================# variance #======================#
  f_ht.base = function(theta, y, do.mix = FALSE, do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    
    out = calc_ht(theta, y)
    return(out)
  }
  # =========================================================#
  
  # ======================# kernel #======================#
  f_kernel.base = function(theta, y, log = TRUE, do.mix = FALSE, do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    lnd = eval_model(theta, y)
    lnd[is.na(lnd) | is.nan(lnd) | is.infinite(lnd)] = -1e+10
    if (!log) 
      lnd = exp(lnd)
    return(lnd)
  }
  # =========================================================
  
  # ===================# predictive density #===================#
  f_pred.base = function(x, theta, y, log = TRUE, do.mix = FALSE, do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    
    out = f.pred.helper(x, theta, y, spec = f.spec(models), log = log)
    return(out)
  }
  # ==============================================================#
  
  # ===========================# PIT #===========================#
  f_pit.base = function(x, theta, y, do.norm = FALSE, do.mix = FALSE, do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    
    out = f.pit(x, theta, y, spec = f.spec(models), do.norm = do.norm)
    return(out)
  }
  # ==============================================================#
  
  # =======================# VaR forecast #========================#
  f_risk.base = function(theta, y, level = 0.95, ES = TRUE, do.mix = FALSE, 
    do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    
    out = list()
    out = f.risk.helper(theta = theta, y = y, spec = f.spec(models), level = level,  ES = ES)
    return(out)
  }
  # ==============================================================#
  
  # =======================# RND #========================#
  f_rnd.base = function(n, theta, K_ = K, y = vector("double", 0), do.state = FALSE, 
    do.mix = FALSE, do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    draws = matrix(data = NA, nrow = nrow(theta), ncol = n)
    state = matrix(data = NA, nrow = nrow(theta), ncol = n)
    for(i in 1:nrow(theta)){
      tmp = f_rnd_Rcpp(n, theta[i,], y)
      if(K_ == 1 ){
        draws[i,] = tmp
        state[i,] = rep(1,n)
      } else{
        draws[i,] = tmp$draws
        state[i,] = tmp$state
      }
      
    }
    out = list()
    out$draws = draws
    if (isTRUE(do.state)) {
      out$state = state
    }
    return(out)
  }
  # ==============================================================#
  
  # =======================# PDF #========================#
  f_pdf.base = function(x, theta, y = vector("double", 0), log = TRUE, do.mix = FALSE, 
    do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    out = matrix(data = NA,nrow = nrow(theta),ncol = length(x))
    for(i in 1:nrow(theta)){
    out[i,] = f_pdf_Rcpp(x, theta[i,], y, log)
    }
    return(out)
  }
  # ==============================================================#
  
  # =======================# CDF #========================#
  f_cdf.base = function(x, theta, y = vector("double", 0), log = TRUE, do.mix = FALSE, 
    do.shape.ind = FALSE) {
    y = as.matrix(y)
    if (isTRUE(do.shape.ind)) {
      theta = f.do.shape.ind(theta)
    }
    
    if (isTRUE(do.mix)) {
      theta = f.do.mix(theta)
    }
    
    if (is.vector(theta)) {
      theta = matrix(theta, nrow = 1)
    }
    out = matrix(data = NA,nrow = nrow(theta),ncol = length(x))
    for(i in 1:nrow(theta)){
      out[i,] = f_cdf_Rcpp(x, theta[i,], y, log)
    }
    
    return(out)
  }
  # ==============================================================#
  
  if (isTRUE(do.mix) && !isTRUE(do.shape.ind)) {
    f.sim = function(n, theta, burnin = 500, do.state = FALSE) {
      return(f.sim.base(n, theta, burnin = burnin, do.state = do.state, 
        do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    f_Pstate = function(theta, y) {
      return(f_Pstate.base(theta, y, do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    f_Plast = function(theta, y) {
      return(f_Plast.base(theta, y, do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    f_unc_vol = function(theta, y = 0) {
      return(f_unc_vol.base(theta, y = y, do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    f_ht = function(theta, y) {
      return(f_ht.base(theta, y, do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    f_kernel = function(theta, y, log = TRUE) {
      return(f_kernel.base(theta, y, log = log, do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    f_pred = function(x, theta, y, log = TRUE) {
      return(f_pred.base(x, theta, y, log = log, do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    f_pit = function(x, theta, y, do.norm = FALSE) {
      return(f_pit.base(x, theta, y, do.norm = do.norm, do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    f_risk = function(theta, y, level = 0.95,  ES = TRUE) {
      return(f_risk.base(theta, y, level = level,  ES = ES, do.mix = TRUE, 
        do.shape.ind = FALSE))
    }
    
    f_rnd = function(n, theta, y = vector("double", 0), do.state = FALSE) {
      return(f_rnd.base(n, theta, y = y, do.state = do.state, do.mix = TRUE, 
        do.shape.ind = FALSE))
    }
    
    f_pdf = function(x, theta, y = vector("double", 0), log = TRUE) {
      return(f_pdf.base(x, theta, y = y, log = log, do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    f_cdf = function(x, theta, y = vector("double", 0), log = TRUE) {
      return(f_cdf.base(x, theta, y = y, log = log, do.mix = TRUE, do.shape.ind = FALSE))
    }
    
    mod$lower = as.vector(f.do.mix.reverse(mod$lower))
    newParamsLength = length(mod$lower)
    mod$upper = as.vector(f.do.mix.reverse(mod$upper))
    mod$theta0 = as.vector(f.do.mix.reverse(mod$theta0))
    mod$Sigma0 = mod$Sigma0[1:newParamsLength]
    mod$label = mod$label[1:newParamsLength]
    
    
    ineq_func = function(theta) {
      theta = as.vector(f.do.mix(theta))
      return(ineq_func.base(theta))
    }
    
  } else if (isTRUE(do.mix) && isTRUE(do.shape.ind)) {
    f.sim = function(n, theta, burnin = 500, do.state = FALSE) {
      return(f.sim.base(n, theta, burnin = burnin, do.state = do.state, 
        do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    f_Pstate = function(theta, y) {
      return(f_Pstate.base(theta, y, do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    f_Plast = function(theta, y) {
      return(f_Plast.base(theta, y, do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    f_unc_vol = function(theta, y = 0) {
      return(f_unc_vol.base(theta, y = y, do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    f_ht = function(theta, y) {
      return(f_ht.base(theta, y, do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    f_kernel = function(theta, y, log = TRUE) {
      return(f_kernel.base(theta, y, log = log, do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    f_pred = function(x, theta, y, log = TRUE) {
      return(f_pred.base(x, theta, y, log = log, do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    f_pit = function(x, theta, y, do.norm = FALSE) {
      return(f_pit.base(x, theta, y, do.norm = do.norm, do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    f_risk = function(theta, y, level = 0.95,  ES = TRUE) {
      return(f_risk.base(theta, y, level = level,  ES = ES, do.mix = TRUE, 
        do.shape.ind = TRUE))
    }
    
    f_rnd = function(n, theta, y = vector("double", 0), do.state = FALSE) {
      return(f_rnd.base(n, theta, y = y, do.state = do.state, do.mix = TRUE, 
        do.shape.ind = TRUE))
    }
    
    f_pdf = function(x, theta, y = vector("double", 0), log = TRUE) {
      return(f_pdf.base(x, theta, y = y, log = log, do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    f_cdf = function(x, theta, y = vector("double", 0), log = TRUE) {
      return(f_cdf.base(x, theta, y = y, log = log, do.mix = TRUE, do.shape.ind = TRUE))
    }
    
    mod$lower = as.vector(f.do.mix.reverse(mod$lower))
    mod$lower = as.vector(f.do.shape.ind.reverse(mod$lower))
    
    newParamsLength = length(mod$lower)
    
    mod$upper = as.vector(f.do.mix.reverse(mod$upper))
    mod$upper = as.vector(f.do.shape.ind.reverse(mod$upper))
    
    mod$theta0 = as.vector(f.do.mix.reverse(mod$theta0))
    mod$theta0 = as.vector(f.do.shape.ind.reverse(mod$theta0))
    
    mod$Sigma0 = mod$Sigma0[1:newParamsLength]
    mod$label = f.do.shape.ind.reverse(mod$label)
    mod$label = mod$label[1:newParamsLength]
    
    ineq_func = function(theta) {
      theta = as.vector(f.do.mix(theta))
      theta = as.vector(f.do.shape.ind(theta))
      return(ineq_func.base(theta))
    }
    
  } else if (!isTRUE(do.mix) && isTRUE(do.shape.ind)) {
    f.sim = function(n, theta, burnin = 500, do.state = FALSE) {
      return(f.sim.base(n, theta, burnin = burnin, do.state = do.state, 
        do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    f_Pstate = function(theta, y) {
      return(f_Pstate.base(theta, y, do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    f_Plast = function(theta, y) {
      return(f_Plast.base(theta, y, do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    f_unc_vol = function(theta, y = 0) {
      return(f_unc_vol.base(theta, y = y, do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    f_ht = function(theta, y) {
      return(f_ht.base(theta, y, do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    f_kernel = function(theta, y, log = TRUE) {
      return(f_kernel.base(theta, y, log = log, do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    f_pred = function(x, theta, y, log = TRUE) {
      return(f_pred.base(x, theta, y, log = log, do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    f_pit = function(x, theta, y, do.norm = FALSE) {
      return(f_pit.base(x, theta, y, do.norm = do.norm, do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    f_risk = function(theta, y, level = 0.95,  ES = TRUE) {
      return(f_risk.base(theta, y, level = level,  ES = ES, do.mix = FALSE, 
        do.shape.ind = TRUE))
    }
    
    f_rnd = function(n, theta, y = vector("double", 0), do.state = FALSE) {
      return(f_rnd.base(n, theta, y = y, do.state = do.state, do.mix = FALSE, 
        do.shape.ind = TRUE))
    }
    
    f_pdf = function(x, theta, y = vector("double", 0), log = TRUE) {
      return(f_pdf.base(x, theta, y = y, log = log, do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    f_cdf = function(x, theta, y = vector("double", 0), log = TRUE) {
      return(f_cdf.base(x, theta, y = y, log = log, do.mix = FALSE, do.shape.ind = TRUE))
    }
    
    mod$lower = as.vector(f.do.shape.ind.reverse(mod$lower))
    
    newParamsLength = length(mod$lower)
    
    mod$upper = as.vector(f.do.shape.ind.reverse(mod$upper))
    
    mod$theta0 = as.vector(f.do.shape.ind.reverse(mod$theta0))
    
    mod$Sigma0 = mod$Sigma0[1:newParamsLength]
    mod$label = f.do.shape.ind.reverse(mod$label)
    mod$label = mod$label[1:newParamsLength]
    
    ineq_func = function(theta) {
      theta = as.vector(f.do.shape.ind(theta))
      return(ineq_func.base(theta))
    }
    
  } else {
    
    f.sim = function(n, theta, burnin = 500, do.state = FALSE) {
      return(f.sim.base(n, theta, burnin = burnin, do.state = do.state, 
        do.mix = FALSE, do.shape.ind = FALSE))
    }
    
    f_Pstate = function(theta, y) {
      return(f_Pstate.base(theta, y, do.mix = FALSE, do.shape.ind = FALSE))
    }
    
    f_Plast = function(theta, y) {
      return(f_Plast.base(theta, y, do.mix = FALSE, do.shape.ind = FALSE))
    }
    
    f_unc_vol = function(theta, y = 0) {
      return(f_unc_vol.base(theta, y = y, do.mix = FALSE, do.shape.ind = FALSE))
    }
    
    f_ht = function(theta, y) {
      return(f_ht.base(theta, y, do.mix = FALSE, do.shape.ind = FALSE))
    }
    
    f_kernel = function(theta, y, log = TRUE) {
      return(f_kernel.base(theta, y, log = log, do.mix = FALSE, do.shape.ind = FALSE))
    }
    
    f_pred = function(x, theta, y, log = TRUE) {
      return(f_pred.base(x, theta, y, log = log, do.mix = FALSE, do.shape.ind = FALSE))
    }
    
    f_pit = function(x, theta, y, do.norm = FALSE) {
      return(f_pit.base(x, theta, y, do.norm = do.norm, do.mix = FALSE, do.shape.ind = FALSE))
    }
    
    f_risk = function(theta, y, level = 0.95,  ES = TRUE) {
      return(f_risk.base(theta, y, level = level,  ES = ES, do.mix = FALSE, 
        do.shape.ind = FALSE))
    }
    
    f_rnd = function(n, theta, y = vector("double", 0), do.state = FALSE) {
      return(f_rnd.base(n, theta, y = y, do.state = do.state, do.mix = FALSE, 
        do.shape.ind = FALSE))
    }
    
    f_pdf = function(x, theta, y = vector("double", 0), log = TRUE) {
      return(f_pdf.base(x, theta, y = y, log = log, do.mix = FALSE, do.shape.ind = FALSE))
    }
    
    f_cdf = function(x, theta, y = vector("double", 0), log = TRUE) {
      return(f_cdf.base(x, theta, y = y, log = log, do.mix = FALSE, do.shape.ind = FALSE))
    }
    ineq_func = ineq_func.base
  }
  names(mod$Sigma0) = mod$label
  out = list(f.sim = f.sim, f.ht = f_ht, f.kernel = f_kernel, f.unc.vol = f_unc_vol, 
    f.pred = f_pred, f.pit = f_pit, f.risk = f_risk, f.rnd = f_rnd, f.pdf = f_pdf, 
    f.cdf = f_cdf, f.ineq = ineq_func, f.Pstate = f_Pstate, f.Plast = f_Plast, 
    theta0 = mod$theta0, is.mix = do.mix, is.shape.ind = do.shape.ind, K = K, sigma0 = diag(mod$Sigma0),
    lower = mod$lower, upper = mod$upper, ineqlb = mod$ineq_lb, inequb = mod$ineq_ub, 
    n.params = n.params, n.params.vol = n.params.vol, do.init = F, 
    label = mod$label, name = mod$name)
  class(out) = "MSGARCH.SPEC"
  return(out)
}
