# creates the spec object from an Rcpp module

# input is now a list of C++ class:       e.g. list(Garch_normal_sym, Garch_student_skew)
# the input can also be a single model:  e.g. Garch_normal_sym
#'@export 
f.spec = function(models, mixture = FALSE, DistRegInd = FALSE) {
  
  ###################################################
# create the relevant model C++ objects
K = length(models) # number of models

if(K == 1){
  mixture = FALSE
  DistRegInd = FALSE
}
options(warn=-1)
if (K > 1) {
  tmp = list()
  for(i in 1:K) {
    tmp[i] = new(models[[i]])
  }
}else{
  tmp = new(models[[1]])
}
options(warn=0)
if(K > 1){
  mod = new(MSgarch, tmp)
}else{
  mod = tmp
}
###################################################
  # detach for speed
  dist = NULL
  name = mod$name
  for (i in 1:length(name)) {
    dist[i] = stringr::str_sub(name[i],start = stringr::str_locate(name,"_")[i,1]+1,nchar(name[i]))
  }
  uniqueDist = unique(dist)
  if(isTRUE(DistRegInd) && length(uniqueDist) > 1){
    stop("The distribution of each regime must be the same if the distribution are not regime dependent")
  }
  calc_ht    = mod$calc_ht
  eval_model = mod$eval_model
  ineq_func.base  = mod$ineq_func
  f_sim      = mod$f_sim
  NbParams   = mod$NbParams
  NbParamsModel   = mod$NbParamsModel
  f_pdf_Rcpp = mod$f_pdf
  f_cdf_Rcpp = mod$f_cdf
  f_rnd_Rcpp = mod$f_rnd
  f_unc_vol_Rcpp = mod$f_unc_vol
  if( K > 1 ){
    f_get_Pstate_Rcpp = mod$f_get_Pstate
  } else {
    f_get_Pstate_Rcpp = function(theta,y,PLast) {
      if (!isTRUE(PLast)) {
        out = matrix(1,nrow = length(y), ncol = 1)
      } else {
        out = matrix(1,nrow = 1, ncol = 1)
      }
      return(out)
    }
  }
  NbtotalParams = sum(NbParams)
  f.mixture = function(theta) {
    return(MSGARCH::f.theta.mixture(K,NbtotalParams,theta))
  }
  
  f.mixture.reverse = function(theta) {
    return(MSGARCH::f.theta.mixture.reverse(K, NbtotalParams, theta))
  }
  
  f.RegIndDist = function(theta) {
    return(MSGARCH::f.theta.RegIndDist(K, NbParams, NbParamsModel, theta))
  }
  
  f.RegIndDist.reverse = function(theta) {
    return(MSGARCH::f.theta.RegIndDist.reverse(K, NbParams, NbParamsModel, theta))
  }
  
  #======================# model simulation #======================#

  f.sim.base = function(n, theta, burnin = 500, outputState = FALSE, mixture = FALSE, DistRegInd = FALSE){
    if(isTRUE(DistRegInd)){
      theta = f.RegIndDist(theta)
    }
    
    if(isTRUE(mixture)){
      theta = f.mixture(theta)
    }
    
    out = f_sim(n, theta, burnin)
    if (K == 1){
      outputState = TRUE
    }
    if (isTRUE(outputState)) {
      return(out)
    } else {
      return(out$value)
    }
  }
 
  #=========================================================#
  f_Pstate.base = function(theta,y, mixture = FALSE, DistRegInd = FALSE){
    
    
    if(isTRUE(DistRegInd)){
      theta = f.RegIndDist(theta)
    }
    
    if(isTRUE(mixture)){
      theta = f.mixture(theta)
    }
    
    f_get_Pstate_Rcpp(theta, y, FALSE)
  }
  
  
  f_Plast.base = function(theta,y, mixture = FALSE, DistRegInd = FALSE){
    
    
    if(isTRUE(DistRegInd)){
      theta = f.RegIndDist(theta)
    }
    
    if(isTRUE(mixture)){
      theta = f.theta.mixture(K, NbParams, theta)
    }
    
    f_get_Pstate_Rcpp(theta, y, TRUE)
  }
  
  #======================# model simulation #======================#
  f_unc_vol.base = function(theta, y = 0, mixture = FALSE, DistRegInd = FALSE){
    
    if(isTRUE(DistRegInd)){
      theta = f.RegIndDist(theta)
    }
    
    if(isTRUE(mixture)){
      theta = f.mixture(theta)
    }
    
    if (is.vector(theta))
      theta = matrix(theta, nrow = 1)
    
  
    
    ht = f_unc_vol_Rcpp(theta, y)
    t(ht)
    return(ht)
  }
  #=========================================================#
  
  #======================# variance #======================#
   f_ht.base = function(theta, y, mixture = FALSE, DistRegInd = FALSE){
     
     if(isTRUE(DistRegInd)){
       theta = f.RegIndDist(theta)
     }
     
     if(isTRUE(mixture)){
       theta = f.mixture(theta)
     }
     
     if (is.vector(theta))
       theta = matrix(theta, nrow = 1)
    
     
     ht = calc_ht(theta, y)
     return(ht)
   }
  #=========================================================#
  
  #======================# kernel #======================#
  f_kernel.base = function(theta, y, log = TRUE, thresh = Inf, mixture = FALSE, DistRegInd = FALSE){
    
    if(isTRUE(DistRegInd)){
      theta = f.RegIndDist(theta)
    }
    
    if(isTRUE(mixture)){
      theta = f.mixture(theta)
    }
    
    if (is.vector(theta))
      theta = matrix(theta, nrow = 1)
    

    lnd   = eval_model(theta, y, thresh)
    lnd[is.na(lnd) | is.nan(lnd) | is.infinite(lnd)] = -1e10
    if (!log)
      lnd = exp(lnd)
    return(lnd)
  } 
  #=========================================================#
  
# =======================# RND  #========================#
   f_rnd.base = function(n, theta, y = vector("double",0), outputState = FALSE, mixture = FALSE, DistRegInd = FALSE){
     
     if(isTRUE(DistRegInd)){
       theta = f.RegIndDist(theta)
     }
     
     if(isTRUE(mixture)){
       theta = f.mixture(theta)
     }
     
     out = list()
     out = f_rnd_Rcpp(n, theta, y)
     if (isTRUE(outputState)) {
       return(out)
     } else {
       return(out$value)
     }
   }
  #==============================================================#
  
  #=======================# PDF  #========================#
   f_pdf.base = function(x, theta, y = vector("double",0), log = TRUE, mixture = FALSE, DistRegInd = FALSE) {
     
     if(isTRUE(DistRegInd)){
       theta = f.RegIndDist(theta)
     }
     
     if(isTRUE(mixture)){
       theta = f.mixture(theta)
     }
     
     out = f_pdf_Rcpp(x, theta, y, log)
     return(out)
   }
  #==============================================================#
  
  #=======================# CDF  #========================#
   f_cdf.base = function(x, theta, y = vector("double",0), log = TRUE, mixture = FALSE, DistRegInd = FALSE) {
     
     if(isTRUE(DistRegInd)){
       theta = f.RegIndDist(theta)
     }
     
     if(isTRUE(mixture)){
       theta = f.mixture(theta)
     }
     
     out = f_cdf_Rcpp(x, theta, y, log)
    return(out)
   }
  #==============================================================#
  
   if(isTRUE(mixture) && !isTRUE(DistRegInd)){
     f.sim = function(n, theta, burnin = 500, outputState = FALSE){
       return (f.sim.base(n, theta, burnin = burnin, outputState = outputState, mixture = TRUE, DistRegInd = FALSE))
     }
     
     f_Pstate = function(theta,y){
      return (f_Pstate.base(theta,y, mixture = TRUE, DistRegInd = FALSE))
     }
     
     f_Plast = function(theta,y){
       return (f_Plast.base(theta,y, mixture = TRUE, DistRegInd = FALSE))
     }
     
     f_unc_vol = function(theta, y = 0){
       return(f_unc_vol.base(theta, y = y, mixture = TRUE, DistRegInd = FALSE))
     }
     
     f_ht = function(theta, y){
       return(f_ht.base(theta, y, mixture = TRUE, DistRegInd = FALSE))
     }
     
     f_kernel = function(theta, y, log = TRUE, thresh = Inf){
       return(f_kernel.base(theta, y, log = log, thresh = thresh, mixture = TRUE, DistRegInd = FALSE))
     }
     
     f_rnd = function(n, theta, y = vector("double",0), outputState = FALSE){
       return(f_rnd.base(n, theta, y = y, outputState = outputState, mixture = TRUE, DistRegInd = FALSE))
     }
     
     f_pdf = function(x, theta, y = vector("double",0), log = TRUE){
       return(f_pdf.base(x, theta, y = y, log = log, mixture = TRUE, DistRegInd = FALSE))
     }
     
     f_cdf = function(x, theta, y = vector("double",0), log = TRUE){
       return(f_cdf.base(x, theta, y = y, log = log, mixture = TRUE, DistRegInd = FALSE))
     }
     
     mod$lower = as.vector(f.mixture.reverse(mod$lower))
     newParamsLength = length(mod$lower)
     mod$upper = as.vector(f.mixture.reverse(mod$upper))
     mod$theta0 =  as.vector(f.mixture.reverse(mod$theta0))
     mod$Sigma0 = mod$Sigma0[1:newParamsLength]
     mod$label = mod$label[1:newParamsLength]
     ineq_func = function(theta){
       theta = as.vector(f.mixture(theta))
       return(mod$ineq_func.base(theta))
     }
     
   } else if(isTRUE(mixture) && isTRUE(DistRegInd)){
     f.sim = function(n, theta, burnin = 500, outputState = FALSE){
       return (f.sim.base(n, theta, burnin = burnin, outputState = outputState, mixture = TRUE, DistRegInd = TRUE))
     }
     
     f_Pstate = function(theta,y){
       return (f_Pstate.base(theta,y, mixture = TRUE, DistRegInd = TRUE))
     }
     
     f_Plast = function(theta,y){
       return (f_Plast.base(theta,y, mixture = TRUE, DistRegInd = TRUE))
     }
     
     f_unc_vol = function(theta, y = 0){
       return(f_unc_vol.base(theta, y = y, mixture = TRUE, DistRegInd = TRUE))
     }
     
     f_ht = function(theta, y){
       return(f_ht.base(theta, y, mixture = TRUE, DistRegInd = TRUE))
     }
     
     f_kernel = function(theta, y, log = TRUE, thresh = Inf){
       return(f_kernel.base(theta, y, log = log, thresh = thresh, mixture = TRUE, DistRegInd = TRUE))
     }
     
     f_rnd = function(n, theta, y = vector("double",0), outputState = FALSE){
       return(f_rnd.base(n, theta, y = y, outputState = outputState, mixture = TRUE, DistRegInd = TRUE))
     }
     
     f_pdf = function(x, theta, y = vector("double",0), log = TRUE){
       return(f_pdf.base(x, theta, y = y, log = log, mixture = TRUE, DistRegInd = TRUE))
     }
     
     f_cdf = function(x, theta, y = vector("double",0), log = TRUE){
       return(f_cdf.base(x, theta, y = y, log = log, mixture = TRUE, DistRegInd = TRUE))
     }
     
     mod$lower = as.vector(f.mixture.reverse(mod$lower))
     mod$lower = as.vector(f.RegIndDist.reverse(mod$lower))
     
     newParamsLength = length(mod$lower)
     
     mod$upper = as.vector(f.mixture.reverse(mod$upper))
     mod$upper = as.vector(f.RegIndDist.reverse(mod$upper))
     
     mod$theta0 =  as.vector(f.mixture.reverse(mod$theta0))
     mod$theta0 =  as.vector(f.RegIndDist.reverse(mod$theta0))
     
     mod$Sigma0 = mod$Sigma0[1:newParamsLength]
     mod$label = f.RegIndDist.reverse(mod$label)
     mod$label = mod$label[1:newParamsLength]
     
     ineq_func = function(theta){
       theta = as.vector(f.mixture(theta))
       theta = as.vector(f.RegIndDist(theta))
       return(mod$ineq_func.base(theta))
     }
     
   } else if (!isTRUE(mixture) && isTRUE(DistRegInd)){
     f.sim = function(n, theta, burnin = 500, outputState = FALSE){
       return (f.sim.base(n, theta, burnin = burnin, outputState = outputState, mixture = FALSE, DistRegInd = TRUE))
     }
     
     f_Pstate = function(theta,y){
       return (f_Pstate.base(theta,y, mixture = FALSE, DistRegInd = TRUE))
     }
     
     f_Plast = function(theta,y){
       return (f_Plast.base(theta,y, mixture = FALSE, DistRegInd = TRUE))
     }
     
     f_unc_vol = function(theta, y = 0){
       return(f_unc_vol.base(theta, y = y, mixture = FALSE, DistRegInd = TRUE))
     }
     
     f_ht = function(theta, y){
       return(f_ht.base(theta, y, mixture = FALSE, DistRegInd = TRUE))
     }
     
     f_kernel = function(theta, y, log = TRUE, thresh = Inf){
       return(f_kernel.base(theta, y, log = log, thresh = thresh, mixture = FALSE, DistRegInd = TRUE))
     }
     
     f_rnd = function(n, theta, y = vector("double",0), outputState = FALSE){
       return(f_rnd.base(n, theta, y = y, outputState = outputState, mixture = FALSE, DistRegInd = TRUE))
     }
     
     f_pdf = function(x, theta, y = vector("double",0), log = TRUE){
       return(f_pdf.base(x, theta, y = y, log = log, mixture = FALSE, DistRegInd = TRUE))
     }
     
     f_cdf = function(x, theta, y = vector("double",0), log = TRUE){
       return(f_cdf.base(x, theta, y = y, log = log, mixture = FALSE, DistRegInd = TRUE))
     }
     
     mod$lower = as.vector(f.RegIndDist.reverse(mod$lower))
     
     newParamsLength = length(mod$lower)
     
     mod$upper = as.vector(f.RegIndDist.reverse(mod$upper))
     
     mod$theta0 =  as.vector(f.RegIndDist.reverse(mod$theta0))
     
     mod$Sigma0 = mod$Sigma0[1:newParamsLength]
     mod$label = f.RegIndDist.reverse(mod$label)
     mod$label = mod$label[1:newParamsLength]
     
     ineq_func = function(theta){
       theta = as.vector(f.RegIndDist(theta))
       return(mod$ineq_func.base(theta))
     }
     
   } else {
     
     f.sim = function(n, theta, burnin = 500, outputState = FALSE){
       return (f.sim.base(n, theta, burnin = burnin, outputState = outputState, mixture = FALSE, DistRegInd = FALSE))
     }
     
     f_Pstate = function(theta,y){
       return (f_Pstate.base(theta,y, mixture = FALSE, DistRegInd = FALSE))
     }
     
     f_Plast = function(theta,y){
       return (f_Plast.base(theta,y, mixture = FALSE, DistRegInd = FALSE))
     }
     
     f_unc_vol = function(theta, y = 0){
       return(f_unc_vol.base(theta, y = y, mixture = FALSE, DistRegInd = FALSE))
     }
     
     f_ht = function(theta, y){
       return(f_ht.base(theta, y, mixture = FALSE, DistRegInd = FALSE))
     }
     
     f_kernel = function(theta, y, log = TRUE, thresh = Inf){
       return(f_kernel.base(theta, y, log = log, thresh = thresh, mixture = FALSE, DistRegInd = FALSE))
     }
     
     f_rnd = function(n, theta, y = vector("double",0), outputState = FALSE){
       return(f_rnd.base(n, theta, y = y, outputState = outputState, mixture = FALSE, DistRegInd = FALSE))
     }
     
     f_pdf = function(x, theta, y = vector("double",0), log = TRUE){
       return(f_pdf.base(x, theta, y = y, log = log, mixture = FALSE, DistRegInd = FALSE))
     }
     
     f_cdf = function(x, theta, y = vector("double",0), log = TRUE){
       return(f_cdf.base(x, theta, y = y, log = log, mixture = FALSE, DistRegInd = FALSE))
     }
     ineq_func = ineq_func.base
   }
   
  out = list( f.sim      = f.sim,                    
              f.ht       = f_ht,                      
              f.kernel   = f_kernel,            
              f.unc.vol  = f_unc_vol,
              f.pred     = f_pred,            
              f.rnd      = f_rnd,             
              f.pdf      = f_pdf,             
              f.cdf      = f_cdf,             
              f.ineq     = ineq_func,
              f.Pstate   = f_Pstate,
              f.Plast    = f_Plast,
              theta0     = mod$theta0,
              mixture    = mixture,
              DistRegInd = DistRegInd,
              K          = K,
              Sigma0     = diag(mod$Sigma0),     
              kSigma     = 1,                   
              lower      = mod$lower,            
              upper      = mod$upper,                         
              ineqlb     = mod$ineq_lb,       
              inequb     = mod$ineq_ub, 
              NbParams   = mod$NbParams,
              NbParamsModel = mod$NbParamsModel,
              do.init    = F,             
              label      = mod$label,         
              name       = mod$name)          
  
  return(out)
}
