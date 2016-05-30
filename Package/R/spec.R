# creates the spec object from an Rcpp module

# input is now a list of C++ class:       e.g. list(Garch_normal_sym, Garch_student_skew)
# the input can also be a single model:  e.g. Garch_normal_sym
#'@export 
f.spec = function(models) {
  
  ###################################################
# create the relevant model C++ objects
K = length(models) # number of models


if (K > 1) {
  tmp = list()
  for(i in 1:K) {
    tmp[i] = new(models[[i]])
  }
}else{
  tmp = new(models[[1]])
}

if(K > 1){
  mod = new(MSgarch, tmp)
}else{
  mod = tmp
}

###################################################
  # detach for speed

  calc_ht    = mod$calc_ht
  eval_model = mod$eval_model
  ineq_func  = mod$ineq_func
  f_sim_Rcpp  = mod$f_sim
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

  
  #======================# model simulation #======================#

  f_sim = function(n, theta, burnin = 500, outputState = FALSE){

    out = f_sim_Rcpp(n, theta, burnin)
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
  f_Pstate = function(theta,y){
    
    f_get_Pstate_Rcpp(theta, y, FALSE)
  }
  
  
  f_Plast = function(theta,y){
    
    f_get_Pstate_Rcpp(theta, y, TRUE)
  }
  
  #======================# model simulation #======================#
  f_unc_vol = function(theta, y = 0){

    
    if (is.vector(theta))
      theta = matrix(theta, nrow = 1)
    
    ht = f_unc_vol_Rcpp(theta, y)
    t(ht)
    return(ht)
  }
  #=========================================================#
  
  #======================# variance #======================#
   f_ht= function(theta, y){
     
     
     if (is.vector(theta))
       theta = matrix(theta, nrow = 1)
    
     
     ht = calc_ht(theta, y)
     return(ht)
   }
  #=========================================================#
  
  #======================# kernel #======================#
  f_kernel = function(theta, y, log = TRUE){
    
    
    if (is.vector(theta))
      theta = matrix(theta, nrow = 1)
    

    lnd   = eval_model(theta, y)
    lnd[is.na(lnd) | is.nan(lnd) | is.infinite(lnd)] = -1e10
    if (!log)
      lnd = exp(lnd)
    return(lnd)
  } 
  #=========================================================#

 

# =======================# RND  #========================#
   f_rnd = function(n, theta, y = vector("double",0), outputState = FALSE){
     
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
   f_pdf = function(x, theta, y = vector("double",0), log = TRUE) {
     

     
     out = f_pdf_Rcpp(x, theta, y, log)
     return(out)
   }
  #==============================================================#
  
  #=======================# CDF  #========================#
   f_cdf = function(x, theta, y = vector("double",0), log = TRUE) {
     
     out = f_cdf_Rcpp(x, theta, y, log)
    return(out)
   }
  #==============================================================#
  
   
  out = list( f.sim      = f_sim,                    
              f.ht       = f_ht,                      
              f.kernel   = f_kernel,            
              f.unc.vol  = f_unc_vol,           
              f.rnd      = f_rnd,             
              f.pdf      = f_pdf,             
              f.cdf      = f_cdf,             
              f.ineq     = ineq_func,
              f.Pstate   = f_Pstate,
              f.Plast    = f_Plast,
              theta0     = mod$theta0,
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
