
#' Default estimation parameters
f.process.ctr = function(ctr = list()) {
  con = list(theta0 = NULL, do.init = FALSE)
  con[names(ctr)] = ctr
  return(con)
}

#' ML estimation of models
#' @param y time series
#' @param spec specification object
#' @param ctr : control parameters
#' @param delta : modify the parameters range
#' @return theta0 and negative log-likelihood
#' @export
f.estim.mle = function(y, spec, ctr = list(), delta = 0){
  
  require("DEoptim")
  ctr.optim   = list(trace = 0, maxit = 50000)
  ctr.deoptim = DEoptim::DEoptim.control(NP = 50 * length(spec$theta0), itermax = 500, trace = 50,initialpop = matrix(spec$theta0,nrow = 50 * length(spec$theta0),ncol = length(spec$theta0)))
  ctr.slsqp   = list(maxeval = 10000, xtol_rel = 1e-8)
  
  ctr = MSGARCH::f.process.ctr(ctr)
  
  f.kernel = function(x, log = TRUE){
      return(spec$f.kernel(x, y = y, log = log))
  }
  
  lower = spec$lower + delta
  upper = spec$upper - delta
  
  
  f.nll = function(x) -f.kernel(x, log = TRUE)

  if (any(ctr$do.init || spec$do.init)) {
      str = "f.find.theta0 -> DEoptim initialization"
      is.ok = tryCatch({ 
        tmp = DEoptim::DEoptim(fn = f.nll, lower = lower, upper = upper, control = ctr.deoptim)
        theta0.init = tmp$optim$bestmem
        TRUE
      }, warning = function(warn){
        MSGARCH::f.error(str)
      }, error = function(err){
        MSGARCH::f.error(str)
      })
    } else {
    theta0.init = spec$theta0
  }
  
  theta0 = MSGARCH::f.find.theta0(f.kernel, theta0 = theta0.init, 
                         lower = lower, upper = upper,
                         f.ineq = spec$f.ineq, ineqlb = spec$ineqlb, inequb = spec$inequb)
  nl_likelihood = f.kernel(theta0)
  
  if(nl_likelihood == -1e+10) {
    tmp = DEoptim::DEoptim(fn = f.nll, lower = lower, upper = upper, control = ctr.deoptim)
    theta0 = tmp$optim$bestmem
    nl_likelihood = f.kernel(theta0)
  }
  
  out = list(theta0 = theta0, nl_likelihood = nl_likelihood)
  return(out)
}
