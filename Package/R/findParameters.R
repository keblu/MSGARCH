f.find.theta0 = function(f.kernel, theta0, lower = NULL, upper = NULL, 
                         f.ineq = NULL, ineqlb = NULL, inequb = NULL, 
                         type = c('robust', 'solnp')){
  
  ctr.optim   = list(trace = 0, maxit = 50000)
  ctr.deoptim = DEoptim::DEoptim.control(NP = 50 * length(theta0), itermax = 500, trace = FALSE)
  ctr.slsqp   = list(maxeval = 10000, xtol_rel = 1e-8)
  
  f.nll = function(x) -f.kernel(x, log = TRUE)
  
  require("nloptr")
  require("DEoptim")
  if (type[1] == 'robust') {
    str = "f.find.theta0 -> Robust estimation sequence - SQP"
    is.ok = tryCatch({ 
      tmp = nloptr::slsqp(x0 = theta0, fn = f.nll, lower = lower, upper = upper, control = ctr.slsqp)
      if (tmp$convergence > 1){
        out = tmp$par
      }
      out = tmp$par
      TRUE
    }, warning = function(warn){
      f.error(str)
    }, error = function(err){
      f.error(str)
    })
    
    if (!isTRUE(is.ok)) {
      str = "f.find.theta0 -> Robust estimation sequence - DEoptim for theta0"
      is.ok = tryCatch({ 
        tmp    = DEoptim::DEoptim(fn = f.nll, lower = lower, upper = upper, control = ctr.deoptim)
        theta0 = tmp$optim$bestmem
        tmp    = nloptr::slsqp(x0 = theta0, fn = f.nll, lower = lower, upper = upper, control = ctr.slsqp)
        if (tmp$convergence > 1){
          out = tmp$par
        }
        out = tmp$par
        TRUE
      }, warning = function(warn){
        f.error(str)
      }, error = function(err){
        f.error(str)
      })
      }
  }
  
  if ( (!isTRUE(is.ok) || type[1] == 'solnp') && !is.null(lower) && !is.null(upper) ) {
    str = "f.find.theta0 -> solnp"
    is.ok = tryCatch({ 
      require("Rsolnp")
      ctr.solnp = list(list(trace = 0))
      out = Rsolnp::solnp(pars = theta0, fun = f.nll, 
                          LB = lower, UB = upper, 
                          ineqfun = f.ineq, ineqLB = ineqlb, ineqUB = inequb,
                          control = ctr.solnp)$pars
      TRUE
    }, warning = function(warn){
      f.error(str)
    }, error = function(err){
      f.error(str)
    })
  }
  
  if (!isTRUE(is.ok)) {
    str = "f.find.theta0 -> starting value"
    f.error(str)
    out = theta0  
  }
  
  return(out)
}