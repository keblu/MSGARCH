optim.default <- function(f.nll, spec, theta0) {
  lower = spec$lower
  upper = spec$upper
  
  tmp = try(dfoptim::nmkb(par = theta0, fn = f.nll, lower = lower,
                          upper = upper,
                          control = list(maxfeval = 100*length(theta0)^2)),
            silent = TRUE)
  if (is(tmp, "try-error")){
    tmp = NULL
    tmp$convergence = 100
  }
  
  if (tmp$convergence > 0){
    str <- "optim.default -> failed optimization, trying with neldermead"
    MSGARCH:::f.error(str)
    tmp = nloptr::neldermead(x0 = theta0, fn = f.nll, lower = lower, upper = upper, control = list(maxeval = 10000))
  } 
  out = tmp$par
  
  return(out)
}