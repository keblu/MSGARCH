f.find.theta0 <- function(f.kernel, theta0, lower = NULL, upper = NULL, f.ineq = NULL,
                          ineqlb = NULL, inequb = NULL) {
  f.nll <- function(x) -f.kernel(x, log = TRUE)
  
  ctr.deoptim <- DEoptim::DEoptim.control(NP = 50 * length(theta0), itermax = 500, trace = FALSE)
  
  tmp = suppressWarnings(dfoptim::nmkb(par = theta0, fn = f.nll, lower = lower, upper = upper,control = list(maxfeval = 100*length(theta0)^2)))
  if (tmp$convergence > 0){
    str <- "f.find.theta0 -> failed optimization, trying with neldermead"
    f.error(str)
    tmp = nloptr::neldermead(x0 = theta0, fn = f.nll, lower = lower, upper = upper, control = list(maxeval = 10000))
    if(tmp$convergence < 0){
      str <- "f.find.theta0 -> failed optimization, trying with DEoptim"
      f.error(str)
      tmp     = DEoptim::DEoptim(fn = f.nll, lower = lower, upper = upper, control = ctr.deoptim)
      tmp$par = tmp$optim$bestmem
    }
  } 
  out = tmp$par
  
  return(out)
}