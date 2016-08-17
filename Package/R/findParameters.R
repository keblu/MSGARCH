f.find.theta0 <- function(f.kernel, theta0, lower = NULL, upper = NULL, f.ineq = NULL,
                         ineqlb = NULL, inequb = NULL, type = c("robust")) {
  ctr.slsqp <- list(maxeval = 10000, xtol_rel = 1e-08)
  f.nll <- function(x) -f.kernel(x, log = TRUE)
  if (type[1] == "robust") {
    str <- "f.find.theta0 -> Robust estimation sequence - SQP"
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
  }
    if (!isTRUE(is.ok)) {
      ctr.deoptim <- DEoptim::DEoptim.control(NP = 50 * length(theta0), itermax = 500,
                                              trace = FALSE)
      str <- "f.find.theta0 -> Robust estimation sequence - DEoptim for theta0"
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
  if (!isTRUE(is.ok)) {
    str <- "f.find.theta0 -> starting value"
    f.error(str)
    out <- theta0
  }
  return(out)
}