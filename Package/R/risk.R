f.risk.helper = function(theta, y, spec, level = c(0.95,0.99), ES = TRUE) {
  if (is.vector(theta)) 
    theta = matrix(theta, nrow = 1)
  p = 1 - level
  n = length(y)
  N = nrow(theta)
  np = length(p)
  xmin = min(y) - sd(y)
  xmax = 0
  itermax = 100
  tol = 1e-05

    tmp.VaR = NULL
    
    f.pdf = function(x) {
      out = spec$f.pred(x, theta, y, log = FALSE)
      return(out)
    }
    
    f.fun = function(x, pi) {
      out = integrate(f.pdf, lower = xmin, upper = x)$value - pi
      return(out)
    }
    
    # gross approximation for VaR
    out = list()
    for (i in 1:np) {
      tmp.VaR[i] = uniroot(f.fun, lower = xmin, upper = xmax, pi = p[i])$root
    }
    
    # add precision by Newton-Raphson
    out$VaR = vector("double", np)
    for (i in 1:np) {
      p_i = p[i]
      level_i = level[i]
      calc.step = function(V) {
        PDF = spec$f.pred(x = V, theta = theta, y = y, log = F)
        CDF = spec$f.pit(x = V, theta = theta, y = y)
        lPDF = log(PDF)
        err = p_i - CDF
        step = err * exp(-lPDF)
        return(list(step = step, err = abs(err)))
      }
      # Starting values
      VaR = tmp.VaR[i]
      # VaR calculation
      newStep = calc.step(VaR)
      delta = newStep$step
      VaR = VaR + delta
      covERR = newStep$err
      for (j in 1:itermax) {
        if (covERR < tol) {
          break
        }
        newStep = calc.step(VaR)
        delta = newStep$step
        VaR = VaR + delta
        covERR = newStep$err
      }
      out$VaR[i] = VaR
    }
  
  if (isTRUE(ES)) {
    out$ES = vector("double", np)
    for (i in 1:np) {
      
      f.condMean = function(x) {
        out = x * spec$f.pred(x, theta, y, log = FALSE)
        return(out)
      }
      out$ES[i] = integrate(f.condMean, lower = -Inf, upper = out$VaR[i], stop.on.error = FALSE)$value/p[i]
    }
  }
  return(out)
}
