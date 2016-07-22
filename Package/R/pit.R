f.pit = function(x, theta, y, spec, do.norm = FALSE) {
  if (is.vector(theta)) 
    theta = matrix(theta, nrow = 1)
  N = nrow(theta)
  nx = length(x)
  
  out = matrix(data = NA, nrow = N, ncol = nx)
  for (i in 1:N) {
    out[i, ] = spec$f.cdf(x, theta = theta[i, ], y = y, log = FALSE)
  }
  out = colMeans(out)
  if (do.norm) {
    out = qnorm(out, mean = 0, sd = 1)
  }
  if (any(is.nan(out))) {
    stop("NaN value in PIT calculation")
  }
  return(out)
}