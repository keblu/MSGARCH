f.pred.helper = function(x, theta, y, spec, log = TRUE) {
  if (is.vector(theta)) 
    theta = matrix(theta, nrow = 1)
  N = nrow(theta)
  nx = length(x)
  
  out = matrix(data = NA, nrow = N, ncol = nx)
  for (i in 1:N) {
    out[i, ] = spec$f.pdf(x, theta = theta[i, ], y = y, log = FALSE)
  }
  out = colMeans(out)
  if (log) {
    out = log(out)
  }
  return(out)
}