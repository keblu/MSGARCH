#' Simulation function at T + 1.
#' @description Method returning random draws at T + 1.
#' @param spec Model specification of class \code{\link{MSGARCH_SPEC}} created with \code{\link{create.spec}}.
#' @param n  Number of random draws to be generated.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param do.state  Boolean indicating if the simulated state are also output. (default: \code{do.state = FALSE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @return A list containing one or two components:
#' \itemize{
#' \item \code{draws}: vector (of size n) or matrix (of size M x n) of simulated draws at T + 1.
#' \item \code{state}: vector (of size n) or matrix (of size M x n) of simulated states at T + 1.
#'  The \code{state} value appear only if \code{do.state = TRUE}.
#' }
#' @examples 
#' spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#' rnd = MSGARCH::rnd(spec = spec, n = 1000, theta = theta, y = sp500ret, do.state = TRUE)
#' @export
rnd <- function(spec, n, theta, y = vector("double", 0), do.state = FALSE)
{
  UseMethod("rnd", spec)
}

#' @export
rnd.MSGARCH_SPEC = function(spec, n, theta, y = vector("double", 0), do.state = FALSE) {
  
  y = as.matrix(y)
  if (isTRUE(spec$is.shape.ind)) {
    theta = spec$func$f.do.shape.ind(theta)
  }
  
  if (isTRUE(spec$is.mix)) {
    theta = spec$func$f.do.mix(theta)
  }
  
  if (is.vector(theta)) {
    theta = matrix(theta, nrow = 1)
  }
  
  draws = matrix(data = NA, nrow = nrow(theta), ncol = n)
  state = matrix(data = NA, nrow = nrow(theta), ncol = n)
  for(i in 1:nrow(theta)){
    tmp = spec$rcpp.func$rnd_Rcpp(n, theta[i,], y)
    if(spec$K == 1){
      draws[i,] = tmp
      state[i,] = rep(1,n)
    } else{
      draws[i,] = tmp$draws
      state[i,] = tmp$state
    }
    
  }
  out = list()
  out$draws = draws
  if (isTRUE(do.state)) {
    out$state = state
  }
  return(out)
}
