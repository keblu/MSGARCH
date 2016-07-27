#' Simulation function.
#' @description  Method returning a simulated process.
#' @param spec Model specification of class \code{\link{MSGARCH_SPEC}} created with \code{\link{create.spec}}.
#' @param n   Simulation length.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param burnin (integer >= 0) Burnin period discarded (first simulation draws). (default: \code{burnin = 500})
#' @param do.state  Boolean  indicating if the simulated state are also output. (default: \code{log = TRUE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#' spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                             do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#' y = MSGARCH::sim(spec = spec, n = 1000, theta = theta, burnin = 500, do.state = TRUE)
#' 
#' @return A list of class \code{\link{MSGARCH_SIM}} containing one or two components.
#' \itemize{
#' \item \code{draws}: vector (of size n) or matrix (of size M x n) of simulated draws.
#' \item \code{state}: vector (of size n) or matrix (of size M x n) of simulated states.
#'  The \code{state} value appear only if \code{do.state = TRUE}.
#' }
#' @export
sim <- function(spec, n, theta, burnin = 500, do.state = FALSE)
{
  UseMethod("sim", spec)
}

#' @export
sim.MSGARCH_SPEC = function(spec, n, theta, burnin = 500, do.state = FALSE) {
  
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
    tmp = spec$rcpp.func$sim(n, theta[i,], burnin)
    if(spec$K == 1 ){
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
  class(out) = "MSGARCH_SIM"
  return(out)
}