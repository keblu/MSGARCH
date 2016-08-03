#' Process simulation method.
#' @description  Method returning a \code{MSGARCH} specification process.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param n   Simulation length.
#' @param m   Number of simulation.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param burnin (integer >= 0) Burnin period discarded (first simulation draws). (Default: \code{burnin = 500})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually and \code{m = nrow(theta)}.
#' @examples 
#'# create model specification
#' spec = MSGARCH::create.spec() 
#'
#'# generate process
#' set.seed(123)
#' sim = MSGARCH::sim(spec = spec, n = 1000, m = 1, theta = spec$theta0, burnin = 500)
#' 
#' plot(sim)
#' @return A list of class \code{MSGARCH_SIM} containing two components.
#' \itemize{
#' \item \code{draws}: Matrix (of size M x n) of simulated draws.
#' \item \code{state}: Matrix (of size M x n) of simulated states.
#' }
#'  The \code{MSGARCH_SIM} class contains the \code{plot} method.
#' @export
sim <- function(spec, n, m, theta, burnin = 500)
{
  UseMethod("sim", spec)
}

#' @export
sim.MSGARCH_SPEC = function(spec, n, m = 1, theta, burnin = 500) {
  
  theta = f.check.theta(spec, theta)
  
  if(nrow(theta) == 1){
    theta = matrix(theta[rep(1,m),], nrow = nrow(theta))
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
  out$state = state

  class(out) = "MSGARCH_SIM"
  return(out)
}