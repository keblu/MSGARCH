#' Process simulation method.
#' @description  Method simulating a \code{MSGARCH} process.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param n   Simulation length. (Default: \code{n = 1000})
#' @param m   Number of simulations. (Default: \code{m = 1})
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param burnin (integer >= 0) Burnin period discarded (first simulation draws). (Default: \code{burnin = 500})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually and \code{m = M}.  The difference between \code{\link{sim}} and \code{\link{simahead}} is that
#' \code{\link{sim}} starts the simulation a t = 0 creating an entire new process while  \code{\link{simahead}} starts the simulation at t = T + 1
#'  taking in consideration all the information available in the original time serie \code{y}.
#' @examples 
#'\dontrun{
#'# create model specification
#' spec = MSGARCH::create.spec() 
#'
#'# generate process
#' set.seed(123)
#' sim = MSGARCH::sim(object = spec, n = 1000, m = 1, theta = spec$theta0, burnin = 500)
#' 
#' plot(sim)
#' }
#' @return A list of class \code{MSGARCH_SIM} containing two components.
#' \itemize{
#' \item \code{draws}: Matrix (of size M x n) of simulated draws.
#' \item \code{state}: Matrix (of size M x n) of simulated states.
#' }
#'  The \code{MSGARCH_SIM} class contains the \code{plot} method.
#' @export
sim <- function(object, n, m, theta, burnin = 500) {
  UseMethod("sim", object)
}

#' @export
sim.MSGARCH_SPEC <- function(object, n = 1000, m = 1, theta = NULL, burnin = 500) {
  theta <- f.check.theta(object, theta)
  if (nrow(theta) == 1) {
    theta <- matrix(theta[rep(1, m), ], ncol = ncol(theta))
  }
  draws <- matrix(data = NA, nrow = nrow(theta), ncol = n)
  state <- matrix(data = NA, nrow = nrow(theta), ncol = n)
  for (i in 1:nrow(theta)) {
    tmp <- object$rcpp.func$sim(n, theta[i, ], burnin)
    if (object$K == 1) {
      draws[i, ] <- tmp
      state[i, ] <- rep(1, n)
    } else {
      draws[i, ] <- tmp$draws
      state[i, ] <- tmp$state
    }
  }
  out <- list()
  out$draws  <- draws
  out$state  <- state
  class(out) <- "MSGARCH_SIM"
  return(out)
}


#' @export
sim.MSGARCH_MLE_FIT <- function(object, n = 1000, m = 1, theta = NULL, burnin = 500) {
  return(MSGARCH::sim(object = object$spec, n = n, m = m, theta = object$theta,
                      burnin = burnin))
}

#' @export
sim.MSGARCH_BAY_FIT <- function(object, n = 1000, m = 1, theta = NULL, burnin = 500) {
  return(MSGARCH::sim(object = object$spec, n = n, m = m, theta = object$theta,
                      burnin = burnin))
}