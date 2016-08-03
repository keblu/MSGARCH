#' Random draws at \code{t = T + 1} simulation method.
#' @description Method returning random draws at \code{t = T + 1}.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param n  Number of random draws to be generated.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @return A list of class \code{MSGARCH_RND} containing  two components:
#' \itemize{
#' \item \code{draws}: vector (of size n) or matrix (of size M x n) of simulated draws at \code{t = T + 1}.
#' \item \code{state}: vector (of size n) or matrix (of size M x n) of simulated states at \code{t = T + 1}.
#' }
#' The \code{MSGARCH_RND} class contains the \code{summary} and \code{plot} method.
#' @usage rnd(spec, n, theta, y)
#' rnd(fit, n)
#' @examples 
#'\dontrun{
#' # load data
#'data("sp500ret")
#'
#'# create model specification
#' spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)
#'  
#'# generate random draws
#'set.seed(123)
#' rnd = MSGARCH::rnd(fit, n = 1000)
#' 
#' plot(rnd)
#' 
#' summary(rnd)
#' }
#' @export
rnd <- function(spec, n, theta, y)
{
  UseMethod("rnd", spec)
}

#' @export
rnd.MSGARCH_SPEC = function(spec, n, theta, y) {
  
  y = f.check.y(y)

  theta = f.check.theta(spec, theta)
  
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
  out$state = state
  
  class(out) = "MSGARCH_RND"
  return(out)
}

#' @export
rnd.MSGARCH_MLE_FIT = function(fit,  n) {
  
  return(MSGARCH::rnd(spec = fit$spec, n = n, theta = fit$theta, y = fit$y))
  
}

#' @export
rnd.MSGARCH_BAY_FIT = function(fit,  n) {
  
  return(MSGARCH::rnd(spec = fit$spec, n = n, theta = fit$theta, y = fit$y))
  
}
