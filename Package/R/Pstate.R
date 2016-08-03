#' Filtered state probabilities.
#' @description Method returning the filtered state probabilities.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#'\dontrun{
#'# load data
#'data("sp500ret")
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)
#'  
#'# compute the filtered state probabilities
#'Pstate  = MSGARCH::Pstate(fit)
#'
#'plot(Pstate)
#'}
#'@return Filtered state probabilities of class \code{MSGARCH_RND} (array of size (T + 1) x M x K).
#'The class \code{MSGARCH_RND} contains the \code{plot} method.
#'@usage Pstate(spec, theta, y)
#'Pstate(fit)
#' @export
Pstate <- function(spec, theta, y)
{
  UseMethod("Pstate", spec)
}

#' @export
Pstate.MSGARCH_SPEC = function(spec, theta, y) {
  
  y = as.matrix(y)
  
  theta = f.check.theta(spec, theta)
  
  out = array(dim = c(nrow(y) + 1, nrow(theta), spec$K))
  for(i in 1:nrow(theta)){
    tmp = spec$rcpp.func$get_Pstate_Rcpp(theta[i,], y, FALSE)
    for(j in 1:spec$K){
      out[,i,j] = tmp[,j]
    }
  }
  class(out) = "MSGARCH_PSTATE"
  return(out)
}

#' @export
rnd.MSGARCH_MLE_FIT = function(fit) {
  
  return(MSGARCH::Pstate(spec = fit$spec,  theta = fit$theta, y = fit$y))
  
}

#' @export
rnd.MSGARCH_BAY_FIT = function(fit) {
  
  return(MSGARCH::Pstate(spec = fit$spec,  theta = fit$theta, y = fit$y))
  
}