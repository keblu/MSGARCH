#' Probability density function.
#' @description Method returning the probability density in-sample or of a vector of points at \code{t = T + 1}.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param x Vector (of size N) of point at \code{t = T + 1} to be evaluated (used when \code{is.its = FALSE}).
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param log  Boolean indicating if the log-density is returned. (Default: \code{log = FALSE})
#' @param is.its  Boolean indicating if the in-sample pdf is returned. (Default: \code{is.its = FALSE})
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually. 
#' If \code{is.its = FALSE}, the points \code{x} are evaluated as \code{t = T + 1} realization and the method uses the variance estimate at \code{t = T + 1}.
#' If \code{is.its = TRUE}, \code{y} is evaluated using their respective variance estimate at each time \code{t}.
#' @usage pdf(spec, theta, y, log = TRUE, is.its = TRUE)
#' pdf(fit, log = FALSE, is.its = TRUE) 
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
#'# run pdf method in-sample
#'pdf.its = MSGARCH::pdf(fit, log = FALSE, is.its = TRUE)
#'    
#'plot(pdf.its)          
#'                                                                                                                                                                                                                              
#'# create mesh
#'x = seq(-3,3,0.01)
#'
#'# run pdf method on mesh at T + 1
#'pdf = MSGARCH::pdf(fit, x = x, log = FALSE, is.its = FALSE)
#'
#'plot(pdf)
#'}
#' @return A list of class \code{MSGARCH_PDF} containing two components:
#' \itemize{
#' \item \code{pdf}:\cr If \code{is.its = FALSE}: (Log-)Probability density of the points \code{x} at \code{t = T + 1} (vector of size N or matrix of size M x N) \cr
#'                   If \code{is.its = TRUE}: In-sample (Log-)Probability density of \code{y} (vector of size T or matrix of size M x T). 
#' \item \code{x}:\cr If \code{is.its = FALSE}: Vector (of size N) of point at \code{t = T + 1} evaluated.\cr
#'                 If \code{is.its = TRUE}: Vector (of size T) of observations.
#' }
#'The class \code{MSGARCH_PDF} contains the \code{plot} method.
#' @usage pdf(spec, theta, y, log = TRUE, is.its = TRUE)
#' pdf(fit, log = TRUE, is.its = TRUE) 
#' @export
pdf <- function(spec, x, theta, y, log = FALSE, is.its = FALSE)
{
  UseMethod("pdf", spec)
}

#' @export
pdf.MSGARCH_SPEC = function(spec, x = NULL, theta, y, log = FALSE, is.its = FALSE) {
  
  y = f.check.y(y)
  
  theta = f.check.theta(spec, theta)
  if(isTRUE(is.its)){
    x = y
    tmp = matrix(data = NA,nrow = nrow(theta), ncol = length(y)-1)
    for(i in 1:nrow(theta)){
     tmp2 = matrix(data = NA,nrow = length(y)-1 , ncol = spec$K)
      if(spec$K == 1){
        tmp2 = spec$rcpp.func$pdf_Rcpp_its(theta[i,] ,y, log)
        tmp[i,] = tmp2
      } else{
        Pstate = MSGARCH::Pstate(spec = spec,theta = theta[i,], y = y)
        Pstate.tmp = matrix(data = NA, nrow = dim(Pstate)[1], ncol = dim(Pstate)[3])
        for(j in 1:dim(Pstate)[3]){
          Pstate.tmp[,j] = Pstate[,,j]
        }
        tmp2 = spec$rcpp.func$pdf_Rcpp_its( theta[i,] ,y, log)
        tmp[i,] = rowSums(tmp2 * Pstate.tmp[2:(nrow(Pstate.tmp)-1),])
      }
    }
    tmp = cbind(NA, tmp)
  } else{
    tmp = matrix(data = NA,nrow = nrow(theta),ncol = length(x))
    for(i in 1:nrow(theta)){
      tmp[i,] = spec$rcpp.func$pdf_Rcpp(x, theta[i,], y, log)
    }
  }
  out = list()
  out$pdf = tmp
  out$x = x
  out$is.its = is.its
  class(out) = "MSGARCH_PDF"
  return(out)
}

#' @export
pdf.MSGARCH_MLE_FIT = function(fit, x = NULL, log = TRUE, is.its = FALSE) {
  
  return(MSGARCH::pdf(spec = fit$spec, x =  x, theta = fit$theta, y = fit$y, log = log, is.its = is.its))
  
}

#' @export
pdf.MSGARCH_BAY_FIT = function(fit, x = NULL, log = TRUE, is.its = FALSE) {
  
  return(MSGARCH::pdf(spec = fit$spec, x =  x, theta = fit$theta, y = fit$y, log = log, is.its = is.its))
  
}