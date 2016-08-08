#' Transition Matrix.
#' @description Method returning the transition matrix.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}}.
#' @param theta Vector (of size d) of parameter estimates (not require when using a fit object).
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
#' fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)
#'
#'# run pdf method in-sample
#'transmat.mle = MSGARCH::transmat(fit)
#'
#'print(transmat.mle)
#'}
#' @return A matrix (of size K x K) in the case of a Markov-Switching model
#'  or a vector (of size K) in the case of a Mixture model.
#' @importFrom stats quantile
#' @export
transmat <- function(object, theta)
{
  UseMethod("transmat", object)
}

transmat.MSGARCH_SPEC = function(object, theta){
  
  if(isTRUE(object$is.shape.ind)){
    theta = object$func$f.do.shape.ind(theta = theta)
  }
  
  Nbparams = object$n.params
  Nmodel = length(Nbparams)
  Nbparams = object$n.params
  paramsLoc = c(0,cumsum(Nbparams))
  
  
  if(!isTRUE(object$is.mix)){
    
    p = matrix(nrow = Nmodel, ncol = Nmodel)
    for (i in 1:(Nmodel-1)){
      p[i,1:Nmodel] = theta[(paramsLoc[Nmodel+1]+1):(paramsLoc[Nmodel+1] + Nmodel)]
    }
    p[Nmodel,] = 1-colSums(matrix(p[1:(Nmodel-1),], ncol = Nmodel ))
  } else {
    p = rep(0,Nmodel)
    for (i in 1:(Nmodel-1)){
      p[i] = theta[(paramsLoc[Nmodel+1]+1)]
    }
    p[Nmodel] = 1-sum(p)
    
  }
  return(p)
}

#' @export
transmat.MSGARCH_MLE_FIT = function(object, theta = NULL) {
  
  return(MSGARCH::transmat(object = object$spec, theta = object$theta))
  
}

