#' Transition Matrix.
#' @description Method returning the transition matrix.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}}.
#' @param theta Vector (of size d) of parameter estimates (not require when using a fit object).
#' @param n Number of steps ahead. (Default: \code{n = 1}
#' @examples 
#'# load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#' fit = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = list(do.init = FALSE))
#'
#'# Extract the transition matrix 10 steps ahead
#'transmat.mle = MSGARCH::transmat(fit, n = 10)
#'
#'print(transmat.mle)
#' @return A matrix (of size K x K) in the case of a Markov-Switching model
#'  or a vector (of size K) in the case of a Mixture model. 
#'  The columns indcates the starting states while the rows indicates the transition states. 
#' @importFrom stats quantile
#' @import expm
#' @export
transmat <- function(object, theta, n) {
  UseMethod(generic = "transmat", object =  object)
}

transmat.MSGARCH_SPEC <- function(object, theta = NULL, n = 1) {
  if (isTRUE(object$is.shape.ind)) {
    theta <- object$func$f.do.shape.ind(theta = theta)
  }
  Nbparams <- object$n.params
  n_model <- length(Nbparams)
  Nbparams <- object$n.params
  params_loc <- c(0, cumsum(Nbparams))
  if (!isTRUE(object$is.mix)) {
    p <- matrix(nrow = n_model, ncol = n_model)
    for (i in 0:(n_model-1)) {

        p[1:(n_model-1),i+1] <- theta[(params_loc[n_model + 1] + n_model*i+1-i):(params_loc[n_model + 1] + n_model*i+n_model-1-i)]
      
    }
    p[n_model, ] <- 1 - colSums(matrix(p[1:(n_model-1), ], ncol = n_model))
  } else {
    p <- matrix(rep(0, n_model), ncol = n_model)
    for (i in 1:(n_model - 1)) {
      p[1,i] <- theta[(params_loc[n_model + 1] + i)]
    }
    p[1,n_model] <- 1 - sum(p)
  }
  if(!object$is.mix){
      p = p %^% n
  }

  if(object$is.mix){
    col_label = paste0("State ", 1:object$K)
    row_label = paste0("Probability")
  } else {     
    col_label = paste0("t = ", 1:object$K)
    row_label = paste0("t + ",n," = ", 1:object$K)
  } 
  rownames(p) = row_label
  colnames(p) = col_label
  return(p)
}

#' @export
transmat.MSGARCH_MLE_FIT <- function(object, theta = NULL, n = 1) {
  return(MSGARCH::transmat(object = object$spec, theta = object$theta, n = n))
}