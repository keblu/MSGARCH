#' Value-at-Risk And Expected-shortfall.
#' @description Method returning the Value-at-Risk and Expected-shortfall in-sample or at \code{t = T + 1} based on the predictive density.
#'@param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not required when using a fit object) where d must have
#'  the same length as the default parameters of the specification.
#' @param y  Vector (of size T) of observations (not required when using a fit object).
#' @param level Vector (of size R) of Value-at-risk and Expected-shortfall levels.\cr
#'  (Default: \code{level = c(0.95,0.99)})
#' @param ES  Boolean indicating if Expected-shortfall is also calculated. (Default: \code{ES = TRUE})
#' @param do.its  Boolean indicating if the in-sample risk estimators are returned.
#'  (Default: \code{do.its = FALSE})
#' @param ctr List of control parameters for risk evaluation.
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian Value-at-Risk and Expected-shortfall are calculated.
#' If \code{do.its = FALSE}, \code{x} the risk estimator at \code{t = T + 1}, the method uses the variance estimated at \code{t = T + 1}.
#' If \code{do.its = TRUE}, the in-sample risk estimator are calculated.
#' @return A list of class \code{MSGARCH_RISK} containing two or three components:
#' \itemize{
#' \item \code{VaR} : \cr If \code{do.its = FALSE}: Value-at-Risk at \code{t = T + 1} at the chosen levels (vector of size R).\cr
#'                    If \code{do.its = TRUE}: In-sample Value-at-Risk at the chosen levels (Matrix of size T x R).
#' \item \code{ES}  :\cr If \code{do.its = FALSE}: Expected-shortfall at \code{t = T + 1} at the chosen levels (vector of size R).\cr
#'                    If \code{do.its = TRUE}: In-sample Expected-shortfall at the chosen levels (Matrix of size T x R).
#' \item \code{y}  : Vector (of size T) of observations.
#' }
#' The \code{MSGARCH_RISK} contains the \code{plot} method. 
#' The Bayesian risk estimator can take long time to calculate depending on the size of the MCMC chain.
#' @examples 
#' require("MSGARCH")
#'# load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = list(do.init = FALSE))
#' 
#'# compute the Value-at-Risk and Expected-shortfall 
#'# Risk estimation in-sample 
#'risk.its = MSGARCH::risk(object = fit, level = 0.95, ES = FALSE, do.its = TRUE)
#'
#'plot(risk.its)                     
#'
#'# Risk estimation at T + 1                     
#'risk = MSGARCH::risk(object = fit, level = 0.95, ES = FALSE, do.its = FALSE)
#' @importFrom stats integrate sd                    
#' @export
risk <- function(object, theta, y, level = c(0.95, 0.99), ES = TRUE, 
                 do.its = FALSE, ctr = list(n.mesh = 1000)) {
  UseMethod("risk", object)
}

#' @export
risk.MSGARCH_SPEC <- function(object, theta, y, level = c(0.95, 0.99), ES = TRUE,
                              do.its = FALSE, ctr = list(n.mesh = 1000)) {
  
  y   <- MSGARCH:::f.check.y(y)
  out <- list()
  ny  <- nrow(y)
  p       <- 1 - level
  np      <- length(p)
  xmin    <- min(y) - sd(y)
  xmax    <- 0
  
  # gross approximation for VaR
  x <- seq(from = xmin, to = xmax, length.out = ctr$n.mesh)
  pdf_x = MSGARCH::pred(object = object, theta = theta, x = x, y = y, do.its = do.its, log = FALSE)$pred
  cumul <- (apply(pdf_x, 2, cumsum)) * (x[2] - x[1])
  out = list()
  out$VaR <- matrix(NA, nrow = ncol(pdf_x), ncol = np)
  for(n in 1:ncol(pdf_x)){
    for (i in 1:np) {
      out$VaR[n,i] <- x[which.min(abs(cumul[,n] - p[i]))]
    }
  }
  
  if(ES == TRUE) {
    out$ES <- matrix(NA, nrow = ncol(pdf_x), ncol = np)
    for(n in 1:ncol(pdf_x)){
      for (i in 1:np) {
        out$ES[n,i] = sum(pdf_x[x <= out$VaR[n,i],n] * (x[2] - x[1])/p[i] * x[x <=  out$VaR[n,i]])
      }
    }
  }
  out$y <- y
  colnames(out$VaR) <- level
  
  if (isTRUE(ES)) {
    colnames(out$ES) <- level
  }
  class(out) <- "MSGARCH_RISK"
  return(out)
}

#' @export
risk.MSGARCH_MLE_FIT <- function(object, theta = NULL, y = NULL, level = c(0.95, 0.99),
                                 ES = TRUE, do.its = FALSE, ctr = list(n.mesh = 1000)) {
  return(risk.MSGARCH_SPEC(object = object$spec, theta = object$theta, y = object$y,
                           level = level, ES = ES, do.its = do.its, ctr = ctr))
}

#' @export
risk.MSGARCH_BAY_FIT <- function(object, theta = NULL, y = NULL, level = c(0.95, 0.99),
                                 ES = TRUE, do.its = FALSE, ctr = list(n.mesh = 1000)) {
  return(risk.MSGARCH_SPEC(object = object$spec, theta = object$theta, y = object$y,
                           level = level, ES = ES, do.its = do.its, ctr = ctr))
}