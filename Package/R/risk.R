#' Value-at-Risk And Expected-shortfall.
#' @description Method returning the Value-at-Risk and Expected-shortfall in-sample or at \code{t = T + 1} based on the predictive density.
#'@param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not require when using a fit object).
#' @param y  Vector (of size T) of observations (not require when using a fit object).
#' @param level Vector (of size R) of Value-at-risk and Expected-shortfall levels.\cr
#'  (Default: \code{level = c(0.95,0.99)})
#' @param ES  Boolean indicating if Expected-shortfall is also calculated. (Default: \code{ES = TRUE})
#' @param do.its  Boolean indicating if the in-sample risk estimator are returned.
#'  (Default: \code{do.its = FALSE})
#' @param ctr List of control parameters for VaR evaluation.
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian Value-at-Risk and Expected-shortfall are calculated.
#' If \code{do.its = FALSE}, \code{x} the risk estimator at \code{t = T + 1}, the method uses the variance estimated at \code{t = T + 1}.
#' If \code{do.its = TRUE}, The in-sample risk estimator are calculated.
#' @return A list containing of class \code{MSGARCH_RISK} containing two or three components:
#' \itemize{
#' \item \code{VaR} : \cr If \code{do.its = FALSE}: Value-at-Risk at \code{t = T + 1} at the choosen levels (vector of size R).\cr
#'                    If \code{do.its = TRUE}: In-sample Value-at-Risk at the choosen levels (Matrix of size T x R).
#' \item \code{ES}  :\cr If \code{do.its = FALSE}: Expected-shortfall at \code{t = T + 1} at the choosen levels (vector of size R).\cr
#'                    If \code{do.its = TRUE}: In-sample Expected-shortfall at the choosen levels (Matrix of size T x R).
#' \item \code{y}  : Vector (of size T) of observations.
#' }
#' The \code{MSGARCH_RISK} contains the \code{plot} method. 
#' The Bayesian risk estimator can take long time to calculate depending on the size of the MCMC chain.
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
                 do.its = FALSE, ctr = list(n.mesh = 500, tol = 1e-04, itermax = 5)) {
  UseMethod("risk", object)
}

#' @export
risk.MSGARCH_SPEC <- function(object, theta, y, level = c(0.95, 0.99), ES = TRUE,
                              do.its = FALSE, ctr = list(n.mesh = 500, tol = 1e-04, itermax = 5)) {
  if (isTRUE(do.its)) {
    y <- c(mean(y), y)
  }
  y   <- f.check.y(y)
  out <- list()
  ny  <- nrow(y)
  if (!isTRUE(do.its)) {
    start <- ny
    step  <- -(ny - 1)
    end   <- ny
  } else {
    start <- 2
    step  <- -1
    end   <- ny - 1
  }
  p       <- 1 - level
  np      <- length(p)
  xmin    <- min(y) - sd(y)
  xmax    <- 0
  itermax <- ctr$itermax
  tol     <- ctr$tol
  tmp.VaR <- NULL
  out$VaR <- matrix(data = NA, nrow = end - start + 1, ncol = np)
  if (isTRUE(ES)) {
    out$ES <- matrix(data = NA, nrow = end - start + 1, ncol = np)
  }
  for (v in start:end) {
    f.pdf <- function(x) {
      out <- MSGARCH::pred(object, x, theta, y[1:v], log = FALSE, do.its = FALSE)$pred
      return(out)
    }
    # gross approximation for VaR
    x <- seq(from = xmin, to = xmax, length.out = ctr$n.mesh)
    cumul <- cumsum(f.pdf(x) * (x[2] - x[1]))
    for (i in 1:np) {
      tmp.VaR[i] <- x[which.min(abs(cumul - p[i]))]
    }
    # add precision by Newton-Raphson
    for (i in 1:np) {
      p_i <- p[i]
      calc.step <- function(V) {
        lPDF <- MSGARCH::pred(object, x = V, theta = theta, y = y[1:v], log = TRUE, do.its = FALSE)$pred
        CDF  <- MSGARCH::pit(object, x = V, theta = theta, y = y[1:v], do.norm = FALSE, do.its = FALSE)$pit
        err  <- p_i - CDF
        step <- err * exp(-lPDF)
        out  <- list(step = step, err = abs(err))
        return(out)
      }
      # Starting values
      VaR <- tmp.VaR[i]
      # VaR calculation
      newStep <- calc.step(VaR)
      delta   <- newStep$step
      VaR     <- VaR + delta
      covERR  <- newStep$err
      for (j in 1:itermax) {
        if (covERR < tol) {
          break
        }
        newStep <- calc.step(VaR)
        delta   <- newStep$step
        VaR     <- VaR + delta
        covERR  <- newStep$err
      }
      out$VaR[v + step, i] <- VaR
    }
    if (isTRUE(ES)) {
      for (i in 1:np) {
        f.condMean <- function(x) {
          out <- x * MSGARCH::pred(object, x, theta, y[1:v], log = FALSE, do.its = FALSE)$pred
          return(out)
        }
        out$ES[v + step, i] <- integrate(f.condMean, lower = -Inf, upper = out$VaR[v + step, i], 
                                         stop.on.error = FALSE)$value / p[i]
      }
    }
  }
  out$y <- y[2:length(y)]
  colnames(out$VaR) <- level
  if (isTRUE(do.its)) {
    out$VaR <- rbind(NA, out$VaR)
    if (isTRUE(ES)) {
      out$ES <- rbind(NA, out$ES)
    }
  }
  if (isTRUE(ES)) {
    colnames(out$ES) <- level
  }
  class(out) <- "MSGARCH_RISK"
  return(out)
}

#' @export
risk.MSGARCH_MLE_FIT <- function(object, theta = NULL, y = NULL, level = c(0.95, 0.99),
                                   ES = TRUE, do.its = FALSE, ctr = list(n.mesh = 500, tol = 1e-04, itermax = 5)) {
  return(risk.MSGARCH_SPEC(object = object$spec, theta = object$theta, y = object$y,
                             level = level, ES = ES, do.its = do.its, ctr = ctr))
}

#' @export
risk.MSGARCH_BAY_FIT <- function(object, theta = NULL, y = NULL, level = c(0.95, 0.99),
                                   ES = TRUE, do.its = FALSE, ctr = list(n.mesh = 500, tol = 1e-04, itermax = 5)) {
  return(risk.MSGARCH_SPEC(object = object$spec, theta = object$theta, y = object$y,
                             level = level, ES = ES, do.its = do.its, ctr = ctr))
}