#' Value-at-Risk and Expected-shortfall.
#' @description Method returning the Value-at-Risk and Expected-shortfall risk measures.
#' @param object Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}} or fit object of type \code{MSGARCH_ML_FIT}
#' created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT}
#' created with \code{\link{FitMCMC}}.
#' @param par Vector (of size d) or matrix (of size \code{nmcmc} x d) of parameter estimates
#' where d must have
#' the same length as the default parameters of the specification.
#' @param data  Vector (of size T) of observations.
#' @param newdata  Vector (of size T*) of new observations. (Default \code{newdata = NULL})
#' @param alpha Vector (of size R) of Value-at-risk and Expected-shortfall levels.\cr
#'  (Default: \code{alpha = c(0.01, 0.05)})
#' @param do.es  Logical indicating if Expected-shortfall is also calculated.
#' (Default: \code{do.es = TRUE})
#' @param do.its  Logical indicating if the in-sample risk estimators are returned.
#'  (Default: \code{do.its = FALSE}).  
#' @param nahead  Scalar indicating the number of step-ahead evaluation. (Default: \code{nahead = 1L}). Not used when
#' \code{do.its = TRUE} as it only return in-sample one-step ahead risk measures.
#' @param do.cumulative  Logical indicating if cumulative risk measure should be return.
#'  (Default: \code{do.cumulative = FALSE}).  
#' @param ctr A list of control parameters:
#'        \itemize{
#'        \item \code{nmesh} (integer >= 0) : Number of points for density
#'         evaluation. (Default: \code{nmesh = 1000L})
#'        \item \code{nsim} (integer >= 0) :
#'        Number indicating the number of simulation done for estimation of the
#'        density at \code{nahead > 1}. (Default: \code{nsim = 10000L})
#'        }
#' @param ... Not used. Other arguments to \code{Risk}.
#' @return A list of class \code{MSGARCH_RISK} with the following elements:
#' \itemize{
#' \item \code{VaR}:\cr
#' If \code{do.its = FALSE}: Value-at-Risk at \code{t = T + T* + 1, ... ,t = T + T* + nahead} at the
#' chosen levels (matrix of size \code{nahead} x R).\cr
#' If \code{do.its = TRUE}: In-sample Value-at-Risk at the chosen levels (Matrix of size (T + T*) x R).
#' \item \code{ES}:\cr
#' If \code{do.its = FALSE}: Expected-shortfall at \code{t = T + T* + 1, ... ,t = T + T* + nahead} at the
#' chosen levels (matrix of size \code{nahead} x R).\cr
#' If \code{do.its = TRUE}: In-sample Expected-shortfall at the chosen levels (Matrix of size (T + T*) x R).
#' }
#' The \code{MSGARCH_RISK} contains the \code{plot} method.
#' Note that the MCMC/Bayesian risk estimator can take long time to calculate
#' depending on the size of the MCMC chain.
#' @details If a matrix of MCMC posterior draws is given, the
#' Bayesian Value-at-Risk and Expected-shortfall are calculated.
#' Two or more step ahead risk measures are estimated via simulation of \code{nsim} paths up to
#' \code{t = T + T* + nahead}.
#' If \code{do.its = FALSE}, the risk estimators at \code{t = T + T* + 1, ... ,t = T + T* + nahead}
#' are computed. \code{do.cumulative = TRUE} indicate the function to compute the risk meausre 
#' over aggregated period up to \code{nahead} period using the \code{cumsum} function on the simulated data. 
#' @examples
#' # load data
#' data("SMI", package = "MSGARCH")
#'
#' # create model specification
#' # MS(2)-GARCH(1,1)-Normal (default)
#' spec <- CreateSpec()
#'
#' # fit the model on the data with ML estimation
#' fit <- FitML(spec = spec, data = SMI)
#'
#' # compute the Value-at-Risk and Expected-shortfall in-sample
#' risk.its <- Risk(object = fit, alpha = 0.05, do.es = FALSE, do.its = TRUE)
#' plot(risk.its)
#'
#' # compute the one-step ahead Value-at-Risk and Expected-shortfall out-of-sample
#' Risk(object = fit, alpha = 0.05, do.es = FALSE, do.its = FALSE, nahead = 1L)
#' @importFrom stats integrate sd
#' @export
Risk <- function(object, ...) {
  UseMethod(generic = "Risk", object)
}

#' @rdname Risk
#' @export
Risk.MSGARCH_SPEC <- function(object, par, data, alpha = c(0.01, 0.05), nahead = 1L, do.es = TRUE,
                              do.its = FALSE, do.cumulative = FALSE, ctr = list(), ...) {
  
  if (is.vector(par)) {
    par <- matrix(par, nrow = 1L)
  }
  if (nrow(par) == 1) {
    ctr   <- f_process_ctr(ctr)
    nsim <- ctr$nsim
  } else {
    if(is.null(ctr$nsim)){
      nsim = 1
    } else {
      nsim = ctr$nsim
    }
  }
  object  <- f_check_spec(object)
  data    <- f_check_y(data)
  ctr     <- f_process_ctr(ctr)
  out     <- list()
  n.alpha <- length(alpha)
  xmin    <- min(data) - sd(data)
  xmax    <- max(data) + sd(data)
  
  x     <- seq(from = xmin, to = xmax, length.out = ctr$nmesh)
  pdf_x <- PredPdf(object = object, par = par, x = x, data = data, do.its = do.its, log = FALSE)
  cumul <- apply(pdf_x, 1L, cumsum) * (x[2L] - x[1L])
  out   <- list()
  draw  <- NULL
  if (do.its == TRUE) {
    out$VaR <- matrix(NA, nrow = nrow(pdf_x), ncol = n.alpha)
    rownames(out$VaR) <-  paste0("t=",1:length(data))
  } else {
    out$VaR <- matrix(NA, nrow = nahead, ncol = n.alpha)
    rownames(out$VaR) <-  paste0("h=",1:nahead)
  }
  for (n in 1:nrow(pdf_x)) {
    for (i in 1:n.alpha) {
      out$VaR[n, i] <- x[which.min(abs(cumul[, n] - alpha[i]))]
    }
  }
  
  if (nahead > 1 & do.its == FALSE) {
    draw <- Sim(object = object, data = data, nahead = nahead, nsim = nsim, par = par)$draw
    if(isTRUE(do.cumulative)){
      draw = apply(draw, 2, cumsum)
    }
    for (j in 2:nahead) {
      out$VaR[j, ] <- quantile(draw[j,], probs = alpha)
    }
  }
  
  if (isTRUE(do.es)) {
    if (do.its == TRUE) {
      out$ES <- matrix(NA, nrow = nrow(pdf_x), ncol = n.alpha)
      rownames(out$ES) <- paste0("t=",1:length(data))
    } else {
      out$ES <- matrix(NA, nrow = nahead, ncol = n.alpha)
      rownames(out$ES) <- paste0("h=",1:nahead)
    }
    for (n in 1:nrow(pdf_x)) {
      for (i in 1:n.alpha) {
        out$ES[n, i] <- sum(pdf_x[n, x <= out$VaR[n, i]] * (x[2L] - x[1L])/alpha[i] * x[x <= out$VaR[n, i]])
      }
    }
    
    if (nahead > 1 & do.its == FALSE) {
      for (i in 1:n.alpha) {
        for (j in 2:nahead) {
          out$ES[j, i] <- mean(draw[j, draw[j, ] <= out$VaR[j, i]])
        }
      }
    }
  }
  colnames(out$VaR) <- alpha
  
  if (isTRUE(do.es)) {
    colnames(out$ES) <- alpha
  }
  class(out) <- "MSGARCH_RISK"
  return(out)
}

#' @rdname Risk
#' @export
Risk.MSGARCH_ML_FIT <- function(object, newdata = NULL, alpha = c(0.01, 0.05),
                                do.es = TRUE, do.its = FALSE, nahead = 1L, ctr = list(), ...) {
  data <- c(object$data, newdata)
  out  <- Risk(object = object$spec, par = object$par, data = data, alpha = alpha,
               do.es = do.es, do.its = do.its, nahead = nahead, ctr = ctr)
  return(out)
}

#' @rdname Risk
#' @export
Risk.MSGARCH_MCMC_FIT <- function(object, newdata = NULL, alpha = c(0.01, 0.05),
                                  do.es = TRUE, do.its = FALSE, nahead = 1L, ctr = list(), ...) {
  data <- c(object$data, newdata)
  out  <- Risk(object = object$spec, par = object$par, data = data, alpha = alpha,
               do.es = do.es, do.its = do.its, nahead = nahead, ctr = ctr)
  return(out)
}
