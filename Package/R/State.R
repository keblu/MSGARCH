#' @title State probabilities.
#' @description Method returning the filtered, predictive, and smoothed probabilities of the states,
#' and the most probable path computed with the Viterbi algorithm.
#' @param object Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}} or fit object of type \code{MSGARCH_ML_FIT}
#' created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT}
#' created with \code{\link{FitMCMC}}.
#' @param par Vector (of size d) or matrix (of size \code{nmcmc} x d) of parameter estimates
#' where d must have
#' the same length as the default parameters of the specification.
#' @param data  Vector (of size T) of observations.
#' @param newdata  Vector (of size T*) of new observations. (Default \code{newdata = NULL})
#' @param ... Not used. Other arguments to \code{State}.
#' @return A list of class \code{MSGARCH_PSTATE} with the following elements:
#' \itemize{
#' \item \code{FiltProb}: Filtered probabilities (array of size (T + T*) x (\code{nmcmc or 1}) x K).
#' \item \code{PredProb}: Predictive probabilities (array of size (T + T* + 1) x (\code{nmcmc or 1}) x K).
#' \item \code{SmoothProb}: Smoothed probabilities (array of size (T + T* + 1) x (\code{nmcmc or 1}) x K).
#' \item \code{Viterbi}:  Most likely path (matrix of size (T + T*) x (\code{nmcmc} or 1)).
#' }
#' The class \code{MSGARCH_PSTATE} contains the \code{plot} method. The plot method contains
#' as input \code{type.prob} which is one of \code{"filtered", "predictive", "smoothed", "viterbi"}.
#' (Default: \code{type.prob = "smoothed"})
#' @details If a matrix of parameter estimates is given, each parameter
#' estimate (each row) is evaluated individually.
#' @examples
#' # create specification
#' spec <- CreateSpec()
#' 
#' # load data
#' data("SMI", package = "MSGARCH")
#' 
#' # state from specification
#' par <- c(0.1, 0.1, 0.8, 0.2, 0.1, 0.8, 0.99, 0.01)
#' state <- State(object = spec, par = par, data = SMI)
#' plot(state, type.prob = "filtered")
#' 
#' # state from ML fit
#' fit <- FitML(spec = spec, data = SMI)
#' state <- State(object = fit)
#' plot(state, type.prob = "smoothed")
#' 
#' \dontrun{
#' # state from MCMC fit
#' set.seed(1234)
#' fit <- FitMCMC(spec = spec, data = SMI)
#' state <- State(object = fit)
#' plot(state, type.prob = "smoothed")
#' }
#' @export
State <- function(object, ...) {
  UseMethod(generic = "State", object)
}

#' @rdname State
#' @export
State.MSGARCH_SPEC <- function(object, par, data, ...) {
  object <- f_check_spec(object)
  par    <- f_check_par(object, par)
  y      <- as.matrix(data)

  if(zoo::is.zoo(data)|| is.ts(data)){
    date = zoo::index(data)
  } else {
    date = 1:nrow(y)
  }
  
  FiltProb <- array(dim = c(nrow(y), nrow(par), object$K),
                    dimnames = list(as.character(date),
                                    paste0("draw #",1:nrow(par)), paste0("k=",1:object$K)))
  PredProb <- array(dim = c(nrow(y) + 1L, nrow(par), object$K),
                    dimnames = list(as.character(c(date, date[length(date)] + 1)),
                                    paste0("draw #",1:nrow(par)), paste0("k=",1:object$K)))
  SmoothProb <- array(dim = c(nrow(y) + 1L, nrow(par), object$K),
                      dimnames = list(as.character(c(date, date[length(date)] + 1)),
                                      paste0("draw #",1:nrow(par)), paste0("k=",1:object$K)))
  
  viterbi <- matrix(data = NA, nrow = nrow(y), ncol = nrow(par),
                    dimnames = list(as.character(date), paste0("draw #",1:nrow(par))))
  
  out <- list(FiltProb = FiltProb, PredProb = PredProb, SmoothProb = SmoothProb, Viterbi = viterbi)
  for (i in 1:nrow(par)) {
    tmp <- object$rcpp.func$get_Pstate_Rcpp(par[i, ], y)
    for (j in 1:object$K) {
      out$FiltProb[, i, j] <- tmp$FiltProb[, j]
      out$PredProb[, i, j] <- tmp$PredProb[, j]
      out$SmoothProb[, i, j] <- tmp$SmoothProb[, j]
    }
    if (object$K > 1) {
      P <- TransMat(object = object, par = par[i, ], nahead = 1)
      if (isTRUE(object$is.mix)) {
        P <- matrix(rep(P, object$K), nrow = object$K, ncol = object$K)
      }
      out$Viterbi[2:length(data), i] <- Viterbi(tmp$LL, P, object$K)
    } else {
      out$Viterbi[2:length(data), i] <- tmp$Viterbi
    }
  }
  out$Viterbi[1,] = out$Viterbi[2, ]
  out$Viterbi =  out$Viterbi + 1
  # missing first value because LL start at time 2 (first observation not included in likelihood
  # evaluation). The + 1 is for starting the state at 1 and not 0
  class(out) <- "MSGARCH_PSTATE"
  return(out)
}

#' @rdname State
#' @export
State.MSGARCH_ML_FIT <- function(object, newdata = NULL, ...) {
  data <- c(object$data, newdata)
  out  <- State(object = object$spec, par = object$par, data = data)
  return(out)
}

#' @rdname State
#' @export
State.MSGARCH_MCMC_FIT <- function(object, newdata = NULL, ...) {
  data <- c(object$data, newdata)
  out  <- State(object = object$spec, par = object$par, data = data)
  return(out)
}
