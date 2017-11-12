#' @title State probabilities.
#' @description Method returning the filtered, predictive, smoothed probabilities of the states,
#' and the most probable path computed with the Viterbi alogirthm.
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
#' # compute the filtered state probabilities
#' state <- State(object = fit)
#' plot(state, type.prob = "smoothed")
#' plot(state, type.prob = "predictive")
#' plot(state, type.prob = "filtered")
#' plot(state, type.prob = "viterbi")
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
  FiltProb <- array(dim = c(nrow(y), nrow(par), object$K),
                    dimnames = list(paste0("t=",1:(nrow(y))),
                                    paste0("draw #",1:nrow(par)), paste0("k=",1:object$K)))
  PredProb <- array(dim = c(nrow(y) + 1L, nrow(par), object$K),
                    dimnames = list(paste0("t=",1:(nrow(y)+1)),
                                    paste0("draw #",1:nrow(par)), paste0("k=",1:object$K)))
  SmoothProb <- array(dim = c(nrow(y) + 1L, nrow(par), object$K),
                      dimnames = list(paste0("t=",1:(nrow(y)+1)),
                                      paste0("draw #",1:nrow(par)), paste0("k=",1:object$K)))
  
  viterbi <- matrix(data = NA, nrow = nrow(y), ncol = nrow(par),
                    dimnames = list(paste0("t=",1:nrow(y)), paste0("draw #",1:nrow(par))))
  
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
