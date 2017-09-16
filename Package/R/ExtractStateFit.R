#' @title Single-regime model extractor.
#' @description Extracts each regime from a fitted multiple regime specificaton
#' and creates a fitted object for each extracted regime.
#' @param object Fit object of type \code{MSGARCH_ML_FIT}
#' created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT}
#' created with \code{\link{FitMCMC}}.
#' @return A list of with \code{K} element where each element is a fit object of type \code{MSGARCH_ML_FIT} or
#' \code{MSGARCH_MCMC_FIT}.
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
#' SR.fit <- ExtractStateFit(fit)
#' @export
ExtractStateFit <- function(object) {
  UseMethod(generic = "ExtractStateFit", object)
}

f_ExtractStateFit <- function(object) {
  K   <- length(object$spec$name)
  out <- vector(mode = "list", length = object$spec$K)
  pos <- c(0,cumsum(object$spec$n.params))
  par <- object$par
  is.vec <- is.vector(par)
  if (isTRUE(is.vec)) {
    par <- t(as.matrix(par))
  }

  for (i in 1:K) {
    out[[i]]$spec <- f_spec(object$spec$name[i], do.mix = FALSE)
    class(out[[i]]$spec) <- "MSGARCH_SPEC"
    out[[i]]$par <- par[, (pos[i] + 1):pos[i + 1]]
    if (is.vec) {
      out[[i]]$par <- as.vector(out[[i]]$par)
      names(out[[i]]$par) <- out[[i]]$spec$label
    } else {
      colnames(out[[i]]$par) <- out[[i]]$spec$label
    }
    out[[i]]$data <- object$data
    out[[i]]$ctr <- object$ctr
    out[[i]]$loglik <- Kernel(object = out[[i]]$spec, par = out[[i]]$par, data = out[[i]]$data, log = TRUE, do.prior = FALSE)
  }
  return(out)
}

#' @rdname ExtractStateFit
#' @export
ExtractStateFit.MSGARCH_ML_FIT <- function(object) {
  out <- f_ExtractStateFit(object)
  K   <- length(object$spec$name)
  for (i in 1:K) {
    class(out[[i]]) <- "MSGARCH_ML_FIT"
  }
  return(out)
}

#' @rdname ExtractStateFit
#' @export
ExtractStateFit.MSGARCH_MCMC_FIT <- function(object) {
  out <- f_ExtractStateFit(object)
  K   <- length(object$spec$name)
  for (i in 1:K) {
    class(out[[i]]) <- "MSGARCH_MCMC_FIT"
  }
  return(out)
}
