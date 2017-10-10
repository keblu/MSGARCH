#' @title Transition matrix.
#' @description Method returning the transition matrix.
#' @param object Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}}
#' or fit object of type \code{MSGARCH_ML_FIT} created with \code{\link{FitML}}.
#' @param par Vector (of size d) of
#' parameter estimates (not required when using a fit object) where d must have
#' the same length as the default parameters of the specification.
#' @param n.ahead Number of steps ahead. (Default: \code{n.ahead = 1L})
#' @param ... Not used. Other arguments to \code{TransMat}.
#' @return A matrix (of size K x K) in the case of a Markov-Switching model
#'  or a vector (of size K) in the case of a Mixture of GARCH model.
#'  The row indicates the starting states while the columns indicates the transition states.
#' @examples
#' # load data
#' data("SMI", package = "MSGARCH")
#'
#' # create model specification
#' # MS(2)-GARCH(1,1)-Normal (default)
#' spec <- CreateSpec()
#'
#' # fit the model on the data by ML
#' fit <- FitML(spec = spec, data = SMI)
#'
#' # Extract the transition matrix 10 steps ahead
#' trans.mat <- TransMat(fit, n.ahead = 10)
#' print(trans.mat)
#' @importFrom stats quantile
#' @import expm
#' @export
TransMat <- function(object, ...) {
  UseMethod(generic = "TransMat", object = object)
}

#' @rdname TransMat
#' @export
TransMat.MSGARCH_SPEC <- function(object, par = NULL, n.ahead = 1L, ...) {
  object <- f_check_spec(object)
  if (isTRUE(object$is.shape.ind)) {
    par <- object$func$f.do.shape.ind(par = par)
  }
  n.params   <- object$n.params
  n.model    <- length(n.params)
  params.loc <- c(0, cumsum(n.params))
  if (!isTRUE(object$is.mix)) {
    p <- matrix(nrow = n.model, ncol = n.model)
    for (i in 0:(n.model - 1L)) {
      p[1:(n.model - 1L), i + 1L] <-
        par[(params.loc[n.model + 1L] + n.model * i + 1L - i):(params.loc[n.model + 1L] + n.model * i + n.model - 1L - i)]
    }
    p[n.model, ] <- 1 - colSums(matrix(p[1:(n.model - 1L), ], ncol = n.model))
    p <- t(p)
  } else {
    p <- matrix(rep(0, n.model), ncol = n.model)
    for (i in 1:(n.model - 1L)) {
      p[1L, i] <- par[(params.loc[n.model + 1L] + i)]
    }
    p[1L, n.model] <- 1 - sum(p)
  }
  if (!object$is.mix) {
    p <- p %^% n.ahead
  }
  if (object$is.mix) {
    col_label <- paste0("State ", 1:object$K)
    row_label <- paste0("Probability")
  } else {
    col_label <- paste0("t+", n.ahead, "|k=", 1:object$K)
    row_label <- paste0("t|k=", 1:object$K)
  }
  rownames(p) <- row_label
  colnames(p) <- col_label
  return(p)
}

#' @rdname TransMat
#' @export
TransMat.MSGARCH_ML_FIT <- function(object, n.ahead = 1L, ...) {
  out <- TransMat(object = object$spec, par = object$par, n.ahead = n.ahead)
  return(out)
}