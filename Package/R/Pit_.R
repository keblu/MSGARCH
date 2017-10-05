#' @title Probability integral transform.
#' @description Method returning the probability integral
#' transform (PIT).
#' @param object Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}} or fit object of type \code{MSGARCH_ML_FIT}
#' created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}.
#' @param x  Vector (of size n). Used when \code{do.its = FALSE}.
#' @param par Vector (of size d) or matrix (of size \code{n.mcmc} x d) of
#' parameter estimates where d must have
#' the same length as the default parameters of the specification.
#' @param data  Vector (of size T) of observations.
#' @param new.data  Vector (of size T*) of new observations. (Default \code{new.data = NULL})
#' @param do.norm  Logical indicating if the PIT values are transformed
#' into standard Normal variate. (Default: \code{do.norm = FALSE})
#' @param do.its  Logical indicating if the in-sample PIT is returned. (Default: \code{do.its = FALSE})
#' @param n.ahead  Scalar indicating the number of step-ahead evaluation.
#'  Valid only when \code{do.its = FALSE}. (Default: \code{n.ahead = 1L})
#' @param ctr A list of control parameters:
#'        \itemize{
#'        \item \code{n.sim} (integer >= 0):
#'        Number indicating the number of simulation done for the
#'        evaluation of the PIT at \code{n.ahead > 1}. (Default: \code{n.sim = 10000L})
#'        }
#' @param ... Not used. Other arguments to \code{PIT}.
#' @return A vector or matrix of class \code{MSGARCH_PIT}. \cr
#' If \code{do.its = FALSE}: Probability integral transform of the
#' points \code{x} at \cr \code{t = T + T* + 1, ... ,t = T + T* + n.ahead} or Normal variate derived from the probability
#' integral transform of \code{x} (matrix of size \code{n.ahead} x n).\cr
#' If \code{do.its = TRUE}: In-sample  probability integral transform or Normal variate
#' derived from the probability integral transform of \code{data} if \code{x = NULL} (vector of
#' size T + T*) or in-sample  probability integral transform or Normal variate
#' derived from the probability integral transform of \code{x} (matrix of size   (T + T*) x n).
#' @details If a matrix of MCMC posterior draws is given, the
#' Bayesian probability integral transform is calculated.
#' Two or more step-ahead probability integral
#' transform are estimated via simulation of \code{n.sim} paths up to \code{t = T + T* + n.ahead}.
#' The empirical probability integral transforms is then inferred from these simulations.\cr
#' If \code{do.its = FALSE}, the vector \code{x} are evaluated as  \code{t = T + T* + 1, ... ,t = T + T* + n.ahead}
#' realizations.\cr
#' If \code{do.its = TRUE}, \code{x} is evaluated
#' at each time \code{t} up to time \code{t = T + T*}.\cr
#' Finally if \code{x = NULL} the vector \code{data} is evaluated for sample evaluation of the PIT.\cr
#' The \code{do.norm} argument transforms the PIT value into Normal variates so that normality test can be done.
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
#' # run PIT method in-sample
#' pit.its <- PIT(object = fit, do.norm = TRUE, do.its = TRUE)
#'
#' # diagnostic of PIT with qqnorm
#' qqnorm(pit.its)
#' qqline(pit.its)
#'
#' # simulate a serie from the model
#' set.seed(123)
#' sim.series <- Sim(object = spec, par = fit$par, n.ahead= 1000L, n.sim = 1L)
#' sim.series <- as.vector(sim.series$draw)
#'
#' # run PIT method on the simualed serie with the true par
#' pit.x <- PIT(object = spec, par = fit$par, data = sim.series, do.norm = TRUE, do.its = TRUE)
#' qqnorm(pit.x)
#' qqline(pit.x)
#' @importFrom stats qnorm
#' @export
PIT <- function(object, ...) {
  UseMethod(generic = "PIT", object)
}

#' @rdname PIT
#' @export
PIT.MSGARCH_SPEC <- function(object, x = NULL, par = NULL, data = NULL,
                             do.norm = FALSE, do.its = FALSE, n.ahead = 1L, ctr = list(), ...) {
  object <- f_check_spec(object)
  data   <- f_check_y(data)
  if (is.vector(par)) {
    par <- matrix(par, nrow = 1L)
  }
  if (nrow(par) == 1) {
    ctr     <- f_process_ctr(ctr)
    n.sim <- ctr$n.sim
  } else {
    if(is.null(ctr$n.sim)){
      n.sim = 1
    } else {
      n.sim = ctr$n.sim
    }
  }
  ctr    <- f_process_ctr(ctr)
  x.is.null  <-  FALSE
  if (is.null(x)) {
    x.is.null <- TRUE
  } 
  draw <- NULL
  par_check <- f_check_par(object, par)
  if (isTRUE(do.its)) {
    if (is.null(x)) {
      x <- matrix(data = data, ncol = length(data))
    } else {
      x <- matrix(x)
      if (ncol(x) == 1L) {
        x <- matrix(x, ncol = length(data), nrow = nrow(x))
      } else {
        stop("x have more than 1 column: x must be a vector, NULL, or a matrix of size n x 1")
      }
    }
    tmp <- matrix(data = 0, nrow = nrow(x), ncol = length(data))
    for (i in 1:nrow(par)) {
      if (object$K == 1) {
        tmp2 <- object$rcpp.func$cdf_Rcpp_its(par_check[i, ], data, x, FALSE)
        tmp <- tmp + tmp2[, , 1L]
      } else {
        Pstate <- State(object = object, par = par[i, ], data = data)$PredProb
        Pstate.tmp <- matrix(data = NA, nrow = dim(Pstate)[1L], ncol = dim(Pstate)[3L])
        for (j in 1:dim(Pstate)[3L]) {
          Pstate.tmp[, j] <- Pstate[, , j]
        }
        tmp2 <- object$rcpp.func$cdf_Rcpp_its(par_check[i, ], data, x, FALSE)
        for (k in 1:object$K) {
          tmp <- tmp + tmp2[, , k] * matrix(Pstate.tmp[1:(nrow(Pstate.tmp) - 1L), k],
                                            ncol = length(data), nrow = nrow(x), byrow = TRUE)
        }
      }
    }
    tmp <- tmp/nrow(par)
    colnames(tmp) =  paste0("t=",1:length(data))
  } else {
    x <- matrix(x)
    if (ncol(x) != 1L) {
      stop("x must be a vector or a matrix of size N x 1")
    }
    tmp <- matrix(data = 0, nrow = nrow(x), ncol = n.ahead)
    for (i in 1:nrow(par)) {
      tmp[, 1] <- tmp[, 1] + object$rcpp.func$cdf_Rcpp(x, par_check[i, ], data, FALSE)
    }
    tmp <- tmp/nrow(par)
    if (n.ahead > 1) {
      draw <- Sim(object = object, data = data, n.ahead = n.ahead, n.sim = n.sim, par = par)$draw
      for (j in 2:n.ahead) {
        tmp[, j] <- f_cdf_empirical(y = draw[j, ], x)
      }
    }
    colnames(tmp) <- paste0("h=",1:n.ahead)
  }
  if (!isTRUE(ctr$do.return.draw)) {
    draw <- NULL
  }
  if (do.norm) {
    tmp2 <- stats::qnorm(tmp, mean = 0, sd = 1)
    colnames(tmp2) <-  colnames(tmp)
    rownames(tmp2) <- rownames(tmp)
    tmp <- tmp2
  }
  if(isTRUE(x.is.null)){
    out <- tmp[1,]
  } else {
    out <- t(tmp)
  }
  class(out) <- "MSGARCH_PIT"
  return(out)
}

#' @rdname PIT
#' @export
PIT.MSGARCH_ML_FIT <- function(object, x = NULL, new.data = NULL,
                               do.norm = TRUE, do.its = FALSE, n.ahead = 1L, ctr = list(), ...) {
  data = c(object$data, new.data)
  out <- PIT(object = object$spec, x = x, par = object$par, data = data,
             do.norm = do.norm, do.its = do.its, n.ahead = n.ahead, ctr = ctr)
  return(out)
}

#' @rdname PIT
#' @export
PIT.MSGARCH_MCMC_FIT <- function(object, x = NULL, new.data = NULL,
                                 do.norm = TRUE, do.its = FALSE, n.ahead = 1L, ctr = list(), ...) {
  data = c(object$data, new.data)
  out <- PIT(object = object$spec, x = x, par = object$par, data = data,
             do.norm = do.norm, do.its = do.its, n.ahead = n.ahead, ctr = ctr)
  return(out)
}
