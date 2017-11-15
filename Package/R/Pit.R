#' @title Probability integral transform.
#' @description Method returning the probability integral
#' transform (PIT).
#' @param object Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}} or fit object of type \code{MSGARCH_ML_FIT}
#' created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}.
#' @param x  Vector (of size n). Used when \code{do.its = FALSE}.
#' @param par Vector (of size d) or matrix (of size \code{nmcmc} x d) of
#' parameter estimates where d must have
#' the same length as the default parameters of the specification.
#' @param data  Vector (of size T) of observations.
#' @param newdata  Vector (of size T*) of new observations. (Default: \code{newdata = NULL})
#' @param do.norm  Logical indicating if the PIT values are transformed
#' into standard Normal variate. (Default: \code{do.norm = FALSE})
#' @param do.its  Logical indicating if the in-sample PIT is returned. (Default: \code{do.its = FALSE})
#' @param nahead  Scalar indicating the number of step-ahead evaluation.
#' Valid only when \code{do.its = FALSE}. (Default: \code{nahead = 1L})
#' @param do.cumulative logical indicating if PIT is computed on the cumulative simulations (typically log-returns, as they can be aggregated).
#'  Only available for \code{do.its = FALSE}. (Default: \code{do.cumulative = FALSE})
#' @param ctr A list of control parameters:
#'        \itemize{
#'        \item \code{nsim} (integer >= 0):
#'        Number indicating the number of simulation done for the
#'        evaluation of the PIT at \code{nahead > 1}. (Default: \code{nsim = 10000L})
#'        }
#' @param ... Not used. Other arguments to \code{PIT}.
#' @return A vector or matrix of class \code{MSGARCH_PIT}. \cr
#' If \code{do.its = FALSE}: Probability integral transform of the
#' points \code{x} at \cr \code{t = T + T* + 1, ... ,t = T + T* + nahead} or Normal variate derived from the probability
#' integral transform of \code{x} (matrix of size \code{nahead} x n).\cr
#' If \code{do.its = TRUE}: In-sample  probability integral transform or Normal variate
#' derived from the probability integral transform of \code{data} if \code{x = NULL} (vector of
#' size T + T*) or in-sample  probability integral transform or Normal variate
#' derived from the probability integral transform of \code{x} (matrix of size   (T + T*) x n).
#' @details If a matrix of MCMC posterior draws is given, the
#' Bayesian probability integral transform is calculated.
#' Two or more step-ahead probability integral
#' transform are estimated via simulation of \code{nsim} paths up to \code{t = T + T* + nahead}.
#' The empirical probability integral transforms is then inferred from these simulations.\cr
#' If \code{do.its = FALSE}, the vector \code{x} are evaluated as  \code{t = T + T* + 1, ... ,t = T + T* + nahead}
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
#' sim.series <- simulate(object = spec, par = fit$par, nahead= 1000L, nsim = 1L)
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
                             do.norm = FALSE, do.its = FALSE, nahead = 1L, do.cumulative = FALSE, ctr = list(), ...) {
  object <- f_check_spec(object)
  data_   <- f_check_y(data)
  if (is.vector(par)) {
    par <- matrix(par, nrow = 1L)
  }
  if (nrow(par) == 1) {
    ctr     <- f_process_ctr(ctr)
    nsim <- ctr$nsim
  } else {
    if(is.null(ctr$nsim)){
      nsim = 1
    } else {
      nsim = ctr$nsim
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
      x <- matrix(data = data_, ncol = length(data_))
    } else {
      x <- matrix(x)
      if (ncol(x) == 1L) {
        x <- matrix(x, ncol = length(data_), nrow = nrow(x))
      } else {
        stop("x have more than 1 column: x must be a vector, NULL, or a matrix of size n x 1")
      }
    }
    tmp <- matrix(data = 0, nrow = length(data_), ncol = nrow(x))
    for (i in 1:nrow(par)) {
      if (object$K == 1) {
        tmp2 <- object$rcpp.func$cdf_Rcpp_its(par_check[i, ], data_, x, FALSE)
        tmp <- tmp + tmp2[, , 1L]
      } else {
        Pstate <- State(object = object, par = par[i, ], data = data_)$PredProb
        Pstate.tmp <- matrix(data = NA, nrow = dim(Pstate)[1L], ncol = dim(Pstate)[3L])
        for (j in 1:dim(Pstate)[3L]) {
          Pstate.tmp[, j] <- Pstate[, , j]
        }
        tmp2 <- object$rcpp.func$cdf_Rcpp_its(par_check[i, ], data_, x, FALSE)
        for (k in 1:object$K) {
          tmp <- tmp + tmp2[, , k] * matrix(Pstate.tmp[1:(nrow(Pstate.tmp) - 1L), k], 
                                            ncol = nrow(x), nrow = length(data))
        }
      }
    }
    tmp <- tmp/nrow(par)
    rownames(tmp) =  paste0("t=",1:length(data_))
    if(zoo::is.zoo(data)){
      tmp = zoo::zooreg(tmp, order.by = zoo::index(data))
    }
    if(is.ts(data)){
      tmp = zoo::zooreg(tmp, order.by = zoo::index(data))
      tmp = as.ts(tmp)
      colnames(tmp) = rep("",ncol(tmp)) 
    }
  } else {
    x <- matrix(x)
    if (ncol(x) != 1L) {
      stop("x must be a vector or a matrix of size N x 1")
    }
    tmp <- matrix(data = 0, nrow = nahead, ncol = nrow(x))
    for (i in 1:nrow(par)) {
      tmp[1,] <- tmp[1, ] + object$rcpp.func$cdf_Rcpp(x, par_check[i, ], data_, FALSE)
    }
    tmp <- tmp/nrow(par)
    if (nahead > 1) {
      draw <- Sim(object = object, data = data_, nahead = nahead, nsim = nsim, par = par)$draw
      if(isTRUE(do.cumulative)){
        draw = apply(draw, 2, cumsum)
      }
      for (j in 2:nahead) {
        tmp[j, ] <- f_cdf_empirical(y = draw[j, ], x)
      }
    }
    rownames(tmp) <- paste0("h=",1:nahead)
    if(zoo::is.zoo(data)){
      tmp = zoo::zooreg(tmp, order.by = zoo::index(data)[length(data)]+(1:nahead))
    }
    if(is.ts(data)){
      tmp = zoo::zooreg(tmp, order.by = zoo::index(data)[length(data)]+(1:nahead))
      tmp = as.ts(tmp)
      colnames(tmp) = rep("",ncol(tmp)) 
    }
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
  if(isTRUE(x.is.null) && !is.ts(data)){
    out <- tmp[,1]
  } else {
    out = tmp
  }
  class(out) <- c("MSGARCH_PIT", class(out))
  return(out)
}

#' @rdname PIT
#' @export
PIT.MSGARCH_ML_FIT <- function(object, x = NULL, newdata = NULL,
                               do.norm = TRUE, do.its = FALSE, nahead = 1L, ctr = list(), ...) {
  data = c(object$data, newdata)
  if(is.ts(object$data)){
    if(is.null(newdata)){
      data = zoo::zooreg(data, order.by =  c(zoo::index(data)))
    } else {
      data = zoo::zooreg(data, order.by =  c(zoo::index(data),zoo::index(data)[length(data)]+(1:length(newdata))))
    }
    data = as.ts(data)
  }
  out <- PIT(object = object$spec, x = x, par = object$par, data = data,
             do.norm = do.norm, do.its = do.its, nahead = nahead, ctr = ctr)
  return(out)
}

#' @rdname PIT
#' @export
PIT.MSGARCH_MCMC_FIT <- function(object, x = NULL, newdata = NULL,
                                 do.norm = TRUE, do.its = FALSE, nahead = 1L, ctr = list(), ...) {
  data = c(object$data, newdata)
  if(is.ts(object$data)){
    if(is.null(newdata)){
      data = zoo::zooreg(data, order.by =  c(zoo::index(data)))
    } else {
      data = zoo::zooreg(data, order.by =  c(zoo::index(data),zoo::index(data)[length(data)]+(1:length(newdata))))
    }
    data = as.ts(data)
  }
  out <- PIT(object = object$spec, x = x, par = object$par, data = data,
             do.norm = do.norm, do.its = do.its, nahead = nahead, ctr = ctr)
  return(out)
}
