#' ML estimation.
#' @description Method that performs Maximum Likelihood estimation of a \code{\link{MSGARCH_SPEC}} object on a set of observations.
#' @param y Vector (of size T) of observations.
#' @param spec Model specification created with \code{\link{create.spec}}.
#' @param ctr List of control parameters.
#'        The control parameters have two components to it:
#'        \itemize{
#'        \item \code{theta0} : Starting parameters (vector of size d). If no starting parameters is provided, the default starting parameters of the specification are used.
#'        \item \code{do.init} : Boolean indicating if there is a pre-optimization with the \R package \code{DEoptim} (Ardia et al., 2011). (default: \code{do.init = FALSE})
#'        }
#' @return A list of class \code{\link{MSGARCH_MLE_FIT}} containing five variables:
#'        \itemize{
#'        \item \code{theta} : Optimal parameters (vector of size d).
#'        \item \code{l_likelihood} : log-likelihood of y given the optimal parameters.
#'        \item \code{spec} : Initial specification.
#'        \item \code{is.init} : Indicating if estimation was made with do.init option.
#'        \item \code{y} : Initial vector (of size T) of observations..
#'        }
#' The \code{\link{MSGARCH_MLE_FIT}} possess these methods:
#' \itemize{
#' \item \code{\link{AIC}} : Compute Akaike information criterion (AIC).
#' \item \code{\link{BIC}} : Compute Bayesian information criterion (BIC).
#' }
#' 
#' @details The Maximum likelihood estimation uses the \R package \code{nloptr} (Johnson, 2014) for main optimizer 
#' while it uses the \R package \code{DEoptim} when \code{do.init = TRUE}.
#' @references Ardia, D.; Mullen, K. M.; Peterson, B. G. & Ulrich, J. (2015). \code{DEoptim}: Differential Evolution in \R. \url{https://cran.r-project.org/web/packages/DEoptim/}.
#' @references Johnson, S. G. (2014). The NLopt Nonlinear-Optimization. \url{https://cran.r-project.org/web/packages/NLopt/}.
#' @examples 
#' data("sp500ret")
#' 
#' spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                               do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#' 
#' fit = MSGARCH::fit.mle(spec = spec, y = sp500ret, ctr = list(do.init = TRUE))
#' @import DEoptim nloptr Rsolnp
#' @export
fit.mle <- function(spec, y, ctr = list())
{
  UseMethod("fit.mle", spec)
}

#' @export
fit.mle.MSGARCH_SPEC = function(spec, y, ctr = list()) {
  
  require("DEoptim")
  ctr.optim = list(trace = 0, maxit = 50000)
  ctr.deoptim = DEoptim::DEoptim.control(NP = 50 * length(spec$theta0), itermax = 500, 
    trace = FALSE, initialpop = matrix(spec$theta0, nrow = 50 * length(spec$theta0), 
      ncol = length(spec$theta0)))
  ctr.slsqp = list(maxeval = 10000, xtol_rel = 1e-08)
  
  ctr = f.process.ctr(ctr)
  
  f.kernel = function(x, log = TRUE) {
    return(MSGARCH::kernel(spec, x, y = y, log = log))
  }
  
  lower = spec$lower
  upper = spec$upper
  
  
  f.nll = function(x) -f.kernel(x, log = TRUE)
  
  if (any(ctr$do.init || spec$do.init)) {
    str = "f.find.theta0 -> DEoptim initialization"
    is.ok = tryCatch({
      tmp = DEoptim::DEoptim(fn = f.nll, lower = lower, upper = upper, control = ctr.deoptim)
      theta0.init = tmp$optim$bestmem
      TRUE
    }, warning = function(warn) {
      f.error(str)
    }, error = function(err) {
     f.error(str)
    })
  } else {
    theta0.init = spec$theta0
  }
  
  theta = f.find.theta0(f.kernel, theta0 = theta0.init, lower = lower, 
    upper = upper, f.ineq = spec$f.ineq, ineqlb = spec$ineqlb, inequb = spec$inequb)
  l_likelihood = f.kernel(theta)
  
  if (l_likelihood == -1e+10) {
    tmp = DEoptim::DEoptim(fn = f.nll, lower = lower, upper = upper, control = ctr.deoptim)
    theta = tmp$optim$bestmem
    l_likelihood = f.kernel(theta)
  }
  
  out = list(theta = theta, l_likelihood = l_likelihood, spec = spec, is.init = any(ctr$do.init || spec$do.init), y = y)
  class(out) = "MSGARCH_MLE_FIT"
  return(out)
}

