#' ML estimation.
#' @description Method that performs Maximum Likelihood estimation of a \code{MSGARCH_SPEC} object on a set of observations.
#' @param y Vector (of size T) of observations.
#' @param spec Model specification created with \code{\link{create.spec}}.
#' @param ctr List of control parameters.
#'        The control parameters have two components to it:
#'        \itemize{
#'        \item \code{theta0} : Starting parameters (vector of size d). If no starting parameters is provided, the default starting parameters of the specification are used.
#'        \item \code{do.init} : Boolean indicating if there is a pre-optimization with the \R package \code{DEoptim} (Ardia et al., 2011). (Default: \code{do.init = FALSE})
#'        \item \code{NP} : Number of parameter vectors in the population in \code{DEoptim} optimization. (Default: \code{NP = 500})
#'        \item \code{itermax} : Maximum iteration (population generation) allowed in \code{DEoptim} optimization. (Default: \code{maxit = 500})
#'        }
#' @return A list of class \code{MSGARCH_MLE_FIT} containing five components:
#'        \itemize{
#'        \item \code{theta} : Optimal parameters (vector of size d).
#'        \item \code{ll_likelihood} : log-likelihood of \code{y} given the optimal parameters.
#'        \item \code{spec} : Specification.
#'        \item \code{is.init} : Indicating if estimation was made with do.init option.
#'        \item \code{y} :  Vector (of size T) of observations..
#'        }
#' The \code{MSGARCH_MLE_FIT} contains these methods:
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
#' # load data
#' data("sp500ret")
#' 
#' 
#' # create model specification
#' spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                               do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#' 
#' set.seed(123)
#' 
#' # fit the model on the data with ML estimation
#' fit = MSGARCH::fit.mle(spec = spec, y = sp500ret, 
#'                        ctr = list(do.init = TRUE, NP = 100, itermax = 100))
#' @import DEoptim nloptr
#' @export
fit.mle <- function(spec, y, ctr = list())
{
  UseMethod("fit.mle", spec)
}

#' @export
fit.mle.MSGARCH_SPEC = function(spec, y, ctr = list()) {
  
  ctr = f.process.ctr(ctr)
  ctr.optim = list(trace = 0, maxit = ctr$maxit)
  ctr.deoptim = DEoptim::DEoptim.control(NP = ctr$NP, itermax = ctr$itermax, 
    trace = FALSE, initialpop = matrix(spec$theta0, nrow = ctr$NP, 
      ncol = length(spec$theta0)))
  ctr.slsqp = list(maxeval = 10000, xtol_rel = 1e-08)
  
  
  
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
  ll_likelihood = f.kernel(theta)
  
  if (ll_likelihood == -1e+10) {
    tmp = DEoptim::DEoptim(fn = f.nll, lower = lower, upper = upper, control = ctr.deoptim)
    theta = tmp$optim$bestmem
    ll_likelihood = f.kernel(theta)
  }
  
  out = list(theta = theta, ll_likelihood = ll_likelihood, spec = spec, is.init = any(ctr$do.init || spec$do.init), y = y)
  class(out) = "MSGARCH_MLE_FIT"
  return(out)
}

