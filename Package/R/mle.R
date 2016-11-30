#' ML estimation.
#' @description Method that performs Maximum Likelihood estimation of a \code{MSGARCH_SPEC} object on a set of observations.
#' @param y Vector (of size T) of observations.
#' @param spec Model specification created with \code{\link{create.spec}}.
#' @param ctr List of control parameters.
#'        The control parameters have two components to it:
#'        \itemize{
#'        \item \code{do.init} : Boolean indicating if there is a pre-optimization with the \R package \code{DEoptim} (Ardia et al., 2011). (Default: \code{do.init = FALSE})
#'        \item \code{NP} : Number of parameter vectors in the population in \code{DEoptim} optimization. (Default: \code{NP = 200})
#'        \item \code{itermax} : Maximum iteration (population generation) allowed in \code{DEoptim} optimization. (Default: \code{maxit = 200})
#'        \item \code{theta0} : Starting value for the chain (if empty the specification default value are used).
#'        \item \code{do.enhance.theta0} : Boolean indicating if the default parameters value are enhance using \code{y} variance. (Default: \code{do.enhance.theta0 = TRUE})
#'        }
#' @return A list of class \code{MSGARCH_MLE_FIT} containing five components:
#'        \itemize{
#'        \item \code{theta} : Optimal parameters (vector of size d).
#'        \item \code{log_kernel} : log-kernel of \code{y} given the optimal parameters.
#'        \item \code{spec} : Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#'        \item \code{is.init} : Indicating if estimation was made with do.init option.
#'        \item \code{y} :  Vector (of size T) of observations.
#'        }
#' The \code{MSGARCH_MLE_FIT} contains these methods:
#' \itemize{
#' \item \code{\link{AIC}} : Compute Akaike information criterion (AIC).
#' \item \code{\link{BIC}} : Compute Bayesian information criterion (BIC).
#' \item \code{\link{ht}}  : Conditional volatility in each regime.
#' \item \code{\link{kernel}} : Kernel method.
#' \item \code{\link{unc.vol}} : Unconditional volatility in each regime.
#' \item \code{\link{pred}} : Predictive method.
#' \item \code{\link{pit}} : Probability Integral Transform.
#' \item \code{\link{risk}} : Value-at-Risk And Expected-Shortfall methods.
#' \item \code{\link{simahead}} : Step ahead simulation method.
#' \item \code{\link{sim}} : Simulation method.
#' \item \code{\link{pdf}} : Probability density function.
#' \item \code{\link{cdf}} : Cumulative function.
#' \item \code{\link{Pstate}} : State probabilities filtering method.
#' \item \code{summary} : Summary of the fit.
#' }
#' 
#' @details The Maximum likelihood estimation uses the \R package \code{nloptr} (Johnson, 2014) for main optimizer 
#' while it uses the \R package \code{DEoptim} when \code{do.init = TRUE} as an initialization for nloptr.
#'  The starting parameters are the specification default parameters.
#'  The argument \code{do.enhance.theta0} uses the volatilities of rolling windows of \code{y} and adjust the starting parameters of
#'  the specification so that the unconditional volatility of each regime
#'  is set to different quantiles of the volatilities of the rolling windows of \code{y}.
#' @references Ardia, D. Boudt, K. Carl, P. Mullen, K. M. & Peterson, B. G. (2011). Differential Evolution with \code{DEoptim}. \emph{R Journal}, 3, pp. 27-34
#' @references Ardia, D. Mullen, K. M. Peterson, B. G. & Ulrich, J. (2015). \code{DEoptim}: Differential Evolution in \R. \url{https://cran.r-project.org/package=DEoptim}
#' @references Mullen, K. M. Ardia, D. Gil, D. L. Windover, D. Cline, J.(2011) \code{DEoptim}: An \R Package for Global Optimization by Differential Evolution. \emph{Journal of Statistical Software}, 40, pp. 1-26, DOI:	\url{http://dx.doi.org/10.18637/jss.v040.i06}
#' @references Johnson, S. G. (2014). The NLopt Nonlinear-Optimization. \url{https://cran.r-project.org/package=nloptr}.
#' @examples 
#' # load data
#' data("sp500")
#' sp500 = sp500[1:1000]
#' 
#' # create model specification
#' spec = MSGARCH::create.spec() 
#' 
#' # fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#' fit = MSGARCH::fit.mle(spec = spec, y = sp500)
#' summary(fit)
#' @importFrom stats runif           
#' @import DEoptim nloptr dfoptim
#' @export
fit.mle <- function(spec, y, ctr = list()) {
  UseMethod("fit.mle", spec)
}

#' @export
fit.mle.MSGARCH_SPEC <- function(spec, y, ctr = list()) {
  y <- f.check.y(y)
  ctr <- f.process.ctr(ctr)
  
  if (is.null(ctr$theta0)) {
    ctr$theta0 <- spec$theta0
  }
  if (isTRUE(ctr$do.enhance.theta0)) {
    ctr$theta0 <- f.enhance.theta(spec = spec, theta = ctr$theta0, y = y)
  }
  theta0.init <- ctr$theta0
  lower <- spec$lower
  upper <- spec$upper
  if (any(ctr$do.init || spec$do.init)) {
    pop        = matrix(runif(length(lower)*10000, min = lower, max = upper), ncol = length(spec$theta0), byrow = TRUE)
    pop_inv    = t(replicate(n = 10000, upper)) +  t(replicate(n = 10000, lower)) - pop
    total_pop  = rbind(pop,pop_inv)
    # DA fix naming to loglik
    likelihood = MSGARCH::kernel(object = spec,theta = total_pop, y = y, log = TRUE)
    ind        = sort(likelihood, decreasing = TRUE, index.return = TRUE)
    initialpop = total_pop[ind$ix[1:ctr$NP], ]
  }
  f.kernel <- function(x, log = TRUE) {
    out = MSGARCH::kernel(spec, x, y = y, log = log)
    return(out)
  }
  f.nll <- function(x) -f.kernel(x, log = TRUE)
  if (any(ctr$do.init || spec$do.init)) {
    ctr.deoptim <- DEoptim::DEoptim.control(NP = ctr$NP, itermax = ctr$itermax, trace = FALSE,
                                            initialpop = initialpop)
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
  }
  theta <- f.find.theta0(f.kernel, theta0 = theta0.init, lower = lower, upper = upper,
                        f.ineq = spec$rcpp.func$ineq_func, ineqlb = spec$ineqlb,
                        inequb = spec$inequb)
  log_kernel <- f.kernel(theta)
  if (log_kernel == -1e+10) {
    str = "f.find.theta0 -> DEoptim initialization"
    ctr.deoptim <- DEoptim::DEoptim.control(NP = ctr$NP, itermax = ctr$itermax, trace = FALSE)
    tmp         <- DEoptim::DEoptim(fn = f.nll, lower = lower, upper = upper, control = ctr.deoptim)
    theta       <- tmp$optim$bestmem
    log_kernel  <- f.kernel(theta)
  }
  theta       <- f.sort.theta(spec = spec, theta)
  theta0.init <- f.sort.theta(spec = spec, theta0.init)
  theta       <- matrix(theta, ncol = length(theta))
  colnames(theta) <- colnames(spec$theta0)
  out <- list( theta = theta, log_kernel = log_kernel, spec = spec,
              is.init = any(ctr$do.init || spec$do.init), y = y,theta0.init = theta0.init)
  class(out) <- "MSGARCH_MLE_FIT"
  return(out)
}