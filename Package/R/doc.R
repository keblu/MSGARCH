#THIS FILE IS FOR DOCUMENTATION PURPOSE ONLY FOR FUNCTION AVAILABLE 
#ONCE A SPECIFICATION HAS BEEN CREATED WITH f.create.spec.

#' Simulation function.
#' @description  Function inside a specification returning a simulated process.
#' @param n   Simulation length.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param burnin (integer >= 0) Burnin period discarded (first simulation draws). (default: \code{burnin = 500})
#' @param do.state  Boolean  indicating if the simulated state are also output. (default: \code{log = TRUE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#' spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#' do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#' y = spec$f.sim(n = 1000, theta = spec$theta0, burnin = 500, do.state = TRUE)
#' 
#' @return A list containing one or two components.
#' \itemize{
#' \item \code{draws}: vector (of size n) or matrix (of size M x n) of simulated draws.
#' \item \code{state}: vector (of size n) or matrix (of size M x n) of simulated states.
#'  The \code{state} value appear only if \code{do.state = TRUE}.
#' }
f.sim = function(n, theta, burnin = 500, do.state = FALSE) {
  return(NULL)
}

#' Conditional variance in each regime.
#' @description Function inside a specification returning the conditional variance of each regime.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'ht = spec$f.ht(theta = spec$theta0, y = sp500ret)
#' @return Condititional variance time serie (array of size T + 1 x M x K) for each regime.
f.ht = function(theta, y) {
  return(NULL)
}

#' Kernel function.
#' @description Function inside a specification returning the kernel value of a vector of observations.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param log  Boolean indicating if the log kernel is returned. (default: \code{log = TRUE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#'  The kernel is a combination of the prior and the likelihood function. 
#'  The kernel is equal to prior(\eqn{\theta}) + L(y|\eqn{\theta}) where L is the likelihood
#'  of y given the parameter \eqn{\theta}. When doing optimization, the goal is to minimize the negative log-kernel.
#'  \itemize{
#'  \item Details on the prior \cr
#'        The prior is different for each specification. It ensures that the \eqn{\theta} makes the conditional variance process stationary, positive,
#'        and that it respect that the sums of the probabilities in the case of a multiple-regime models are all equal to 1. If any of these three conditions is not respected the prior return \code{-1e10}, meaning that the optimizer or sampler
#'        will know that \eqn{\theta} is not a good candidate.
#'   }
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'kernel = spec$f.kernel(theta = spec$theta0, y = sp500ret, log = TRUE)
#' @references Hamilton, J. D. (1989) A New Approach to the Economic Analysis of Nonstationary Time Series and the Business Cycle. \emph{Econometrica}, 57, pp.357-38
#' @return Kernel or log-kernel value (scalar or vector of size M) of the vector of observations. 
f.kernel = function(theta, y, log = TRUE) {
  return(NULL)
}

#' Unconditional variance of each regime.
#' @description Function inside a specification returning the unconditional variance of the process in each state.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'unc.vol = spec$f.unc.vol(theta = spec$theta0)
#' @return Unconditional variance (vector of size K or matrix of size M x K) of each regime. 
f.unc.vol = function(theta) {
  return(NULL)
}

#' Predictive density function.
#' @description Function inside a specification returning the predictive probability density of a vector of points.
#' @param x  Vector (of size N) of point to be evaluated
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param log  Boolean indicating if the log-density is returned. (default: \code{log = TRUE})
#' @details  If a matrix of MCMC posterior draws estimates is given, the Bayesian predictive density is calculated.
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'                              
#'set.seed(123)
#'
#'x = rnorm(100)
#'
#'pred = spec$f.pred(x = x, theta = spec$theta0, y = sp500ret, log = TRUE)
#' @return Predictive density or log-density of \code{x} (vector of size N).
f.pred = function(x, theta, y, log = TRUE) {
  return(NULL)
}

#'Probability Integral Transform at T + 1.
#' @description Function inside a specification returning the predictive Probability integral transform (PIT).
#' @param x Vector (of size N) of point to be evaluated
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param do.norm  Boolean indicating if the PIT value are transforms into standard Normal variate. (\code{do.norm = FALSE})
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian Probability integral transform is calculated.
#'          The \code{do.norm} argument transforms the PIT value into Normal variate so that normality test can be done.
#' @return Probability integral transform of the points \code{x} or Normal variate derived from the Probability integral transform of \code{x} (vector of size N).
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'                              
#'set.seed(123)
#'
#'x = rnorm(100)
#'
#'pit = spec$f.pit(x = x theta = spec$theta0, y = sp500ret, do.norm = FALSE)
f.pit = function(x, theta, y, do.norm = FALSE) {
  return(NULL)
}

#' Value-at-Risk And Expected-shortfall functions.
#' @description Function inside a specification returning the Value-at-risk and Expected-shortfall.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param level Vector (of size A) of Value-at-risk and Expected-shortfall levels. (default: \code{level = c(0.95,0.99)})
#' @param ES  Boolean indicating if Expected-shortfall is also calculated. (default: \code{ES = TRUE})
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian Value-at-Risk and Expected-shortfall are calculated.
#' @return A list containing one or two components:
#' \itemize{
#' \item VaR : Value-at-Risk at the choosen levels (vector of size A).
#' \item ES  : Expected-shortfall at the choosen levels (vector of size A).
#' }
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'risk = spec$f.risk(theta = spec$theta0, y = sp500ret,level = c(0.95,0.99), ES = TRUE)
f.risk = function(theta, y, level = c(0.95,0.99), ES = TRUE) {
  return(NULL)
}

#' Simulation function at T + 1.
#' @description Function inside a specification returning random draws at T + 1.
#' @param n  Number of random draws to be generated.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param do.state  Boolean indicating if the simulated state are also output. (default: \code{do.state = FALSE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @return A list containing one or two components:
#' \itemize{
#' \item \code{draws}: vector (of size n) or matrix (of size M x n) of simulated draws at T + 1.
#' \item \code{state}: vector (of size n) or matrix (of size M x n) of simulated states at T + 1.
#'  The \code{state} value appear only if \code{do.state = TRUE}.
#' }
#' @examples 
#' spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#' do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#' rnd = spec$f.rnd(n = 1000, theta = theta, y = sp500ret, do.state = TRUE)
f.rnd = function(n, theta, y, do.state = FALSE) {
  return(NULL)
}

#' Probability density function at T + 1.
#' @description Function inside a specification returning the probability density of a vector of points.
#' @param x Vector (of size N) of point to be evaluated
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param log  Boolean indicating if the log-density is returned. (default: \code{log = TRUE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#'  The \code{f.pdf} function uses the last variance estimate by filtering.
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'                              
#'set.seed(123)
#'
#'x = rnorm(100)
#'
#'pdf = spec$f.pdf(x = x, theta = spec$theta0, y = sp500ret, log = TRUE)
#' @return Probability density or log-density of the points \code{x} (vector of size N or matrix of size M x N).
f.pdf = function(x, theta, y, log = TRUE) {
  return(NULL)
}

#' Cumulative density function at T + 1.
#' @description Function inside a specification returning the cumulative density of a vector of points.
#' @param x Vector (of size N) of point to be evaluated.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param log  Boolean indicating if the log cumulative is returned. (default: \code{log = TRUE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' The \code{f.cdf} function uses the last variance estimate by filtering.
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'                              
#'set.seed(123)
#'
#'x = rnorm(100)
#'
#'cdf = spec$f.cdf(x = x, theta = spec$theta0, y = sp500ret, log = FALSE)
#' @return Cumulative density or log-density of the points \code{x} (vector of size N or matrix of size M x N).
f.cdf = function(x, theta, y , log = TRUE) {
  return(NULL)
}


#' State probabilities filtering function.
#' @description Function inside a specification returning the filtered state probabilities.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'Pstate  = spec$f.Pstate(theta = spec$theta0, y = sp500ret)
#'@return Filtered state probabilities (array of size T x M x K).
f.Pstate = function(theta, y) {
  return(NULL)
}

#' State probabilities at T + 1.
#' @description Function inside a specification returning the state probabilities at  T + 1.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @return State probabilities at T + 1 (matrix of size M x K).
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'Plast = spec$f.Plast(theta = spec$theta0, y = sp500ret) 
f.Plast = function(theta, y) {
  return(NULL)
}
