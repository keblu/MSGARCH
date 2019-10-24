#' @title Model specification.
#' @description Creates a model specification before fitting and
#' using the \pkg{MSGARCH} functionalities.
#' @param variance.spec \code{list} with element \code{model}.
#' \code{model} is a \code{character} vector (of size K, number of regimes)
#' with the variance model specifications. Valid models  are  \code{"sARCH"}, \code{"sGARCH"},
#' \code{"eGARCH"}, \code{"gjrGARCH"}, and \code{"tGARCH"} (see *Details*).
#' Default: \code{model = c("sGARCH", "sGARCH")}.
#' @param distribution.spec \code{list} with element \code{distribution}.
#' \code{distribution} is a \code{character} vector (of size K)
#' of conditional distributions. Valid distributions are \code{"norm"}, \code{"snorm"},
#' \code{"std"}, \code{"sstd"}, \code{"ged"}, and \code{"sged"}  (see *Details*).
#' The vector must be of the same length as the models' vector in \code{variance.spec}.\cr
#' Default: \code{distribution = c("norm", "norm")}.
#' @param switch.spec \code{list} with element \code{do.mix} and \code{K}.\cr
#'  \code{do.mix} is a \code{logical} indicating if the specification is a mixture type.
#' If \code{do.mix = TRUE}, a Mixture of GARCH is created, while if \code{do.mix = FALSE},
#' a Markov-Switching GARCH is created (see *Details*). Default: \code{do.mix = FALSE}.\cr
#' \code{K} is a optional \code{numeric} scalar indicating the number of regime.
#'  In the case where a single regime is specified in \code{variance.spec} and \code{distribution.spec},
#'  this parameter allows to automatically expand this single regime to \code{K} similar
#'  regimes without the need to explicitly define them in \code{variance.spec}
#'  and \code{distribution.spec} (see *Examples*).
#' @param constraint.spec  \code{list} with element \code{fixed} and \code{regime.const}.
#' Only one of \code{fixed} and \code{regime.const}
#'  can be set by the user as it is not allowed to set both at the same time. \cr
#'  \code{fixed} is a  \code{list} with \code{numeric} entries and named elements.
#'  This argument controls for fixed parameters defined by the user.  The names of the entries
#' in the \code{list}  have to coincide with the names of the model parameters.\cr
#'  For instance, if \code{contraint.spec  = list(fixed = list(beta_1 = 0))}, \code{beta_1}
#'  will be fixed to \code{0} during optimization. \cr
#'  \code{regime.const} is a \code{character} vector.
#'  This argument controls for the parameters which are set equal across regimes.
#'  The names of the entries in the \code{list}  have to coincide with the names of the model parameters
#'  minus the regime indicator.\cr
#'  For instance, if \code{contraint.spec  = list(regime.const = c("beta"))},
#'  all the parameters named \code{beta} will be the same in all regimes during optimization.\cr
#' @param prior \code{list} with element \code{mean} and \code{sd}. The element \code{mean} and \code{sd}
#' are \code{list}  with \code{numeric} and named elements that allow to adjust the prior mean and
#' standard deviation of the truncated Normal prior.
#' The names of the entries in the lists have to coincide with the names of the model parameters.\cr
#' For instance, if \code{prior = list(mean = list(beta_1 = 0.7), sd = list(beta_1 = 0.1))},
#' the prior mean of \code{beta_1} will be set to \code{0.7} while the prior standard
#' deviation will set to \code{0.1}.
#' @return A list of class \code{MSGARCH_SPEC} with the following elements:\cr
#' \itemize{
#' \item \code{par0}: Vector (of size d) of default parameters.
#' \item \code{is.mix}: Logical indicating if the specification is a mixture.
#' \item \code{K}: Number of regimes.
#' \item \code{lower}: Vector (of size d) of lower parameters' bounds.
#' \item \code{upper}: Vector (of size d) of upper parameters' bounds.
#' \item \code{n.params}:  Vector (of size K) of the total number of parameters by
#' regime including distributions' parameters.
#' \item \code{n.params.vol}:  Vector (of size K) of the total number of parameters
#' by regime excluding distributions' parameters.
#' \item \code{label}: Vector (of size d) of parameters' labels.
#' \item \code{name}: Vector (of size K) of model specifications' names.
#' \item \code{func}: List of internally used \R functions.
#' \item \code{rcpp.func}: List of internally used \code{Rcpp} functions.
#' \item \code{fixed.pars}: List of user inputed fixed parameters.
#' \item \code{regime.const.pars}: Vector of user imputed parameter set equal across regimes.
#' \item \code{regime.fixed.pars}: Logical indicating if there is any fixed parameteter set by the user.
#' \item \code{regime.const.pars.bool}: Logical indicating if there is any parameteters
#'  equal across regime set by the user.
#' }
#' The \code{MSGARCH_SPEC} class has the following methods:
#' \itemize{
#' \item \code{simulate}: Simulation.
#' \item \code{\link{Volatility}}: In-sample conditional volatility.
#' \item \code{predict}: Forecast of the conditional volatility (and predictive distribution).
#' \item \code{\link{UncVol}}: Unconditional volatility.
#' \item \code{\link{PredPdf}}: Predictive density (pdf).
#' \item \code{\link{PIT}}: Probability Integral Transform.
#' \item \code{\link{Risk}}: Value-at-Risk and Expected-Shortfall.
#' \item \code{\link{State}}: State probabilities (smoothed, filtered, predictive, Viterbi).
#' \item \code{\link{FitML}}: Maximum Likelihood estimation.
#' \item \code{\link{FitMCMC}}: Bayesian estimation.
#' \item \code{print} and \code{summary}: Summary of the created specification.
#' }
#' @details The Markov-Switching specification is based on the
#' Haas et al. (2004a) MSGARCH specification. It is a MSGARCH model that is separated
#' in K single-regime specifications  which are updated in parallel. Under the Haas et al. (2004a)
#' specification, the conditional variance is a function of past data and the current state.
#' The Mixture of GARCH option (\code{do.mix = TRUE}) is based on Haas et al. (2004b). A Mixture of GARCH is a mixture of distributions
#' where the variance process of each distribution is a single-regime process.\cr
#'  For the models, \code{"sARCH"} is the ARCH(1) model (Engle, 1982), \code{"sGARCH"} the GARCH(1,1) model
#' (Bollerslev, 1986), \code{"eGARCH"} the EGARCH(1,1) model (Nelson, 1991), \code{"gjrGARCH"}
#' the GJR(1,1) model (Glosten et al., 1993), and \code{"tGARCH"} the TGARCH(1,1) model (Zakoian, 1994).\cr
#' For the distributions, \code{"norm"} is the Normal distribution, \code{"std"} the
#' Student-t distribution, and \code{"ged"} the GED distribution.
#' Their skewed version, implemented via the Fernandez and & Steel (1998) transformation,
#' are \code{"snorm"}, \code{"sstd"} and \code{"sged"}.
#' Please see Ardia et al. (2019) for more details on the models and distributions.\cr
#' The user must choose between \code{fixed} or \code{regime.const} in \code{contraint.spec}
#' as both cannot be set at the same time. The \code{list} \code{fixed.pars}
#' will ensure that the chosen fixed parameters will be fixed during optimization according to
#' the values set by the user.
#' Thus only the non-fixed parameters are optimized. The vector \code{regime.const} will
#' ensure that the chosen parameters will be the same across regime during optimization.\cr
#' The \code{list} \code{mean} and \code{sd} in \code{prior} will adjust the prior mean and
#' prior standard deviation of the truncated Normal prior for MCMC
#' estimation via \code{\link{FitMCMC}} according to the inputed prior mean and standard deviation.
#' Those prior means and standard deviations that are not set will take on preset default values (a mean 
#' of zero and a variance of 1,000).
#' 
#' @references Ardia, D. Bluteau, K. Boudt, K. Catania, L. Trottier, D.-A. (2019).
#' Markov-switching GARCH models in \R: The \pkg{MSGARCH} package.
#' \emph{Journal of Statistical Software}, 91(4), 1-38.
#' \url{http://doi.org/10.18637/jss.v091.i04}
#' 
#' @references Engle, R. (1982).
#' Autoregressive conditional heteroscedasticity with estimates of the variance of United Kingdom inflation
#' \emph{Econometrica}, 50, 987-1008.
#' 
#' @references Bollerslev, T. (1986).
#' Generalized autoregressive conditional heteroskedasticity.
#' \emph{Journal of Econometrics}, 31, 307-327.
#' \url{http://doi.org/10.1016/0304-4076(86)90063-1}
#' 
#' @references Fernandez, C. & Steel, M. F. (1998).
#' On Bayesian modeling of fat tails and skewness.
#' \emph{Journal of the American Statistical Association}, 93, 359-371.
#' \url{http://doi.org/10.1080/01621459.1998.10474117}
#' 
#' @references Glosten, L. R. Jagannathan, R. & Runkle, D. E. (1993).
#' On the relation between the expected value and the volatility of the nominal excess return on stocks.
#' \emph{Journal of Finance}, 48, 1779-1801.
#' \url{http://doi.org/10.1111/j.1540-6261.1993.tb05128.x}
#' 
#' @references Haas, M. Mittnik, S. & Paolella, M. S. (2004a).
#' A new approach to Markov-switching GARCH models.
#' \emph{Journal of Financial Econometrics}, 2, 493-530.
#' \url{http://doi.org/10.1093/jjfinec/nbh020}
#' 
#' @references Haas, M. Mittnik, S. & Paolella, M. S. (2004b).
#' Mixed normal conditional heteroskedasticity.
#' \emph{Journal of Financial Econometrics}, 2, 211-250.
#' \url{http://doi.org/10.1093/jjfinec/nbh009}
#' 
#' @references Nelson, D. B. (1991).
#' Conditional heteroskedasticity in asset returns: A new approach.
#' \emph{Econometrica}, 59, 347-370.
#' 
#' @references Zakoian, J.-M. (1994).
#' Threshold heteroskedastic models.
#' \emph{Journal of Economic Dynamics and Control}, 18, 931-955.
#' \url{http://doi.org/10.1016/0165-1889(94)90039-6}
#' 
#' @examples
#' # create a Markov-switching specification
#' # MS-GARCH(1,1)-GJR(1,1)-Student
#' spec <- CreateSpec(variance.spec = list(model = c("sGARCH","gjrGARCH")),
#'                    distribution.spec = list(distribution = c("std","std")),
#'                    switch.spec = list(do.mix = FALSE))
#' print(spec)
#'
#' # create a 3-regime Markov-switching specification with the help of variable K
#' # MS(3)-GARCH(1,1)- Student
#' spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
#'                    distribution.spec = list(distribution = c("std")),
#'                    switch.spec = list(do.mix = FALSE, K = 3))
#' print(spec)
#'
#' # create a mixture specification
#' # MIX-GARCH(1,1)-GJR(1,1)-Student
#' spec <- CreateSpec(variance.spec = list(model = c("sGARCH","gjrGARCH")),
#'                    distribution.spec = list(distribution = c("std","std")),
#'                    switch.spec = list(do.mix = TRUE))
#' print(spec)
#'
#' # setting fixed parameter for the sGARCH beta parameter
#' # MS-GARCH(1,1)-GJR(1,1)-Student with beta_1 fixed to 0
#' spec <- CreateSpec(variance.spec = list(model = c("sGARCH","gjrGARCH")),
#'                    distribution.spec = list(distribution = c("std","std")),
#'                    switch.spec = list(do.mix = FALSE),
#'                    constraint.spec = list(fixed = list(beta_1 = 0)))
#' print(spec)
#'
#' # setting restriction for the shape parameter of the Student-t across regimes
#' # MS-GARCH(1,1)-GJR(1,1)-Student with shape parameter constraint across regime
#' spec <- CreateSpec(variance.spec = list(model = c("sGARCH","gjrGARCH")),
#'                    distribution.spec = list(distribution = c("std","std")),
#'                    switch.spec = list(do.mix = FALSE),
#'                    constraint.spec = list(regime.const = c("nu")))
#' print(spec)
#'
#' # setting custom parameter priors for the beta parameters
#' # MS-GARCH(1,1)-GJR(1,1)-Student with prior modification
#' spec <- CreateSpec(variance.spec = list(model = c("sGARCH","gjrGARCH")),
#'                    distribution.spec = list(distribution = c("std","std")),
#'                    switch.spec = list(do.mix = FALSE),
#'                    prior = list(mean = list(beta_1 = 0.9, beta_2 = 0.3),
#'                                 sd = list(beta_1 = 0.05, beta_2 = 0.01)))
#' print(spec)
#' @import Rcpp
#' @export
CreateSpec <- function(variance.spec = list(model = c("sGARCH", "sGARCH")),
                       distribution.spec = list(distribution = c("norm", "norm")),
                       switch.spec = list(do.mix = FALSE, K = NULL),
                       constraint.spec  = list(fixed = list(), regime.const = NULL),
                       prior = list(mean = list(), sd = list())) {
  
  ## check
  variance.spec     <- f_check_variance_spec(variance.spec)
  distribution.spec <- f_check_distribution_spec(distribution.spec, length(variance.spec$model))
  switch.spec       <- f_check_markov_spec(switch.spec)
  fixed.pars        <- constraint.spec$fixed
  regime.const.pars <- constraint.spec$regime.const
  prior.mean        <- prior$mean
  prior.sd          <- prior$sd
  
  if (!is.null(switch.spec$K)) {
    if (length(variance.spec$model) > 1 | length(distribution.spec$model) > 1) {
      stop("you can only use the variable K if you specified one
           regime in variance.spec and distribution.spec")
    } else {
      variance.spec$model = rep(variance.spec$model, switch.spec$K)
      distribution.spec$distribution = rep(distribution.spec$distribution, switch.spec$K)
    }
  }
  
  model  <- variance.spec$model
  do.mix <- switch.spec$do.mix
  distribution <- distribution.spec$distribution
  
  distribution <- distribution[1:length(model)]
  valid.distribution <- c("norm", "std", "ged", "snorm", "sstd", "sged")
  
  valid.model <- c("sARCH", "sGARCH", "eGARCH", "gjrGARCH", "tGARCH")
  
  if (length(distribution) != length(model)) {
    stop("\nCreateSpec-->error: model vector and distribution
         vector must be of the same length")
  }
  for (i in 1:length(distribution)) {
    if (is.null(distribution[i])) {
      distribution[i] <- "norm"
    }
    if (!is.character(distribution[1L])) {
      stop(paste0("\nCreateSpec-->error: The distribution #", i, "
                  argument must be a character"))
    }
    if (!any(distribution[i] == valid.distribution)) {
      stop(paste0("\nCreateSpec-->error: The distribution #", i, "
                  does not appear to be a valid choice."))
    }
    }
  
  
  if (length(distribution) == 1 & do.mix == TRUE) {
    do.mix = FALSE
    warning("Single regime model, automatically setting do.mix = FALSE")
  }
  
  for (i in 1:length(model)) {
    if (!any(model[i] == valid.model)) {
      stop(paste0("\nCreateSpec-->error: Model #", i, "
                  does not appear to be a valid choice.\n",
                  call. = FALSE))
    }
  }
  
  if (is.null(do.mix)) {
    do.mix <- FALSE
  }
  
  if (isTRUE(do.mix) || !isTRUE(do.mix)) {
  } else {
    stop("\nCreateSpec-->error: do.mix must be a TRUE or FALSE\n", call. = FALSE)
  }
  
  
  models.merge <- paste0(model, "_", distribution)
  
  out <- f_spec(models = models.merge, do.mix = do.mix)
  
  ## rixed par
  fixed.pars     <- f_check_parameterConstraints(fixed.pars, out$label)
  out$fixed.pars <- fixed.pars
  
  if (length(fixed.pars) > 0) {
    out$fixed.pars.bool = TRUE
  } else {
    out$fixed.pars.bool = FALSE
  }
  
  ## regime.const.pars
  if (!is.null(regime.const.pars)) {
    regime.const.pars <- f_check_regime_const_pars(regime.const.pars, length(variance.spec$model), out$label,
                                                   variance.spec$model, distribution.spec$distribution)
  }
  
  if (!is.null(regime.const.pars) && length(fixed.pars) > 0) {
    stop("Only regime.const.pars or fixed.pars can be defined.")
  }
  out$regime.const.pars <- regime.const.pars
  
  if (length(regime.const.pars) > 0) {
    out$regime.const.pars.bool = TRUE
  } else {
    out$regime.const.pars.bool = FALSE
  }
  
  ## prior Mean
  if (length(prior.mean) >= 1) {
    prior.mean <- f_check_parameterPriorMean(prior.mean, out$label)
    out$prior.mean <- f_substitute_fixedpar(out$prior.mean, prior.mean)
    out$rcpp.func$set_mean(out$prior.mean)
  }
  ## prior Sd
  if (length(prior.sd) >= 1) {
    prior.sd <- f_check_parameterPriorMean(prior.sd, out$label)
    out$prior.sd <- f_substitute_fixedpar(out$prior.sd, prior.sd)
    out$rcpp.func$set_sd(out$prior.sd)
  }
  
  class(out) <- "MSGARCH_SPEC"
  return(out)
  }

# Function used to build the specification
#' @importFrom methods new
f_spec <- function(models, do.mix = FALSE) {
  
  models.list <- NULL
  for (j in 1:length(models)) {
    models.list[[j]] <- get(models[j])
  }
  
  K <- length(models)
  if (K == 1L) {
    do.mix <- FALSE
  }
  if (K > 1L) {
    tmp <- list()
    for (i in 1:K) {
      tmp[[i]] <- methods::new(models.list[[i]])
    }
  } else {
    tmp <- methods::new(models.list[[1L]])
  }
  if (K > 1L) {
    mod <- methods::new(MSgarch, tmp)
  } else {
    mod <- tmp
  }
  mod$name               <- models
  n.params               <- mod$NbParams
  n.params.vol           <- mod$NbParamsModel
  rcpp.func              <- list()
  rcpp.func$calc_ht      <- mod$calc_ht
  rcpp.func$eval_model   <- mod$eval_model
  rcpp.func$sim          <- mod$f_sim
  rcpp.func$pdf_Rcpp     <- mod$f_pdf
  rcpp.func$cdf_Rcpp     <- mod$f_cdf
  rcpp.func$rnd_Rcpp     <- mod$f_rnd
  rcpp.func$pdf_Rcpp_its <- mod$f_pdf_its
  rcpp.func$simahead     <- mod$f_simAhead
  rcpp.func$cdf_Rcpp_its <- mod$f_cdf_its
  rcpp.func$unc_vol_Rcpp <- mod$f_unc_vol
  rcpp.func$get_sd       <- mod$f_get_sd
  set_sd.base            <- mod$f_set_sd
  rcpp.func$set_sd       <- mod$f_set_sd
  rcpp.func$get_mean     <- mod$f_get_mean
  set_mean.base          <- mod$f_set_mean
  rcpp.func$set_mean     <- mod$f_set_mean
  prior.mean             <- mod$f_get_mean()
  prior.sd               <- mod$f_get_sd()
  
  if (K > 1L) {
    rcpp.func$get_Pstate_Rcpp <- mod$f_get_Pstate
  } else {
    rcpp.func$get_Pstate_Rcpp <- function(par, y) {
      FiltProb   <- matrix(data = 1, nrow = nrow(y), ncol = 1)
      PredProb   <- matrix(data = 1, nrow = nrow(y) + 1L, ncol = 1)
      SmoothProb <- matrix(data = 1, nrow = nrow(y) + 1L, ncol = 1)
      viterbi    <- rep(x = 0, length.out = nrow(y) - 1)
      out        <- list(FiltProb = FiltProb, PredProb = PredProb, SmoothProb = SmoothProb, Viterbi = viterbi)
      return(out)
    }
  }
  
  nb_total_params <- sum(n.params)
  func <- list()
  func$f.do.mix <- function(par) {
    out <- f_par_mixture(K, nb_total_params, par)
    return(out)
  }
  func$f.do.mix.reverse <- function(par) {
    out <- f_par_mixture_reverse(K, nb_total_params, par)
    return(out)
  }
  loc <- c(0, cumsum(n.params))
  for (i in 1:K) {
    mod$label[(loc[i] + 1L):loc[i + 1L]] <- paste0(mod$label[(loc[i] + 1L):loc[i + 1L]], "_", i)
  }
  if (K > 1) {
    for (i in 0:(K - 1L)) {
      mod$label[(loc[K + 1L] + K * i + 1L - i):(loc[K + 1L] + K * i + K - 1L - i)] <-
        paste0(mod$label[(loc[K + 1L] + K * i + 1L - i):(loc[K + 1L] + K * i + K - 1L - i)], "_", i + 1, "_", 1:(K - 1))
    }
  }
  if (isTRUE(do.mix)) {
    mod$lower       <- as.vector(func$f.do.mix.reverse(mod$lower))
    newParamsLength <- length(mod$lower)
    mod$upper       <- as.vector(func$f.do.mix.reverse(mod$upper))
    mod$theta0      <- as.vector(func$f.do.mix.reverse(mod$theta0))
    mod$label       <- mod$label[1:newParamsLength]
    mod$label[(nb_total_params + 1):length(mod$label)] <- paste0("P_", 1:(K - 1))
  }
  mod$theta0 <- matrix(mod$theta0, ncol = length(mod$theta0))
  names(prior.mean) = names(prior.sd) = mod$label[1:length(prior.mean)]
  colnames(mod$theta0) <- mod$label
  out <- list(par0 = mod$theta0[1L,], is.mix = do.mix, K = K, lower = mod$lower, upper = mod$upper,
              n.params = n.params, n.params.vol = n.params.vol, label = mod$label, name = mod$name,
              prior.mean = prior.mean, prior.sd = prior.sd, rcpp.func = rcpp.func, func = func)
  return(out)
}
