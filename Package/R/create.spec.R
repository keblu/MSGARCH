#' Model specification
#' @description Function for creating a model specification before fitting and using the \R package \code{MSGARCH} functionalities.
#' @param model  Vector (of size K) containing the variance model specifications.
#'                       Valid models  are  \code{"sGARCH"},  \code{"eGARCH"},
#'                      \code{"gjrGARCH"}, \code{"tGARCH"}, and  \code{"GAS"}. \cr (Default: \code{model = c("sGARCH", "sGARCH"}))
#' @param distribution  Vector (of size K) of conditional densities. Valid
#'                      distribution are  \code{"norm"}, \code{"std"}, and  \code{"ged"}. The vector must be of the same length as the models vector. \cr (Default: \code{distribution = c("norm", "norm"}))
#' @param do.skew Vector (of size K) of boolean indicating if the conditional density is skewed. The vector must be of the same length as the distributions vector. \cr (Default: \code{do.skew = c(FALSE, FALSE}))
#' @param do.mix  Boolean indicating if the specification is a mixture type.  If the argument is \code{TRUE}, a Mixture of GARCH is created,  
#'                while if the argument is \code{FALSE}, a Markov-Switching GARCH is created (see details). (Default: \code{do.mix = FALSE}) 
#' @param do.shape.ind  Boolean indicating if the distribution are Regime-Independent. If the argument is \code{TRUE}, all distributions are
#'                           the same and the distribution parameters does not dependent on the regime in which the distribution is attributed to.
#'                           If the argument is \code{TRUE}, all distributions in the distribution argument and all skew argument must be the same. (Default: \code{do.shape.ind = FALSE}) 
#' @return A list  of class \code{MSGARCH_SPEC} containing variables related to the created specification. \cr
#' The list contains:\cr
#' 
#' \itemize{
#' \item \code{theta0} : Vector (of size d) of default parameters.
#' \item \code{is.mix} : Boolean indicating if the specification is a mixture.
#' \item \code{is.shape.ind} : Boolean indicating if the distribution parameters are regime-independent.
#' \item \code{K} : Number of regimes.
#' \item \code{sigma0} : Default variance-covariance matrix (of size K x K) used for the Bayesian esimation.
#' \item \code{lower} : Vector (of size d) of lower parameters bound.
#' \item \code{upper} : Vector (of size d) of upper parameters bound.
#' \item \code{ineqlb} : Vector (of size d) of lower inequality bound.
#' \item \code{inequb} :  Vector (of size d) of upper inequality bound.
#' \item \code{n.params} :  Vector (of size K) of the total number of parameters by regime including distribution parameters.
#' \item \code{n.params.vol} :  Vector (of size K) of the total number of parameters by regime excluding distribuion parameters.
#' \item \code{do.init} : Boolean indicating the default \code{do.init}  argument.
#' \item \code{label} : Vector (of size d) of parameters label.
#' \item \code{name} : Vector (of size K) of model specification name.
#' \item \code{func} : List of \R functions internaly used.
#' \item \code{rcpp.func} : List of \code{Rcpp} functions internaly used.
#' }
#' The \code{MSGARCH_SPEC} class possesses these methods:
#' \itemize{
#' \item \code{\link{sim}} : Simulation method.
#' \item \code{\link{simahead}} : Step ahead simulation method.
#' \item \code{\link{ht}}  : Conditional volatility in each regime.
#' \item \code{\link{kernel}} : Kernel method.
#' \item \code{\link{unc.vol}} : Unconditional volatility in each regime.
#' \item \code{\link{pred}} : Predictive method.
#' \item \code{\link{pit}} : Probability Integral Transform.
#' \item \code{\link{risk}} : Value-at-Risk And Expected-Shortfall methods.
#' \item \code{\link{pdf}} : Probability density function.
#' \item \code{\link{cdf}} : Cumulative function.
#' \item \code{\link{Pstate}} : State probabilities filtering method.
#' \item \code{\link{fit.mle}} : Maximum Likelihood estimation.
#' \item \code{\link{fit.bayes}} : Bayesian estimation.
#' \item \code{print} and \code{summary} : Summary of the created specification.
#' }
#' @useDynLib MSGARCH
#' @details The Markov-Switching specification created is based on the Haas et al. (2004a) MSGARCH specification.
#'   It is a MSGARCH model that is separated in K single-regimes specifications  which are updated in parallel.
#'   Under this specification, the conditional variance is a function of the past data and the current state.
#'   The Mixture of GARCH option is based on the Haas et al. (2004b). A Mixture of GARCH is a mixture of distribution 
#'   where the variance process of each distribution is a single-regime process. Every single-regime specification is a one-lag process (e.g., GARCH(1,1))
#'   since it has proved to be sufficient in financial econometrics and it reduces models complexity which can become a problem during the optimization procedure
#' @references Bollerslev, T. (1986). Generalized Autoregressive Conditional Heteroskedasticity. \emph{Journal of Econometrics}, 31, pp. 307-327.
#' @references Creal, D. Koopman, S. J. & Lucas, A. (2013). Generalized Autoregressive Score Models with Applications. \emph{Journal of Applied Econometrics}, 28, pp. 777-795.
#' @references Fernandez, C. & Steel, M. F. (1998). On Bayesian Modeling of Fat Tails and Skewness. \emph{Journal of the American Statistical Association}, 93, pp. 359-371.
#' @references Glosten, L. R. Jagannathan, R. & Runkle, D. E. (1993). On the Relation Between the Expected Value and the Volatility of the Nominal Excess Return on Stocks. \emph{Journal of Finance}, 48, pp. 1779-1801.
#' @references Haas, M. Mittnik, S. & Paolella, M. S. (2004a). A New Approach to Markov-Switching GARCH Models. \emph{Journal of Financial Econometrics}, 2, pp. 493-530.
#' @references Haas, M. Mittnik, S. & Paolella, M. S. (2004b). Mixed Normal Conditional Heteroskedasticity. \emph{Journal of Financial Econometrics}, 2, pp. 211-250.
#' @references Nelson, D. B. (1991). Conditional Heteroskedasticity in Asset Returns: A New Approach. \emph{Econometrica}, 59, pp. 347-370.
#' @references Zakoian, J.-M. (1994). Threshold Heteroskedastic Models. \emph{Journal of Economic Dynamics and Control}, 18, pp. 931-955.
#' @examples 
#' # create model specification
#' spec = MSGARCH::create.spec(model = c("sGARCH","gjrGARCH"), distribution = c("norm","std"),
#'                              do.skew = c(TRUE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#' print(spec)
#' @import Rcpp
#' @export
create.spec <- function(model = c("sGARCH", "sGARCH"),
                        distribution = c("norm", "norm"),
                        do.skew = c(FALSE, FALSE),
                        do.mix = FALSE,
                        do.shape.ind = FALSE) {
  require("MSGARCH")
  distribution <- distribution[1:length(model)]
  do.skew <- do.skew[1:length(model)]
  valid.distribution <- c("norm", "std", "ged")
  valid.model <- c("sGARCH", "eGARCH", "gjrGARCH", "tGARCH", "GAS")
  if (length(distribution) != length(model)) {
    stop("\ncreate.spec-->error: model vector and distribution
         vector must be of the same length")
  }
  for (i in 1:length(distribution)) {
    if (is.null(distribution[i])) {
      distribution[i] <- "norm"
    }
    if (!is.character(distribution[1])) {
      stop(paste0("\ncreate.spec-->error: The distribution #", i, "
                  argument must be a character"))
    }
    if (!any(distribution[i] == valid.distribution)) {
      stop(paste0("\ncreate.spec-->error: The distribution #", i, "
                  does not appear to be a valid choice."))
    }
  }
  skew_tag <- do.skew
  for (i in 1:length(do.skew)) {
    if (is.null(do.skew[i]) || !isTRUE(do.skew[i])) {
      skew_tag[i] <- "sym"
    } else if (isTRUE(do.skew[i])) {
      skew_tag[i] <- "skew"
    } else {
      stop(paste0("\ncreate.spec-->error: do.skew #", i, "
                  argument must be a boolean"))
    }
  }
  dist.merge <- paste0(distribution, "_", skew_tag)
  for (i in 1:length(model)) {
    if (is.null(model[i])) {
      model[i] <- "sGARCH"
    } else {
      if (!is.character(model[i])) {
        stop(paste0("\ncreate.spec-->error: Model #", i, "
                    argument must be a character.\n",
          call. = FALSE))
      }
      if (!any(model[i] == valid.model)) {
        stop(paste0("\ncreate.spec-->error: Model #", i, "
                    does not appear to be a valid choice.\n",
          call. = FALSE))
      }
    }
  }
  if (is.null(do.mix)) {
    do.mix <- FALSE
  }
  if (isTRUE(do.mix) || !isTRUE(do.mix)) {
  } else {
    stop("\ncreate.spec-->error: do.mix must be a TRUE or FALSE\n",
         call. = FALSE)
  }
  if (is.null(do.shape.ind)) {
    do.shape.ind <- FALSE
  }
  if (isTRUE(do.shape.ind) || !isTRUE(do.shape.ind)) {
  } else {
    stop("\ncreate.spec-->error: do.shape.ind must be a TRUE or FALSE\n",
         call. = FALSE)
  }
  models.merge <- paste0(model, "_", dist.merge)
  models.list <- NULL
  for (j in 1:length(models.merge)) {
    models.list[[j]] <- get(models.merge[j])
  }
  # KB: use loop to ensure that spec is correctly created
  uncvol = NA
  while (any(is.na(uncvol))){
    out <- suppressWarnings(expr = f.spec(models = models.list,
                                          do.mix = do.mix,
                                          do.shape.ind = do.shape.ind))
    
    class(out) <- "MSGARCH_SPEC"
    uncvol = MSGARCH::unc.vol(out, out$theta0)
  }
  
  test.pred = NA
  test.pit  = NA
  test.uc   = NA
  test.all  = is.na(c(test.pred, test.pit, test.uc))
  test.k    = 1
  max.k     = 10
  y0 = c(-0.46, 0.42, 0.85, 0.83, 2.10, -0.47, 0.26, 0.52, -0.76, -1.51)
  while (any(test.all)) {
    out <- suppressWarnings(expr = f.spec(models = models.list,
                                          do.mix = do.mix,
                                          do.shape.ind = do.shape.ind))
    class(out) <- "MSGARCH_SPEC"
    test.pred = MSGARCH::pred(object = out, x = 0, theta = out$theta0, y = y0, log = FALSE, do.its = FALSE)
    test.pit  = MSGARCH::pit(object = out, x = 0, theta = out$theta0, y = y0, do.norm = FALSE, do.its = FALSE)
    test.uc   = MSGARCH::unc.vol(object = out, theta = out$theta0)
    test.all  = is.na(c(test.pred, test.pit, test.uc))
    test.k    = test.k + 1 
    if (test.k > max.k) {
      break
    }
  }
  
  return(out)
}