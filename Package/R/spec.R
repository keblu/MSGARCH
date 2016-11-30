#Function used to build the specification
#' @importFrom methods new
#' @import stringr
f.spec <- function(models, do.mix = FALSE, do.shape.ind = FALSE) {
  K <- length(models)  # number of models
  if (K == 1) {
    do.mix <- FALSE
    do.shape.ind <- FALSE
  }
  options(warn = -1)
  if (K > 1) {
    tmp <- list()
    for (i in 1:K) {
      tmp[i] <- new(models[[i]])
    }
  } else {
    tmp <- new(models[[1]])
  }
  options(warn = 0)
  if (K > 1) {
    mod <- new(MSgarch, tmp)
  } else {
    mod <- tmp
  }
  dist <- NULL
  name <- mod$name
  for (i in 1:length(name)) {
    # DA consider replacing with simple code without "stringr"
    dist[i] <- stringr::str_sub(name[i], start = stringr::str_locate(name, "_")[i, 1] + 1, nchar(name[i]))
  }
  uniqueDist <- unique(dist)
  if (isTRUE(do.shape.ind) && length(uniqueDist) > 1) {
    stop("The distribution of each regime must be the same if the distribution are not regime dependent")
  }
  n.params     <- mod$NbParams
  n.params.vol <- mod$NbParamsModel
  rcpp.func <- list()
  rcpp.func$calc_ht      <- mod$calc_ht
  rcpp.func$eval_model   <- mod$eval_model
  ineq_func.base         <- mod$ineq_func
  rcpp.func$sim          <- mod$f_sim
  rcpp.func$pdf_Rcpp     <- mod$f_pdf
  rcpp.func$cdf_Rcpp     <- mod$f_cdf
  rcpp.func$rnd_Rcpp     <- mod$f_rnd
  rcpp.func$pdf_Rcpp_its <- mod$f_pdf_its
  rcpp.func$ineq_func    <- mod$ineq_func
  rcpp.func$cdf_Rcpp_its <- mod$f_cdf_its
  rcpp.func$unc_vol_Rcpp <- mod$f_unc_vol
  if (K > 1) {
    rcpp.func$get_Pstate_Rcpp <- mod$f_get_Pstate
  } else {
    rcpp.func$get_Pstate_Rcpp <- function(theta, y, PLast) {
      if (!isTRUE(PLast)) {
        out <- matrix(1, nrow = length(y) + 1, ncol = 1)
      } else {
        out <- matrix(1, nrow = 1, ncol = 1)
      }
      return(out)
    }
  }
  nb_total_params <- sum(n.params)
  func <- list()
  func$f.do.mix <- function(theta) {
    return(f.theta.mixture(K, nb_total_params, theta))
  }
  func$f.do.mix.reverse <- function(theta) {
    return(f.theta.mixture.reverse(K, nb_total_params, theta))
  }
  func$f.do.shape.ind <- function(theta) {
    return(f.theta.RegIndDist(K, n.params, n.params.vol, theta))
  }
  func$f.do.shape.ind.reverse <- function(theta) {
    return(f.theta.RegIndDist.reverse(K, n.params, n.params.vol, theta))
  }
  loc <- c(0, cumsum(n.params))
  for (i in 1:K) {
    mod$label[(loc[i] + 1):loc[i + 1]] <- paste0(mod$label[(loc[i] + 1):loc[i + 1]], "_", i)
  }
  if (isTRUE(do.mix) && !isTRUE(do.shape.ind)) {
    mod$lower <- as.vector(func$f.do.mix.reverse(mod$lower))
    newParamsLength <- length(mod$lower)
    mod$upper  <- as.vector(func$f.do.mix.reverse(mod$upper))
    mod$theta0 <- as.vector(func$f.do.mix.reverse(mod$theta0))
    mod$Sigma0 <- mod$Sigma0[1:newParamsLength]
    mod$label  <- mod$label[1:newParamsLength]
    rcpp.func$ineq_func <- function(theta) {
      theta <- as.vector(f.theta.mixture(K, nb_total_params, theta))
      return(ineq_func.base(theta))
    }
  } else if (isTRUE(do.mix) && isTRUE(do.shape.ind)) {
    mod$lower <- as.vector(func$f.do.mix.reverse(mod$lower))
    mod$lower <- as.vector(func$f.do.shape.ind.reverse(mod$lower))
    newParamsLength <- length(mod$lower)
    mod$upper  <- as.vector(func$f.do.mix.reverse(mod$upper))
    mod$upper  <- as.vector(func$f.do.shape.ind.reverse(mod$upper))
    mod$theta0 <- as.vector(func$f.do.mix.reverse(mod$theta0))
    mod$theta0 <- as.vector(func$f.do.shape.ind.reverse(mod$theta0))
    mod$Sigma0 <- mod$Sigma0[1:newParamsLength]
    mod$label  <- func$f.do.shape.ind.reverse(mod$label)
    mod$label  <- mod$label[1:newParamsLength]
    rcpp.func$ineq_func <- function(theta) {
      theta <- as.vector(f.theta.mixture(K, nb_total_params, theta))
      theta <- as.vector(f.theta.RegIndDist(K, n.params, n.params.vol, theta))
      return(ineq_func.base(theta))
    }
  } else if (!isTRUE(do.mix) && isTRUE(do.shape.ind)) {
    mod$lower <- as.vector(func$f.do.shape.ind.reverse(mod$lower))
    newParamsLength <- length(mod$lower)
    mod$upper  <- as.vector(func$f.do.shape.ind.reverse(mod$upper))
    mod$theta0 <- as.vector(func$f.do.shape.ind.reverse(mod$theta0))
    mod$Sigma0 <- mod$Sigma0[1:newParamsLength]
    mod$label  <- func$f.do.shape.ind.reverse(mod$label)
    mod$label  <- mod$label[1:newParamsLength]
    rcpp.func$ineq_func <- function(theta) {
      theta <- as.vector(f.theta.RegIndDist(K, n.params, n.params.vol, theta))
      return(ineq_func.base(theta))
    }
  }
  mod$theta0 <- matrix(mod$theta0, ncol = length(mod$theta0))
  colnames(mod$theta0) <- mod$label
  out <- list(theta0 = mod$theta0, is.mix = do.mix, is.shape.ind = do.shape.ind,
    K = K, sigma0 = diag(mod$Sigma0), lower = mod$lower, upper = mod$upper, ineqlb = mod$ineq_lb,
    inequb = mod$ineq_ub, n.params = n.params, n.params.vol = n.params.vol, do.init = F,
    label = mod$label, name = mod$name, rcpp.func = rcpp.func, func = func)
  return(out)
}