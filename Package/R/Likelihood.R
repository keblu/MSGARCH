f_nll <- function(vPw, data, spec, do.plm) {

  if (is.null(names(vPw))) {
    vPw <- f_rename_par(vPw, spec)
  }

  vPn <- f_mapPar(vPw, spec, do.plm)

  if (isTRUE(spec$fixed.pars.bool)) {
    vPn <- f_add_fixedpar(vPn, spec$fixed.pars)
    vPn <- vPn[spec$label]
  }

  if (isTRUE(spec$regime.const.pars.bool)) {
    vPn <- f_add_regimeconstpar(vPn, spec$K, spec$label)
  }

  # the working-to-natural map can overflow for extreme trial values; treat that
  # as an infeasible point rather than letting the strict parameter check throw
  if (anyNA(vPn) || any(!is.finite(vPn))) {
    return(1e+10)
  }

  dLLK <- Kernel(spec, vPn, data, log = TRUE, do.prior = FALSE)

  if (!is.finite(dLLK)) {
    dLLK <- -1e+10
  }

  return(-dLLK)
}

#' @importFrom stats logLik
#' @export
logLik.MSGARCH_ML_FIT <- function(object, ...){
  out = structure(object$loglik, df = dofMSGARCH(object), 
                  nobs = length(object$data))
  return(out)
}