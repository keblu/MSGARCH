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

  dLLK <- Kernel(spec, vPn, data, log = TRUE, do.prior = FALSE)

  if (!is.finite(dLLK)) {
    dLLK <- -1e+10
  }

  return(-dLLK)
}


logLik.MSGARCH_ML_FIT <- function(object, ...){
  out = structure(object$loglik, df = dofMSGARCH(object), 
                  nobs = length(object$data))
  return(out)
}