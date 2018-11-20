f_posterior <- function(vPw, data, spec, PriorFun) {
  if (is.null(names(vPw))) {
    vPw <- f_rename_par(vPw, spec)
  }
  vPn <- f_mapPar(vPw, spec, TRUE)      
  mJacob <- f_mapJacob(vPw, spec)

  if (isTRUE(spec$fixed.pars.bool)) {
    vPn <- f_add_fixedpar(vPn, spec$fixed.pars)
    vPn <- vPn[names(spec$par0)]
  }

  if (isTRUE(spec$regime.const.pars.bool)) {
    vPn <- f_add_regimeconstpar(vPn, spec$K, spec$label)
  }

  dLLK <- Kernel(spec, vPn, data, log = TRUE, do.prior = TRUE) + sum(log(diag(abs(mJacob))))

  if (!is.finite(dLLK)) {
    dLLK <- -1e10
  }

  return(dLLK)
}

