#' @importFrom stats pnorm
f_InferenceFun <- function(vPw, data, spec, do.plm, mNegHessian = NULL) {
  spec <- f_check_spec(spec)
  out  <- matrix(data = NA, nrow = length(vPw), ncol = 4L, dimnames = list(names(vPw), c("Estimate", "Std. Error",
                                                                                         "t value", "Pr(>|t|)")))
  #correction if the parameter is on the boundary
  vPn <- f_mapPar(vPw, spec, do.plm)
  vPn_mod = vPn
  vLower = spec$lower
  vUpper = spec$upper
  
  names(vLower) = names(spec$par0)
  names(vUpper) = names(spec$par0)
  
  vPn_mod[vPn >= vUpper[names(vPn)]] = vPn[vPn >= vUpper[names(vPn)]] - 0.000001
  vPn_mod[vPn <= vLower[names(vPn)]] = vPn[vPn <= vLower[names(vPn)]] + 0.000001
  vPw_mod = f_unmapPar(vPn_mod, spec, do.plm)
  ###
  
  if (is.null(mNegHessian)) {
    mNegHessian <- stats::optim(par = vPw_mod, fn = f_nll, data = data,
                 spec = spec, do.plm = do.plm, method = "BFGS",
                 control = list(maxit = 1L),
                 hessian = TRUE)$hessian
    
  }
  
  dMinEigen <- min(eigen(mNegHessian)$values)
  
  if (dMinEigen < .Machine$double.eps) {
    
    mNegHessian <- mNegHessian + abs(dMinEigen) * 1.1 * diag(length(vPw_mod))
    
  }
  
  mJacob      <- numDeriv::jacobian(f_mapPar, vPw_mod, spec = spec, do.plm = do.plm)
  mInvHessian <- MASS::ginv(mNegHessian)
  mSandwitch  <- t(mJacob) %*% mInvHessian %*% mJacob
  
  vSE   <- sqrt(diag(mSandwitch))
  vTest <- vPn/vSE
  
  vPvalues <- 1 - pnorm(abs(vTest))
  
  out[, "Estimate"]   <- vPn
  out[, "Std. Error"] <- vSE
  out[, "t value"]    <- vTest
  out[, "Pr(>|t|)"]   <- vPvalues
  
  out <- list(MatCoef = out, Hessian = mNegHessian)
  return(out)
}
