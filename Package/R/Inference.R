f_InferenceFun <- function(vPw, data, spec, do.plm, mNegHessian = NULL) {
  spec <- f_check_spec(spec)
  out  <- matrix(data = NA, nrow = length(vPw), ncol = 4L, dimnames = list(names(vPw), c("Estimate", "Std. Error",
                                                                                         "t value", "Pr(>|t|)")))
  
  vPn <- f_mapPar(vPw, spec, do.plm)
  
  if (is.null(mNegHessian)) {
    
    mNegHessian <- numDeriv::hessian(f_nll, vPw, data = data, spec = spec, do.plm = do.plm, method.args = list(d = 1e-05))
    
  }
  
  dMinEigen <- min(eigen(mNegHessian)$values)
  
  if (dMinEigen < .Machine$double.eps) {
    
    mNegHessian <- stats::optim(par = vPw, fn = f_nll, data = data, spec = spec, do.plm = do.plm, control = list(maxit = 1L),
                                hessian = TRUE)$hessian
    
  }
  
  dMinEigen <- min(eigen(mNegHessian)$values)
  
  if (dMinEigen < .Machine$double.eps) {
    
    mNegHessian <- mNegHessian + abs(dMinEigen) * 1.1 * diag(length(vPw))
    
  }
  
  mJacob      <- numDeriv::jacobian(f_mapPar, vPw, spec = spec, do.plm = do.plm)
  mInvHessian <- MASS::ginv(mNegHessian)
  mSandwitch  <- t(mJacob) %*% mInvHessian %*% mJacob
  
  vSE   <- sqrt(diag(mSandwitch))
  vTest <- vPn/(vSE/sqrt(length(data)))
  
  vPvalues <- 1 - pnorm(abs(vTest))
  
  out[, "Estimate"]   <- vPn
  out[, "Std. Error"] <- vSE/sqrt(length(data))
  out[, "t value"]    <- vTest
  out[, "Pr(>|t|)"]   <- vPvalues
  
  out <- list(MatCoef = out, Hessian = mNegHessian)
  return(out)
}
