#################################################### fixed.pars ####

f_check_parameterConstraints <- function(fixed.pars, vParNames) {
  
  if (any(!names(fixed.pars) %in% vParNames)) {
    vWrongPars <- names(fixed.pars)[!names(fixed.pars) %in% vParNames]
    stop(cat(paste("Wrong name in fixed.pars:", vWrongPars)))
  }
  
  return(fixed.pars)
  
}

f_remove_fixedpar <- function(vPar, fixed.pars) {
  vPar <- vPar[!names(vPar) %in% names(fixed.pars)]
  return(vPar)
}

f_substitute_fixedpar <- function(vPar, fixed.pars) {
  vFixed = names(vPar) %in% names(fixed.pars)
  if (any(vFixed)) {
    vPar[vFixed] <- unlist(fixed.pars[names(vPar)[vFixed]])
  }
  return(vPar)
}

f_add_fixedpar <- function(vPar, fixed.pars) {
  
  vPar <- c(vPar, unlist(fixed.pars))
  
  return(vPar)
}

f_recover_fixedpar_SR <- function(spec) {
  K <- spec$K
  lFixed_SR <- list()
  
  if (isTRUE(spec$fixed.pars.bool)) {
    vFixdPar = names(spec$fixed.pars)
    for (k in 1:K) {
      vFixdPar_SR <- which(gsub(paste("_", k, sep = ""), "", vFixdPar) != vFixdPar)
      if (length(vFixdPar_SR) > 0) {
        lFixed_SR[[k]] <- spec$fixed.pars[vFixdPar_SR]
        names(lFixed_SR[[k]]) <- gsub(paste("_", k, sep = ""), "_1", names(lFixed_SR[[k]]), fixed = TRUE)
      } else {
        lFixed_SR[[k]] <- list()
      }
    }
  } else {
    for (k in 1:K) {
      lFixed_SR[[k]] <- list()
    }
  }
  
  return(lFixed_SR)
}

f_check_fixedpars <- function(fixed.pars, spec) {
  
  vLower <- spec$lower
  vUpper <- spec$upper
  
  names(vLower) <- names(spec$par0)
  names(vUpper) <- names(spec$par0)
  
  vLower = vLower[names(fixed.pars)]
  vUpper = vUpper[names(fixed.pars)]
  
  vFixedPar <- unlist(fixed.pars)
  
  vLowerCheck <- vFixedPar < vLower
  vUpperCheck <- vFixedPar > vUpper
  
  if (any(vLowerCheck)) {
    stop(paste("The following fixed parameters are below their lower bound:\n", vFixedPar[vLowerCheck], "\nThe lower bounds are:", vLower[vLowerCheck]))
  }
  if (any(vUpperCheck)) {
    stop(paste("The following fixed parameters are above their upper bound:\n", vFixedPar[vUpperCheck], "\nThe lower bounds are:", vUpper[vUpperCheck]))
  }
  
}

#################################################### regime.const.pars ####

f_check_regime_const_pars <- function(regime.const.pars, K, vParNames, model, distribution) {
  
  if (K > 1L) {
    
    regime.const.pars_original <- regime.const.pars
    distribution_original <- distribution
    
    regime.const.pars <- sapply(regime.const.pars, function(x) unlist(strsplit(x, split = "_"))[1L])
    names(regime.const.pars) <- NULL
    
    if ((any(regime.const.pars %in% c("alpha0", "alpha1", "beta"))) && (length(unique(model)) != 1)) {
      vTest = regime.const.pars_original[regime.const.pars %in% c("alpha0", "alpha1", "beta")]
      warning(paste("The following regime constant parameters are defined:\n", vTest, "\nBut the variance specification is different across the regimes."), call. = FALSE)
    }
    
    if ((any(regime.const.pars %in% c("nu", "xi"))) && (length(unique(distribution)) != 1)) {
      vTest = regime.const.pars_original[regime.const.pars %in% c("nu", "xi")]
      warning(paste("The following regime constant parameters are defined:\n", vTest, "\nBut the conditional distribution specification is different across the regimes."), 
              call. = FALSE)
    }
    
    vParNames_Unique <- unique(sapply(vParNames, function(x) unlist(strsplit(x, split = "_"))[1L]))
    names(vParNames_Unique) <- NULL
    
    vTest = !regime.const.pars %in% vParNames_Unique
    
    if (any(vTest)) {
      stop(paste("The following regime constant parameters are wrongly defined:\n", regime.const.pars_original[vTest]))
    }
    
  } else {
    if (!is.null(regime.const.pars)) {
      stop("The argument regime.const.pars is not NULL, but the model is not Markov Switching.")
    }
    regime.const.pars <- NULL
  }
  return(regime.const.pars)
}

f_remove_regimeconstpar <- function(vPar, regime.const.pars, K, for.se = FALSE) {
  
  if(isTRUE(for.se)){
    FUN <- function(x){return(x)}
  } else {
    FUN <- function(x){return(abs(x))}
  }
  
  for (k in 1:(K - 1)) {
    for (sFixParName in regime.const.pars) {
      
      vTest_k <- names(vPar) == paste(sFixParName, k, sep = "_")
      vTest_kp12K <- names(vPar) %in% paste(sFixParName, (k + 1):(K), sep = "_")
      
      if (any(vTest_k) && any(vTest_kp12K)) {
        vPar[vTest_k] <- min(FUN(vPar[c(which(vTest_k), which(vTest_kp12K))]))  # for eGARCH we put abs
        vPar = vPar[!vTest_kp12K]
      }
    }
  }
  
  return(vPar)
}

f_add_regimeconstpar <- function(vPar, K, vFullParNames) {
  
  vPar <- vPar[vFullParNames]
  names(vPar) <- vFullParNames
  
  for (i in 1:length(vPar)) {
    if (is.na(vPar[i])) {
      foo <- unlist(strsplit(vFullParNames[i], split = "_"))
      sParName <- foo[1]
      k_foo <- as.numeric(foo[2])
      l = 1
      while (is.na(vPar[i])) {
        vPar[i] = vPar[paste(sParName, k_foo - l, sep = "_")]
        l = l + 1
      }
    }
  }
  
  return(vPar)
}

f_add_regimeconstpar_matrix <- function(mPar, K, vFullParNames) {
  
  vAddNames <- vFullParNames[!vFullParNames %in% colnames(mPar)]
  mAddPar <- matrix(rep(NA, length(vAddNames) * nrow(mPar)), nrow = nrow(mPar))
  colnames(mAddPar) <- vAddNames
  
  mPar <- cbind(mPar, mAddPar)
  mPar <- mPar[, vFullParNames]
  
  for (i in 1:ncol(mPar)) {
    if (is.na(mPar[1L, i])) {
      foo <- unlist(strsplit(vFullParNames[i], split = "_"))
      sParName <- foo[1]
      k_foo <- as.numeric(foo[2])
      l = 1
      while (is.na(mPar[1L, i])) {
        sName_foo = paste(sParName, k_foo - l, sep = "_")
        if (sName_foo %in% vFullParNames) {
          mPar[, i] = mPar[, paste(sParName, k_foo - l, sep = "_")]
        }
        l = l + 1
      }
    }
  }
  
  return(mPar)
}