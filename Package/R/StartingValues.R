f_DistParNames <- function(sDist, bSkew) {

  if (bSkew) {
    if (sDist == "norm") {
      vNames <- c("xi")
    }
    if (sDist == "std") {
      vNames <- c("nu", "xi")
    }
    if (sDist == "ged") {
      vNames <- c("nu", "xi")
    }
  } else {
    if (sDist == "std") {
      vNames <- c("nu")
    }
    if (sDist == "ged") {
      vNames <- c("nu")
    }
  }
  return(vNames)
}

f_ModelParNames <- function(sModel) {

  if (sModel == "sARCH") {
    vNames <- c("alpha0", "alpha1")
  }
  if (sModel == "sGARCH") {
    vNames <- c("alpha0", "alpha1", "beta")
  }
  if (sModel %in% c("eGARCH", "gjrGARCH", "tGARCH")) {
    vNames <- c("alpha0", "alpha1", "alpha2", "beta")
  }

  return(vNames)
}

f_SwitchDistLabels <- function(vLabels) {

  vSwithLabels <- vLabels

  for (i in 1:length(vLabels)) {
    sLabel <- vLabels[i]
    if (sLabel == "normal") {
      vSwithLabels[i] <- "norm"
    }
    if (sLabel == "student") {
      vSwithLabels[i] <- "std"
    }
  }

  return(vSwithLabels)
}

f_StaticStarting_Uni <- function(sDist, bSkew) {

  if (bSkew) {
    if (sDist == "norm") {
      vpar <- c(xi = 1)
    } else {
      vpar <- c(nu = 7, xi = 1)
    }
  } else {
    vpar <- c(nu = 7)
  }

  vpar_tilde <- as.numeric(UnmapParameters_univ(vpar, sDist, bSkew))

  return(vpar_tilde)
}

f_Fit_StaticDist <- function(vZ, sDist, bSkew) {

  vpar_tilde <- f_StaticStarting_Uni(sDist, bSkew)

  iT <- length(vZ)

  vNames <- f_DistParNames(sDist, bSkew)

  optimizer <- optim(vpar_tilde, function(vpar_tilde, vZ, iT, sDist, bSkew, vNames) {

    vpar <- as.numeric(MapParameters_univ(vpar_tilde, sDist, bSkew))
    names(vpar) <- vNames

    dLLK <- dUnivLike(vZ, sDist, bSkew, dXi = vpar["xi"], dNu = vpar["nu"])

    return(-dLLK)

  }, vZ = vZ, iT = iT, sDist = sDist, bSkew = bSkew, vNames = vNames, method = "BFGS")

  vpar <- as.numeric(MapParameters_univ(optimizer$par, sDist, bSkew))

  names(vpar) <- f_DistParNames(sDist, bSkew)

  if (sDist == "std") {
    if (vpar["nu"] > 10) {
      vpar["nu"] = 10
    }
  }

  return(vpar)

}

f_VarianceTargeting <- function(dSigma2, sModel, vpar) {

  if (sModel == "sARCH") {
    dAlpha0 <- dSigma2 * (1 - vpar["alpha1_1"])
  }
  if (sModel == "sGARCH") {
    dAlpha0 <- dSigma2 * (1 - vpar["alpha1_1"] - vpar["beta_1"])
  }

  if (sModel == "gjrGARCH") {
    dAlpha0 <- dSigma2 * (1 - vpar["alpha1_1"] - 0.5 * vpar["alpha2_1"] - vpar["beta_1"])
  }

  if (sModel == "eGARCH") {
    dAlpha0 <- log(dSigma2) * (1 - vpar["beta_1"])
  }

  if (sModel == "tGARCH") {
    dAlpha0 <- dSigma2 * (1 + (vpar["alpha1_1"] + vpar["alpha2_1"]) * 0.5 - vpar["beta_1"])
  }

  return(dAlpha0)

}

f_StartingValueMSGARCH <- function(y, spec, ctr = NULL) {
  spec = f_check_spec(spec)
  K <- spec$K
  vSpec <- spec$name
  vModel <- sapply(vSpec, function(x) unlist(strsplit(x, split = "_"))[1L])
  vDist.or <- sapply(vSpec, function(x) unlist(strsplit(x, split = "_"))[2L])
  vSkew <- sapply(vDist.or, FUN = function(x) any(c("snorm", "sstd", "sged") == x))
  vDist <- vDist.or
  vDist[vSkew] <- substring(vDist[vSkew], 2)
  names(vSpec) <- NULL
  names(vModel) <- NULL
  names(vDist) <- NULL
  names(vSkew) <- NULL

  do.mix <- spec$is.mix

  ## Do EM
  if (do.mix) {
    EM_Fit <- EM_MM(y, K, constraintZero = TRUE)
    vDecoding <- EM_Fit$vDecoding + 1L
  } else {
    EM_Fit <- EM_HMM(y, K, constraintZero = TRUE)
    vDecoding <- EM_Fit$vDecoding + 1L
  }

  ## local decoding
  lY <- list()
  for (k in 1:K) {
    vDec_foo <- which(vDecoding == k)
    if (length(vDec_foo) > 100L) {
      lY[[k]] <- y[vDec_foo]
    } else {
      lY[[k]] <- y
    }
  }

  ## initialize transition probability matrix
  if (do.mix) {
    vP <- EM_Fit$vP[1:(K - 1L)]
    names(vP) <- paste("P", 1:(K - 1), sep = "_")
  } else {
    StartingGamma <- EM_Fit$mGamma
    vP <- c(t(StartingGamma[, -K]))
    names(vP) <- tail(spec$label, K * (K - 1))
  }


  ## Initialize unconditional volatilties
  vSigma2 <- EM_Fit$vSigma2

  ## initiale shape parameters
  lShape <- list()


  for (k in 1:K) {

    lShape[[k]] <- NULL

    if (vDist[k] != "norm") {

      vZ <- (lY[[k]] - mean(lY[[k]]))/sd(lY[[k]])
      lShape[[k]] <- f_Fit_StaticDist(vZ, vDist[k], FALSE)
    }
  }

  ## initialize shape and skew
  lSingleRegimeSpec <- list()
  lSingleRegimeCoef <- list()

  ## Fixed Parameters
  lFixed_SR <- f_recover_fixedpar_SR(spec)

  for (k in 1:K) {
    lSingleRegimeSpec[[k]] <- CreateSpec(variance.spec = list(model = vModel[k]),
                                         distribution.spec = list(distribution = vDist.or[k]),
                                         constraint.spec = list(fixed = lFixed_SR[[k]]))

    if (vDist[k] != "norm") {
      lSingleRegimeSpec[[k]]$par0["nu_1"] <- lShape[[k]]
    }

    dAlpha0 <- f_VarianceTargeting(vSigma2[k], vModel[k], lSingleRegimeSpec[[k]]$par0)

    lSingleRegimeSpec[[k]]$par0["alpha0_1"] <- dAlpha0
    Fit <- MSGARCH::FitML(spec = lSingleRegimeSpec[[k]], data = lY[[k]], ctr = ctr)

    lSingleRegimeCoef[[k]] <- Fit$par

    names(lSingleRegimeCoef[[k]]) <- sapply(names(lSingleRegimeCoef[[k]]), function(x) {
      unlist(strsplit(x, split = "_"))[1L]
    })

    names(lSingleRegimeCoef[[k]]) <- paste(names(lSingleRegimeCoef[[k]]), k, sep = "_")

  }

  vpar0 <- do.call(c, lSingleRegimeCoef)

  vpar0 <- c(vpar0, vP)

  return(vpar0)
}

f_StargingValues <- function(y, spec, ctr = NULL) {
  spec <- f_check_spec(spec)
  K <- spec$K
  if (K > 1L) {

    vPn <- f_StartingValueMSGARCH(y, spec, ctr)
    if (is(vPn, "try-error")) {
      vPn <- spec$par0
    }

    ### Regime Constant Parameters ###
    if (isTRUE(spec$regime.const.pars.bool)) {
      vPn <- f_remove_regimeconstpar(vPn, spec$regime.const.pars, K)
    }

    if (isTRUE(spec$fixed.pars.bool)) {
      vPn <- f_remove_fixedpar(vPn, spec$fixed.pars)
    }
    
    vPw <- f_unmapPar(vPn, spec, ctr$do.plm)

    dLLK <- f_nll(vPw, y, spec, ctr$do.plm)

    if (dLLK == 1e+10) {

      vPn <- spec$par0

      ### Fixed Parameters ###
      if (isTRUE(spec$fixed.pars.bool)) {
        vPn <- f_substitute_fixedpar(vPn, spec$fixed.pars)
      }

      ### Regime Constant Parameters ###
      if (isTRUE(spec$regime.const.pars.bool)) {
        vPn <- f_remove_regimeconstpar(vPn, spec$regime.const.pars, K)
      }
      
      if (isTRUE(spec$fixed.pars.bool)) {
        vPn <- f_remove_fixedpar(vPn, spec$fixed.pars)
      }
      
      vPw <- f_unmapPar(vPn, spec, ctr$do.plm)

    }

  } else {

    vPn <- spec$par0

    ### Fixed Parameters ###
    if (isTRUE(spec$fixed.pars.bool)) {
      vPn <- f_substitute_fixedpar(vPn, spec$fixed.pars)
    }

    sModel <- unlist(strsplit(spec$name, split = "_"))[1]

    vPn["alpha0_1"] <- f_VarianceTargeting(var(c(y)), sModel, vPn)

    if (isTRUE(spec$fixed.pars.bool)) {
      vPn <- f_remove_fixedpar(vPn, spec$fixed.pars)
    }
    vPw <- f_unmapPar(vPn, spec, ctr$do.plm)

  }

  return(vPw)

}
