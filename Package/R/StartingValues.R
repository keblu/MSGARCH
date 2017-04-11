
DistParNames <- function(sDist, bSkew) {

  if (bSkew) {
    if (sDist == "norm") {
      vNames = c("xi")
    }
    if (sDist == "std") {
      vNames = c("nu", "xi")
    }
    if (sDist == "ged") {
      vNames = c("nu", "xi")
    }
  } else {
    if (sDist == "std") {
      vNames = c("nu")
    }
    if (sDist == "ged") {
      vNames = c("nu")
    }
  }
  return(vNames)
}

SwitchDistLabels <- function(vLabels) {

  vSwithLabels = vLabels

  for (i in 1:length(vLabels)) {
    sLabel = vLabels[i]
    if (sLabel == "normal") {
      vSwithLabels[i] = "norm"
    }
    if (sLabel == "student") {
      vSwithLabels[i] = "std"
    }
  }

  return(vSwithLabels)
}

StaticStarting_Uni <- function(sDist, bSkew) {

  if (bSkew) {
    if (sDist == "norm") {
      vTheta = c(xi = 1.0)
    } else {
      vTheta = c(nu = 7.0, xi = 1.0)
    }
  } else {
    vTheta = c(nu = 7.0)
  }

  vTheta_tilde = as.numeric(UnmapParameters_univ(vTheta, sDist, bSkew))

  return(vTheta_tilde)
}

ddist_ThetaParam <- function(dZ, vTheta, sDist, bSkew, bLog) {

  dPDF = 0.0

  if (!bSkew) {
    if (sDist == "norm") {
      dPDF = ddist(dZ, "norm", log = bLog)
    } else {
      dPDF = ddist(dZ, sDist, shape = vTheta[1], log = bLog)
    }
  } else {
    if (sDist == "norm") {
      dPDF = ddist(dZ, "norm", skew = vTheta[1], log = bLog)
    } else {
      dPDF = ddist(dZ, sDist, shape = vTheta[1], skew = vTheta[2], log = bLog)
    }
  }

  if (bLog) {
    if (dPDF < -1e50) {
      dPDF = -1e50
    }
  } else {
    if (dPDF < 1e-50) {
      dPDF = 1e-50
    }
  }

  return(dPDF)

}

pdist_ThetaParam <- function(dQ, vTheta, sDist, bSkew, bLog) {

  dP = 0.0

  if (!bSkew) {
    if (sDist == "norm") {
      dP = pdist(dQ, "norm")
    } else {
      dP = pdist(dQ, sDist, shape = vTheta[1])
    }
  } else {
    if (sDist == "norm") {
      dP = pdist(dQ, "norm", skew = vTheta[1])
    } else {
      dP = pdist(dQ, sDist, shape = vTheta[1], skew = vTheta[2])
    }
  }

    if (dP < 1e-50) {
      dP = 1e-50
    }

  return(dP)

}

Fit_StaticDist <- function(vZ, sDist, bSkew) {

  vTheta_tilde = StaticStarting_Uni(sDist, bSkew)

  iT = length(vZ)

  vNames = DistParNames(sDist, bSkew)

  optimizer = optim(vTheta_tilde, function(vTheta_tilde, vZ, iT, sDist, bSkew, vNames){

    vTheta = as.numeric(MapParameters_univ(vTheta_tilde, sDist, bSkew))
    names(vTheta) = vNames

    dLLK = dUnivLike(vZ, sDist, bSkew, dXi = vTheta["xi"], dNu = vTheta["nu"])

    return(-dLLK)

  }, vZ = vZ, iT = iT, sDist = sDist, bSkew = bSkew, vNames = vNames, method = "BFGS")

  vTheta =  as.numeric(MapParameters_univ(optimizer$par, sDist, bSkew))

  names(vTheta) = DistParNames(sDist, bSkew)

  return(vTheta)

}

VarianceTargeting <- function(dSigma2, sModel, vTheta) {

  if (sModel == "sGARCH") {
    dAlpha0 = dSigma2 * (1.0 - vTheta[1, "alpha1_1"] - vTheta[1, "beta_1"])
  }

  if (sModel == "gjrGARCH") {
    dAlpha0 = dSigma2 * (1.0 - vTheta[1, "alpha1_1"] - 0.5 * vTheta[1, "alpha2_1"] - vTheta[1, "beta_1"])
  }

  if (sModel == "eGARCH") {
    dAlpha0 = log(dSigma2) * (1.0 - vTheta[1, "beta_1"])
  }

  if (sModel == "tGARCH") {
    dAlpha0 = dSigma2 * (1.0 + (vTheta[1, "alpha1_1"] + vTheta[1, "alpha2_1"]) * 0.5 - vTheta[1, "beta_1"])
  }

  return(dAlpha0)

}

StartingValueMSGARCH <- function(y, spec, optim.fun = NULL) {

  K      = spec$K
  vSpec  = spec$name
  vModel = sapply(vSpec, function(x) unlist(strsplit(x, split = "_"))[1])
  vDist  = SwitchDistLabels(sapply(vSpec, function(x) unlist(strsplit(x, split = "_"))[2]))
  vSkew   = sapply(vSpec, function(x) unlist(strsplit(x, split = "_"))[3]) == "skew"

  names(vSpec)  = NULL
  names(vModel) = NULL
  names(vDist)  = NULL
  names(vSkew)  = NULL

  do.mix = spec$is.mix
  do.shape.ind = spec$is.shape.ind

  ## Do EM
  if (do.mix) {
    EM_Fit = EM_MM(y, K, constraintZero = TRUE)
    vDecoding = EM_Fit$vDecoding + 1
  } else {
    EM_Fit = EM_HMM(y, K, constraintZero = TRUE)
    vDecoding = EM_Fit$vDecoding + 1
  }

  ## local decoding
  lY = list()
  for (k in 1:K) {
    vDec_foo = which(vDecoding == k)
    if (length(vDec_foo) > 300) {
      lY[[k]] = y[vDec_foo]
    } else {
      lY[[k]] = y[vDec_foo]
    }
  }

  ## initialize transition probability matrix
  if (do.mix) {
    vP = matrix(EM_Fit$vP[1:(K - 1)], nrow = 1)
  } else {
    StartingGamma = t(EM_Fit$mGamma)
    vP = matrix(c(StartingGamma[-K, ]), nrow = 1)
  }
  colnames(vP) = rep("P", length(vP))

  ## Initialize unconditional volatilties
  vSigma2 = EM_Fit$vSigma2

  ## initiale shape parameters
  lShape = list()

  if (do.shape.ind) {

    vZ = (y - mean(y))/sd(y)

    lShape = Fit_StaticDist(vZ, vDist[1], FALSE)

  } else {

    for (k in 1:K) {

      lShape[[k]] = NULL

      if (vDist[k] != "norm") {

        vZ = (lY[[k]] - mean(lY[[k]]))/sd(lY[[k]])
        lShape[[k]] = Fit_StaticDist(vZ, vDist[k], FALSE)

      }
    }
  }

  ## initialize shape and skew

  lSingleRegimeSpec = list()
  lSingleRegimeCoef = list()

  for (k in 1:K) {

    if (do.shape.ind) {

      lSingleRegimeSpec[[k]] = create.spec(model = vModel[k], distribution = "norm",
                                           do.skew = FALSE, do.mix = FALSE,
                                           do.shape.ind = FALSE)

    } else {

      lSingleRegimeSpec[[k]] = create.spec(model = vModel[k], distribution = vDist[k],
                                           do.skew = vSkew[k], do.mix = FALSE,
                                           do.shape.ind = FALSE)

      if (vDist[k] != "norm") {
        lSingleRegimeSpec[[k]]$theta0[1, "nu_1"] = lShape[[k]]
      }

    }

    dAlpha0 = VarianceTargeting(vSigma2[k], vModel[k], lSingleRegimeSpec[[k]]$theta0)

    lSingleRegimeSpec[[k]]$theta0[1, "alpha0_1"] = dAlpha0

    Fit = MSGARCH::fit.mle(optim.fun = optim.fun, spec = lSingleRegimeSpec[[k]], y = lY[[k]])
    lSingleRegimeCoef[[k]] = Fit$theta

    colnames(lSingleRegimeCoef[[k]]) = sapply(colnames(lSingleRegimeCoef[[k]]), function(x) {
      unlist(strsplit(x, split = "_"))[1]
    })

    colnames(lSingleRegimeCoef[[k]]) = paste(colnames(lSingleRegimeCoef[[k]]), k, sep = "_")

  }

  vTheta0 = do.call(cbind, lSingleRegimeCoef)

  if (do.shape.ind) {

    if (vDist[1] != "norm") {
      vTheta0 = cbind(vTheta0, nu_1 = lShape)
      rownames(vTheta0) = NULL
    }

    if (vSkew[1]) {
      vTheta0 = cbind(vTheta0, xi_1 = 1.0)
      rownames(vTheta0) = NULL
    }

  }

  vTheta0 = cbind(vTheta0, vP)

  return(vTheta0)

}


