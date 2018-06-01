###################################### MAPPING Par ####

f_mapPar <- function(vPw, spec, do.plm = FALSE) {
  if (isTRUE(do.plm)) {
    vLower = spec$lower
    vUpper = spec$upper

    names(vLower) = names(spec$par0)
    names(vUpper) = names(spec$par0)

    vPn <- f_map(vPw, vLower[names(vPw)], vUpper[names(vPw)])[1, ]

    names(vPn) = names(vPw)

  } else {
    K = spec$K

    if (K > 1) {
      vPn <- f_mapMS(vPw, spec)
    } else {
      vPn <- f_mapSR(vPw, spec)
    }
  }

  return(vPn)
}

f_mapJacob <- function(vPw, spec) {

  vLower <- spec$lower
  vUpper <- spec$upper

  names(vLower) <- spec$label
  names(vUpper) <- spec$label

  vLower <- vLower[names(vPw)]
  vUpper <- vUpper[names(vPw)]

  vJ <- -vPw + log(vUpper - vLower) - 2 * log(1 + exp(-vPw))
  vJ <- exp(vJ)

  mJ <- diag(as.vector(vJ))

  return(mJ)
}

f_unmapPar <- function(par, spec, do.plm = FALSE) {
  if (isTRUE(do.plm)) {
    vLower = spec$lower
    vUpper = spec$upper

    names(vLower) = names(spec$par0)
    names(vUpper) = names(spec$par0)

    vPw <- f_unmap(par, vLower[names(par)], vUpper[names(par)])[1, ]

    names(vPw) = names(par)

  } else {
    K <- spec$K

    if (K > 1) {
      vPw <- f_unmapMS(par, spec)
    } else {
      vPw <- f_unmapSR(par, spec)
    }
  }

  return(vPw)
}

###################################### MAPPING SR ####

f_pw2pn_SR <- function(vPw, k, sModel, skew = NA, shape = NA, sDist = NULL) {

  if (is.na(skew)) {
    skew <- NULL
  }

  if (is.na(shape)) {
    shape <- NULL
  }

  vPn <- switch(sModel, sARCH = f_pw2pn_sARCH_SR(vPw, k), sGARCH = f_pw2pn_sGARCH_SR(vPw, k), eGARCH = f_pw2pn_eGARCH_SR(vPw, k), gjrGARCH = f_pw2pn_gjrGARCH_SR(vPw, k, skew, shape,
                                                                                                                                                                 sDist), tGARCH = f_pw2pn_tGARCH_SR(vPw, k, skew, shape, sDist))
  return(vPn)

}

f_pn2pw_SR <- function(vPn, k, sModel, skew = NA, shape = NA, sDist = NULL) {

  if (is.na(skew)) {
    skew <- NULL
  }

  if (is.na(shape)) {
    shape <- NULL
  }

  vPw <- switch(sModel, sARCH = f_pn2pw_sARCH_SR(vPn, k), sGARCH = f_pn2pw_sGARCH_SR(vPn, k), eGARCH = f_pn2pw_eGARCH_SR(vPn, k), gjrGARCH = f_pn2pw_gjrGARCH_SR(vPn, k, skew, shape,
                                                                                                                                                                 sDist), tGARCH = f_pn2pw_tGARCH_SR(vPn, k, skew, shape, sDist))
  return(vPw)

}

f_unmapSR <- function(vPn, spec) {
  sModel <- f_getModel(spec)
  sDist <- f_getDist(spec)

  lower <- spec$lower
  upper <- spec$upper

  names(lower) <- names(vPn)
  names(upper) <- names(vPn)

  sSkew <- vPn[paste("xi_1")]  #nas for non skew cond dist
  sShape <- vPn[paste("nu_1")]  #nas for gauss cond dist

  sSkew_tilde <- f_unmap(sSkew, lower["xi_1"], upper["xi_1"])[1, ]
  names(sSkew_tilde) <- names(sSkew)

  sShape_tilde <- f_unmap(sShape, lower["nu_1"], upper["nu_1"])[1, ]
  names(sShape_tilde) <- names(sShape)

  vPw <- vPn
  vPw[] <- NA

  vPw_foo <- f_pn2pw_SR(vPn, 1L, sModel, sSkew, sShape, sDist)

  vPw[names(vPw_foo)] <- vPw_foo

  if (!is.na(sSkew)) {
    vPw["xi_1"] <- sSkew_tilde
  }
  if (!is.na(sShape)) {
    vPw["nu_1"] <- sShape_tilde
  }

  return(vPw)

}

f_mapSR <- function(vPw, spec) {

  sModel <- f_getModel(spec)
  sDist <- f_getDist(spec)
  lower <- spec$lower
  upper <- spec$upper

  names(lower) <- names(vPw)
  names(upper) <- names(vPw)

  sSkew_tilde <- vPw[paste("xi_1")]  #nas for non skew cond dist
  sShape_tilde <- vPw[paste("nu_1")]  #nas for gauss cond dist

  sSkew <- f_map(sSkew_tilde, lower["xi_1"], upper["xi_1"])[1, ]

  sShape <- f_map(sShape_tilde, lower["nu_1"], upper["nu_1"])[1, ]

  vPn <- vPw
  vPn[] <- NA

  vPn_foo <- f_pw2pn_SR(vPw, 1L, sModel, sSkew, sShape, sDist)

  vPn[names(vPn_foo)] <- vPn_foo

  if (!is.na(sSkew)) {
    vPn["xi_1"] <- sSkew
  }
  if (!is.na(sShape)) {
    vPn["nu_1"] <- sShape
  }

  return(vPn)

}

###################################### MAPPING MS ####

f_unmapGamma <- function(vGamma_pn, K) {

  mGamma <- matrix(NA, K, K)

  mGamma[-K, ] <- vGamma_pn
  mGamma       <- t(mGamma)
  mGamma[, K]  <- 1 - apply(mGamma[, -K, drop = FALSE], 1, sum)

  foo          <- log(mGamma/diag(mGamma))
  vGamma_tilde <- as.vector(foo[!diag(K)])

  names(vGamma_tilde) <- names(vGamma_pn)

  return(vGamma_tilde)

}

f_mapGamma <- function(vGamma_tilde, K) {

  vFoo <- exp(vGamma_tilde)

  vFoo[vFoo < 1e-10] <- 1e-10
  vFoo[vFoo > 1e+10] <- 1e+10

  mGamma <- diag(K)
  mGamma[!mGamma] <- vFoo
  mGamma <- mGamma/apply(mGamma, 1, sum)

  vGamma_pn <- c(t(mGamma[, -K]))
  names(vGamma_pn) <- names(vGamma_tilde)

  return(vGamma_pn)

}

f_unmapMS <- function(vPn, spec) {

  vModel <- f_getModel(spec)
  vDist <- f_getDist(spec)
  is.mix <- spec$is.mix

  K <- spec$K

  if (is.mix) {
    vP_pn <- tail(vPn, (K - 1))
    vP_tilde <- as.numeric(SimplexUnmapping(vP_pn, K))
  } else {
    vP_pn <- tail(vPn, K * (K - 1))
    vP_tilde <- f_unmapGamma(vP_pn, K)
  }
  vSkew <- vPn[paste("xi", 1:K, sep = "_")]  #nas for non skew cond dist
  vShape <- vPn[paste("nu", 1:K, sep = "_")]  #nas for gauss cond dist

  vPw <- vPn
  vPw[] <- NA

  for (k in 1:K) {
    vParNames_Model <- paste(f_ModelParNames(vModel[k]), k, sep = "_")
    vPn_foo <- vPn[vParNames_Model]

    dSkew_foo <- vSkew[k]
    dShape_foo <- vShape[k]


    vPw[vParNames_Model] <- f_pn2pw_SR(vPn_foo, k, vModel[k], dSkew_foo, dShape_foo, vDist[k])

  }

  lower <- spec$lower
  upper <- spec$upper

  names(lower) <- names(spec$par0)
  names(upper) <- names(spec$par0)

  vSkew <- vSkew[!is.na(vSkew)]
  vShape <- vShape[!is.na(vShape)]

  vSkew_tilde <- f_unmap(vSkew, lower[names(vSkew)], upper[names(vSkew)])[1, ]
  names(vSkew_tilde) <- names(vSkew)

  vPw[names(vSkew_tilde)] <- vSkew_tilde

  vShape_tilde <- f_unmap(vShape, lower[names(vShape)], upper[names(vShape)])[1, ]
  names(vShape_tilde) <- names(vShape)

  vPw[names(vShape_tilde)] <- vShape_tilde

  L <- length(vPw)

  if (is.mix) {
    vPw[(L - (K - 1) + 1):L] <- vP_tilde
  } else {
    vPw[(L - K * (K - 1) + 1):L] <- vP_tilde
  }

  return(vPw)
}

#' @importFrom stats na.omit
f_mapMS <- function(vPw, spec) {
  vModel <- f_getModel(spec)
  vDist <- f_getDist(spec)
  is.mix <- spec$is.mix

  K <- spec$K

  if (is.mix) {
    vP_tilde <- tail(vPw, K - 1)
    vP_pn <- as.numeric(SimplexMapping(vP_tilde, K))
  } else {
    vP_tilde <- tail(vPw, K * (K - 1))
    vP_pn <- f_mapGamma(vP_tilde, K)
  }

  vPn <- vPw
  vPn[] <- NA

  vSkew_tilde <- vPw[paste("xi", 1:K, sep = "_")]  #nas for non skew cond dist
  vShape_tilde <- vPw[paste("nu", 1:K, sep = "_")]  #nas for gauss cond dist

  lower <- spec$lower
  upper <- spec$upper

  names(lower) <- names(spec$par0)
  names(upper) <- names(spec$par0)

  vSkew <- f_map(vSkew_tilde, lower[names(vSkew_tilde)], upper[names(vSkew_tilde)])[1, ]
  names(vSkew) <- names(vSkew_tilde)

  vShape <- f_map(vShape_tilde, lower[names(vShape_tilde)], upper[names(vShape_tilde)])[1, ]
  names(vShape) <- names(vShape_tilde)

  for (k in 1:K) {

    vParNames_Model <- paste(f_ModelParNames(vModel[k]), k, sep = "_")
    vPw_foo <- vPw[vParNames_Model]

    dSkew_foo <- vSkew[k]
    dShape_foo <- vShape[k]


    vPn[vParNames_Model] <- f_pw2pn_SR(vPw_foo, k, vModel[k], dSkew_foo, dShape_foo, vDist[k])

  }

  vSkew <- na.omit(vSkew)
  vShape <- na.omit(vShape)

  vPn[names(vShape)] <- vShape
  vPn[names(vSkew)] <- vSkew

  L <- length(vPn)

  if (is.mix) {
    vPn[(L - (K - 1) + 1):L] <- vP_pn
  } else {
    vPn[(L - K * (K - 1) + 1):L] <- vP_pn
  }
  return(vPn)
}

###################################### sARCH ####

f_pw2pn_sARCH_SR <- function(vPw, k) {

  dAlpha0_tilde <- vPw[paste("alpha0", k, sep = "_")]
  dAlpha1_tilde <- vPw[paste("alpha1", k, sep = "_")]

  dAlpha0 <- exp(dAlpha0_tilde)
  dAlpha1 <- f_map(dAlpha1_tilde, 1e-10, 0.9999)[1, 1]

  names(dAlpha1) <- paste("alpha1", k, sep = "_")

  vPn <- c(dAlpha0, dAlpha1)

  return(vPn)
}

f_pn2pw_sARCH_SR <- function(vPn, k) {

  dAlpha0 <- vPn[paste("alpha0", k, sep = "_")]
  dAlpha1 <- vPn[paste("alpha1", k, sep = "_")]

  dAlpha0_tilde <- log(dAlpha0)
  dAlpha1_tilde <- f_unmap(dAlpha1, 1e-10, 0.9999)[1, 1]

  names(dAlpha1_tilde) <- paste("alpha1", k, sep = "_")

  vPw <- c(dAlpha0_tilde, dAlpha1_tilde)

  return(vPw)
}


###################################### sGARCH ####

f_pw2pn_sGARCH_SR <- function(vPw, k) {

  dAlpha0_tilde <- vPw[paste("alpha0", k, sep = "_")]
  dAlpha1_tilde <- vPw[paste("alpha1", k, sep = "_")]
  dBeta_tilde <- vPw[paste("beta", k, sep = "_")]

  dAlpha0 <- exp(dAlpha0_tilde)
  dAlpha1 <- f_map(dAlpha1_tilde, 1e-10, 0.9999)[1, 1]
  dBeta <- f_map(dBeta_tilde, 1e-10, 0.9999 - dAlpha1)[1, 1]

  names(dAlpha1) <- paste("alpha1", k, sep = "_")
  names(dBeta) <- paste("beta", k, sep = "_")

  vPn <- c(dAlpha0, dAlpha1, dBeta)

  return(vPn)
}

f_pn2pw_sGARCH_SR <- function(vPn, k) {

  dAlpha0 <- vPn[paste("alpha0", k, sep = "_")]
  dAlpha1 <- vPn[paste("alpha1", k, sep = "_")]
  dBeta <- vPn[paste("beta", k, sep = "_")]

  dAlpha0_tilde <- log(dAlpha0)
  dAlpha1_tilde <- f_unmap(dAlpha1, 1e-10, 0.9999)[1, 1]
  dBeta_tilde <- f_unmap(dBeta, 1e-10, 0.9999 - dAlpha1)[1, 1]

  names(dAlpha1_tilde) <- paste("alpha1", k, sep = "_")
  names(dBeta_tilde) <- paste("beta", k, sep = "_")

  vPw <- c(dAlpha0_tilde, dAlpha1_tilde, dBeta_tilde)

  return(vPw)
}

###################################### eGARCH ####

f_pw2pn_eGARCH_SR <- function(vPw, k) {

  dAlpha0_tilde <- vPw[paste("alpha0", k, sep = "_")]
  dAlpha1_tilde <- vPw[paste("alpha1", k, sep = "_")]
  dAlpha2_tilde <- vPw[paste("alpha2", k, sep = "_")]
  dBeta_tilde <- vPw[paste("beta", k, sep = "_")]

  dBeta <- f_map(dBeta_tilde, -0.9999, 0.9999)[1, 1]

  names(dBeta) <- paste("beta", k, sep = "_")
  dAlpha0 <- dAlpha0_tilde
  dAlpha1 <- dAlpha1_tilde
  dAlpha2 <- dAlpha2_tilde

  vPn <- c(dAlpha0, dAlpha1, dAlpha2, dBeta)

  return(vPn)
}

f_pn2pw_eGARCH_SR <- function(vPn, k) {

  dAlpha0 <- vPn[paste("alpha0", k, sep = "_")]
  dAlpha1 <- vPn[paste("alpha1", k, sep = "_")]
  dAlpha2 <- vPn[paste("alpha2", k, sep = "_")]
  dBeta <- vPn[paste("beta", k, sep = "_")]

  dBeta_tilde <- f_unmap(dBeta, -0.9999, 0.9999)[1, 1]

  names(dBeta_tilde) <- paste("beta", k, sep = "_")
  dAlpha0_tilde <- dAlpha0
  dAlpha1_tilde <- dAlpha1
  dAlpha2_tilde <- dAlpha2

  vPw <- c(dAlpha0_tilde, dAlpha1_tilde, dAlpha2_tilde, dBeta_tilde)

  return(vPw)
}

###################################### gjrGARCH ####

f_pw2pn_gjrGARCH_SR <- function(vPw, k, skew = NULL, shape = NULL, sDist) {

  if (is.null(skew)) {
    skew = 1
  }
  if (is.null(shape)) {
    shape = 100
  }

  dAlpha0_tilde = vPw[paste("alpha0", k, sep = "_")]
  dAlpha1_tilde = vPw[paste("alpha1", k, sep = "_")]
  dAlpha2_tilde = vPw[paste("alpha2", k, sep = "_")]
  dBeta_tilde = vPw[paste("beta", k, sep = "_")]

  dCdf = pdist(0, dist = sDist, shape = shape, skew = skew)

  dAlpha0 = exp(dAlpha0_tilde)
  dAlpha1 = f_map(dAlpha1_tilde, 1e-10, 0.9999)[1, 1]
  dAlpha2 = f_map(dAlpha2_tilde, 1e-10, 2.0 * (0.9999 - dAlpha1))[1, 1]

  dBeta = f_map(dBeta_tilde, 1e-10, 0.9999 - dAlpha1 - dAlpha2 * dCdf)[1, 1]

  names(dAlpha1) = paste("alpha1", k, sep = "_")
  names(dAlpha2) = paste("alpha2", k, sep = "_")
  names(dBeta) = paste("beta", k, sep = "_")

  vPn = c(dAlpha0, dAlpha1, dAlpha2, dBeta)

  return(vPn)
}

f_pn2pw_gjrGARCH_SR <- function(vPn, k, skew = NULL, shape = NULL, sDist) {

  if (is.null(skew)) {
    skew = 1
  }
  if (is.null(shape)) {
    shape = 100
  }

  dAlpha0 = vPn[paste("alpha0", k, sep = "_")]
  dAlpha1 = vPn[paste("alpha1", k, sep = "_")]
  dAlpha2 = vPn[paste("alpha2", k, sep = "_")]
  dBeta = vPn[paste("beta", k, sep = "_")]

  dCdf = pdist(0, dist = sDist, shape = shape, skew = skew)

  dAlpha0_tilde = log(dAlpha0)
  dAlpha1_tilde = f_unmap(dAlpha1, 1e-10, 0.9999)[1, 1]
  dAlpha2_tilde = f_unmap(dAlpha2, 1e-10, 2.0 * (0.9999 - dAlpha1))[1, 1]

  dBeta_tilde = f_unmap(dBeta, 1e-10, 0.9999 - dAlpha1 - dAlpha2 * dCdf)[1, 1]

  names(dAlpha1_tilde) = paste("alpha1", k, sep = "_")
  names(dAlpha2_tilde) = paste("alpha2", k, sep = "_")
  names(dBeta_tilde) = paste("beta", k, sep = "_")

  vPw = c(dAlpha0_tilde, dAlpha1_tilde, dAlpha2_tilde, dBeta_tilde)

  return(vPw)
}

###################################### tGARCH ####

f_pw2pn_tGARCH_SR <- function(vPw, k, skew = NULL, shape = NULL, sDist) {

  if (is.null(skew)) {
    skew = 1
  }
  if (is.null(shape)) {
    shape = 100
  }

  dAlpha0_tilde = vPw[paste("alpha0", k, sep = "_")]
  dAlpha1_tilde = vPw[paste("alpha1", k, sep = "_")]
  dAlpha2_tilde = vPw[paste("alpha2", k, sep = "_")]
  dBeta_tilde = vPw[paste("beta", k, sep = "_")]

  Ez = f_EzDist(shape, skew, sDist)

  EzIneg = Ez["EzIneg"]
  Ez2Ineg = Ez["Ez2Ineg"]

  dAlpha0 = exp(dAlpha0_tilde)
  dAlpha1 = f_map(dAlpha1_tilde, 1e-10, 0.9999)[1, 1]
  dAlpha2 = f_map(dAlpha2_tilde, 1e-10, 0.9999 - dAlpha1)[1, 1]

  dDelta = (dAlpha1 + dAlpha2)^2 * EzIneg^2 + Ez2Ineg * (dAlpha1^2 - dAlpha2^2) + 1 - dAlpha1

  vRange = (dAlpha1 + dAlpha2) * EzIneg + c(-1, 1) * sqrt(dDelta)

  if (vRange[1] < 0) {
    vRange[1] = 0
  }

  dBeta = f_map(dBeta_tilde, vRange[1], vRange[2])[1, 1]

  names(dAlpha1) = paste("alpha1", k, sep = "_")
  names(dAlpha2) = paste("alpha2", k, sep = "_")
  names(dBeta) = paste("beta", k, sep = "_")

  vPn = c(dAlpha0, dAlpha1, dAlpha2, dBeta)

  return(vPn)
}

f_pn2pw_tGARCH_SR <- function(vPn, k, skew = NULL, shape = NULL, sDist) {

  if (is.null(skew)) {
    skew <- 1
  }
  if (is.null(shape)) {
    shape <- 100
  }

  dAlpha0 = vPn[paste("alpha0", k, sep = "_")]
  dAlpha1 = vPn[paste("alpha1", k, sep = "_")]
  dAlpha2 = vPn[paste("alpha2", k, sep = "_")]
  dBeta = vPn[paste("beta", k, sep = "_")]

  Ez = f_EzDist(shape, skew, sDist)

  EzIneg = Ez["EzIneg"]
  Ez2Ineg = Ez["Ez2Ineg"]

  dDelta = (dAlpha1 + dAlpha2)^2 * EzIneg^2 + Ez2Ineg * (dAlpha1^2 - dAlpha2^2) + 1 - dAlpha1

  vRange = (dAlpha1 + dAlpha2) * EzIneg + c(-1, 1) * sqrt(dDelta)

  if (vRange[1] < 0) {
    vRange[1] = 0
  }

  dAlpha0_tilde = log(dAlpha0)
  dAlpha1_tilde = f_unmap(dAlpha1, 1e-10, 0.9999)[1, 1]
  dAlpha2_tilde = f_unmap(dAlpha2, 1e-10, 0.9999 - dAlpha1)[1, 1]
  dBeta_tilde = f_unmap(dBeta, vRange[1], vRange[2])[1, 1]

  names(dAlpha1_tilde) = paste("alpha1", k, sep = "_")
  names(dAlpha2_tilde) = paste("alpha2", k, sep = "_")
  names(dBeta_tilde) = paste("beta", k, sep = "_")

  vPw = c(dAlpha0_tilde, dAlpha1_tilde, dAlpha2_tilde, dBeta_tilde)

  return(vPw)
}
