#Error function
f.error <- function(message) {
  cat(paste0("try/catch ", message, "\n"))
  return(FALSE)
}

#Default parameters
f.process.ctr <- function(ctr = list()) {
  con <- list(theta0 = NULL, do.init = TRUE, N.mcmc = 10000, N.burn = 5000, N.thin = 10,
    NP = 500, itermax = 500, do.enhance.theta0 = TRUE)
  con[names(ctr)] <- ctr
  return(con)
}

#Function that checks if the passed y is one of the good format
f.check.y <- function(y) {
  if (is.null(y)) {
    stop("y is NULL")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (all(is.nan(y))) {
    stop("nan dectected in y")
  }
  if (!is.null(dim(y))) {
    if (any(dim(y) == 1)) {
      y <- as.vector(as.matrix(y))
    } else {
      stop("y is not a vector")
    }
  }
  y <- as.vector(y)
  y <- as.matrix(y)
  return(y)
}

#Function that checks if the passed theta are of the good format
f.check.theta <- function(spec, theta) {
  if (is.null(theta)) {
    stop("theta is NULL")
  }
  if (!is.numeric(theta)) {
    stop("theta must be a numeric")
  }
  if (all(is.nan(theta))) {
    stop("nan dectected in theta")
  }
  len.theta <- length(spec$theta0)
  if (is.vector(theta)) {
    theta <- matrix(theta, nrow = 1)
  }
  if (is.data.frame(theta)) {
    theta <- data.matrix(theta)
  }
  if (dim(theta)[2] != len.theta) {
    stop(paste0("Each parameter estimate in theta must be of length ", len.theta))
  }
  if (isTRUE(spec$is.shape.ind)) {
    theta <- spec$func$f.do.shape.ind(theta)
  }
  if (isTRUE(spec$is.mix)) {
    theta <- spec$func$f.do.mix(theta)
  }
  return(theta)
}

#Function that sorts the theta according to the unconditional variance (Used for Bayesian estimation)
f.sort.theta <- function(spec, theta) {
  thetaUncVol <- theta
  if (isTRUE(spec$is.shape.ind)) {
    theta <- spec$func$f.do.shape.ind(theta = theta)
  }
  Nbparams <- spec$n.params
  Nmodel <- length(Nbparams)
  if (Nmodel == 1) {
    return(theta)
  }
  name <- spec$name
  unique.spec <- unique(name, FALSE)
  params_loc <- c(0, cumsum(Nbparams))
  tmp <- theta
  pos <- 1:Nmodel
  for (i in 1:length(unique.spec)) {
    postmp <- pos
    idx <- name == unique.spec[i]
    posidx <- pos[idx]
    Nmodelidx <- length(posidx)
    idx_loc <- params_loc[c(idx)]
    idx_params <- spec$n.params[c(idx)][1]
    unc.vol <- MSGARCH::unc.vol(object = spec, thetaUncVol)
    unc.vol.idx <- unc.vol[idx]
    unc.vol.sort <- sort(unc.vol.idx, index.return = TRUE)
    new.pos.index = 1
    for (j in 1:Nmodelidx) {
      new.pos <- which(unc.vol == unc.vol.sort$x[j])
      if(length(new.pos) > 1){
         new.pos <- new.pos[new.pos.index]
         new.pos.index <- new.pos.index + 1
      }
      postmp[posidx[j]] <- new.pos
    }
    for (j in 1:Nmodelidx) {
      ind <- unc.vol.sort$ix[j]
      tmp[(idx_loc[j] + 1):(idx_loc[j] + idx_params)] <- theta[(idx_loc[ind] + 1):(idx_loc[ind] + idx_params)]
    }
  }
  if (!isTRUE(spec$is.mix)) {
    p <- matrix(nrow = Nmodel, ncol = Nmodel)
    for (i in 1:(Nmodel - 1)) {
      p[i, 1:Nmodel] <- tmp[(params_loc[Nmodel + 1] + 1):(params_loc[Nmodel + 1] + Nmodel)]
    }
    p[Nmodel, ] <- 1 - colSums(matrix(p[1:(Nmodel - 1), ], ncol = Nmodel))
    tmpp <- p
    for (i in 1:(Nmodel)) {
      for (j in 1:(Nmodel)) {
        tmpp[i, j] <- p[postmp[i], postmp[j]]
      }
    }
    new.p <- as.vector(tmpp[1:(Nmodel - 1), ])
    tmp[(params_loc[Nmodel + 1] + 1):length(tmp)] <- new.p
  } else {
    p <- rep(0, Nmodel)
    for (i in 1:(Nmodel - 1)) {
      p[i] <- tmp[(params_loc[Nmodel + 1] + 1)]
    }
    p[Nmodel] <- 1 - sum(p)
    tmpp <- p
    for (j in 1:(Nmodel)) {
      tmpp[j] <- p[postmp[j]]
    }
    new.p <- tmpp[1:(Nmodel - 1)]
    tmp[(params_loc[Nmodel + 1] + 1):length(tmp)] <- new.p
  }
  if (isTRUE(spec$is.shape.ind)) {
    tmp <- spec$func$f.do.shape.ind.reverse(tmp)
  }
  return(tmp)
}

#Function that enhance theta0 according to the variance of moving windows of y 
f.enhance.theta <- function(spec, theta, y) {
  K <- spec$K
  l.y <- length(y)
  sep <- seq(from = 1, to = l.y, length.out = 11)
  vol <- NULL
  if (spec$K == 1) {
    vol <- sqrt(var(y))
  }
  for (i in 1:(length(sep) - 1)) {
    vol[i] <- sqrt(var(y[sep[i]:sep[i + 1]]))
  }
  vol.goal <- quantile(vol, prob = seq(0.1, 0.9, length.out = K))
  pos <- c(1, cumsum(spec$n.params) + 1)
  if (isTRUE(spec$is.shape.ind)) {
    theta <- spec$func$f.do.shape.ind(theta)
  }
  if (isTRUE(spec$is.mix)) {
    theta <- spec$func$f.do.mix(theta)
  }
  for (i in 1:K) {
    f.fun <- function(x) {
      theta.try <- theta
      theta.try[, pos[i]] <- x
      if (isTRUE(spec$is.shape.ind)) {
        theta.try <- spec$func$f.do.shape.ind.reverse(theta.try)
      }
      if (isTRUE(spec$is.mix)) {
        theta.try <- spec$func$f.do.mix.reverse(theta.try)
      }
      unc.vol <- MSGARCH::unc.vol(spec, theta = theta.try)
      if (isTRUE(spec$is.shape.ind)) {
        theta.try <- spec$func$f.do.shape.ind(theta.try)
      }
      if (isTRUE(spec$is.mix)) {
        theta.try <- spec$func$f.do.mix(theta.try)
      }
      return(unc.vol[i] - vol.goal[i])
    }
    theta[, pos[i]] <- uniroot(f.fun, lower = spec$lower[pos[i]], upper = spec$upper[pos[i]])$root
  }
  if (spec$K > 1) {
    if (!spec$is.mix) {
      pos <- pos[length(pos)]
      theta[pos:length(theta)] <- (1 - 0.8) / (K - 1)
      for (i in 1:(K - 1)) {
        theta[pos] <- 0.8
        pos <- pos + K + 1
      }
    }
  }
  if (isTRUE(spec$is.shape.ind)) {
    theta <- spec$func$f.do.shape.ind.reverse(theta)
  }
  if (isTRUE(spec$is.mix)) {
    theta <- spec$func$f.do.mix.reverse(theta)
  }
  return(theta)
}