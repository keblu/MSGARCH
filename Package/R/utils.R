#Error function
f.error <- function(message) {
  cat(paste0("try/catch ", message, "\n"))
  return(FALSE)
}

#Default parameters
f.process.ctr <- function(ctr = list()) {
  
  con <- list(theta0 = NULL, do.init = FALSE, N.mcmc = 1000, N.burn = 500, N.thin = 1,
              NP = 500, itermax = 500, do.enhance.theta0 = FALSE)
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
# 
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
  params_loc  <- c(0, cumsum(Nbparams))
  tmp <- theta
  pos <- 1:Nmodel
  for (i in 1:length(unique.spec)) {
    postmp       <- pos
    idx          <- name == unique.spec[i]
    posidx       <- pos[idx]
    Nmodelidx    <- length(posidx)
    idx_loc      <- params_loc[c(idx)]
    idx_params   <- spec$n.params[c(idx)][1]
    unc.vol      <- MSGARCH::unc.vol(object = spec, thetaUncVol)
    unc.vol.idx  <- unc.vol[idx]
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
    for (i in 0:(Nmodel-1)) {
      
      p[1:(Nmodel-1),i+1] <- theta[(params_loc[Nmodel + 1] + Nmodel*i+1-i):(params_loc[Nmodel + 1] + Nmodel*i+Nmodel-1-i)]
      
    }
    p[Nmodel, ] <- 1 - colSums(matrix(p[1:(Nmodel-1), ], ncol = Nmodel))
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
      p[i] <- tmp[(params_loc[Nmodel + 1] + i)]
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
#' @importFrom stats optim 
f.enhance.theta <- function(spec, theta, y) {
  
  if (is.vector(theta)) {
    # DA needed for the Rcpp framework
    theta = matrix(theta, nrow = 1)
  }
  
  K   <- spec$K
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
  if (isTRUE(spec$is.shape.ind)) {
    pos <- c(1, cumsum(spec$n.params.vol) + 1)
    pos[length(pos)] = pos[length(pos)] + 1
  } else {
    pos <- c(1, cumsum(spec$n.params) + 1)
  }
  is.ok = NULL
  for (i in 1:K) {
    f.fun <- function(x) {
      theta.try <- theta
      theta.try[, pos[i]] <- x
      unc.vol <- MSGARCH::unc.vol(object = spec, theta.try)
      #out = unc.vol[i] - vol.goal[i]
      out = (unc.vol[i] - vol.goal[i])^2
      return(out)
    }
    # DA replace the unitroot which fails sometimes in cluster by optimization 
    #theta[, pos[i]] <- uniroot(f.fun, lower = spec$lower[pos[i]], upper = spec$upper[pos[i]])$root
    theta.hat = spec$theta0[pos[i]]
    
    str = "f.enhance.theta -> fail in optimization"
    is.ok[i] = tryCatch({
      tmp = stats::optim(par = theta.hat, fn = f.fun, method = "Brent",
                  lower = spec$lower[pos[i]]+0.001, upper = spec$upper[pos[i]]-0.001)
      if (tmp$convergence == 0) {
        theta.hat = tmp$par
      }
      TRUE
    }, warning = function(warn) {
      f.error(str)
    }, error = function(err) {
      f.error(str)
    })
  }
   if (any(!isTRUE(is.ok))){
     theta = spec$theta0
   }
  
  if (spec$K > 1) {
    if (!spec$is.mix) {
      pos <- pos[length(pos)]
      theta[pos:length(theta)] <- (1 - 0.8) / (K - 1)
      for (i in 1:(K - 1)) {
        theta[pos] <- 0.8
        pos <- pos + K 
      }
    }
  }
  
  return(theta)
}