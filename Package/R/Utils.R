# Get Gamma par names
f_GammaParNames <- function(K) {
  vNames = character(K * (K - 1))
  iC = 1
  for (i in 1:(K - 1)) {
    for (j in 1:K) {
      vNames[iC] = paste("P", j, i, sep = "_")
      iC = iC + 1
    }
  }
  return(vNames)
}

# Get models from spec
f_getModel <- function(spec) {
  name = spec$name
  vModels = sapply(name, function(x) unlist(strsplit(x, split = "_"))[1])
  names(vModels) = NULL
  return(vModels)
}

# Get conditional distributions from spec
f_getDist <- function(spec) {
  name = spec$name
  vDist <- sapply(name, function(x) unlist(strsplit(x, split = "_"))[2L])
  vSkew <- sapply(vDist, FUN = function(x) any(c("snorm", "sstd", "sged") == x))
  vDist[vSkew] <- substring(vDist[vSkew], 2)
  names(vDist) = NULL
  return(vDist)
}

# Check variance spec
f_check_variance_spec <- function(variance.spec) {

  if (!is(variance.spec, "list")) {
    stop("variance.spec has to be a list.")
  }
  if (is.null(variance.spec$model)) {
    variance.spec$model = c("sGARCH", "sGARCH")
  } else {
    if (!is(variance.spec$model, "character")) {
      stop("variance.spec$model has to be a character vector
           with the variance specifications.")
    }
  }
  if (is.null(variance.spec$do.mix)) {
    variance.spec$do.mix = FALSE
  } else {
    if (!is(variance.spec$do.mix, "logical")) {
      stop("variance.spec$do.mix has to be a logical.")
    }
  }
  return(variance.spec)
}

# Check variance spec
f_check_distribution_spec <- function(distribution.spec, K) {
  if (!is(distribution.spec, "list")) {
    stop("distribution.spec has to be a list.")
  }
  if (is.null(distribution.spec$distribution)) {
    distribution.spec$distribution = rep("norm", K)
  } else {
    if (!is(distribution.spec$distribution, "character")) {
      stop("distribution.spec$distribution has to be a character
           vector with the distribution specifications.")
    }
  }
  if (is.null(distribution.spec$do.skew)) {
    distribution.spec$do.skew = rep(FALSE, K)
  } else {
    if (!is(distribution.spec$do.skew, "logical")) {
      stop("distribution.spec$do.skew has to be a vector of logicals.")
    }
  }
  if (is.null(distribution.spec$do.shape.ind)) {
    distribution.spec$do.shape.ind = FALSE
  } else {
    if (!is(distribution.spec$do.shape.ind, "logical")) {
      stop("distribution.spec$do.shape.ind has to be a logical.")
    }
  }
  return(distribution.spec)
}

# Check markov spec
f_check_markov_spec <- function(markov.spec) {
  if (!is(markov.spec, "list")) {
    stop("markov.spec has to be a list.")
  }
  if (is.null(markov.spec$do.mix)) {
    markov.spec$do.mix = FALSE
  } else {
    if (!is(markov.spec$do.mix, "logical")) {
      stop("markov.spec$do.mix has to be a logical.")
    }
  }
  return(markov.spec)
}

# Error function
f_error <- function(message) {
  cat(paste0("try/catch ", message, "\n"))
  return(FALSE)
}

# Inverse generalize logit map to constraint
f_map <- function(x, lower, upper) {
  if (is.vector(x)) {
    x <- matrix(data = x, nrow = 1L, ncol = length(x), dimnames = list(NULL, names(x)))
  }
  x_map <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))

  if (is.matrix(x)) {
    colnames(x_map) <- colnames(x)
  }
  for (i in 1:nrow(x)) {
    x_map[i, ] <- lower + (upper - lower)/(1 + exp(-x[i, ]))
  }
  return(x_map)
}

# Generalize logit map to real line
f_unmap <- function(x, lower, upper) {
  if (is.vector(x)) {
    x <- matrix(data = x, nrow = 1L, ncol = length(x), dimnames = list(NULL, names(x)))
  }
  x_unmap <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))

  if (is.matrix(x)) {
    colnames(x_unmap) <- colnames(x)
  }
  for (i in 1:nrow(x)) {
    x_unmap[i, ] <- log((x[i, ] - lower)/(upper - x[i, ]))
  }
  return(x_unmap)
}

# Inverse generalize logit gradient
f_map_deriv <- function(x, lower, upper) {
  if (is.vector(x)) {
    x <- matrix(data = x, nrow = 1L, ncol = length(x))
  }
  x_map_deriv <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:nrow(x)) {
    x_map_deriv[i, ] <- -x[i, ] + log(upper - lower) - 2 * log(1 + exp(-x[i, ]))
  }
  return(exp(x_map_deriv))
}

# Default parameters
f_process_ctr <- function(ctr = list(), type = 1) {
  if (type == 1) {
    con <- list(par0 = NULL, n.mcmc = 10000L,
                OptimFUN = f_OptimFUNDefault,
                SamplerFUN = f_SamplerFUNDefault,
                n.burn = 5000L, n.thin = 10L,  do.se = TRUE, do.plm = FALSE,
                n.sim = 10000L, n.mesh = 1000L)
  } else if (type == 2) {
    con <- list(n.sim = 250L, n.burn = 5000L, n.ahead = 1000L)
  }
  con[names(ctr)] <- ctr
  return(con)
}

# Function that checks if the passed y is one of the good format
f_check_y <- function(y) {
  if (zoo::is.zoo(y)) {
    y = zoo::coredata(y)
  }
  if (is.null(y)) {
    stop("y is NULL")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (all(is.nan(y))) {
    stop("nan dectected in y")
  }
  return(y)
}

# Function that checks if the passed par are of the good format
f_check_par <- function(spec, par) {
  if (is.null(par)) {
    stop("par is NULL")
  }
  if (!is.numeric(par)) {
    stop("par must be a numeric")
  }
  if (all(is.nan(par))) {
    stop("nan dectected in par")
  }
  len.par <- length(spec$par0)
  if (is.vector(par)) {
    par <- matrix(par, nrow = 1L)
  }
  if (is.data.frame(par)) {
    par <- data.matrix(par)
  }
  if (dim(par)[2L] != len.par) {
    stop(paste0("Each parameter estimate in par must be of length ", len.par))
  }
  if (isTRUE(spec$is.mix)) {
    par <- spec$func$f.do.mix(par)
  }
  colnames(par) = spec$label[1:ncol(par)]
  return(par)
}

# Function that sorts the par according to the unconditional variance (Used for Bayesian estimation)
f_sort_par <- function(spec, par) {
  parUncVol <- par
  Nbparams <- spec$n.params
  Nmodel <- length(Nbparams)
  if (Nmodel == 1L) {
    return(par)
  }
  name <- spec$name
  unique.spec <- unique(name, FALSE)
  params_loc  <- c(0, cumsum(Nbparams))
  tmp <- par
  pos <- 1:Nmodel
  unc.vol.all <- spec$rcpp.func$unc_vol_Rcpp(parUncVol)
  for (f in 1:nrow(par)) {
    for (i in 1:length(unique.spec)) {
      postmp <- pos
      idx    <- name == unique.spec[i]
      posidx       <- pos[idx]
      Nmodelidx    <- length(posidx)
      idx_loc      <- params_loc[c(idx)]
      idx_params   <- spec$n.params[c(idx)][1]
      unc.vol      <- unc.vol.all[f,]
      unc.vol.idx  <- unc.vol[idx]
      unc.vol.sort <- sort(unc.vol.idx, index.return = TRUE)
      if (all(unc.vol.idx == unc.vol.sort$x)) {
        next()
      }
      new.pos.index = 1
      for (j in 1:Nmodelidx) {
        new.pos <- which(unc.vol[j] == unc.vol.sort$x)
        if (length(new.pos) > 1L) {
          new.pos <- new.pos[new.pos.index]
          new.pos.index <- new.pos.index + 1L
        }

        postmp[posidx[j]] <- new.pos
      }
      for (j in 1:Nmodelidx) {
        ind <- unc.vol.sort$ix[j]
        tmp[f, (idx_loc[j] + 1L):(idx_loc[j] + idx_params)] <- par[f,(idx_loc[ind] + 1L):(idx_loc[ind] + idx_params)]
      }
    }
    if (all(pos == postmp)) {
      next()
    }
    if (!isTRUE(spec$is.mix)) {
      p <- matrix(nrow = Nmodel, ncol = Nmodel)
      for (i in 0:(Nmodel - 1L)) {
        p[1:(Nmodel - 1), i + 1L] <- par[f,(params_loc[Nmodel + 1L] + Nmodel * i + 1 - i):(params_loc[Nmodel + 1L] + Nmodel * i + Nmodel - 1L - i)]

      }
      p[Nmodel, ] <- 1 - colSums(matrix(p[1:(Nmodel - 1L), ], ncol = Nmodel))
      tmpp <- p
      for (i in 1:(Nmodel)) {
        for (j in 1:(Nmodel)) {
          tmpp[i, j] <- p[postmp[i], postmp[j]]
        }
      }
      new.p <- as.vector(tmpp[1:(Nmodel - 1L), ])
      tmp[f, (params_loc[Nmodel + 1] + 1):ncol(tmp)] <- new.p
    } else {
      p <- rep(0, Nmodel)
      for (i in 1:(Nmodel - 1L)) {
        p[i] <- tmp[f, (params_loc[Nmodel + 1L] + i)]
      }
      p[Nmodel] <- 1 - sum(p)
      tmpp <- p
      for (j in 1:(Nmodel)) {
        tmpp[j] <- p[postmp[j]]
      }
      new.p <- tmpp[1:(Nmodel - 1L)]
      tmp[f, (params_loc[Nmodel + 1L] + 1L):ncol(tmp)] <- new.p
    }
  }
  return(tmp)
}

# Get gamma !!! DA we should put this function elsewhere
f_getGamma <- function(par, K) {
  vP <- c(tail(par[1L, ], K * (K - 1L)), rep(0, K))
  mGamma <- matrix(vP, K, K)
  mGamma[, K] <- 1 - apply(mGamma[, -K, drop = FALSE], 1L, sum)
  return(mGamma)
}


f_par_mixture <- function(K, nb_params, par) {
  if (is.vector(par)) {
    par <- matrix(par, nrow = 1L)
  }
  n_total_params <- ncol(par)
  n_par <- nrow(par)
  new_par <- matrix(nrow = n_par, ncol = nb_params + K * (K - 1L))
  for (i in 1:nrow(par)) {
    p_vector <- par[i, (nb_params + 1L):n_total_params]
    new_p <- rep(p_vector, K)
    new_par[i, ] <- c(par[i, 1:nb_params], new_p)
  }
  return(new_par)
}

f_par_mixture_reverse <- function(K, nb_params, par) {
  if (is.vector(par)) {
    par <- matrix(par, nrow = 1L)
  }
  n_total_params <- ncol(par)
  n_par <- nrow(par)
  new_par <- matrix(nrow = n_par, ncol = nb_params + K - 1L)
  for (i in 1:n_par) {
    p_vector <- par[i, (nb_params + 1L):n_total_params]
    new_par[i, 1:nb_params] <- par[i, 1:nb_params]
    idx <- 1L
    for (j in 1:(K - 1L)) {
      new_par[i, (nb_params + j)] <- p_vector[idx]
      idx <- idx + K
    }
  }
  return(new_par)
}

f_match <- function(x, target) {
  return(toupper(substr(x, 1, 3)) == target)
}

f_check_spec <- function(spec) {
  is.OK = tryCatch({
    spec$rcpp.func$get_sd()
    TRUE
  }, error = function(e) {
    FALSE
  })
  if (!isTRUE(is.OK)) {
    spec.new = f_spec(models = spec$name, do.mix = spec$is.mix)
    prior.mean = spec$rcpp.func$get_mean()
    prior.sd = spec$rcpp.func$get_sd()
    names(prior.mean) = names(prior.sd) = spec$label[1:length(prior.mean)]
    prior.mean[names(spec$prior.mean)] = spec$prior.mean
    prior.sd[names(spec$prior.sd)] = spec$prior.sd
    spec$rcpp.func = spec.new$rcpp.func
    spec$rcpp.func$set_mean(spec$prior.mean)
    spec$rcpp.func$set_sd(spec$prior.sd)
  }
  return(spec)
}

f_rename_par <- function(vPar, spec) {
  vNames <- spec$label

  if (isTRUE(spec$regime.const.pars.bool)) {
    for (i in 1:length(spec$regime.const.pars)) {
      vNames <- vNames[vNames != paste0(spec$regime.const.pars[i], "_" ,2:spec$K)]
    }
  }

  if (isTRUE(spec$fixed.pars.bool)) {
    names(vPar) <- vNames[!vNames %in% names(spec$fixed.pars)]
  } else {
    names(vPar) <- vNames
  }
  return(vPar)
}

#' @importFrom stats ecdf
f_cdf_empirical = function(y, x) {
  return((stats::ecdf(y))(x))
}

#' @importFrom stats density
f_pdf_kernel = function(y, x) {
  d <- density(x = y)
  h = d$bw
  x.length = length(x)
  out = vector(mode = "numeric", length = x.length)
  for (i in 1:x.length) {
    out[i] = sum(dnorm(x[i], mean = y, sd = h))/length(y)
  }
  return(out)
}

f_rbindrep = function(mat, n) {
  matFull = NULL
  for (i in 1:n) {
    matFull = rbind(matFull, mat)
  }
  return(matFull)
}

f_check_parameterPriorMean <- function(prior.mean, vParNames) {
  if (any(!names(prior.mean) %in% vParNames)) {
    vWrongPars <- names(prior.mean)[!names(prior.mean) %in% vParNames]
    stop(cat(paste("Wrong name in prior.mean:", vWrongPars)))
  }
  return(prior.mean)
}

f_check_parameterPriorSd <- function(prior.sd, vParNames) {
  if (any(!names(prior.sd) %in% vParNames)) {
    vWrongPars <- names(prior.sd)[!names(prior.sd) %in% vParNames]
    stop(cat(paste("Wrong name in prior.sd:", vWrongPars)))
  }
  return(prior.sd)
}
