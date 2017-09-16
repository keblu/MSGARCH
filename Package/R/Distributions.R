ddist <- function(x, dist = "norm", shape = 100, skew = 1, log = FALSE) {
  require("MSGARCH")
  valid.distribution <- c("norm", "std", "ged")
  if (!any(dist == valid.distribution)) {
    stop(paste0("error: The distribution ", dist, "
                does not appear to be a valid choice."))
  }
  if (dist == "norm") {
    ind <- 1L
  } else {
    ind <- 0L
  }
  if (skew == 1) {
    skew_lab <- "sym"
  } else {
    skew_lab <- "skew"
  }
  par <- c(shape, skew)
  dist_func <- get(x = paste0(".", dist, "_", skew_lab, "_created"))
  dist_func$load_theta(par, ind)
  out <- vector(mode = "numeric", length = length(x))
  for (i in 1:length(x)) {
    out[i] <- dist_func$f_pdf(x[i])
  }
  if (isTRUE(log)) {
    out <- log(out)
  }
  return(out)
  }

rdist <- function(n, dist = "norm", shape = 100, skew = 1) {
  require("MSGARCH")
  valid.distribution <- c("norm", "std", "ged")
  if (!any(dist == valid.distribution)) {
    stop(paste0("error: The distribution ", dist, "
                does not appear to be a valid choice."))
  }
  if (dist == "norm") {
    ind <- 1L
  } else {
    ind <- 0L
  }
  if (skew == 1) {
    skew_lab <- "sym"
  } else {
    skew_lab <- "skew"
  }
  par <- c(shape, skew)
  dist_func <- get(x = paste0(".", dist, "_", skew_lab, "_created"))
  dist_func$load_theta(par, ind)
  out <- dist_func$f_rnd(n)
  return(out)
  }

qdist <- function(p, dist = "norm", shape = 100, skew = 1) {
  require("MSGARCH")
  valid.distribution <- c("norm", "std", "ged")
  if (!any(dist == valid.distribution)) {
    stop(paste0("error: The distribution ", dist, "
                does not appear to be a valid choice."))
  }
  if (dist == "norm") {
    ind <- 1L
  } else {
    ind <- 0L
  }
  if (skew == 1) {
    skew_lab <- "sym"
  } else {
    skew_lab <- "skew"
  }
  par <- c(shape, skew)
  dist_func <- get(x = paste0(".", dist, "_", skew_lab, "_created"))
  dist_func$load_theta(par, ind)
  out <- vector(mode = "numeric", length = length(p))
  for (i in 1:length(p)) {
    out[i] <- dist_func$f_invsample(p[i])
  }
  return(out)
  }

pdist <- function(q, dist = "norm", shape = 100, skew = 1) {
  require("MSGARCH")
  valid.distribution <- c("norm", "std", "ged")
  if (!any(dist == valid.distribution)) {
    stop(paste0("error: The distribution ", dist, "
                does not appear to be a valid choice."))
  }
  if (dist == "norm") {
    ind <- 1L
  } else {
    ind <- 0L
  }
  if (skew == 1) {
    skew_lab <- "sym"
  } else {
    skew_lab <- "skew"
  }
  par <- c(shape, skew)
  dist_func <- get(x = paste0(".", dist, "_", skew_lab, "_created"))
  dist_func$load_theta(par, ind)
  out <- vector(mode = "numeric", length = length(q))
  for (i in 1:length(q)) {
    out[i] <- dist_func$f_cdf(q[i])
  }
  return(out)
  }

f_EzDist <- function(shape = 100, skew = 1, dist = "norm") {
  require("MSGARCH")
  valid.distribution <- c("norm", "std", "ged")
  if (!any(dist == valid.distribution)) {
    stop(paste0("error: The distribution ", dist, "
                does not appear to be a valid choice."))
  }
  if (dist == "norm") {
    ind <- 1L
  } else {
    ind <- 0L
  }
  if (skew == 1) {
    skew_lab <- "sym"
  } else {
    skew_lab <- "skew"
  }
  par <- c(shape, skew)
  
  dist_func <- get(x = paste0(".", dist, "_", skew_lab, "_created"))
  dist_func$load_theta(par, ind)
  dist_func$set_Eabsz()
  dist_func$set_EzIpos()
  dist_func$set_EzIneg()
  dist_func$set_Ez2Ineg()
  out <- c(dist_func$Eabsz, dist_func$EzIpos, dist_func$EzIneg, dist_func$Ez2Ineg)
  names(out) <- c("Eabsz", "EzIpos", "EzIneg", "Ez2Ineg")
  return(out)
  }