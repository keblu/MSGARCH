f_CondVol <- function(object, par, data, do.its = FALSE, n.ahead = 1L, ctr = list(), ...) {
  object    <- f_check_spec(object)
  data      <- f_check_y(data)
  par.check <- f_check_par(object, par)
  if (nrow(par.check) == 1) {
    ctr     <- f_process_ctr(ctr)
    n.sim <- ctr$n.sim
  } else {
    if(is.null(ctr$n.sim)){
      n.sim = 1
    } else {
      n.sim = ctr$n.sim
    }
  }
  ctr       <- f_process_ctr(ctr)
  variance  <- object$rcpp.func$calc_ht(par.check, data)
  P         <- State(object = object,par = par, data = data)
  PredProb  <- P$PredProb
  vol <- matrix(NA, nrow = dim(PredProb)[1], ncol = nrow(par.check))
  if (object$K == 1) {
    for (i in 1:nrow(par.check)) {
      vol[,i] <- variance[,i]
    }
  } else {
    for (i in 1:nrow(par.check)) {
      vol[,i] <- rowSums(PredProb[,i,] * variance[,i,])
    }
  }
  vol <- sqrt(vol)
  draw <- NULL
  if (!isTRUE(do.its)) {
    tmp    <- mean(vol[dim(PredProb)[1]])
    vol    <- vector(mode = "numeric", length = n.ahead)
    vol[1] <- tmp
    if (n.ahead > 1) {
      draw <- Sim(object = object, data = data, n.ahead = n.ahead,
                  n.sim = n.sim, par = par)$draw
      vol[2:n.ahead] = apply(draw[2:n.ahead,, drop = FALSE], 1, sd)
    }
    names(vol) <- paste0("h=", 1:n.ahead)
  } else {
    draw <- NULL
    if (nrow(par.check) > 1) {
      vol <- rowMeans(vol[1:length(data),])
      vol <- as.vector(vol)
    } else {
      vol <- as.vector(vol)
      vol <- vol[1:length(data)]
    }
    names(vol) <- paste0("t=", 1:(length(data)))
  }
  out = list()
  class(vol) <- "MSGARCH_CONDVOL"
  out$vol = vol
  out$draw = draw
  return(out)
}


