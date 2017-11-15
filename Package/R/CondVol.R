f_CondVol <- function(object, par, data, do.its = FALSE, nahead = 1L, do.cumulative = FALSE, ctr = list(), ...) {
  object    <- f_check_spec(object)
  data_      <- f_check_y(data)
  par.check <- f_check_par(object, par)
  if (nrow(par.check) == 1) {
    ctr     <- f_process_ctr(ctr)
    nsim <- ctr$nsim
  } else {
    if(is.null(ctr$nsim)){
      nsim = 1
    } else {
      nsim = ctr$nsim
    }
  }
  ctr       <- f_process_ctr(ctr)
  variance  <- object$rcpp.func$calc_ht(par.check, data_)
  P         <- State(object = object,par = par, data = data_)
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
    vol    <- vector(mode = "numeric", length = nahead)
    vol[1] <- tmp
    if (nahead > 1) {
      draw <- Sim(object = object, data = data, nahead = nahead,
                  nsim = nsim, par = par)$draw
      if(isTRUE(do.cumulative)){
        draw = apply(draw, 2, cumsum)
      }
      vol[2:nahead] = apply(draw[2:nahead,, drop = FALSE], 1, sd)
    }
    names(vol) <- paste0("h=", 1:nahead)
    if(zoo::is.zoo(data)){
      vol = zoo::zooreg(vol, order.by = zoo::index(data)[length(data)]+(1:nahead))
    }
    if(is.ts(data)){
      vol = zoo::zooreg(vol, order.by = zoo::index(data)[length(data)]+(1:nahead))
      vol = as.ts(vol)
    }
  } else {
    draw <- NULL
    if (nrow(par.check) > 1) {
      vol <- rowMeans(vol[1:length(data_),])
      vol <- as.vector(vol)
    } else {
      vol <- as.vector(vol)
      vol <- vol[1:length(data_)]
    }
    names(vol) <- paste0("t=", 1:(length(data_)))
    if(zoo::is.zoo(data)){
      vol = zoo::zooreg(vol, order.by = zoo::index(data))
    }
    if(is.ts(data)){
      vol = zoo::zooreg(vol, order.by = zoo::index(data))
      vol = as.ts(vol)
    }
  }
  out = list()
  class(vol) <- c("MSGARCH_CONDVOL",class(vol))
  out$vol = vol
  out$draw = draw
  return(out)
}


