#' @title Simulation of MSGARCH processes.
#' @description  Method for simulating \code{MSGARCH} processes.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{CreateSpec}}
#' or fit object of type \code{MSGARCH_ML_FIT} created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT}
#' created with \code{\link{FitMCMC}}.
#' @param data  Vector (of size T) of observations for filtering.
#' @param new.data  Vector (of size T*) of new observations.. (Default \code{new.data = NULL})
#' @param n.ahead Simulation length. (Default: \code{n.ahead = 1L})
#' @param n.sim Number of simulations. (Default: \code{n.sim = 1L})
#' @param par Vector (of size d) or matrix (of size \code{n.ahead} x d) of parameter
#' estimates where d must have
#' the same length as the default parameters of the specification.
#' @param n.burnin Burnin period discarded (first simulation draws).
#' Not used when \code{data} is provided. (Default: \code{n.burnin = 500L})
#' @param ... Not used. Other arguments to \code{Sim}.
#' @return A list of class \code{MSGARCH_SIM} with the following elements:.
#' \itemize{
#' \item \code{draw}: Matrix (of size \code{n.ahead} x \code{n.sim}) of simulated draws.
#' \item \code{state}: Matrix (of size \code{n.ahead} x \code{n.sim}) of simulated states.
#' \item \code{CondVol}: Array (of size \code{n.ahead} x \code{n.sim} x K) of simulated conditional volatility.  
#' }
#' The \code{MSGARCH_SIM} class contains the \code{plot} method.
#' @details If a matrix of parameters estimates is provided, \code{n.sim} simuations will be done for each row..
#' When \code{data} is provided, the conditional variance and state probability are update up to time \code{T + T* + 1}
#' before beginning the simulations.
#' If \code{data = NULL}, new simulations will start from scratch where the
#' first \code{n.burnin} simulation will be discarded.
#' Passing a \code{MSGARCH_ML_FIT} or \code{MSGARCH_MCMC_FIT} object
#' will automatically filter the corresponding
#' \code{data} in the object and start the simulation ahead of \code{data}.
#' @examples
#' # create model specification
#' # MS(2)-GARCH(1,1)-Normal (default)
#' spec <- CreateSpec()
#'
#' # generate process
#' par.sim <- c(0.1,0.6,0.2,0.2,0.8,0.1,0.99,0.01)
#' set.seed(123)
#' sim <- Sim(object = spec, n.ahead = 1000L, n.sim = 1L, par = par.sim, n.burnin = 500L)
#' plot(sim)
#'
#' # generate process after filtering for fitted model
#' # load data
#' data("SMI", package = "MSGARCH")
#'
#' # fit the model on the data with ML estimation
#' fit <- FitML(spec = spec, data = SMI)
#'
#' set.seed(123)
#' sim <- Sim(fit, n.ahead = 30L, n.sim = 1000L)
#' plot(sim)
#' @export
Sim <- function(object, ...) {
  UseMethod(generic = "Sim", object)
}

#' @rdname Sim
#' @export
Sim.MSGARCH_SPEC <- function(object, data = NULL, n.ahead = 1L,
                             n.sim = 1L, par = NULL, n.burnin = 500L, ...) {
  object <- f_check_spec(object)
  if (is.vector(par)) {
    par <- matrix(par, nrow = 1L)
  }
  # New simulation
  if (is.null(data)) {
    par   <- f_check_par(object, par)
    start <- 1
    end   <- n.sim
    draw  <- matrix(data = NA, nrow = n.ahead + n.burnin, ncol = n.sim * nrow(par) )
    state <- matrix(data = NA, nrow = n.ahead + n.burnin, ncol = n.sim * nrow(par) )
    CondVol <- array(data = NA, dim = c(n.ahead + n.burnin, n.sim * nrow(par), object$K),
                     dimnames =  list(paste0("t=",1:(n.ahead+n.burnin)),
                                      paste0("Sim #",1:(n.sim * nrow(par))),paste0("k=",1:object$K)))
    for (i in 1:nrow(par)) {
      tmp <- object$rcpp.func$sim(n.ahead  + n.burnin, n.sim, par[i, ])
      if (object$K == 1L) {
        draw[,start:end]  <- t(tmp$draws)
        state[,start:end] <- matrix(0, nrow = n.ahead + n.burnin, ncol = n.sim)
        CondVol[,start:end,] <- t(tmp$CondVol)
      } else {
        draw[,start:end]  <- t(tmp$draws)
        state[,start:end] <- t(tmp$state)
        CondVol[,start:end,] <- aperm(tmp$CondVol,perm = c(2,1,3))
      }
      start <- start + n.sim
      end   <- end + n.sim
    }
    draw  <- draw[-(1:n.burnin),,drop = FALSE]
    state <- state[-(1:n.burnin),,drop = FALSE]
    CondVol <- CondVol[-(1:n.burnin),,,drop = FALSE]
    rownames(draw) = rownames(state) = paste0("t=",1:n.ahead)
    colnames(draw) = colnames(state) =  paste0("Sim #",1:(n.sim * nrow(par)))
    dimnames(CondVol)[[1]] = paste0("t=",1:n.ahead)
  } else {
    # Simulation ahead of data
    data  <- f_check_y(data)
    P_0   <- matrix(State(object, par = par, data = data)$PredProb[(length(data) + 1L), ,], ncol = object$K)
    par   <- f_check_par(object, par)
    start <- 1
    end   <- n.sim
    draw  <- matrix(data = NA, nrow = n.ahead, ncol =  n.sim * nrow(par))
    state <- matrix(data = NA, nrow = n.ahead, ncol =  n.sim * nrow(par))
    CondVol <- array(data = NA, dim = c(n.ahead, n.sim * nrow(par), object$K),
                     dimnames =  list(paste0("h=",1:(n.ahead)),
                                      paste0("Sim #",1:(n.sim * nrow(par))),paste0("k=",1:object$K)))
    for (i in 1:nrow(par)) {
      tmp <- object$rcpp.func$simahead(y = data, n = n.ahead, m = n.sim, par = par[i, ], P_0[i, ])
      if (object$K == 1L) {
        draw[,start:end]  <- t(tmp$draws)
        state[,start:end] <- matrix(0, nrow = n.ahead, ncol = n.sim)
        CondVol[,start:end,1] <- t(tmp$CondVol)
      } else {
        draw[,start:end]  <- t(tmp$draws)
        state[,start:end] <- t(tmp$state)
        CondVol[,start:end,] <- aperm(tmp$CondVol,perm = c(2,1,3))
      }
      start <- start + n.sim
      end   <- end + n.sim
    }
    rownames(draw) = rownames(state) = paste0("h=",1:n.ahead)
    colnames(draw) = colnames(state) =  paste0("Sim #",1:(n.sim * nrow(par)))
  }
  out <- list()
  out$draw <- draw
  out$state <- state + 1
  out$CondVol <- CondVol
  class(out) <- "MSGARCH_SIM"
  return(out)
}

#' @rdname Sim
#' @export
Sim.MSGARCH_ML_FIT <- function(object, new.data = NULL, n.ahead = 1L,
                               n.sim = 1L,  n.burnin = 500L, ...) {
  data <- c(object$data, new.data)
  out  <- Sim(object = object$spec, data = data, n.ahead = n.ahead,
              n.sim = n.sim, par = object$par, n.burnin = n.burnin)
  return(out)
}

#' @rdname Sim
#' @export
Sim.MSGARCH_MCMC_FIT <- function(object, new.data = NULL, n.ahead = 1L,
                                 n.sim = 1L, n.burnin = 500L, ...) {
  data <- c(object$data, new.data)
  out  <- Sim(object = object$spec, data = data, n.ahead = n.ahead,
              n.sim = n.sim, par = object$par, n.burnin = n.burnin)
  return(out)
}
