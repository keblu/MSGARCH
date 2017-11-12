#' @title Simulation of MSGARCH processes.
#' @description  Method for simulating \code{MSGARCH} processes.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{CreateSpec}}
#' or fit object of type \code{MSGARCH_ML_FIT} created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT}
#' created with \code{\link{FitMCMC}}.
#' @param nsim Number of simulations. (Default: \code{nsim = 1L})
#' @param nahead Simulation length. (Default: \code{nahead = 1L})
#' @param nburn Burnin period discarded (first simulation draws).
#' @param par Vector (of size d) or matrix (of size \code{nahead} x d) of parameter
#' @param seed 	Integer indicateing if and how the random number generator should be initialized. 
#' If \code{seed = NULL},the state of the random generator will not change. (Default: \code{seed = NULL})
#' estimates where d must have the same length as the default parameters of the specification.
#' @param ... Not used. Other arguments to \code{simulate}.
#' @return A list of class \code{MSGARCH_SIM} with the following elements:.
#' \itemize{
#' \item \code{draw}: Matrix (of size \code{nahead} x \code{nsim}) of simulated draws.
#' \item \code{state}: Matrix (of size \code{nahead} x \code{nsim}) of simulated states.
#' \item \code{CondVol}: Array (of size \code{nahead} x \code{nsim} x K) of simulated conditional volatility.  
#' }
#' The \code{MSGARCH_SIM} class contains the \code{plot} method.
#' @details If a matrix of parameters estimates is provided, \code{nsim} simuations will be done for each row.
#' @examples
#' # create model specification
#' # MS(2)-GARCH(1,1)-Normal (default)
#' spec <- CreateSpec()
#'
#' # generate process
#' par.sim <- c(0.1,0.6,0.2,0.2,0.8,0.1,0.99,0.01)
#' set.seed(123)
#' sim <- simulate(object = spec, nsim = 1L, nahead = 1000L, nburnin = 500L, par = par.sim)
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
#' sim <- simulate(fit, nsim = 1000L, nahead = 30L)
#' plot(sim)
#' @rdname simulate

#' @rdname simulate
#' @export
simulate.MSGARCH_SPEC  <- function(object, nsim = 1L, seed = NULL, nahead = 1L,
                                   par = NULL, nburn = 500L, ...) {
  out <- Sim(object = object, data = NULL, nahead = nahead, 
            nsim = nsim, par = par, nburn = nburn, seed = seed)
  return(out)
}

#' @rdname simulate
#' @export
simulate.MSGARCH_ML_FIT <- function(object, nsim = 1L, seed = NULL, nahead = 1L,
                                nburn = 500L, ...) {
  out <- Sim(object = object$spec, data = NULL, nahead = nahead,
              nsim = nsim, par = object$par, nburn = nburn, seed = seed)
  return(out)
}

#' @rdname simulate
#' @export
simulate.MSGARCH_MCMC_FIT <- function(object, nsim = 1L, seed = NULL, nahead = 1L,
                                      nburn = 500L, ...) {
  out <- Sim(object = object$spec, data = NULL, nahead = nahead,
              nsim = nsim, par = object$par, nburn = nburn, seed = seed)
  return(out)
}


#For internal use and simulate function
Sim <- function(object, data = NULL, nahead = 1L,
                nsim = 1L, par = NULL, nburn = 500L, seed = NULL, ...) {
  UseMethod(generic = "Sim", object)
}

Sim.MSGARCH_SPEC <- function(object, data = NULL, nahead = 1L,
                             nsim = 1L, par = NULL, nburn = 500L, seed = NULL, ...) {
  
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)
  }
  if (is.null(seed)){ 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  object <- f_check_spec(object)
  if (is.vector(par)) {
    par <- matrix(par, nrow = 1L)
  }
  # New simulation
  if (is.null(data)) {
    par   <- f_check_par(object, par)
    start <- 1
    end   <- nsim
    draw  <- matrix(data = NA, nrow = nahead + nburn, ncol = nsim * nrow(par) )
    state <- matrix(data = NA, nrow = nahead + nburn, ncol = nsim * nrow(par) )
    CondVol <- array(data = NA, dim = c(nahead + nburn, nsim * nrow(par), object$K),
                     dimnames =  list(paste0("t=",1:(nahead+nburn)),
                                      paste0("Sim #",1:(nsim * nrow(par))),paste0("k=",1:object$K)))
    for (i in 1:nrow(par)) {
      tmp <- object$rcpp.func$sim(nahead  + nburn, nsim, par[i, ])
      if (object$K == 1L) {
        draw[,start:end]  <- t(tmp$draws)
        state[,start:end] <- matrix(0, nrow = nahead + nburn, ncol = nsim)
        CondVol[,start:end,] <- t(tmp$CondVol)
      } else {
        draw[,start:end]  <- t(tmp$draws)
        state[,start:end] <- t(tmp$state)
        CondVol[,start:end,] <- aperm(tmp$CondVol,perm = c(2,1,3))
      }
      start <- start + nsim
      end   <- end + nsim
    }
    draw  <- draw[-(1:nburn),,drop = FALSE]
    state <- state[-(1:nburn),,drop = FALSE]
    CondVol <- CondVol[-(1:nburn),,,drop = FALSE]
    rownames(draw) = rownames(state) = paste0("t=",1:nahead)
    colnames(draw) = colnames(state) =  paste0("Sim #",1:(nsim * nrow(par)))
    dimnames(CondVol)[[1]] = paste0("t=",1:nahead)
  } else {
    # Simulation ahead of data
    data  <- f_check_y(data)
    P_0   <- matrix(State(object, par = par, data = data)$PredProb[(length(data) + 1L), ,], ncol = object$K)
    par   <- f_check_par(object, par)
    start <- 1
    end   <- nsim
    draw  <- matrix(data = NA, nrow = nahead, ncol =  nsim * nrow(par))
    state <- matrix(data = NA, nrow = nahead, ncol =  nsim * nrow(par))
    CondVol <- array(data = NA, dim = c(nahead, nsim * nrow(par), object$K),
                     dimnames =  list(paste0("h=",1:(nahead)),
                                      paste0("Sim #",1:(nsim * nrow(par))),paste0("k=",1:object$K)))
    for (i in 1:nrow(par)) {
      tmp <- object$rcpp.func$simahead(y = data, n = nahead, m = nsim, par = par[i, ], P_0[i, ])
      if (object$K == 1L) {
        draw[,start:end]  <- t(tmp$draws)
        state[,start:end] <- matrix(0, nrow = nahead, ncol = nsim)
        CondVol[,start:end,1] <- t(tmp$CondVol)
      } else {
        draw[,start:end]  <- t(tmp$draws)
        state[,start:end] <- t(tmp$state)
        CondVol[,start:end,] <- aperm(tmp$CondVol,perm = c(2,1,3))
      }
      start <- start + nsim
      end   <- end + nsim
    }
    rownames(draw) = rownames(state) = paste0("h=",1:nahead)
    colnames(draw) = colnames(state) =  paste0("Sim #",1:(nsim * nrow(par)))
  }
  out <- list()
  out$draw <- draw
  out$state <- state + 1
  out$CondVol <- CondVol
  class(out) <- "MSGARCH_SIM"
  return(out)
}

Sim.MSGARCH_ML_FIT <- function(object, newdata = NULL, nahead = 1L,
                               nsim = 1L,  nburn = 500L, seed = NULL, ...) {
  data <- c(object$data, newdata)
  out  <- Sim(object = object$spec, data = data, nahead = nahead,
              nsim = nsim, par = object$par, nburn = nburn, seed = seed)
  return(out)
}

Sim.MSGARCH_MCMC_FIT <- function(object, newdata = NULL, nahead = 1L,
                                 nsim = 1L, nburn = 500L, seed = NULL, ...) {
  data <- c(object$data, newdata)
  out  <- Sim(object = object$spec, data = data, nahead = nahead,
              nsim = nsim, par = object$par, nburn = nburn, seed = seed)
  return(out)
}
