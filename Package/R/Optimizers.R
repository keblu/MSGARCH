#' @importFrom stats optim
f_OptimFUNDefault <- function(vPw, f_nll, spec, data, do.plm) {
  out <- try(stats::optim(vPw, f_nll, spec = spec, data = data, do.plm = do.plm, hessian = TRUE, method = "BFGS"), silent = TRUE)
  return(out)
}

#' #' @import Rsolnp
#' f_solnp <- function(vPw, f_nll, spec, y, do.plm) {
#'
#'   tmp <- try(Rsolnp::solnp(vPw, f_nll, spec = spec, y = y, do.plm = do.plm), silent = TRUE)
#'
#'   out = list(par = tmp$pars, value = tail(tmp$values, 1), hessian = tmp$hessian)
#'
#'   return(out)
#' }
#'
#' #' @import DEoptim
#' #' @importFrom stats rnorm
#' f_DEoptim <- function(vPw, f_nll, spec, y, do.plm) {
#'
#'   NP = 15 * length(vPw)
#'
#'   mInitialPop = matrix(rep(vPw, NP) * rnorm(length(vPw) * NP) * 0.1, nrow = NP, byrow = TRUE)
#'
#'   tmp <- try(DEoptim::DEoptim(f_nll, lower = rep(-10, length(vPw)), upper = rep(10, length(vPw)),
#'   control = DEoptim.control(initialpop = mInitialPop, NP = NP, trace = TRUE, itermax = 500),
#'                               spec = spec, y = y, do.plm = do.plm), silent = TRUE)
#'
#'   vPw_optim <- tmp$optim$bestmem
#'   dnllk <- tmp$optim$bestval
#'
#'   names(vPw_optim) = names(vPw)
#'
#'   out = list(par = vPw_optim, value = dnllk)
#'
#'   return(out)
#' }
