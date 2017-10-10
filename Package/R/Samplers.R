#' @import adaptMCMC
f_SamplerFUNDefault <- function(f_posterior, data, spec, par0, ctr) {
  p.log <- function(vPw) {
    return(f_posterior(vPw = vPw, data = data, spec = spec, PriorFun = TRUE))
  }
  draw <- f_RCPP_adaptMCMC(theta0 = par0, acc_rate = 0.25, 
                           sigma = diag(length(par0)), 
                           func = p.log, n_mcmc = ctr$n.burn + ctr$n.mcmc)
  colnames(draw) = names(par0)
  return(draw)
}
