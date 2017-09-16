#' @import adaptMCMC
f_SamplerFUNDefault <- function(f_posterior, data, spec, par0, ctr) {
  out <- list()
  out <- adaptMCMC::MCMC(p = f_posterior, data = data, spec = spec, 
                         n = ctr$n.burn + ctr$n.mcmc, init = par0, adapt = TRUE,
                         acc.rate = 0.23)$samples
  return(out)
}
