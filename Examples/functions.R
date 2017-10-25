library("MSGARCH")

f_mc <- function(IN) {
  library("MSGARCH")
  data("SMI", package = "MSGARCH")
  spec <- MSGARCH::CreateSpec(variance.spec = list(model = "gjrGARCH"),
                              distribution.spec = list(distribution = "std"),
                              switch.spec = list(K = 2),
                              constraint.spec = list(regime.const = "nu"))
  
  ## function for custom optimization with DEoptim 
  f_DEoptim <- function(vPw, f_nll, spec, data, do.plm) {
    lower <- rep(-10, length(vPw)) # unmap space
    upper <- rep(10, length(vPw))  # unmap space
    
    tmp <- DEoptim::DEoptim(f_nll, lower = lower, upper = upper,
                            control = DEoptim::DEoptim.control(initialpop = NULL,
                                                               NP = 110, trace = TRUE, itermax = 500),
                            spec = spec, data = data, do.plm = do.plm)
    
    vPw_optim <- tmp$optim$bestmem
    dnllk     <- tmp$optim$bestval
    
    names(vPw_optim) <- names(vPw)
    out <- list(par = vPw_optim, value = dnllk)
    return(out)
  }
  
  set.seed(IN$seed)
  fit <- MSGARCH::FitML(spec, data = SMI, ctr = list(OptimFUN = f_DEoptim))
  out <- MSGARCH::State(fit)$FiltProb[,,2]
  return(out)
}