f.error = function(message) {
  cat(paste0("try/catch ", message, "\n"))
  return(FALSE)
}


f.process.ctr = function(ctr = list()) {
  con = list(theta0 = NULL, do.init = FALSE, N.mcmc = 5000, N.burn = 1000, N.thin = 10)
  con[names(ctr)] = ctr
  return(con)
}
