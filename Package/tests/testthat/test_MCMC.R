testthat::context("Test MCMC Estimation")

set.seed(1234)
data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")),
                   switch.spec = list(do.mix = FALSE, K = 2))
fit <- MSGARCH::FitMCMC(spec, data = SMI, ctr = list(n.burn = 500, n.mcmc = 500, n.thin = 1))

testthat::test_that("MCMC Estimation MS GARCH NORMAL", {
  
  tol <- 0.05
  est.par <- colMeans(fit$par)
  exp.par <- c(0.02483046, 0.05832523, 0.90097897, 0.46547319, 
               0.02984767, 0.89541786, 0.97862471, 0.24272290)
  
  testthat::expect_true(max(abs(est.par - exp.par)) < tol)
  
})

testthat::test_that("DIC", {
  
  tol <- 0.1
  est.par <- DIC(fit)$DIC
  exp.par <- 6825.8205675690924
  
  testthat::expect_true(abs(est.par - exp.par) < tol)
  
})