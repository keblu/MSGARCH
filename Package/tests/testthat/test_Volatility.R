testthat::context("Test Volatility")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                            distribution.spec = list(distribution = c("norm")),
                            switch.spec = list(do.mix = FALSE, K = 2))
par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
         0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)

testthat::test_that("Forecast", {
  
  tol <- 0.05
  set.seed(1234)
  est.forecast <- predict(object = spec, par = par, newdata = SMI,nahead = 2)$vol
  exp.forecast <- c(1.0304257211510406, 1.0340222685323162)
  
  testthat::expect_true(max(abs(est.forecast - exp.forecast)) < tol)
  
})

testthat::test_that("Forecast from an MCMC fit averages over the posterior draws", {

  y <- SMI[1:500]
  set.seed(1234)
  fit  <- MSGARCH::FitMCMC(spec, data = y,
                           ctr = list(nburn = 100L, nmcmc = 200L, nthin = 20L))
  mPar <- as.matrix(fit$par)

  # one-step-ahead volatility draw by draw, through the public interface; each of
  # these calls carries a single parameter vector, so it cannot depend on how the
  # draws are pooled afterwards
  vVol <- vapply(seq_len(nrow(mPar)), function(i) {
    as.numeric(predict(object = spec, par = mPar[i, ], newdata = y, nahead = 1L)$vol)
  }, FUN.VALUE = numeric(1))

  # guard: the posterior mean has to be distinguishable from the first draw,
  # otherwise this test proves nothing
  testthat::expect_true(abs(mean(vVol) - vVol[1]) > 1e-6)

  est.forecast <- as.numeric(predict(object = fit, nahead = 1L)$vol)
  testthat::expect_true(abs(est.forecast - mean(vVol)) < 1e-10)

})

testthat::test_that("Conditional Vol", {

  tol <- 0.05
  est.Vol <- Volatility(object = spec, par = par, data = SMI)[2000]
  exp.Vol <- c(2.1321725800180471)

  testthat::expect_true(max(abs(est.Vol - exp.Vol)) < tol)

})