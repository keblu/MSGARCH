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
  est.forecast <- MSGARCH::Forecast(object = spec,par = par,data = SMI,n.ahead = 2)$vol
  exp.forecast <- c(1.0304257211510406, 1.0340222685323162)
  
  testthat::expect_true(max(abs(est.forecast - exp.forecast)) < tol)
  
})

testthat::test_that("Conditional Vol", {
  
  tol <- 0.05
  est.Vol <- MSGARCH::Volatility(object = spec, par = par, data = SMI)[2000]
  exp.Vol <- c(2.1321725800180471)
  
  testthat::expect_true(max(abs(est.Vol - exp.Vol)) < tol)
  
})