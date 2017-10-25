testthat::context("Test MLE Estimation")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")),
                   switch.spec = list(do.mix = FALSE, K = 2))
fit <- MSGARCH::FitML(spec, data = SMI)
summary(fit)

testthat::test_that("Estimation MS GARCH NORMAL", {
  
  tol <- 0.05
  est.par <- fit$par
  exp.par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
               0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894 )
  
  testthat::expect_true(max(abs(est.par - exp.par)) < tol)
  
})

testthat::test_that("Estimation BIC", {
  
  tol <- 0.1
  est.BIC <- MSGARCH::BIC(fit)
  exp.BIC <- 6841.1848542696416
  
  testthat::expect_true(abs(exp.BIC - exp.BIC) < tol)
  
})


testthat::test_that("Estimation AIC", {
  
  tol <- 0.1
  est.AIC <- MSGARCH::AIC(fit)
  exp.AIC <- 6794.5924861827916
  
  testthat::expect_true(abs(est.AIC - exp.AIC) < tol)
  
})

testthat::test_that("Estimation SR ARCH NORMAL", {
  
  spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- MSGARCH::FitML(spec, data = SMI)
  
  tol <- 0.05
  est.par <- fit$par
  exp.par <- c(0.80504578435066854, 0.31371193670568287)
  
  testthat::expect_true(max(abs(est.par - exp.par)) < tol)
  
})


testthat::test_that("Estimation SR GARCH NORMAL", {
  
  spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- MSGARCH::FitML(spec, data = SMI)
  
  tol <- 0.05
  est.par <- fit$par
  exp.par <- c(0.085912610592369154, 0.124325097982214264, 0.796421589649365935)
  
  testthat::expect_true(max(abs(est.par - exp.par)) < tol)
  
})

testthat::test_that("Estimation SR GJR NORMAL", {
  
  spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- MSGARCH::FitML(spec, data = SMI)
  
  tol <- 0.05
  est.par <- fit$par
  exp.par <- c(0.096745854166231729, 0.042005697086443979, 
               0.174648268027545489, 0.786944835620325578)
  
  testthat::expect_true(max(abs(est.par - exp.par)) < tol)
  
})


testthat::test_that("Estimation SR GJR STUDENT", {
  
  spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                     distribution.spec = list(distribution = c("std")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- MSGARCH::FitML(spec, data = SMI)
  
  tol <- 0.05
  est.par <- fit$par
  exp.par <- c(0.042068160604000746, 0.041341103857468385, 
               0.123041130026077569, 0.862000796059737562, 
               8.408571456679910128)
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})

