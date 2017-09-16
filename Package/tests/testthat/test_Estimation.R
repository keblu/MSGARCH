testthat::context("Test Estimation")

tol <- 1e-4

testthat::test_that("Estimation MSGARCH", {

  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = FALSE, K = 2))
  fit <- FitML(spec = spec, data = SMI)
  bic <- BIC(fit)
  tmp <- abs(bic - 6841.18485426964)
  test1 <- tmp < tol
  aic <- AIC(fit)
  tmp <- abs(aic - 6794.59248618279)
  test2 <- tmp < tol
  testthat::expect_true(test1 & test2)
})

testthat::test_that("Estimation MIXGARCH", {
  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = TRUE, K = 2))
  fit <- FitML(spec = spec, data = SMI)
  bic <- BIC(fit)
  tmp <- abs(bic - 6838.9210075872)
  test1 <- tmp < tol
  aic <- AIC(fit)
  tmp <- abs(aic - 6798.15268498273)
  test2 <- tmp < tol
  testthat::expect_true(test1 & test2)
})
