testthat::context("Test Extract State Fit")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                            distribution.spec = list(distribution = c("norm")),
                            switch.spec = list(do.mix = FALSE, K = 2))
fit <- MSGARCH::FitML(spec, data = SMI)

testthat::test_that("ExtractStateFit", {
 
  est.len = length(MSGARCH::ExtractStateFit(fit))
  exp.len = 2
  
  testthat::expect_true(est.len == exp.len)
})
