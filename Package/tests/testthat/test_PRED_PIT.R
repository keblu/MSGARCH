testthat::context("Test Pred and PIT")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")),
                   switch.spec = list(do.mix = FALSE, K = 2))

par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
         0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)

testthat::test_that("Pred", {
  tol <- 0.05
  set.seed(1234)
  est.pred <- as.vector(MSGARCH::Pred(object = spec, x = 0, data = SMI, par = par, n.ahead = 2))
  exp.pred <- c(0.39670008854976474, 0.39453144146303032)

  testthat::expect_true(max(abs(est.pred - exp.pred)) < tol)
  
})

testthat::test_that("PIT", {
  tol <- 0.05
  set.seed(1234)
  est.PIT <- as.vector(MSGARCH::PIT(object = spec, x = 1, data = SMI, par = par, n.ahead = 2))
  exp.PIT <- c(0.83930591484143302, 0.83970000000000000)
  
  testthat::expect_true(max(abs(est.PIT - exp.PIT)) < tol)
})
