testthat::context("Test State")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")),
                   switch.spec = list(do.mix = FALSE, K = 2))

par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
         0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)

est.state <- MSGARCH::State(object = spec, x = 0, data = SMI, par = par, n.ahead = 2)

testthat::test_that("Preditive Prob", {
 
  tol <- 0.05
  est.Pred <- est.state$PredProb[2000,,]
  exp.Pred <- c(0.97874614696735351, 0.02125385303264642)

  testthat::expect_true(max(abs(est.Pred - exp.Pred)) < tol)
  
})

testthat::test_that("Filtered Prob", {
  
  tol <- 0.05
  est.Filt <- est.state$FiltProb[2000,,]
  exp.Filt <- c(0.979703177331200159, 0.020296822668799768)
  
  testthat::expect_true(max(abs(est.Filt - exp.Filt)) < tol)
  
})

testthat::test_that("Smoothed Prob", {
  
  tol <- 0.05
  est.Smoo <- est.state$SmoothProb[2000,,]
  exp.Smoo <- c(0.979635245161457835, 0.020364754838541596 )
  
  testthat::expect_true(max(abs(est.Smoo - exp.Smoo)) < tol)
  
})

testthat::test_that("Viterbi", {
  
  est.Vit <- est.state$Viterbi[2000]
  exp.Vit <- c(1)
  
  testthat::expect_true(all(est.Vit == exp.Vit))
  
})