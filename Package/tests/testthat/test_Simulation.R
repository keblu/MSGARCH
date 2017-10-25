testthat::context("Test Simulation")

spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                            distribution.spec = list(distribution = c("norm")),
                            switch.spec = list(do.mix = FALSE, K = 2))

par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
         0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)

set.seed(1234)
sim <- MSGARCH::Sim(object = spec,n.ahead = 20,n.sim = 2,par = par,n.burnin = 10)

testthat::test_that("Simulation Draw", {
  
  tol <- 0.05
  est.draw <- sim$draw[20,]
  exp.draw <- c(1.08979324613423567, 0.75420478439175198)
  
  testthat::expect_true(max(abs(est.draw - exp.draw)) < tol)
})

testthat::test_that("Simulation State", {
  est.state <- sim$state[20,]
  exp.state <- c(1, 1)
  
  testthat::expect_true(all(est.state == exp.state))
})

testthat::test_that("Simulation CondVol", {
  
  tol <- 0.05
  est.condvol <- sim$CondVol[20,,]
  exp.condvol <- matrix(data = c(1.05850076010763927, 0.90701776672387635, 
                                 5.6188545586716181, 5.6137137219119113),
                        ncol = 2)
  
  testthat::expect_true(max(est.condvol - exp.condvol) < tol)
})
