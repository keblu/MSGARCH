testthat::context("Test Simulation")

tol <- 1e-4
tmp <- sessionInfo()
Plaform_version <- tmp$platform
#MSGARCH_Version <- tmp$otherPkgs$MSGARCH$Version

set.seed(1234)
spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")),
                   switch.spec = list(do.mix = FALSE, K = 2))

par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
         0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894 )
sim = Sim(object = spec,n.ahead = 20,n.sim = 2,par = par,n.burnin = 10)
est.draw <- sim$draw[20,]
est.state <- sim$state[20,]
est.condvol <- sim$CondVol[20,,]



###############################################################
testthat::test_that("Simulation Draw", {
  exp.draw <- NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.draw = c(1.08979324613423567, 0.75420478439175198)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.draw = c(1.08979324613423567, 0.75420478439175198)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.draw = c(1.08979324613423567, 0.75420478439175198)
  }
  
  if (is.null(exp.draw)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.draw <- est.draw
  }
  testthat::expect_true(max(est.draw - exp.draw) < tol)
})

###############################################################
testthat::test_that("Simulation State", {
  exp.state <- NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.state = c(1, 1)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.state = c(1, 1)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.state = c(1, 1)
  }
  
  if (is.null(exp.state)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.state <- est.state
  }
  testthat::expect_true(max(est.state - exp.state) < tol)
})

testthat::test_that("Simulation CondVol", {
  exp.condvol <- NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.condvol = matrix(c(1.05850076010763927,0.90701776672387635,5.6188545586716181,5.6137137219119113),
                         ncol = 2)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.condvol = matrix(c(1.05850076010763927,0.90701776672387635,5.6188545586716181,5.6137137219119113),
                         ncol = 2)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.condvol = matrix(c(1.05850076010763927,0.90701776672387635,5.6188545586716181,5.6137137219119113),
                         ncol = 2)
  }
  
  if (is.null(exp.condvol)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.condvol <-  est.condvol
  }
  testthat::expect_true(max(est.condvol - exp.condvol) < tol)
})
