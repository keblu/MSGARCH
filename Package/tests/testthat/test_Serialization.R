testthat::context("Test Serialization (saveRDS / readRDS round trip)")

# A MSGARCH_SPEC carries Rcpp module objects. R serializes external pointers as NULL,
# so any spec or fit written with saveRDS comes back with dead pointers and has to be
# rebuilt by f_check_spec. That happens in one session too, which is what these tests
# exploit; it is the same failure a user hits when reloading an overnight fit or when
# shipping a spec to a parallel worker.

data("SMI", package = "MSGARCH")

f_roundtrip <- function(object) {
  sFile <- tempfile(fileext = ".rds")
  on.exit(unlink(sFile))
  saveRDS(object, file = sFile)
  return(readRDS(sFile))
}

spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                            distribution.spec = list(distribution = c("norm")),
                            switch.spec = list(do.mix = FALSE, K = 2))
par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566,
         0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)

testthat::test_that("The round trip really does invalidate the Rcpp pointers", {

  # if this ever stops holding the tests below become vacuous
  spec.rt <- f_roundtrip(spec)
  testthat::expect_error(spec.rt$rcpp.func$get_sd())

})

testthat::test_that("A reloaded spec is usable and gives identical results", {

  spec.rt <- f_roundtrip(spec)

  exp.vol <- Volatility(object = spec,    par = par, data = SMI)
  est.vol <- Volatility(object = spec.rt, par = par, data = SMI)
  testthat::expect_true(max(abs(as.numeric(est.vol) - as.numeric(exp.vol))) < 1e-12)

  exp.llk <- MSGARCH:::Kernel(spec,    par, SMI, log = TRUE, do.prior = FALSE)
  est.llk <- MSGARCH:::Kernel(spec.rt, par, SMI, log = TRUE, do.prior = FALSE)
  testthat::expect_true(abs(est.llk - exp.llk) < 1e-12)

  exp.state <- State(object = spec,    par = par, data = SMI)$SmoothProb
  est.state <- State(object = spec.rt, par = par, data = SMI)$SmoothProb
  testthat::expect_true(max(abs(est.state - exp.state)) < 1e-12)

})

testthat::test_that("A reloaded spec keeps its user-supplied priors", {

  spec.prior <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                    distribution.spec = list(distribution = c("norm")),
                                    switch.spec = list(do.mix = FALSE, K = 2),
                                    prior = list(mean = list(beta_1 = 0.7),
                                                 sd   = list(beta_1 = 0.1)))
  exp.mean <- spec.prior$rcpp.func$get_mean()
  exp.sd   <- spec.prior$rcpp.func$get_sd()

  spec.rt <- MSGARCH:::f_check_spec(f_roundtrip(spec.prior))

  testthat::expect_true(max(abs(spec.rt$rcpp.func$get_mean() - exp.mean)) < 1e-12)
  testthat::expect_true(max(abs(spec.rt$rcpp.func$get_sd()   - exp.sd))   < 1e-12)

})

testthat::test_that("A reloaded MSGARCH_ML_FIT is usable", {

  fit    <- MSGARCH::FitML(spec, data = SMI[1:500])
  fit.rt <- f_roundtrip(fit)

  testthat::expect_true(max(abs(as.numeric(Volatility(fit.rt)) -
                                as.numeric(Volatility(fit)))) < 1e-12)

  set.seed(1234)
  exp.pred <- predict(object = fit,    nahead = 1L)$vol
  set.seed(1234)
  est.pred <- predict(object = fit.rt, nahead = 1L)$vol
  testthat::expect_true(abs(as.numeric(est.pred) - as.numeric(exp.pred)) < 1e-12)

  testthat::expect_true(abs(AIC(fit.rt) - AIC(fit)) < 1e-12)
  testthat::expect_silent(invisible(capture.output(summary(fit.rt))))

})

testthat::test_that("A reloaded MSGARCH_MCMC_FIT is usable", {

  set.seed(1234)
  fit    <- MSGARCH::FitMCMC(spec, data = SMI[1:500],
                             ctr = list(nburn = 100L, nmcmc = 100L, nthin = 1L))
  fit.rt <- f_roundtrip(fit)

  testthat::expect_true(max(abs(as.numeric(Volatility(fit.rt)) -
                                as.numeric(Volatility(fit)))) < 1e-12)
  testthat::expect_true(abs(DIC(fit.rt)$DIC - DIC(fit)$DIC) < 1e-12)

})
