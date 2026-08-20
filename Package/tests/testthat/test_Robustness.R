testthat::context("Test numerical robustness and input validation")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")),
                            distribution.spec = list(distribution = c("norm", "norm")))

testthat::test_that("The filter survives extreme between-regime likelihood ratios", {

  # One tight regime and one loose one: at y = 1 their log densities differ by
  # far more than the ~1400 that used to overflow the exponential and turn the
  # filtered probabilities into NaN.
  y   <- c(0, 1)
  par <- c(1e-5, 1e-6, 0.001, 1, 1e-6, 0.001, 0.9, 0.1)
  h   <- c(1e-5, 1) / (1 - 1e-6 - 0.001)
  testthat::expect_true(abs(diff(stats::dnorm(1, 0, sqrt(h), log = TRUE))) > 1000)

  dLLK <- MSGARCH:::Kernel(spec, par, y, log = TRUE, do.prior = FALSE)
  testthat::expect_true(is.finite(dLLK))

  # regime 1 is impossible at this observation, so the likelihood collapses onto
  # regime 2 weighted by its predictive probability (here the ergodic 1/2)
  testthat::expect_true(abs(dLLK - (stats::dnorm(1, 0, sqrt(h[2]), log = TRUE) + log(0.5))) < 1e-8)

  mProb <- State(object = spec, par = par, data = y)$FiltProb
  testthat::expect_true(all(is.finite(mProb)))
  testthat::expect_true(max(abs(apply(mProb, 1L, sum) - 1)) < 1e-12)

  testthat::expect_true(is.finite(as.numeric(PredPdf(object = spec, par = par, x = 0, data = y))))

})

testthat::test_that("A singular transition matrix is rejected, not thrown from C++", {

  # the plain mapping bounds each free transition probability separately, so a
  # row can leave the simplex and make the stationary-distribution solve singular
  spec3 <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                               switch.spec = list(do.mix = FALSE, K = 3))
  par <- spec3$par0
  par[grep("^P_", names(par))] <- c(1, 0.999973, 0, 1, 0, 0)

  testthat::expect_true(is.finite(MSGARCH:::Kernel(spec3, par, SMI[1:200], log = TRUE,
                                                   do.prior = FALSE)))
  testthat::expect_equal(MSGARCH:::Kernel(spec3, par, SMI[1:200], log = TRUE,
                                          do.prior = TRUE), -1e10)

  # a constrained three-regime fit used to die on such a point mid-optimization
  spec3c <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                distribution.spec = list(distribution = c("std")),
                                switch.spec = list(do.mix = FALSE, K = 3),
                                constraint.spec = list(regime.const = c("nu")))
  set.seed(1234)
  fit <- MSGARCH::FitML(spec3c, data = SMI, ctr = list(do.se = FALSE))
  testthat::expect_true(is.finite(fit$loglik))

})

testthat::test_that("Data and parameters must be finite and long enough", {

  y <- SMI[1:100]
  for (bad in list(replace(y, 10L, NA), replace(y, 10L, NaN), replace(y, 10L, Inf))) {
    testthat::expect_error(Volatility(object = spec, par = spec$par0, data = bad))
    testthat::expect_error(MSGARCH::FitML(spec, data = bad))
  }
  testthat::expect_error(Volatility(object = spec, par = spec$par0, data = y[1L]))
  testthat::expect_error(Volatility(object = spec, par = replace(spec$par0, 1L, NA), data = y))

  # a well-formed call must still work
  testthat::expect_true(all(is.finite(as.numeric(Volatility(object = spec, par = spec$par0,
                                                            data = y)))))

})

testthat::test_that("simulate() accepts a zero burn-in", {

  set.seed(1234)
  sim0 <- simulate(object = spec, nsim = 2L, nahead = 3L, par = spec$par0, nburn = 0L)
  testthat::expect_equal(dim(sim0$draw), c(3L, 2L))

  set.seed(1234)
  sim5 <- simulate(object = spec, nsim = 2L, nahead = 3L, par = spec$par0, nburn = 5L)
  testthat::expect_equal(dim(sim5$draw), c(3L, 2L))

})

testthat::test_that("UncVol averages the horizons after the burn-in", {

  set.seed(1234)
  fit <- MSGARCH::FitML(spec, data = SMI, ctr = list(do.se = FALSE))
  ctr <- list(nsim = 200L, nburn = 20L, nahead = 5L)

  set.seed(7)
  est <- UncVol(object = fit, ctr = ctr)
  set.seed(7)
  vol <- MSGARCH:::f_CondVol(fit$spec,
                             matrix(fit$par, nrow = 1L, dimnames = list(NULL, names(fit$par))),
                             data = c(1, 1), do.its = FALSE, nahead = 25L,
                             ctr = list(nsim = 200L))$vol

  testthat::expect_true(abs(est - mean(vol[21:25])) < 1e-12)
  testthat::expect_true(abs(est - mean(vol[20:5]))  > 1e-8)

})

testthat::test_that("The mixture Viterbi path is the per-observation MAP decoding", {

  # For a mixture every row of the transition matrix is the same weight vector,
  # so Viterbi decoding collapses to maximising log(w_k) + log f_k(y_t) at each
  # observation. (The matrix State() builds used to be filled column-major and
  # was therefore not row-stochastic; that is fixed, but note it never changed
  # the decoded path, because the misplaced factor is constant in the index
  # being maximised over and so cancels.)
  spec.mix <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")),
                                  switch.spec = list(do.mix = TRUE))
  par <- spec.mix$par0
  par["P_1"]      <- 0.8
  par["alpha0_2"] <- 2
  y <- SMI[1:300]

  vPath <- as.numeric(State(object = spec.mix, par = par, data = y)$Viterbi)

  mPar <- spec.mix$func$f.do.mix(matrix(par, nrow = 1L))
  aHt  <- spec.mix$rcpp.func$calc_ht(mPar, y)
  vW   <- c(par["P_1"], 1 - par["P_1"])
  mLL  <- sapply(seq_len(spec.mix$K), function(k) {
    log(vW[k]) + stats::dnorm(y[-1L], 0, sqrt(aHt[2:length(y), 1L, k]), log = TRUE)
  })
  vRef <- apply(mLL, 1L, which.max)
  vRef <- c(vRef[1L], vRef)   # State() copies the first state, see its comment

  testthat::expect_equal(vPath, as.numeric(vRef))
  testthat::expect_true(all(vPath %in% seq_len(spec.mix$K)))

})

testthat::test_that("CreateSpec validates K and the regime it expands", {

  testthat::expect_error(MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                             distribution.spec = list(distribution = c("norm", "std")),
                                             switch.spec = list(K = 3)))
  testthat::expect_error(MSGARCH::CreateSpec(switch.spec = list(K = 2.5)))
  testthat::expect_error(MSGARCH::CreateSpec(switch.spec = list(K = 0)))

  # the documented expansion must keep working
  spec3 <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                               distribution.spec = list(distribution = c("std")),
                               switch.spec = list(do.mix = FALSE, K = 3))
  testthat::expect_equal(spec3$name, rep("sGARCH_std", 3L))

})

testthat::test_that("Risk validates its arguments and warns on an inadequate grid", {

  set.seed(1234)
  fit <- MSGARCH::FitML(spec, data = SMI[1:400], ctr = list(do.se = FALSE))
  testthat::expect_error(Risk(fit, alpha = 0))
  testthat::expect_error(Risk(fit, alpha = 1.5))
  testthat::expect_error(Risk(fit, nahead = 0L))
  testthat::expect_error(Risk(fit, ctr = list(nmesh = 1L)))

  # the grid spans the observed range, which here holds a sliver of the
  # predictive mass, so the tail cannot be resolved
  spec.sr <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                 distribution.spec = list(distribution = c("norm")),
                                 switch.spec = list(do.mix = FALSE, K = 1))
  testthat::expect_warning(Risk(object = spec.sr, par = c(1, 0, 0),
                                data = c(-0.1, 0.1, -0.05, 0.05), alpha = 0.01))

  # and no warning when the grid does cover the distribution
  testthat::expect_silent(Risk(fit, alpha = 0.01))

})

testthat::test_that("Prior means and standard deviations are validated", {

  testthat::expect_error(MSGARCH::CreateSpec(prior = list(sd = list(beta_1 = 0))))
  testthat::expect_error(MSGARCH::CreateSpec(prior = list(sd = list(beta_1 = -1))))
  testthat::expect_error(MSGARCH::CreateSpec(prior = list(sd = list(wrong_name = 1))))
  testthat::expect_error(MSGARCH::CreateSpec(prior = list(mean = list(beta_1 = Inf))))

  spec.p <- MSGARCH::CreateSpec(prior = list(mean = list(beta_1 = 0.7),
                                             sd = list(beta_1 = 0.1)))
  testthat::expect_equal(unname(spec.p$rcpp.func$get_sd()[3L]), 0.1)

})

testthat::test_that("Fixed parameters are honoured and transition probabilities refused", {

  testthat::expect_error(MSGARCH::CreateSpec(constraint.spec = list(fixed = list(P_1_1 = 0.99))))
  testthat::expect_error(MSGARCH::CreateSpec(constraint.spec = list(fixed = list(beta_1 = NA))))

  spec.f <- MSGARCH::CreateSpec(constraint.spec = list(fixed = list(beta_1 = 0.8)))

  set.seed(1234)
  fit <- MSGARCH::FitML(spec.f, data = SMI[1:400], ctr = list(do.se = FALSE))
  testthat::expect_true(abs(fit$par["beta_1"] - 0.8) < 1e-12)

  # a user-supplied par0 must not shift the sampler's parameter vector, and the
  # identification sort must not relabel the regime the constraint refers to
  set.seed(1234)
  mcmc <- suppressMessages(MSGARCH::FitMCMC(spec.f, data = SMI[1:400],
                                            ctr = list(par0 = spec.f$par0, nburn = 50L,
                                                       nmcmc = 100L, nthin = 1L)))
  mPar <- as.matrix(mcmc$par)
  testthat::expect_equal(ncol(mPar), length(spec.f$label))
  testthat::expect_true(max(abs(mPar[, "beta_1"] - 0.8)) < 1e-12)

})
