testthat::context("Test Inference (standard errors, p-values, degrees of freedom)")

data("SMI", package = "MSGARCH")

# Single-regime GARCH(1,1)-Normal is the cleanest probe for the orientation of the
# delta method: the working -> natural map is triangular (the upper bound on beta is
# 0.9999 - alpha1, so d beta / d alpha1_tilde != 0 while d alpha1 / d beta_tilde == 0)
# and no parameter sits near a bound, so the observed information in the natural
# parameterisation is well conditioned and can be used as an independent reference.
spec.sr <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                               distribution.spec = list(distribution = c("norm")),
                               switch.spec = list(do.mix = FALSE, K = 1))
fit.sr <- MSGARCH::FitML(spec.sr, data = SMI)

f_jacob <- function(fit) {
  vPw <- MSGARCH:::f_unmapPar(fit$par, fit$spec, fit$ctr$do.plm)
  numDeriv::jacobian(MSGARCH:::f_mapPar, vPw, spec = fit$spec, do.plm = fit$ctr$do.plm)
}

testthat::test_that("Standard errors use J V J', not its transpose", {

  J  <- f_jacob(fit.sr)
  V  <- MASS::ginv(fit.sr$Inference$Hessian)
  se <- unname(fit.sr$Inference$MatCoef[, "Std. Error"])

  se.delta      <- sqrt(diag(J %*% V %*% t(J)))
  se.transposed <- sqrt(diag(t(J) %*% V %*% J))

  # guard against a vacuous test: the two orientations must actually differ here
  testthat::expect_true(max(abs(se.delta - se.transposed)) > 1e-4)

  testthat::expect_true(max(abs(se - se.delta)) < 1e-8)

})

testthat::test_that("Standard errors match the observed information in the natural scale", {

  f_nll_natural <- function(vPn) {
    -MSGARCH:::Kernel(fit.sr$spec, vPn, SMI, log = TRUE, do.prior = FALSE)
  }

  # central-difference Hessian with a small step, to stay away from the constraint
  # boundaries at which the kernel is floored to -1e10
  vPn  <- fit.sr$par
  d    <- length(vPn)
  step <- pmax(abs(vPn), 1) * 1e-4
  mH   <- matrix(data = 0, nrow = d, ncol = d)
  for (i in 1:d) {
    for (j in i:d) {
      ei <- rep(0, d); ei[i] <- step[i]
      ej <- rep(0, d); ej[j] <- step[j]
      mH[i, j] <- mH[j, i] <- (f_nll_natural(vPn + ei + ej) - f_nll_natural(vPn + ei - ej)
                               - f_nll_natural(vPn - ei + ej) + f_nll_natural(vPn - ei - ej)) /
                              (4 * step[i] * step[j])
    }
  }
  se.natural <- sqrt(diag(solve(mH)))
  se         <- unname(fit.sr$Inference$MatCoef[, "Std. Error"])

  # agreement is ~1e-4 in relative terms; the transposed sandwich is off by 44% and 85%
  testthat::expect_true(max(abs(se / se.natural - 1)) < 0.02)

})

testthat::test_that("Transition-probability standard errors are not swapped across regimes", {

  # The P-block of the Jacobian is anti-diagonal (f_mapGamma enumerates the
  # off-diagonal entries column-major but carries the row-major parameter names),
  # so transposing the sandwich exchanges the two transition probabilities' errors.
  spec.ms <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                 distribution.spec = list(distribution = c("norm")),
                                 switch.spec = list(do.mix = FALSE, K = 2))
  par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566,
           0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)
  names(par) <- spec.ms$label

  vPw <- MSGARCH:::f_unmapPar(par, spec.ms, FALSE)
  inf <- MSGARCH:::f_InferenceFun(vPw, SMI, spec.ms, do.plm = FALSE)

  J  <- numDeriv::jacobian(MSGARCH:::f_mapPar, vPw, spec = spec.ms, do.plm = FALSE)
  V  <- MASS::ginv(inf$Hessian)
  se <- unname(inf$MatCoef[, "Std. Error"])

  se.delta      <- sqrt(diag(J %*% V %*% t(J)))
  se.transposed <- sqrt(diag(t(J) %*% V %*% J))
  iP            <- match(c("P_1_1", "P_2_1"), spec.ms$label)

  # the two orientations must disagree on the P block, otherwise this proves nothing
  testthat::expect_true(max(abs(se.delta[iP] - se.transposed[iP])) > 1e-4)

  testthat::expect_true(max(abs(se - se.delta)) < 1e-8)
  testthat::expect_true(all(is.finite(se)) && all(se > 0))

})

testthat::test_that("Pr(>|t|) is a two-sided p-value", {

  mCoef <- fit.sr$Inference$MatCoef
  testthat::expect_true(max(abs(mCoef[, "Pr(>|t|)"] -
                                2 * (1 - stats::pnorm(abs(mCoef[, "t value"]))))) < 1e-12)
  testthat::expect_true(all(mCoef[, "Pr(>|t|)"] >= 0 & mCoef[, "Pr(>|t|)"] <= 1))

})

testthat::test_that("Degrees of freedom drop K - 1 values per regime-constant parameter", {

  # a regime-constant parameter leaves one free value where there were K, so it
  # removes K - 1 degrees of freedom, not one. K = 2 cannot tell the two apart.
  for (K in 2:4) {
    spec.rc <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                   distribution.spec = list(distribution = c("std")),
                                   switch.spec = list(do.mix = FALSE, K = K),
                                   constraint.spec = list(regime.const = c("nu")))
    exp.dof <- length(spec.rc$label) - (K - 1L)
    testthat::expect_equal(MSGARCH:::dofMSGARCH(list(spec = spec.rc)), exp.dof)
  }

  # unconstrained and fixed-parameter specifications must be unaffected
  spec.free <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                   distribution.spec = list(distribution = c("norm")),
                                   switch.spec = list(do.mix = FALSE, K = 2))
  testthat::expect_equal(MSGARCH:::dofMSGARCH(list(spec = spec.free)),
                         length(spec.free$label))

  spec.fix <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                  distribution.spec = list(distribution = c("norm")),
                                  switch.spec = list(do.mix = FALSE, K = 2),
                                  constraint.spec = list(fixed = list(beta_1 = 0.8)))
  testthat::expect_equal(MSGARCH:::dofMSGARCH(list(spec = spec.fix)),
                         length(spec.fix$label) - 1L)

})

testthat::test_that("AIC and BIC use the free-parameter count of a constrained fit", {

  # end-to-end on a fitted K = 2 model: logLik()'s df attribute and the AIC/BIC
  # arithmetic must both follow dofMSGARCH
  spec.rc <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                 distribution.spec = list(distribution = c("std")),
                                 switch.spec = list(do.mix = FALSE, K = 2),
                                 constraint.spec = list(regime.const = c("nu")))
  fit.rc  <- MSGARCH::FitML(spec.rc, data = SMI, ctr = list(do.se = FALSE))
  exp.dof <- length(spec.rc$label) - 1L

  testthat::expect_equal(as.integer(attr(stats::logLik(fit.rc), "df")), exp.dof)
  testthat::expect_true(abs(AIC(fit.rc) - (-2 * fit.rc$loglik + 2 * exp.dof)) < 1e-8)
  testthat::expect_true(abs(BIC(fit.rc) -
                            (-2 * fit.rc$loglik + log(length(SMI)) * exp.dof)) < 1e-8)

  # K = 2 cannot separate "one" from "K - 1", so repeat for K = 3 and K = 4 on a
  # fit-shaped object: logLik.MSGARCH_ML_FIT only reads $loglik, $data and $spec
  for (K in 3:4) {
    spec.K <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                                  distribution.spec = list(distribution = c("std")),
                                  switch.spec = list(do.mix = FALSE, K = K),
                                  constraint.spec = list(regime.const = c("nu")))
    fit.K  <- structure(list(loglik = -3000, data = SMI, spec = spec.K),
                        class = "MSGARCH_ML_FIT")
    exp.K  <- length(spec.K$label) - (K - 1L)

    testthat::expect_equal(as.integer(attr(stats::logLik(fit.K), "df")), exp.K)
    testthat::expect_true(abs(AIC(fit.K) - (-2 * fit.K$loglik + 2 * exp.K)) < 1e-8)
    testthat::expect_true(abs(BIC(fit.K) -
                              (-2 * fit.K$loglik + log(length(SMI)) * exp.K)) < 1e-8)
  }

})
