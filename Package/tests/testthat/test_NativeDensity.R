testthat::context("Test the native mixture density and CDF entry points")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                            distribution.spec = list(distribution = c("norm")),
                            switch.spec = list(do.mix = FALSE, K = 2))
par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566,
         0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)
y <- SMI[1:300]

testthat::test_that("In-sample CDF fills the first time slice over the whole grid", {

  # At t = 1 the conditional variance of each regime is its unconditional variance
  # and the predictive state distribution is the ergodic one, so for a Normal
  # specification the in-sample PIT has a closed form to check against.
  x   <- c(-2, -1, 0, 1, 2)
  h1  <- as.numeric(spec$rcpp.func$unc_vol_Rcpp(matrix(par, nrow = 1L)))
  P0  <- State(object = spec, par = par, data = y)$PredProb[1L, 1L, ]
  exp.pit <- vapply(x, function(z) sum(P0 * stats::pnorm(z / sqrt(h1))),
                    FUN.VALUE = numeric(1))

  est.pit <- as.numeric(PIT(object = spec, x = x, par = par, data = y,
                            do.its = TRUE)[1L, ])

  testthat::expect_true(max(abs(est.pit - exp.pit)) < 1e-12)

  # the whole first row must be a genuine CDF, not a single value padded with zeros
  testthat::expect_true(all(est.pit > 0 & est.pit < 1))
  testthat::expect_true(!is.unsorted(est.pit, strictly = TRUE))

})

testthat::test_that("In-sample CDF accepts a grid longer than the sample", {

  # the cube is (length(data), length(x), K); writing it transposed ran off the end
  testthat::expect_silent(
    PIT(object = spec, x = seq(from = -5, to = 5, length.out = 100L), par = par,
        data = y[1:2], do.its = TRUE)
  )

})
