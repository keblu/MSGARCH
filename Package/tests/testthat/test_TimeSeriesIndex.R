testthat::context("Test ts and zoo handling of fitted-object methods")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                            distribution.spec = list(distribution = c("norm")),
                            switch.spec = list(do.mix = FALSE, K = 2))

y.num <- as.numeric(SMI[1:200])
y.ts  <- stats::ts(y.num, start = c(2000, 1), frequency = 12)
newdata <- c(0.1, 0.2)

set.seed(1234)
fit.ts  <- MSGARCH::FitML(spec, data = y.ts,  ctr = list(do.se = FALSE))
set.seed(1234)
fit.num <- MSGARCH::FitML(spec, data = y.num, ctr = list(do.se = FALSE))

testthat::test_that("Combining a fit's data with newdata keeps index and values aligned", {

  # The index used to be derived from the already-concatenated series and then
  # extended by another length(newdata) points, so zooreg() recycled
  # observations and the model conditioned on values that are not in the sample.
  combined <- MSGARCH:::f_combine_data(y.ts, newdata)
  testthat::expect_equal(length(combined), length(y.num) + length(newdata))
  testthat::expect_equal(as.numeric(combined), c(y.num, newdata))
  testthat::expect_equal(stats::frequency(combined), stats::frequency(y.ts))
  testthat::expect_equal(stats::start(combined), stats::start(y.ts))

  # without newdata the series must come back untouched, index included
  testthat::expect_equal(stats::tsp(MSGARCH:::f_combine_data(y.ts, NULL)), stats::tsp(y.ts))

  # and a plain numeric series is just concatenated
  testthat::expect_equal(MSGARCH:::f_combine_data(y.num, newdata), c(y.num, newdata))

})

testthat::test_that("ts input gives the same numbers as the equivalent numeric input", {

  testthat::expect_equal(as.numeric(Volatility(fit.ts, newdata = newdata)),
                         as.numeric(Volatility(fit.num, newdata = newdata)))
  testthat::expect_equal(length(Volatility(fit.ts, newdata = newdata)),
                         length(y.num) + length(newdata))

  testthat::expect_equal(as.numeric(predict(fit.ts,  nahead = 1L)$vol),
                         as.numeric(predict(fit.num, nahead = 1L)$vol))
  testthat::expect_equal(as.numeric(PIT(fit.ts,  do.its = TRUE)),
                         as.numeric(PIT(fit.num, do.its = TRUE)))
  testthat::expect_equal(as.numeric(Risk(fit.ts,  alpha = 0.05)$VaR),
                         as.numeric(Risk(fit.num, alpha = 0.05)$VaR))

})

testthat::test_that("Forecast indexes advance by the series' own time step", {

  pred <- predict(fit.ts, nahead = 3L)$vol
  testthat::expect_equal(stats::frequency(pred), stats::frequency(y.ts))
  # three monthly steps beyond the end of the sample, not three whole years
  testthat::expect_equal(as.numeric(stats::time(pred)),
                         stats::tsp(y.ts)[2L] + (1:3) / stats::frequency(y.ts))

  vol <- Volatility(fit.ts)
  testthat::expect_equal(stats::tsp(vol)[1:2], stats::tsp(y.ts)[1:2])

})

testthat::test_that("zoo input is handled the same way", {

  y.zoo <- zoo::zoo(y.num, order.by = seq_len(length(y.num)))
  set.seed(1234)
  fit.zoo <- MSGARCH::FitML(spec, data = y.zoo, ctr = list(do.se = FALSE))

  testthat::expect_equal(as.numeric(Volatility(fit.zoo, newdata = newdata)),
                         as.numeric(Volatility(fit.num, newdata = newdata)))
  testthat::expect_equal(as.numeric(zoo::index(predict(fit.zoo, nahead = 2L)$vol)),
                         length(y.num) + 1:2)

})
