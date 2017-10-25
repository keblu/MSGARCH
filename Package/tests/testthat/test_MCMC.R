testthat::context("Test Estimation")

tol <- 1e-4
tmp <- sessionInfo()
Plaform_version <- tmp$platform
#MSGARCH_Version <- tmp$otherPkgs$MSGARCH$Version
set.seed(1234)
data("SMI", package = "MSGARCH")
spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")),
                   switch.spec = list(do.mix = FALSE, K = 2))
fit <- FitMCMC(spec, data = SMI, ctr = list(n.burn = 500, n.mcmc = 500, n.thin = 1))

testthat::test_that("MCMC Estimation GARCH NORMAL", {
  est.par <- colMeans(fit$par)
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.024830248370305032, 0.058325465217369206, 0.900978983936512878, 0.465471511033367147, 0.029848278434822500,
                 0.895415983876721944, 0.978624434569418389, 0.242716760045187135)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.024830248370305032, 0.058325465217369206, 0.900978983936512878, 0.465471511033367147, 0.029848278434822500,
                 0.895415983876721944, 0.978624434569418389, 0.242716760045187135)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.024830248370305032, 0.058325465217369206, 0.900978983936512878, 0.465471511033367147, 0.029848278434822500,
                 0.895415983876721944, 0.978624434569418389, 0.242716760045187135)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})

testthat::test_that("DIC", {
  est.par <- DIC(fit)$DIC
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(6825.8205675690924)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(6825.8205675690924)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(6825.8205675690924)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})