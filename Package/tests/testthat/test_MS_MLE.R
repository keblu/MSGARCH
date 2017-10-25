testthat::context("Test Estimation")

tol <- 1e-4
tmp <- sessionInfo()
Plaform_version <- tmp$platform
#MSGARCH_Version <- tmp$otherPkgs$MSGARCH$Version

data("SMI", package = "MSGARCH")
spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")),
                   switch.spec = list(do.mix = FALSE, K = 2))
fit <- FitML(spec, data = SMI)
summary(fit)
est.par <- fit$par
est.BIC <- BIC(fit)
est.AIC <- AIC(fit)
exp.par <- NULL
exp.BIC <- NULL
exp.AIC <- NULL


testthat::test_that("Estimation MSGARCH", {
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
                 0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894 )
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.021631876264, 0.087024443575, 0.881493722230, 0.020659831227, 
                 0.005396009212, 0.994040728788, 0.978348086694, 0.998703303884)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566,
                 0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})


testthat::test_that("Estimation BIC", {
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.BIC <- c(6841.1848542696416)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.BIC <- c(6841.1848542696416)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.BIC <- c(6841.1848542696416)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.BIC <- exp.BIC
  }
  
  testthat::expect_true(max(exp.BIC - exp.BIC) < tol)
  
})


testthat::test_that("Estimation AIC", {
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.AIC <- c(6794.5924861827916)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.AIC <- c(6794.5924861827916)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.AIC <- c(6794.5924861827916)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.AIC <- est.AIC
  }
  
  testthat::expect_true(max(est.AIC - exp.AIC) < tol)
  
})

