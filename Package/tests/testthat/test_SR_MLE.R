testthat::context("Test Estimation")

tol <- 1e-4
tmp <- sessionInfo()
Plaform_version <- tmp$platform
#MSGARCH_Version <- tmp$otherPkgs$MSGARCH$Version


testthat::test_that("Estimation ARCH NORMAL", {
  
  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("sARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- FitML(spec, data = SMI)
  summary(fit)
  est.par <- fit$par
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.80504578435066854, 0.31371193670568287)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.80504578435066854, 0.31371193670568287)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.80504578435066854, 0.31371193670568287)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})

testthat::test_that("Estimation GARCH NORMAL", {
  
  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- FitML(spec, data = SMI)
  
  est.par <- fit$par
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.085912610592369154, 0.124325097982214264, 0.796421589649365935)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.085912610592369154, 0.124325097982214264, 0.796421589649365935)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.085912610592369154, 0.124325097982214264, 0.796421589649365935)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})

testthat::test_that("Estimation GJR NORMAL", {
  
  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- FitML(spec, data = SMI)
  
  est.par <- fit$par
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.096745854166231729, 0.042005697086443979, 0.174648268027545489, 0.786944835620325578)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.096745854166231729, 0.042005697086443979, 0.174648268027545489, 0.786944835620325578)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.096745854166231729, 0.042005697086443979, 0.174648268027545489, 0.786944835620325578)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})

testthat::test_that("Estimation EGARCH NORMAL", {
  
  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("eGARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- FitML(spec, data = SMI)
  
  est.par <- fit$par
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.0089803682227889703,  0.1892548543632634195, -0.0983816411668451546,  0.9242572117676925991)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.0089803682227889703,  0.1892548543632634195, -0.0983816411668451546,  0.9242572117676925991)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.0089803682227889703,  0.1892548543632634195, -0.0983816411668451546,  0.9242572117676925991)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})

testthat::test_that("Estimation TGARCH NORMAL", {
  
  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("tGARCH")),
                     distribution.spec = list(distribution = c("norm")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- FitML(spec, data = SMI)
  
  est.par <- fit$par
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.086895299072832388, 0.057689890206626467, 0.182831991103102276, 0.824269321819041778)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.086895299072832388, 0.057689890206626467, 0.182831991103102276, 0.824269321819041778)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.086895299072832388, 0.057689890206626467, 0.182831991103102276, 0.824269321819041778)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})

testthat::test_that("Estimation GJR STUDENT", {
  
  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                     distribution.spec = list(distribution = c("std")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- FitML(spec, data = SMI)
  
  est.par <- fit$par
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.042068160604000746, 0.041341103857468385, 0.123041130026077569, 0.862000796059737562, 8.408571456679910128)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.042068160604000746, 0.041341103857468385, 0.123041130026077569, 0.862000796059737562, 8.408571456679910128)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.042068160604000746, 0.041341103857468385, 0.123041130026077569, 0.862000796059737562, 8.408571456679910128)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})

testthat::test_that("Estimation GJR GED", {
  
  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                     distribution.spec = list(distribution = c("ged")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- FitML(spec, data = SMI)
  
  est.par <- fit$par
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.065574607143466801, 0.043629701610720460, 0.148869843737306995, 0.826662611816671999, 1.438472937139093855)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.065574607143466801, 0.043629701610720460, 0.148869843737306995, 0.826662611816671999, 1.438472937139093855)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.065574607143466801, 0.043629701610720460, 0.148869843737306995, 0.826662611816671999, 1.438472937139093855)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})

testthat::test_that("Estimation GJR STUDENT SKEWED", {
  
  data("SMI", package = "MSGARCH")
  spec <- CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                     distribution.spec = list(distribution = c("sstd")),
                     switch.spec = list(do.mix = FALSE, K = 1))
  fit <- FitML(spec, data = SMI)
  
  est.par <- fit$par
  exp.par = NULL
  if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
    exp.par <- c(0.039325587246528841, 0.043003103795721226, 0.114555021691679665, 0.870109422052372516, 8.113661279136989535, 0.855334642627564756)
  }
  
  if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
    exp.par <- c(0.039325587246528841, 0.043003103795721226, 0.114555021691679665, 0.870109422052372516, 8.113661279136989535, 0.855334642627564756)
  }
  
  if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
    exp.par <- c(0.039325587246528841, 0.043003103795721226, 0.114555021691679665, 0.870109422052372516, 8.113661279136989535, 0.855334642627564756)
  }
  
  if (is.null(exp.par)) {
    mess = paste0("Current platform never tested\n", 
                  "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
    warning(mess)
    exp.par <- est.par
  }
  
  testthat::expect_true(max(est.par - exp.par) < tol)
  
})
