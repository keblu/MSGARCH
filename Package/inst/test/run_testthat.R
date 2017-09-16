## TestThat comparison
## Obtained with R 3.4.1 & MSGARCH 0.19.7
## 20170729

rm(list = ls())
library("MSGARCH")
library("testthat")
options(digits = 10, max.print = 40, prompt = "R> ", warn = 1)
tol <- 1e-5 # tolerance tests

tmp <- sessionInfo()
Plaform_version <- tmp$platform
MSGARCH_Version <- tmp$otherPkgs$MSGARCH$Version

#TestedPlaform_version = c("x86_64-w64-mingw32/x64 (64-bit)",
#                          "x86_64-apple-darwin15.6.0 (64-bit)",
#                          "x86_64-pc-linux-gnu (64-bit)")

## x86_64-w64-mingw32/x64 (64-bit) results
exp.par <- c(0.002719246, 0.003884714, 0.987470360, 0.202271584,
             0.065907799, 0.844972475, 0.965680895, 0.0355244463)

#if (!isTRUE(Plaform_version %in% TestedPlaform_version)) {
#  mess = paste0("Current platform never tested\n", "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
#  warning(mess)
#}

#if (Plaform_version == "x86_64-apple-darwin15.6.0 (64-bit)") {
#}

#if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
#}

## asset SMI
data("SMI", package = "MSGARCH")
spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")),
                   switch.spec = list(do.mix = FALSE, K = 2))
fit <- FitML(spec, data = SMI)

## testthat
cat("running tests\n")

est.par <- fit$par
testthat::expect_true(max(est.par - exp.par) < tol)
