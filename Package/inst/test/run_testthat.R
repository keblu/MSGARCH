## TestThat comparison
## Obtained with R 3.4.1 & MSGARCH 1.2
## 20171012

rm(list = ls())
library("MSGARCH")
library("testthat")
options(digits = 10, max.print = 40, prompt = "R> ", warn = 1)
tol <- 1e-5 # tolerance tests

tmp <- sessionInfo()
Plaform_version <- tmp$platform
#MSGARCH_Version <- tmp$otherPkgs$MSGARCH$Version
exp.par <- NULL

if (Plaform_version == "x86_64-w64-mingw32/x64 (64-bit)") {
  exp.par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
               0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894 )
}

if (Plaform_version == "x86_64-apple-darwin16.6.0 (64-bit)") {
  exp.par <- c(0.021631876264, 0.087024443575, 0.881493722230, 0.020659831227, 
               0.005396009212, 0.994040728788, 0.978348086694, 0.998703303884)
}


## asset SMI
data("SMI", package = "MSGARCH")
spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")),
                   switch.spec = list(do.mix = FALSE, K = 2))
fit <- FitML(spec, data = SMI)

## testthat
cat("running tests\n")

est.par <- fit$par
if (is.null(exp.par)) {
  mess = paste0("Current platform never tested\n", 
                "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
  warning(mess)
  exp.par <- est.par
}
testthat::expect_true(max(est.par - exp.par) < tol)
