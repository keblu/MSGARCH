#################################################################################
### DESCRIPTION

### This code is used in the illustrations of
### regime-switching breaks in the volatility of GARCH processes

# require("devtools")
# devtools::install_github("keblu/MSGARCH", subdir = "Package")

rm(list = ls())
library("MSGARCH")
spec.MS  = MSGARCH::CreateSpec()
spec.SR  = MSGARCH::CreateSpec(variance.spec = list(model = "sGARCH"))
par.SR.1 = c(0.02, 0.1, 0.8)
par.SR.2 = c(0.04, 0.1, 0.8)

set.seed(26)
SR.sim.1 = simulate(spec.SR, nsim = 1,nahead = 250, nburn = 250, par = par.SR.1)
SR.sim.2 = simulate(spec.SR, nsim = 1,nahead = 250, nburn = 1, par = par.SR.2)
y = c(SR.sim.1$draw, SR.sim.2$draw)

par(mfrow = c(1,1))
plot(1:250, SR.sim.1$draw, col = "blue", type = "l", xlim = c(0, 500), ylim = c(-1.5, 1.5))
lines(251:500, SR.sim.2$draw, col = "red", type = "l", xlim = c(0, 500), ylim = c(-1.5, 1.5))
plot(y, type = 'l', ylim = c(-1.5, 1.5))

## Single-regime fit
fit.SR = FitML(spec = spec.SR, data =  y) # covariance stationary
fit.SR$par[2] + fit.SR$par[3]

## not covariance stationary (using rugarch)
spec.SR.ru = rugarch::ugarchspec(variance.model = list(model = "sGARCH"),
                                 mean.model = list(armaOrder = c(0,0), include.mean = FALSE))

fit.SR.ru = rugarch::ugarchfit(spec.SR.ru, data = y,fit.control = list(stationarity = 0))
fit.SR.ru@fit$coef[2] + fit.SR.ru@fit$coef[3]

## MSGARCH fit
fit.MS = FitML(spec = spec.MS, data =  y)
summary(fit.MS)
state.MS = State(fit.MS)
plot(state.MS, type =  "viterbi")
plot(state.MS, type =  "smoothed")

# compare the unconditional variances in each state
fit.MS$par["alpha0_1"] / (1 - fit.MS$par["alpha1_1"]- fit.MS$par["beta_1"])
par.SR.1[1] / (1 - par.SR.1[2]-par.SR.1[3])

fit.MS$par["alpha0_2"] / (1 - fit.MS$par["alpha1_2"]- fit.MS$par["beta_2"])
par.SR.2[1] / (1 - par.SR.2[2]-par.SR.2[3])
