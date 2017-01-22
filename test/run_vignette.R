# http://kbroman.org/pkg_primer/pages/github.html
# install.packages("devtools")
# require("devtools")
# devtools::install_github("keblu/MSGARCH", subdir = "Package")
# install.packages("MSGARCH")

rm(list = ls())
options(digits = 4, max.print = 80, prompt = "R> ", width = 80)
require("MSGARCH")
print(sessionInfo())

####################################################################
## Model specification
####################################################################

spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                            distribution = c("norm", "norm"),
                            do.skew = c(FALSE, FALSE), 
                            do.mix = FALSE, do.shape.ind = FALSE) 
spec

## Example 1: A single-regime model
spec = MSGARCH::create.spec(model = c("sGARCH"),
                            distribution = "norm",
                            do.skew = FALSE) 
spec

## Example 2: A model with heterogeneous regimes
spec = MSGARCH::create.spec(model = c("sGARCH", "tGARCH", "eGARCH"),
                            distribution = c("norm", "std", "ged"),
                            do.skew = c(TRUE, FALSE, TRUE),
                            do.mix = FALSE, do.shape.ind = FALSE)

## Summary of the specification object
spec = MSGARCH::create.spec()
summary(spec)

## Example: Estimation of a Markov{switching model
data("sp500")
spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                            distribution = c("norm", "norm"),
                            do.skew = c(FALSE, FALSE), 
                            do.mix = FALSE, do.shape.ind = FALSE) 
out.mle = MSGARCH::fit.mle(spec = spec, y = sp500)
summary(out.mle)

### Mixture of GARCH
spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                            distribution = c("norm", "norm"),
                            do.skew = c(FALSE, FALSE), 
                            do.mix = TRUE, do.shape.ind = FALSE)
out.mle = MSGARCH::fit.mle(spec = spec, y = sp500)
summary(out.mle)

## Regime-independent shape parameters
spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                            distribution = c("std", "std"),
                            do.skew = c(FALSE, FALSE), 
                            do.mix = FALSE, do.shape.ind = TRUE) 
out.mle = MSGARCH::fit.mle(spec = spec, y = sp500)
summary(out.mle)

####################################################################
## Estimation
####################################################################

## Maximum likelihood estimation
data("sp500")
spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                            distribution = c("std", "std"),
                            do.skew = c(FALSE, FALSE), 
                            do.mix = FALSE, do.shape.ind = FALSE) 
out.mle = MSGARCH::fit.mle(spec = spec, y = sp500)
summary(out.mle)

# with DEoptim
#set.seed(123)
#ctr.mle = list(do.init = TRUE, NP = 10 * length(spec$theta0), itermax = 500)
#out.mle = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = ctr.mle)
#summary(out.mle)

## Bayesian estimation
ctr.bay = list(N.burn = 5000, N.mcmc = 10000, N.thin = 10)
set.seed(123)
out.bay = MSGARCH::fit.bayes(spec = spec, y = sp500, ctr = ctr.bay)
summary(out.bay)
out.bay$theta

####################################################################
## Empirical illustration
####################################################################

rm(list = ls())
options(digits = 4, max.print = 80, prompt = "R> ", width = 80)
require("MSGARCH") 
require("coda")
require("DEoptim") # used for SMI data
data("SMI") # in DEoptim

postscript(file = "fig_smi.eps")
par(mfrow = c(1,1))
plot(y, xlab = "Date", ylab = "Log-return")
dev.off()

SMI  = as.matrix(y)
date = as.Date(rownames(SMI))
date = c(date, date[length(date)] + 1)

## single regime specification
spec.1 = MSGARCH::create.spec(model = "gjrGARCH", distribution = "std",
                              do.skew = TRUE, do.mix = FALSE, do.shape.ind = FALSE) 
set.seed(123)
out.mle.1 = MSGARCH::fit.mle(spec = spec.1, y = SMI, 
                             ctr = list(do.init = TRUE))
summary(out.mle.1)

postscript(file = "fig_vol.eps")
par(mfrow = c(1,1))
ht = sqrt(250) * MSGARCH::ht(out.mle.1)
plot(ht, date = date) # annual vol
dev.off()

## MS specification
spec.2 = MSGARCH::create.spec(model = c("gjrGARCH", "gjrGARCH"),
                              distribution = c("std", "std"),
                              do.skew = c(TRUE, TRUE), 
                              do.mix = FALSE, do.shape.ind = FALSE) 
## use DEoptim initialization to reproduce the results in Mullen et al.
set.seed(123)
out.mle.2 = MSGARCH::fit.mle(spec = spec.2, y = SMI, 
                             ctr = list(do.init = TRUE))
summary(out.mle.2)

state = MSGARCH::Pstate(out.mle.2)
# postscript(file = "fig_prob_mle.eps")
par(mfrow = c(1,1))
plot(state, date = date)
# dev.off()

# Bayesian estimation
ctr.bay.1 = list(N.burn = 5000, N.mcmc = 10000, N.thin = 10, theta0 = out.mle.1$theta)
set.seed(123)
out.bay.1 = MSGARCH::fit.bayes(spec = spec.1, y = SMI, ctr = ctr.bay.1)
summary(out.bay.1)

ctr.bay.2 = list(N.burn = 5000, N.mcmc = 10000, N.thin = 10, theta0 = out.mle.2$theta)
set.seed(123)
out.bay.2 = MSGARCH::fit.bayes(spec = spec.2, y = SMI, ctr = ctr.bay.2)
summary(out.bay.2)

postscript(file = "fig_mcmc1.eps")
par(mfrow = c(3,2))
coda::traceplot(out.bay.2$theta[,1:6])
dev.off()

postscript(file = "fig_mcmc2.eps")
pairs(x = as.matrix(out.bay.2$theta[,c(1,3,4,7,9,10)]), pch = 20, cex = 0.8)
dev.off()

state = MSGARCH::Pstate(out.bay.2)
# postscript(file = "fig_prob_bay.eps")
par(mfrow = c(1,1))
plot(state, date = date)
# dev.off()

out.bay.2.diff.leverage = out.bay.2$theta[,9] - out.bay.2$theta[,3]
postscript(file = "fig_leverage.eps")
par(mfrow = c(1,1))
hist(out.bay.2.diff.leverage, breaks = 100, main = "Difference in leverage effect", xlab = "Difference")
dev.off()

# Information criteria
c(MSGARCH::AIC(out.mle.1), MSGARCH::AIC(out.mle.2))
c(MSGARCH::BIC(out.mle.1), MSGARCH::BIC(out.mle.2))
c(MSGARCH::DIC(out.bay.1)$DIC, MSGARCH::DIC(out.bay.2)$DIC)

# Risk analysis
# !!! SUPER SLOW !!!
risk.mle.1 = MSGARCH::risk(out.mle.1, level = 0.95, ES = FALSE, do.its = TRUE)
risk.mle.2 = MSGARCH::risk(out.mle.2, level = 0.95, ES = FALSE, do.its = TRUE)
risk.bay.1 = MSGARCH::risk(out.bay.1, level = 0.95, ES = FALSE, do.its = TRUE)
risk.bay.2 = MSGARCH::risk(out.bay.2, level = 0.95, ES = FALSE, do.its = TRUE)
risk = cbind(risk.mle.1$VaR, risk.mle.2$VaR, risk.bay.1$VaR, risk.bay.2$VaR)
tsRainbow = rainbow(ncol(risk), alpha = 0.8)
colnames(risk) = c("GJR mle", "MSGARCH GJR mle", "GJR bay", "MSGARCH GJR bay")
postscript(file = "fig_var.eps")
par(mfrow = c(1,1))
plot(zoo::zoo(risk, order.by = date),plot.type = "single", col = tsRainbow, ylab = "VaR",xlab = "Date")
legend("bottomright",legend =  colnames(risk), lty = 1, col = tsRainbow)
dev.off()


































rm(list = ls())
options(digits = 4, max.print = 80, prompt = "R> ", width = 80)
require("MSGARCH")
print(sessionInfo())

####################################################################
spec = MSGARCH::create.spec(model = c("sGARCH"),
                            distribution = "norm",
                            do.skew = FALSE) 
spec

spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                            distribution = c("norm", "norm"),
                            do.skew = c(FALSE, FALSE), 
                            do.mix = FALSE, do.shape.ind = FALSE) 

spec

spec = MSGARCH::create.spec(model = c("sGARCH", "tGARCH", "eGARCH"),
                            distribution = c("norm", "std", "ged"),
                            do.skew = c(TRUE, FALSE, TRUE),
                            do.mix = FALSE, do.shape.ind = FALSE)

spec = MSGARCH::create.spec()
summary(spec)

####################################################################
## SPECIFICATION : Markov-switchingsingle-regime
data("sp500")
spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                            distribution = c("norm", "norm"),
                            do.skew = c(FALSE, FALSE), 
                            do.mix = FALSE, do.shape.ind = FALSE) 
out.mle = MSGARCH::fit.mle(spec = spec, y = sp500)
summary(out.mle)

####################################################################
## SPECIFICATION : Mixture
spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                            distribution = c("norm", "norm"),
                            do.skew = c(FALSE, FALSE), 
                            do.mix = TRUE, do.shape.ind = FALSE)
out.mle = MSGARCH::fit.mle(spec = spec, y = sp500)
summary(out.mle)

####################################################################
## SPECIFICATION : Shape independent parameters
spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                   distribution = c("std", "std"),
                   do.skew = c(FALSE, FALSE), 
                   do.mix = FALSE, do.shape.ind = TRUE) 
out.mle = MSGARCH::fit.mle(spec = spec, y = sp500)
summary(out.mle)

####################################################################
## ESTIMATION : MLE 
data("sp500")
spec = MSGARCH::create.spec(model = c("sGARCH", "sGARCH"),
                   distribution = c("std", "std"),
                   do.skew = c(FALSE, FALSE), 
                   do.mix = FALSE, do.shape.ind = FALSE) 
out.mle = MSGARCH::fit.mle(spec = spec, y = sp500)
summary(out.mle)

# with DEoptim
set.seed(123)
ctr.mle = list(do.init = TRUE, NP = 10 * length(spec$theta0), itermax = 500)
out.mle = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = ctr.mle)
summary(out.mle)

####################################################################
## ESTIMATION : Bayesian
ctr.bay = list(N.burn = 5000, N.mcmc = 10000, N.thin = 10)
set.seed(123)
out.bay = MSGARCH::fit.bayes(spec = spec, y = sp500, ctr = ctr.bay)
summary(out.bay)

####################################################################
## EMPIRICAL APPLICATION
rm(list = ls())
options(digits = 4, max.print = 80, prompt = "R> ", width = 80)
require("MSGARCH") 
require("DEoptim") # used for SMI data
data("SMI") # in DEoptim

postscript(file = "fig_smi.eps")
par(mfrow = c(1,1))
plot(y, xlab = "Date", ylab = "Log-return")
dev.off()

SMI  = as.matrix(y)
date = as.Date(rownames(SMI))
date = c(date, date[length(date)] + 1)

## single regime specification
spec.1 = MSGARCH::create.spec(model = "gjrGARCH", distribution = "std",
                     do.skew = TRUE, do.mix = FALSE, do.shape.ind = FALSE) 
out.mle.1 = MSGARCH::fit.mle(spec = spec.1, y = SMI)
summary(out.mle.1)

postscript(file = "fig_vol.eps")
par(mfrow = c(1,1))
ht = sqrt(250) * MSGARCH::ht(out.mle.1)
plot(ht, date = date) # annual vol
dev.off()

## MS specification
spec.2 = MSGARCH::create.spec(model = c("gjrGARCH", "gjrGARCH"),
                     distribution = c("std", "std"),
                     do.skew = c(TRUE, TRUE), 
                     do.mix = FALSE, do.shape.ind = FALSE) 
## use DEoptim initialization to reproduce the results in Mullen et al.
set.seed(123)
out.mle.2 = MSGARCH::fit.mle(spec = spec.2, y = SMI, ctr = list(do.init = TRUE, do.enhanced = TRUE))
summary(out.mle.2)

state = MSGARCH::Pstate(out.mle.2)
postscript(file = "fig_prob_mle.eps")
par(mfrow = c(1,1))
plot(state, date = date)
dev.off()

# Bayesian estimation
ctr.bay.1 = list(N.burn = 5000, N.mcmc = 10000, N.thin = 10, theta0 = out.mle.1$theta)
set.seed(123)
out.bay.1 = MSGARCH::fit.bayes(spec = spec.1, y = SMI, ctr = ctr.bay.1)
summary(out.bay.1)

ctr.bay.2 = list(N.burn = 5000, N.mcmc = 10000, N.thin = 10, theta0 = out.mle.2$theta)
set.seed(123)
out.bay.2 = MSGARCH::fit.bayes(spec = spec.2, y = SMI, ctr = ctr.bay.2)
summary(out.bay.2)

postscript(file = "fig_mcmc1.eps")
par(mfrow = c(3,2))
coda::traceplot(out.bay.2$theta[,1:6])
dev.off()

postscript(file = "fig_mcmc2.eps")
pairs(x = as.matrix(out.bay.2$theta[,c(1,3,4,7,9,10)]), pch = 20, cex = 0.8)
dev.off()

state = MSGARCH::Pstate(out.bay.2)
postscript(file = "fig_prob_bay.eps")
par(mfrow = c(1,1))
plot(state, date = date)
dev.off()

out.bay.2.diff.leverage = out.bay.2$theta[,9] - out.bay.2$theta[,3]
postscript(file = "fig_leverage.eps")
par(mfrow = c(1,1))
hist(out.bay.2.diff.leverage, breaks = 100, main = "Difference in leverage effect", xlab = "Difference")
dev.off()

# Information criterion
c(MSGARCH::AIC(out.mle.1), MSGARCH::AIC(out.mle.2))
c(MSGARCH::BIC(out.mle.1), MSGARCH::BIC(out.mle.2))
c(MSGARCH::DIC(out.bay.1)$DIC, MSGARCH::DIC(out.bay.2)$DIC)

# Risk analysis
risk.mle.1 = MSGARCH::risk(out.mle.1, level = 0.95, ES = FALSE, do.its = TRUE)
risk.mle.2 = MSGARCH::risk(out.mle.2, level = 0.95, ES = FALSE, do.its = TRUE)
risk.bay.1 = MSGARCH::risk(out.bay.1, level = 0.95, ES = FALSE, do.its = TRUE)
risk.bay.2 = MSGARCH::risk(out.bay.2, level = 0.95, ES = FALSE, do.its = TRUE)
risk = cbind(risk.mle.1$VaR, risk.mle.2$VaR, risk.bay.1$VaR, risk.bay.2$VaR)
tsRainbow = rainbow(ncol(risk), alpha = 0.8)
colnames(risk) = c("GJR mle", "MSGARCH GJR mle", "GJR bay", "MSGARCH GJR bay")
postscript(file = "fig_var.eps")
par(mfrow = c(1,1))
plot(zoo::zoo(risk, order.by = date),plot.type = "single", col = tsRainbow, ylab = "VaR",xlab = "Date")
legend("bottomright",legend =  colnames(risk), lty = 1, col = tsRainbow)
dev.off()




