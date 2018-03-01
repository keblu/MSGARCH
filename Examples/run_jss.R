#################################################################################
### DESCRIPTION

### This code is used in the illustrations of
### Ardia, Bluteau, Boudt, Catania & Trottier (2017)
### 'Markov-Switching GARCH Models in R: The MSGARCH Package'.

### !!! Results of the paper were obtained with the following setup:
### !!! R version 3.4.3 (2017-11-30)
### !!! Platform: x86_64-w64-mingw32/x64 (64-bit)
### !!! Running under: Windows 7 x64 (build 7601) Service Pack 1

### !!! RESULTS ARE PLATFORM DEPENDENT (but similar up to the 8th digits) !!!
### !!! ALSO SET THE SEED PROPERLY !!!

#################################################################################
### LOAD THE PACKAGE

# Load the package from CRAN or from the tar
# install.packages("MSGARCH") # !!! install version 2.1
# install.packages("MSGARCH_2.1.tar.gz", repos = NULL)

#################################################################################
### PACKAGE MSGARCH
### Reports the code used to generate full results in the paper
#################################################################################

rm(list = ls())
options(prompt = "R> ", continue = "+  ", width = 70,
        digits = 4, max.print = 80, useFancyQuotes = FALSE)
tmp <- sessionInfo()
nam <- paste0("PART_I_R", tmp$R.version$major, ".",
              tmp$R.version$minor, "_", tmp$R.version$os)
sink(file = paste0("sink_", nam, ".txt"), append = FALSE, split = TRUE) # output printed in txt
print(tmp)

######################################################
cat("\n\n")
cat("SECTION 3.1\n")
cat("-----------\n\n")

library("MSGARCH")
spec <- CreateSpec()
summary(spec)

## Example 1: A single-regime model
spec <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                   distribution.spec = list(distribution = c("norm")))

## Example 2: A model with heterogeneous regimes
spec <- CreateSpec(variance.spec = list(model = c("sGARCH", "tGARCH", "eGARCH")),
                   distribution.spec = list(distribution = c("snorm", "sstd", "sged")))

## Example 3: A model with non-switching shape parameters
spec <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")),
                   distribution.spec = list(distribution = c("sstd", "sstd")),
                   constraint.spec = list(regime.const = c("nu", "xi")))

######################################################
cat("\n\n")
cat("SECTION 3.2\n")
cat("-----------\n\n")

data("dem2gbp", package = "MSGARCH")
ms2.garch.n <- CreateSpec(variance.spec = list(model = c("sGARCH")),
                          distribution.spec = list(distribution = c("norm")),
                          switch.spec = list(K = 2))

fit.ml <- FitML(spec = ms2.garch.n, data = dem2gbp)
summary(fit.ml)

set.seed(1234)
fit.mcmc <- FitMCMC(spec = ms2.garch.n, data = dem2gbp)
summary(fit.mcmc)

######################################################
cat("\n\n")
cat("SECTION 3.3\n")
cat("-----------\n\n")

pred <- predict(fit.ml, nahead = 5, do.return.draw = TRUE)
pred$vol

pred$draw[, 1:4]

######################################################
cat("\n\n")
cat("SECTION 3.4\n")
cat("-----------\n\n")

risk <- Risk(fit.ml, alpha = c(0.01, 0.05), nahead = 5)
risk$VaR
risk$ES

######################################################
cat("\n\n")
cat("SECTION 3.5\n")
cat("-----------\n\n")

simulate(fit.ml, nsim = 2, nahead = 4, nburin = 500)

BIC(fit.ml)

sr.fit <- ExtractStateFit(fit.ml)
risk1 <- Risk(sr.fit[[1]], alpha = 0.05, nahead = 5)
risk2 <- Risk(sr.fit[[2]], alpha = 0.05, nahead = 5)
VaR <- cbind(risk1$VaR, risk2$VaR)
colnames(VaR) <- c("State 1", "State 2")
VaR

sink()

#################################################################################
### EMPIRICAL ILLUSTRATIONS
### Reports the code used to generate full results in the paper
#################################################################################

rm(list = ls())
library("MSGARCH")
options(prompt = "R> ", continue = "+  ", width = 70,
        digits = 4, max.print = 80, useFancyQuotes = FALSE)
tmp <- sessionInfo()
nam <- paste0("PART_II_R", tmp$R.version$major, ".",
              tmp$R.version$minor, "_", tmp$R.version$os)
sink(file = paste0("sink_", nam, ".txt"), append = FALSE, split = TRUE) # output printed in txt
print(tmp)

## load SMI
data("SMI", package = "MSGARCH")

## Create MS(2)-GJR-std specification (Ardia 2008 and Ardia et al. 2008)
ms2.gjr.s <- CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                        distribution.spec = list(distribution = c("std")),
                        switch.spec = list(K = 2),
                        constraint.spec = list(regime.const = "nu"))

## ML estimation
fit.ml <- FitML(ms2.gjr.s, data = SMI)

## BIC
BIC(fit.ml)

## Summary
summary(fit.ml)

## Unconditional vol
set.seed(1234) 
sqrt(250) * sapply(ExtractStateFit(fit.ml), UncVol)

## Smoothed probabilities in regime 2 and volatility
smoothed.prob <- State(fit.ml)$SmoothProb[, 1, 2, drop = TRUE]
vol <- sqrt(250) * Volatility(fit.ml)

## MCMC estimation
nmcmc <- 12500
nburn <- 5000
nthin <- 5
ctr <- list(nmcmc = nmcmc, nburn = nburn,
            nthin = nthin, par0 = fit.ml$par)
fit.mcmc <- FitMCMC(ms2.gjr.s, data = SMI, ctr = ctr)
summary(fit.mcmc)

## Convergence of the chain
#par(mfrow = c(3, 4))
#coda::traceplot(fit.MCMC$par)
#coda::heidel.diag(fit.MCMC$par)
#coda::acfplot(fit.MCMC$par)

## Posterior draws
draws <- as.matrix(fit.mcmc$par)

## This function computes the unconditional volatility
## for a GJR model with symmeetric disturbances
f_ucvol <- function(par) {
  if (is.vector(par)) {
    par <- matrix(data = par, nrow = 1, dimnames = list(1, names(par)))
  }
  ucvol_1 <- sqrt(250) * par[,"alpha0_1"] / (1 - (par[,"alpha1_1"] + 0.5 * par[,"alpha2_1"] + par[,"beta_1"]))
  ucvol_2 <- sqrt(250) * par[,"alpha0_2"] / (1 - (par[,"alpha1_2"] + 0.5 * par[,"alpha2_2"] + par[,"beta_2"]))
  out <- list(ucvol_1 = ucvol_1, ucvol_2 = ucvol_2)
  return(out)
}

## Compute unconditional volatility
ucvol.draws <- f_ucvol(draws)
ucvol.bay   <- lapply(ucvol.draws, mean)
ucvol.mle   <- f_ucvol(fit.ml$par)

## Posterior mean
unlist(ucvol.bay)

## Quantiles of unconditional volatility
sapply(ucvol.draws, quantile, probs = c(0.025, 0.975))

## Impact of paramter uncertainty in pred
n.mesh <- 1000
x <- seq(from = -5, to = 0, length.out = n.mesh)
pred.mle <- as.vector(PredPdf(fit.ml, x = x, nahead = 1))
pred.bay <- as.vector(PredPdf(fit.mcmc, x = x, nahead = 1))

pred.draws <- matrix(data = NA, nrow = nrow(draws), ncol = n.mesh)
for (i in 1:nrow(draws)) {
  tmp <- PredPdf(ms2.gjr.s, par = draws[i,], x = x, data = SMI, nahead = 1)
  pred.draws[i,] <- as.vector(tmp)
}

## Backtesting

## Create GJR-std specification for comparison
gjr.s <- CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                    distribution.spec = list(distribution = c("std")),
                    switch.spec = list(K = 1))

models <- list(gjr.s, ms2.gjr.s)

n.ots    <- 1000 # number of out-of-sample evaluation 
n.its    <- 1500 # fit sample size
alpha    <- 0.05 # risk Level 
k.update <- 100  # estimation frequency

## Initialization 
VaR   <- matrix(NA, nrow = n.ots, ncol = length(models))
y.ots <- matrix(NA, nrow = n.ots, ncol = 1)
model.fit <- vector(mode = "list", length = length(models))

# iterate over out-of-sample time
for (i in 1:n.ots) { 
  cat("Backtest - Iteration: ", i, "\n")
  y.its    <- SMI[i:(n.its + i - 1)] # in-sample data 
  y.ots[i] <- SMI[n.its + i]         # out-of-sample data
  
  # iterate over models
  for (j in 1:length(models)) { 
    
    # do we update the model estimation
    if (k.update == 1 || i %% k.update == 1) {
      cat("Model", j, "is reestimated\n")
      model.fit[[j]] <- FitML(spec = models[[j]], data = y.its, 
                              ctr = list(do.se = FALSE)) 
    }
    
    # calculate VaR 1-step ahead
    VaR[i,j] <- Risk(model.fit[[j]]$spec, par = model.fit[[j]]$par,
                     data = y.its,
                     n.ahead = 1,
                     alpha   = alpha,
                     do.es   = FALSE,
                     do.its  = FALSE)$VaR
  }                                
}

## Test the VaR
# install.packages("GAS")
library("GAS")
CC.pval <- DQ.pval <- vector("double", length(models))
for (j in 1:length(models)) { 
  test <- GAS::BacktestVaR(data  = y.ots,
                           VaR   = VaR[,j],
                           alpha = alpha)
  
  CC.pval[j] <- test$LRcc[2]      
  DQ.pval[j] <- test$DQ$pvalue  
}
names(CC.pval) <- names(DQ.pval) <- c("GJR-std", "MS2-GJR-std")

print(CC.pval)
print(DQ.pval)

sink()

############################################
################# FIGURES ##################
############################################

######################
##     FIGURE 1     ##
######################

pdf(file = "figure1.pdf", height = 13, width = 16, compress = TRUE)
par(mfrow = c(1, 1), mar = c(5,3,2,2) + 0.1)
plot(SMI, type = 'l', las = 1, lwd = 2, xlab = "Date (year)",
     ylab = "", col = "black", cex.axis = 1.5, cex.lab = 1.5)
title("SMI log-returns (%)", cex.main = 1.5)
dev.off()

######################
##     FIGURE 2     ##
######################

pdf(file = "figure2.pdf", height = 13, width = 16, compress = TRUE)
op <- par(mfrow = c(2,1),
          oma = c(1,1,0,0) + 0.0,
          mar = c(2,2,2,2) + 0.0)
plot(as.vector(SMI), las = 1, type = 'p', pch = 20, col = 'black',
     cex = 1.5, axes = FALSE, ann = FALSE)
par(new = TRUE)
ylabel <- expression(paste("Pr(", s[t], " = 2 | ", hat(psi), ", ", I[t], ")"))
plot(zoo::zoo(smoothed.prob, order.by = zoo::index(SMI)), lty = 1, plot.type = "single",
     col = "red", las = 1, ylab = "", xlab = "Date", lwd = 3, cex.axis = 1.5, cex.lab = 1.5)
title(main = "Smoothed probabilities", cex.main = 1.5)
plot(zoo::zoo(vol, order.by = zoo::index(SMI)), lty = 1, plot.type = "single",
     col = "black", las = 1, ylab = "", xlab = "Date", lwd = 3, cex.axis = 1.5, cex.lab = 1.5)
title(main = "Volatility (%)", cex.main = 1.5)
par(op)
dev.off()

######################
##     FIGURE 3     ##
######################

sel <- c("alpha1_1", "alpha2_1")
tmp <- draws[, sel]
par.mle <- fit.ml$par[sel]
par.bay <- apply(tmp, 2, mean)
xlim <- range(c(tmp[,1], par.mle[1]))
ylim <- range(c(tmp[,2], par.mle[2]))
pdf(file = "figure3.pdf", height = 13, width = 13, compress = TRUE)
par(mfrow = c(1, 1))
plot(tmp, pch = 20, las = 1, lwd = 2, cex = 2, xlim = xlim, ylim = ylim,
     col = "lightsteelblue", cex.axis = 1.5, cex.lab = 1.5)
grid()
par(new = TRUE)
points(par.bay[1], par.bay[2], cex = 4, lwd = 2,
       pch = 15, col = "blue", xlim = xlim, ylim = ylim)
points(par.mle[1], par.mle[2], cex = 4, lwd = 2,
       pch = 17, col = "red", xlim = xlim, ylim = ylim)
dev.off()

######################
##     FIGURE 4     ##
######################

n <- length(ucvol.draws$ucvol_1)
pdf(file = "figure4.pdf", height = 13, width = 16, compress = TRUE)
op <- par(mar = c(2, 2, 2, 2),
          mfrow = c(1, 2),
          oma = c(2, 2, 0.2, 0.2))
hist(ucvol.draws$ucvol_1, nclass = round(10 * log(n)), prob = TRUE,
     col = "lightsteelblue", las = 1, xlab = "Volatility (%)",
     ylab = "", cex.lab = 1.5, cex.axis = 1.5, main = "")
title(main = "Regime 1", cex.main = 1.5)
lines(density(ucvol.draws$ucvol_1), col = "black", lwd = 2)
rug(ucvol.draws$ucvol_1); box()
points(ucvol.bay$ucvol_1, 0, pch = 15, col = "blue", lwd = 2, cex = 3)
points(ucvol.mle$ucvol_1, 0, pch = 17, col = "red", lwd = 2, cex = 3)
hist(ucvol.draws$ucvol_2, nclass = round(10 * log(n)), prob = TRUE,
     col = "lightsteelblue", las = 1, xlab = "Volatility (%)",
     ylab = "", cex.lab = 1.5, cex.axis = 1.5, main = "")
rug(ucvol.draws$ucvol_2); box()
points(ucvol.bay$ucvol_2, 0, pch = 15, col = "blue", lwd = 2, cex = 3)
points(ucvol.mle$ucvol_2, 0, pch = 17, col = "red", lwd = 2, cex = 3)
title(main = "Regime 2", cex.main = 1.5)
lines(density(ucvol.draws$ucvol_2), col = "black", lwd = 2)
par(op)
dev.off()

######################
##     FIGURE 5     ##
######################

pdf(file = "figure5.pdf", height = 13, width = 16, compress = TRUE)
#png(filename = "figure5.png", height = 13, width = 16, units = "in", res = 600)
xlim <- c(-4, -1.2)
ylim <- c(0, 0.1)
par(mfrow = c(1, 1))
matplot(x, t(pred.draws), xlim = xlim, ylim = ylim,
        type = "l", col = "lightsteelblue",
        xlab = "Return (%)", ylab = "Predictives",
        lty = 1.5, las = 1, cex.axis = 1.5, cex.lab = 1.5)
title(main = "Left-tail forecast of SMI index return", cex.main = 1.5)
lines(x, pred.bay, xlim = xlim, ylim = ylim,
      type = "l", lty = "solid", col = "blue", lwd = 3)
lines(x, pred.mle, xlim = xlim, ylim = ylim,
      type = "l", pch = "o", lty = "dashed", col = "red", lwd = 3)
legend("topleft", c("MCMC draws", "Bayesian","ML"),
       col = c("lightsteelblue", "blue", "red"), lwd = 3,
       lty = c(1, 1, 2), bty = "n", cex = 2)
box()
dev.off()

######################
##     FIGURE 6     ##
######################

library("zoo")
time.index <- zoo::index(SMI)[(n.its + 1):(n.ots + n.its)]
y_ots <- zoo::zoo(y.ots, order.by = time.index)
VaR   <- zoo::zoo(VaR, order.by = time.index)

pdf(file = "figure6.pdf", height = 13, width = 16, compress = TRUE)
par(mfrow = c(1, 1))
plot(y_ots, type = 'p', las = 1, lwd = 1, xlab = "Date (year)",
     ylab = "", col = "black", cex.axis = 1.5, cex.lab = 1.5, pch = 19)
lines(VaR[,1], type = 'l', col = "red", lwd = 3, lty = "dashed")
lines(VaR[,2], type = 'l', col = "blue", lwd = 3)
legend("topleft", legend = c("VaR 5% - GJR-std", "VaR 5% - MS2-GJR-std"), 
       col = c("red", "blue"), lwd = 3, cex = 1.5, lty = c("dashed", "solid"))
abline(h = 0)
title("Backtesing VaR at 5% risk level", cex.main = 1.5)
dev.off()

