#################################################################################
### DESCRIPTION

### This code aims at illustrating how MSGARCH can be used to estimate
### MSGARCH-type models published in JSS
### - Ardia et al. (2009). "Adaptive Mixture of Student-t Distributions
###   as a Flexible Candidate Distribution for Effiient Simulation: The R Package AdMit"
###   Journal of Statistical Software, 29(3)
###   doi:10.18637/jss.v029.i03
### - Mullen et al. (2011). "DEoptim: An R Package for Global Optimization by Differential Evolution"
###   Journal of Statistical Software, 40(6), 1-26
###   doi:10.18637/jss.v040.i06

#################################################################################
### PART I
### Ardia et al. (2009)

## Posterior mode (optimization)
##    alpha0_1 0.0350
##    alpha0_2 0.2782
##    alpha1_1 0.2129
##    P_1      0.5826

##  Posterior mean (MCMC - AdMit)
##    alpha0_1 0.0457
##    alpha0_2 0.3519
##    alpha1_1 0.2316
##    P_1      0.6389

rm(list = ls())
options(prompt = "R> ", continue = "+  ", width = 70, digits = 4,
        max.print = 80, useFancyQuotes = FALSE)
library("MSGARCH")
library("AdMit")
library("MitISEM")
data("dem2gbp", package = "MSGARCH")
dem2gbp <- dem2gbp[1:250]

par(mfrow = c(1,1))
plot(dem2gbp, type = "l", las = 1, cex.axis = 1.1, cex.lab = 1.2,
     ylab = "log-returns (in-percent)", xlab = "Time index")

## specification
# spec <- CreateSpec(variance.spec = list(model = "sARCH"),
#                    distribution.spec = list(distribution = "norm"),
#                    switch.spec = list(do.mix = TRUE, K = 2),
#                    constraint.spec = list(regime.const = "alpha1"))

spec <- CreateSpec(variance.spec = list(model = "sARCH"),
                   distribution.spec = list(distribution = "norm"),
                   switch.spec = list(do.mix = TRUE, K = 2),
                   constraint.spec = list(regime.const = "alpha1"),
                   prior = list(mean = c(alpha0_1 = 0, alpha0_2 = 0, alpha1_1 = 0.2),
                                sd = c(alpha0_1 = 2, alpha0_2 = 2, alpha1_1 = 0.5)))
summary(spec)

## estimation with MSGARCH setup
fit.ML <- FitML(spec, data = dem2gbp)
summary(fit.ML)

set.seed(1234)
fit.MCMC <- FitMCMC(spec, data = dem2gbp,
                    ctr = list(n.burn = 1000L, n.mcmc = 50000L, n.thin = 1))
summary(fit.MCMC)
draws <- fit.MCMC$par

par(mfrow = c(2,2))
par(cex.axis = 1.1, cex.lab = 1.2)
acf(draws[,"alpha0_1"], lag.max = 30, las = 1, main = expression(alpha["0,1"]))
acf(draws[,"alpha0_2"], lag.max = 30, las = 1, main = expression(alpha["0,2"]))
acf(draws[,"alpha1_1"], lag.max = 30, las = 1,main = expression(alpha[1]))
acf(draws[,"P_1"], lag.max = 30, las = 1, main = expression(P[11]))

## draws from the marginal distribution (omega_2,p)'
draws <- as.matrix(draws)
par(mfrow = c(1,1))
plot(draws[, c("alpha0_2", "P_1")], pch = 19, cex = .7, las = 1,
     cex.axis = 1.1, cex.lab = 1.2,
     xlab = expression(alpha["0,2"]), ylab = expression(P[11]),
     xlim = c(0.0, 1.5), ylim = c(0.0, 1.0), axes = FALSE)
axis(side = 1, at = seq(from = 0.0, to = 1.6, by = 0.2))
axis(side = 2, at = seq(from = 0.0, to = 1.0, by = 0.2), las = 1)
box()
abline(v = c(0.8, 1.0, 1.2), lwd = 2, lty = "dotted")
abline(h = c(0.8, 0.9), lwd = 2, lty = "dotted")

  f_nelder_mead <- function(vPw, f_nll, spec, data, do.plm){
  out <- stats::optim(vPw, f_nll, spec = spec, data = data,
                      do.plm = do.plm, method = "Nelder-Mead",
                      control = list(trace = 1, maxit = 5000))
  return(out)
}

f_MH <- function(f_posterior, data, spec, par0, ctr){
  KERNEL <- function(vPw, data, log = TRUE) {
    if (is.matrix(vPw)) {
      out <- vector(mode = "numeric", length = nrow(vPw))
      for (i in 1:nrow(vPw)) {
        out[i] <- f_posterior(vPw = vPw[i,],
                              data = data,
                              spec = spec)
      }
    } else {
      out <- f_posterior(vPw = vPw,
                         data = data,
                         spec = spec)
    }
    return(out)
  }
  cat("Construct proposal\n")
  out.mit <- MitISEM::MitISEM(KERNEL = KERNEL, mu0 = par0, data = data)

  cat("MH sampler\n")
  out.MH <- AdMit::AdMitMH(N = ctr$n.mcmc + ctr$n.burn, KERNEL = KERNEL,
                           mit = out.mit$mit, data = data)

  colnames(out.MH$draws) <- names(par0)
  return(out.MH$draws)
}

par.start <- c(alpha0_1 = 0.1, alpha1_1 = 0.5, alpha0_2 = 0.1, P_1 = 0.5)
fit.ML <- FitML(spec, data = dem2gbp, ctr = list(OptimFUN = f_nelder_mead, par0 = par.start))
summary(fit.ML)

set.seed(1234)
fit.MCMC <- FitMCMC(spec, data = dem2gbp, ctr = list(SamplerFUN = f_MH,
                                                     n.burn = 1000L, n.mcmc = 50000L, n.thin = 1,
                                                     par0 = fit.ML$par))
summary(fit.MCMC)
draws <- fit.MCMC$par

par(mfrow = c(2,2))
par(cex.axis = 1.1, cex.lab = 1.2)
acf(draws[,"alpha0_1"], lag.max = 30, las = 1, main = expression(alpha["0,1"]))
acf(draws[,"alpha0_2"], lag.max = 30, las = 1, main = expression(alpha["0,2"]))
acf(draws[,"alpha1_1"], lag.max = 30, las = 1,main = expression(alpha[1]))
acf(draws[,"P_1"], lag.max = 30, las = 1, main = expression(P[11]))

draws <- as.matrix(draws)
par(mfrow = c(1,1))
plot(draws[, c("alpha0_2", "P_1")], pch = 19, cex = .7, las = 1,
     cex.axis = 1.1, cex.lab = 1.2,
     xlab = expression(alpha["0,2"]), ylab = expression(P[11]),
     xlim = c(0.0, 1.5), ylim = c(0.0, 1.0), axes = FALSE)
axis(side = 1, at = seq(from = 0.0, to = 1.6, by = 0.2))
axis(side = 2, at = seq(from = 0.0, to = 1.0, by = 0.2), las = 1)
box()
abline(v = c(0.8, 1.0, 1.2), lwd = 2, lty = "dotted")
abline(h = c(0.8, 0.9), lwd = 2, lty = "dotted")

#################################################################################
### PART II
### Mullen et al. (2011)

rm(list = ls())
options(prompt = "R> ", continue = "+  ", width = 70, digits = 4,
        max.print = 80, useFancyQuotes = FALSE)
library("MSGARCH")
data("SMI", package = "MSGARCH")

## Parameters at the global optimum (NLL = 3350.6979) obtained by
## solnp and DEoptim (after a longer run with itermax = 2500)
##   alpha0_1 0.2062
##   alpha1_1 0.0000
##   alpha2_1 0.2123
##   beta_1   0.5295
##   nu_1     9.2480
##   alpha0_2 0.0930
##   alpha1_2 0.0043
##   alpha2_2 0.1566
##   beta_2   0.8717
##   P_1_1    0.9981
##   P_2_1    0.0031

0.2062 / (1 - 0.5 * 0.0000 - 0.5 * 0.2123 - 0.5295)
0.0930 / (1 - 0.5 * 0.0043 - 0.5 * 0.1566 - 0.8717)

## specification
spec <- CreateSpec(variance.spec = list(model = "gjrGARCH"),
                   distribution.spec = list(distribution = "std"),
                   switch.spec = list(K = 2),
                   constraint.spec = list(regime.const = "nu"))
summary(spec)

## estimation with MSGARCH setup
fit <- FitML(spec, data = SMI)
summary(fit)

f_optim <- function(vPw, f_nll, spec, data, do.plm){
  tmp <- Rsolnp::solnp(vPw, f_nll, spec = spec, data = data,
                      do.plm = do.plm)
  
  out <- list(par = tmp$pars, value = tmp$values[length(tmp$values)])
  return(out)
}
fit <- FitML(spec, data = SMI, ctr = list(OptimFUN = f_optim))
summary(fit)

source("functions.R")

n.mc <- 50
IN.list = vector('list', n.mc)
for (i in 1:n.mc) {
  IN.list[[i]] = list(seed = i)
}
# f_mc(IN.list[[1]])

# parallel
library("snowfall")
nCore <- min(n.mc, 4)
snowfall::sfSetMaxCPUs(number = nCore)
snowfall::sfInit(parallel = TRUE, cpus = nCore, type = "SOCK")
snowfall::sfSource("functions.R")
OUT.list <- snowfall::sfClusterApplyLB(IN.list, fun = f_mc)
snowfall::sfStop()
save(OUT.list, file = paste0("DE.rda"), compress = TRUE)
# load("DE.rda"); n.mc <- length(OUT.list)

prob   <- matrix(data = unlist(OUT.list), nrow = 2500, ncol = n.mc, byrow = FALSE)
prob.q <- t(apply(prob, 1, quantile, probs = c(0.25, 0.75)))

## compute 50% area
par(mfrow = c(1, 1))
ylabel <- expression(paste("Pr(", s[t], " = 2 | ", hat(theta), ", ", I[t], ")"))
matplot(prob.q, type = 'l', lty = 'dashed', col = 'blue', ylim = c(0, 1),
        ylab = ylabel, xlab = "Date (year)")
par(new = TRUE)
plot(as.vector(SMI), type = 'p', col = 'black', cex = 0.5, axes = FALSE, ann = FALSE)




