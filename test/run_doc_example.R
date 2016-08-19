
#TOFIX: summary method for fit objects (MLE AND BAYESIAN)
#TOFIX: error with DIC computation

# http://kbroman.org/pkg_primer/pages/github.html
# install.packages("devtools")
# require("devtools")
# devtools::install_github("keblu/MSGARCH")
# install.packages("./MSGARCH_0.16.tar.gz", repos = NULL, type = "source")

rm(list = ls())
require("MSGARCH")
require("coda")

# load data
data("sp500")
sp500 = sp500[1:1000]
plot(sp500)

# create model specification
spec = MSGARCH::create.spec() 
#spec = MSGARCH::create.spec(do.mix = TRUE)
# fit the model on the data with ML estimation without DEoptim intialization
fit = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = list(do.init = FALSE))
# fit the model on the data with ML estimation using DEoptim intialization
#set.seed(123)
#fit = MSGARCH::fit.mle(spec = spec, y = sp500, 
#                       ctr = list(do.init = TRUE, NP = 500, itermax = 500))
#TODO summary method for fit object? # KB OK, DA does not work
summary(fit)

# compute the unconditional volatility in each regime
unc.vol = MSGARCH::unc.vol(fit)
unc.vol

# run pdf method in-sample
transmat.mle = MSGARCH::transmat(fit)
transmat.mle

# simulate a MSGARCH
set.seed(123)
sim = MSGARCH::sim(object = spec, n = 1000, m = 1, theta = fit$theta, burnin = 500)
plot(sim)

# generate random draws
set.seed(123)
simahead = MSGARCH::simahead(object = fit, n = 15, m = 100)
plot(simahead)
matplot(t(simahead$draws), pch = 20, type = "l")
#TOFIX why "average simulated state" in title? # answer: Because we see 1 and 0 everywhere and the plot makes no sense at all so the average simulated state seems more adequat since we have a better view to what are the states that comes next and the "trend"
#TODO impove speed answer: have to do C++ code so version
#TOFIX why plot of states and no simulation? answer: it is there but in a different plot

# compute the Value-at-Risk and Expected-shortfall 
risk.its = MSGARCH::risk(object = fit, level = c(0.95,0.99), ES = TRUE, do.its = TRUE)
plot(risk.its)  

# Risk estimation at T + 1                     
risk.ots = MSGARCH::risk(object = fit, level = c(0.95,0.99), ES = TRUE, do.its = FALSE)
risk.ots$VaR
risk.ots$ES

# compute the filtered state probabilities
Pstate = MSGARCH::Pstate(object = fit)
plot(Pstate, ylim = c(0, 1))
#TODO we should be abot to pass arguments, ... # OK but for what parameters

# run pred method in-sample     
pred.its = MSGARCH::pred(object = fit, log = TRUE, do.its = TRUE)  
sum(pred.its$pred, na.rm = TRUE)

# run pred method on mesh at T + 1
x = seq(from = min(fit$y), to = max(fit$y), length.out = 100)
pred = MSGARCH::pred(object = fit, x = x, log = FALSE, do.its = FALSE)
plot(pred)

# run pit method in-sample              
pit.its = MSGARCH::pit(object = fit, do.norm = FALSE, do.its = TRUE)                              
plot(pit.its)  
                                                                          
# run pdf method in-sample
pdf.its = MSGARCH::pdf(object = fit, log = FALSE, do.its = TRUE)
                  
# run pdf method on mesh at T + 1
x = seq(from = min(fit$y), to = max(fit$y), length.out = 100)
pdf = MSGARCH::pdf(object = fit, x = x, log = FALSE, do.its = FALSE)
plot(pdf)

# run cdf method on mesh at T + 1
x = seq(from = min(fit$y), to = max(fit$y), length.out = 100)
cdf = MSGARCH::cdf(object = fit, x = x, log = FALSE, do.its = FALSE)
plot(cdf)

# compute the kernel
kernel = MSGARCH::kernel(fit, log = TRUE)
kernel

# compute AIC
AIC = MSGARCH::AIC(fit) 
AIC

# compute BIC
BIC = MSGARCH::BIC(fit) 
BIC

# compute the conditional volatility
ht = MSGARCH::ht(object = fit)
plot(ht)
#y <- f.check.y(y)
#theta <- f.check.theta(object, theta)
#out <- object$rcpp.func$calc_ht(theta, y)
#class(out) <- "MSGARCH_HT"
#out <- sqrt(out)
                                                             
# fit the model by Bayesian estimation 
set.seed(123)                                                           
fit = MSGARCH::fit.bayes(spec = spec, y = sp500, ctr = list(N.burn = 500, N.mcmc = 1000, N.thin = 1))
summary(fit)
par(mfrow = c(3,3))
coda::traceplot(fit$theta)
coda::autocorr.plot(fit$theta)
pairs(x = as.matrix(fit$theta[,1:5]), pch = 20, cex = 0.8)
par(mfrow = c(1,1))
# run pred method on mesh at T + 1
x = seq(from = min(fit$y), to = max(fit$y), length.out = 100)
pred = MSGARCH::pred(object = fit, x = x, log = FALSE, do.its = FALSE)
plot(pred)

# compute DIC
DIC = MSGARCH::DIC(fit)
DIC
#TOFIX WE still have a problem with DIC output => pD is missing
