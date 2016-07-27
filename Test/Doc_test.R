spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

unc.vol = MSGARCH::unc.vol(spec = spec, theta = spec$theta0)
##################################################################
spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

y = MSGARCH::sim(spec = spec, n = 1000, theta = spec$theta0, burnin = 500, do.state = TRUE)
##################################################################
data("sp500ret")
spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                             do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

rnd = MSGARCH::rnd(spec = spec, n = 1000, theta = spec$theta0, y = sp500ret, do.state = TRUE) 
##################################################################
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                             do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

risk = MSGARCH::risk(spec = spec, theta = spec$theta0, y = sp500ret, level = c(0.95,0.99), ES = TRUE)

##################################################################
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                             do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

Pstate  = MSGARCH::Pstate(spec = spec, theta = spec$theta0, y = sp500ret)
##################################################################
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                              
set.seed(123)

x = rnorm(100)

pred = MSGARCH::pred(spec = spec, x = x, theta = spec$theta0, y = sp500ret, log = TRUE)
##################################################################

data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                             do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

Plast = MSGARCH::Plast(spec = spec, theta = spec$theta0, y = sp500ret) 

##################################################################
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                              
set.seed(123)

x = rnorm(100)

pit = MSGARCH::pit(spec = spec, x = x, theta = spec$theta0, y = sp500ret, do.norm = FALSE)
##################################################################
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                             
set.seed(123)

x = rnorm(100)

pdf = MSGARCH::pdf(spec = spec, x = x, theta = spec$theta0, y = sp500ret, log = FALSE)
##################################################################
 data("sp500ret")
 
 spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                               do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
 
 fit = MSGARCH::fit.mle(spec = spec, y = sp500ret, ctr = list(do.init = FALSE))
 ################################################################## 
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

kernel = MSGARCH::kernel(spec = spec, theta = spec$theta0, y = sp500ret, log = TRUE)
##################################################################
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

ht = MSGARCH::ht(spec = spec,theta = spec$theta0, y = sp500ret)
##################################################################
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                              
fit = MSGARCH::fit.bayes(spec = spec, y = sp500ret)

DIC = MSGARCH::DIC(fit = fit)

 spec = MSGARCH::create.spec(model = c("sGARCH","gjrGARCH"), distribution = c("norm","std"),
                              do.skew = c(TRUE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
 ##################################################################                              
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                              
set.seed(123)

x = rnorm(100)

cdf = MSGARCH::cdf(spec = spec, x = x, theta = spec$theta0, y = sp500ret, log = FALSE)
##################################################################
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                              
fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)

BIC = MSGARCH::BIC(fit = fit)
##################################################################
 data("sp500ret")
 
 spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                               do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
 
 fit = MSGARCH::fit.bayes(spec = spec, y = sp500ret, 
                          ctr = list(N.burn = 100,N.mcmc = 1000, N.thin = 1))
 ##################################################################                          
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                              
fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)

AIC = MSGARCH::AIC(fit = fit)                                                                                                                                                                                                                                        