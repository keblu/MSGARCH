require(MSGARCH)
testthat::context("Documentation Test")

test_that("unconditional variance",{
spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

unc.vol = MSGARCH::unc.vol(spec = spec, theta = spec$theta0)
expect_that(unc.vol[2],equals(1))
}
)
##################################################################
test_that("simulation",{
spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
set.seed(123)

sim = MSGARCH::sim(spec = spec, n = 1000, theta = spec$theta0, burnin = 500, do.state = TRUE)

plot(sim)
expect_that(round(sim$draws[1],6),equals(0.213346))
}
)
##################################################################
test_that("rnd",{
data("sp500ret")
spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                             do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

set.seed(123)

rnd = MSGARCH::rnd(spec = spec, n = 1000, theta = spec$theta0, y = sp500ret, do.state = TRUE)

summary(rnd)

plot(rnd)

expect_that(round(rnd$draws[1],6),equals(0.566268))
}
)
##################################################################
test_that("risk",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                             do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

set.seed(123)

risk = MSGARCH::risk(spec = spec, theta = spec$theta0, y = sp500ret, level = c(0.95,0.99), ES = TRUE)
expect_that(round(risk$VaR[1],6),equals(-1.163479))
expect_that(round(risk$ES[1],6),equals(-1.459050))
expect_that(round(risk$VaR[2],6),equals(-1.645531))
expect_that(round(risk$ES[2],6),equals(-1.885226))
}
)
##################################################################
test_that("Pstate",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                             do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

Pstate  = MSGARCH::Pstate(spec = spec, theta = spec$theta0, y = sp500ret)
expect_that(round(Pstate[5,1,1],6),equals(0.5))
expect_that(round(Pstate[5,1,2],6),equals(0.5))
}
)
##################################################################
test_that("pred",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                              
set.seed(123)

x = rnorm(100)

pred = MSGARCH::pred(spec = spec, x = x, theta = spec$theta0, y = sp500ret, log = TRUE)
expect_that(round(pred[1],6),equals(-0.886623))
}
)
##################################################################
test_that("Plast",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                             do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

Plast = MSGARCH::Plast(spec = spec, theta = spec$theta0, y = sp500ret) 
expect_that(round(Plast[1,1],6),equals(0.5))
}
)
##################################################################
test_that("pit",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                              
set.seed(123)

x = rnorm(100)

pit = MSGARCH::pit(spec = spec, x = x, theta = spec$theta0, y = sp500ret, do.norm = FALSE)
expect_that(round(pit[1],6),equals(0.214074))
}
)
##################################################################
test_that("pdf",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                             
set.seed(123)

x = rnorm(100)

pdf = MSGARCH::pdf(spec = spec, x = x, theta = spec$theta0, y = sp500ret, log = FALSE)
expect_that(round(pdf[1],6),equals(0.412045))
}
)
##################################################################
test_that("mle",{
 data("sp500ret")
 
 spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                               do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
 
 set.seed(123)
 
 fit = MSGARCH::fit.mle(spec = spec, y = sp500ret, ctr = list(do.init = TRUE, NP = 100, itermax = 100))
 expect_that(round(fit$theta[1],6),equals(0.009084))
 expect_that(round(fit$ll_likelihood[1],2),equals(17336.25))
 expect_that(fit$is.init,is_true())
}
)
 ################################################################## 
test_that("kernel",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

kernel = MSGARCH::kernel(spec = spec, theta = spec$theta0, y = sp500ret, log = TRUE)
expect_that(round(kernel,3),equals(-3186.016))
}
)
##################################################################
test_that("ht",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

ht = MSGARCH::ht(spec = spec,theta = spec$theta0, y = sp500ret)

plot(ht)

expect_that(round(ht[5,1,1],6),equals(0.839528))
expect_that(round(ht[5,1,2],6),equals(0.839528))
}
)
##################################################################
test_that("DIC",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

set.seed(123)

fit = MSGARCH::fit.bayes(spec = spec, y = sp500ret)

DIC = MSGARCH::DIC(fit = fit)
expect_that(round(DIC$DIC,2),equals(-26861.06))
expect_that(round(DIC$IC,2),equals(-26841.27))
expect_that(round(DIC$pD,2),equals(19.79))
expect_that(round(DIC$pV,1),equals(202807.3))
expect_that(round(DIC$D.bar,2),equals(-26880.85))
expect_that(round(DIC$D.hat,2),equals(-26900.63))
}
)
 ##################################################################
test_that("cdf",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
                              
set.seed(123)

x = rnorm(100)

cdf = MSGARCH::cdf(spec = spec, x = x, theta = spec$theta0, y = sp500ret, log = FALSE)
expect_that(round(cdf[,1],6),equals(0.214074))
}
)
##################################################################
test_that("BIC",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

set.seed(123)

fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)

BIC = MSGARCH::BIC(fit = fit)
expect_that(round(BIC,2),equals(28830.82))
}
)
##################################################################
test_that("bay",{
 data("sp500ret")
 
 spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                               do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
 
 set.seed(123)
 
 fit = MSGARCH::fit.bayes(spec = spec, y = sp500ret, 
                          ctr = list(N.burn = 100,N.mcmc = 1000, N.thin = 1))
 expect_that(round(fit$accept,3),equals(0.016))
 expect_that(unname(round(fit$theta[300,1],6)),equals(0.000104))
}
)
 ##################################################################
test_that("AIC",{
data("sp500ret")

spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 

set.seed(123)

fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)

AIC = MSGARCH::AIC(fit = fit) 
expect_that(round(AIC,2),equals(28883.76))
}
)
##################################################################  
