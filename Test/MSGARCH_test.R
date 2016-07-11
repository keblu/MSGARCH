rm(list = ls())

# WARNING : MLE STARTING NEEDS A BETTER STARTING VALUE => 
# - DEOPTIM
# - MCMC ADAPTIVE WITH NAIVE STARTING POINT

require("MSGARCH")     # <--- YOU MUST CREATE THIS PACKAGE FIRST

# create spec

models = list(MSGARCH::Garch_normal_sym,MSGARCH::Garch_normal_sym)
models = list(MSGARCH::Garch_normal_sym,MSGARCH::Garch_normal_sym,MSGARCH::Garch_normal_sym)
models = list(MSGARCH::Gjr_normal_sym,MSGARCH::Gjr_normal_sym)
models = list(MSGARCH::Gjr_normal_sym,MSGARCH::Gjr_normal_sym,MSGARCH::Gjr_normal_sym)
models = list(MSGARCH::Egarch_normal_sym,MSGARCH::Egarch_normal_sym)
models = list(MSGARCH::Egarch_normal_sym,MSGARCH::Egarch_normal_sym,MSGARCH::Egarch_normal_sym)
models = list(MSGARCH::Tgarch_normal_sym,MSGARCH::Tgarch_normal_sym)
models = list(MSGARCH::Tgarch_normal_sym,MSGARCH::Tgarch_normal_sym,MSGARCH::Tgarch_normal_sym)
models = list(MSGARCH::Gas_normal_sym,MSGARCH::Gas_normal_sym)
models = list(MSGARCH::Gas_normal_sym,MSGARCH::Gas_normal_sym,MSGARCH::Gas_normal_sym)

models = list(MSGARCH::Garch_student_sym,MSGARCH::Garch_student_sym)
models = list(MSGARCH::Garch_student_sym,MSGARCH::Garch_student_sym,MSGARCH::Garch_student_sym)
models = list(MSGARCH::Gjr_student_sym,MSGARCH::Gjr_student_sym)
models = list(MSGARCH::Gjr_student_sym,MSGARCH::Gjr_student_sym,MSGARCH::Gjr_student_sym)
models = list(MSGARCH::Egarch_student_sym,MSGARCH::Egarch_student_sym)
models = list(MSGARCH::Egarch_student_sym,MSGARCH::Egarch_student_sym,MSGARCH::Egarch_student_sym)
models = list(MSGARCH::Tgarch_student_sym,MSGARCH::Tgarch_student_sym)
models = list(MSGARCH::Tgarch_student_sym,MSGARCH::Tgarch_student_sym,MSGARCH::Tgarch_student_sym)
models = list(MSGARCH::Gas_student_sym,MSGARCH::Gas_student_sym)
models = list(MSGARCH::Gas_student_sym,MSGARCH::Gas_student_sym,MSGARCH::Gas_student_sym)

models = list(MSGARCH::Garch_ged_sym,MSGARCH::Garch_ged_sym)
models = list(MSGARCH::Garch_ged_sym,MSGARCH::Garch_ged_sym,MSGARCH::Garch_ged_sym)
models = list(MSGARCH::Gjr_ged_sym,MSGARCH::Gjr_ged_sym)
models = list(MSGARCH::Gjr_ged_sym,MSGARCH::Gjr_ged_sym,MSGARCH::Gjr_ged_sym)
models = list(MSGARCH::Egarch_ged_sym,MSGARCH::Egarch_ged_sym)
models = list(MSGARCH::Egarch_ged_sym,MSGARCH::Egarch_ged_sym,MSGARCH::Egarch_ged_sym)
models = list(MSGARCH::Tgarch_ged_sym,MSGARCH::Tgarch_ged_sym)
models = list(MSGARCH::Tgarch_ged_sym,MSGARCH::Tgarch_ged_sym,MSGARCH::Tgarch_ged_sym)
models = list(MSGARCH::Gas_ged_sym,MSGARCH::Gas_ged_sym)
models = list(MSGARCH::Gas_ged_sym,MSGARCH::Gas_ged_sym,MSGARCH::Gas_ged_sym)

models = list(MSGARCH::Garch_normal_skew,MSGARCH::Garch_normal_skew)
models = list(MSGARCH::Garch_normal_skew,MSGARCH::Garch_normal_skew,MSGARCH::Garch_normal_skew)
models = list(MSGARCH::Gjr_normal_skew,MSGARCH::Gjr_normal_skew)
models = list(MSGARCH::Gjr_normal_skew,MSGARCH::Gjr_normal_skew,MSGARCH::Gjr_normal_skew)
models = list(MSGARCH::Egarch_normal_skew,MSGARCH::Egarch_normal_skew)
models = list(MSGARCH::Egarch_normal_skew,MSGARCH::Egarch_normal_skew,MSGARCH::Egarch_normal_skew)
models = list(MSGARCH::Tgarch_normal_skew,MSGARCH::Tgarch_normal_skew)
models = list(MSGARCH::Tgarch_normal_skew,MSGARCH::Tgarch_normal_skew,MSGARCH::Tgarch_normal_skew)
models = list(MSGARCH::Gas_normal_skew,MSGARCH::Gas_normal_skew)
models = list(MSGARCH::Gas_normal_skew,MSGARCH::Gas_normal_skew,MSGARCH::Gas_normal_skew)

models = list(MSGARCH::Garch_student_skew,MSGARCH::Garch_student_skew)
models = list(MSGARCH::Garch_student_skew,MSGARCH::Garch_student_skew,MSGARCH::Garch_student_skew)
models = list(MSGARCH::Gjr_student_skew,MSGARCH::Gjr_student_skew)
models = list(MSGARCH::Gjr_student_skew,MSGARCH::Gjr_student_skew,MSGARCH::Gjr_student_skew)
models = list(MSGARCH::Egarch_student_skew,MSGARCH::Egarch_student_skew)
models = list(MSGARCH::Egarch_student_skew,MSGARCH::Egarch_student_skew,MSGARCH::Egarch_student_skew)
models = list(MSGARCH::Tgarch_student_skew,MSGARCH::Tgarch_student_skew)
models = list(MSGARCH::Tgarch_student_skew,MSGARCH::Tgarch_student_skew,MSGARCH::Tgarch_student_skew)
models = list(MSGARCH::Gas_student_skew,MSGARCH::Gas_student_skew)
models = list(MSGARCH::Gas_student_skew,MSGARCH::Gas_student_skew,MSGARCH::Gas_student_skew)

models = list(MSGARCH::Garch_ged_skew,MSGARCH::Garch_ged_skew)
models = list(MSGARCH::Garch_ged_skew,MSGARCH::Garch_ged_skew,MSGARCH::Garch_ged_skew)
models = list(MSGARCH::Gjr_ged_skew,MSGARCH::Gjr_ged_skew)
models = list(MSGARCH::Gjr_ged_skew,MSGARCH::Gjr_ged_skew,MSGARCH::Gjr_ged_skew)
models = list(MSGARCH::Egarch_ged_skew,MSGARCH::Egarch_ged_skew)
models = list(MSGARCH::Egarch_ged_skew,MSGARCH::Egarch_ged_skew,MSGARCH::Egarch_ged_skew)
models = list(MSGARCH::Tgarch_ged_skew,MSGARCH::Tgarch_ged_skew)
models = list(MSGARCH::Tgarch_ged_skew,MSGARCH::Tgarch_ged_skew,MSGARCH::Tgarch_ged_skew)
models = list(MSGARCH::Gas_ged_skew,MSGARCH::Gas_ged_skew)
models = list(MSGARCH::Gas_ged_skew,MSGARCH::Gas_ged_skew,MSGARCH::Gas_ged_skew)

spec = MSGARCH::f.spec(models,mixture = FALSE, DistRegInd = FALSE)
spec = MSGARCH::f.spec(models,mixture = TRUE, DistRegInd = FALSE)
spec = MSGARCH::f.spec(models,mixture = FALSE, DistRegInd = TRUE)
spec = MSGARCH::f.spec(models,mixture = TRUE, DistRegInd = TRUE)

spec$K
spec$name
spec$label
spec$NbParamsModel
spec$NbParams
spec$inequb
spec$ineqlb
spec$upper
spec$lower
spec$kSigma
spec$Sigma0
spec$theta0
spec$mixture
spec$DistRegInd

################## SIM TEST ########################################################################
thetaSim = spec$theta0
nobs = 5000
y = spec$f.sim(n = nobs, theta = thetaSim, burnin = 500, outputState = TRUE)$value

################## Real Data TEST ########################################################################
quantmod::getSymbols('AAPL',src='yahoo')
y = quantmod::periodReturn(AAPL,period='daily',subset='2007::')
y = as.vector(y)*100
################## MLE ESTIMATION TEST ######################################################################## 
theta_mle = MSGARCH::f.estim.mle(y = y, spec = spec)
theta_mle_init = MSGARCH::f.estim.mle(y = y, spec = spec, ctr = list(do.init = TRUE))

################## ht and Unconditional Volatility function TEST ###################################################### 
spec$f.ht(theta = theta_mle$theta0, y = y)
spec$f.unc.vol(theta = theta_mle$theta0)

################## PDF, CDF and RND TEST ########################################################################################## 

x = seq(from = -5, to = 5, length.out = 100)
pdf_test = spec$f.pdf(x = x, theta = theta_mle$theta0, y = y, log = FALSE)
rnd_test = spec$f.rnd(n = length(y), theta = theta_mle$theta0, y = y, outputState = TRUE)$value

f.g = function(x){
  out = x ^ 2 * spec$f.pdf(x = x, theta = theta_mle$theta0, y = y, log = FALSE)
  return(out)
}

integrate(f = f.g, lower = -10, upper = 10)
mean(rnd_test ^ 2)

hist(rnd_test, nclass = round(10 * log(length(y))), freq = FALSE)
lines(x, pdf_test)

############################################################################################################ 

