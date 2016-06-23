rm(list = ls())

# WARNING : MLE STARTING NEEDS A BETTER STARTING VALUE => 
# - DEOPTIM
# - MCMC ADAPTIVE WITH NAIVE STARTING POINT

require("MSGARCH")     # <--- YOU MUST CREATE THIS PACKAGE FIRST



# create spec

models = list(MSGARCH::Garch_normal_sym)
models = list(MSGARCH::Gjr_normal_sym)
models = list(MSGARCH::Egarch_normal_sym)
models = list(MSGARCH::Tgarch_normal_sym)

models = list(MSGARCH::Garch_student_sym)
models = list(MSGARCH::Gjr_student_sym)
models = list(MSGARCH::Egarch_student_sym)
models = list(MSGARCH::Tgarch_student_sym)

models = list(MSGARCH::Garch_ged_sym)
models = list(MSGARCH::Gjr_ged_sym)
models = list(MSGARCH::Egarch_ged_sym)
models = list(MSGARCH::Tgarch_ged_sym)

# y = Data[,1]
spec = MSGARCH::f.spec(models)

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
################## SIM TEST ########################################################################
thetaSim = spec$theta0
nobs = 5000
y = spec$f.sim(n = nobs, theta = thetaSim, burnin = 500, outputState = FALSE)

################## MLE ESTIMATION TEST ######################################################################## 
theta_mle = MSGARCH::f.estim.mle(y = y,spec = spec)
theta_mle_init = MSGARCH::f.estim.mle(y = y,spec = spec,ctr = list(do.init = TRUE))

################## ht and Unconditional Volatility function TEST ###################################################### 
spec$f.ht(theta = theta_mle$theta0, y = y)
spec$f.unc.vol(theta = theta_mle$theta0)

################## PDF, CDF and RND TEST ########################################################################################## 

x = seq(from = -5, to = 5, length.out = 100)
pdf_test = spec$f.pdf(x = x, theta = theta_mle$theta0, y = y, log = FALSE)
rnd_test = spec$f.rnd(n = length(y), theta = theta_mle$theta0, y = y, outputState = TRUE)

f.g = function(x){
  out = x ^ 2 * spec$f.pdf(x = x, theta = theta_mle$theta0, y = y, log = FALSE)
  return(out)
}

integrate(f = f.g, lower = -10, upper = 10)
mean(rnd_test ^ 2)

hist(rnd_test, nclass = round(10 * log(length(y))), freq = FALSE)
lines(x, pdf_test)

############################################################################################################ 

  