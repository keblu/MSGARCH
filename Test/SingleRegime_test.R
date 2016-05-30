rm(list = ls())

require("MSGARCH")  

nobs = 5000

# create spec

models = list(MSGARCH::Garch_normal_sym)
thetaSim = c(0.01,  0.15,  0.8)

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
y = spec$f.sim(n = nobs, theta = thetaSim, burnin = 500, outputState = FALSE)

################## MLE ESTIMATION TEST ######################################################################## 
MSGARCH::f.estim.mle(y = y,spec = spec)
MSGARCH::f.estim.mle(y = y,spec = spec,ctr = list(do.init = TRUE))

################## ht and Unconditional Volatility function TEST ###################################################### 
ht = spec$f.ht(theta = thetaSim, y = y)
spec$f.unc.vol(theta = thetaSim)

################## PDF, CDF and RND TEST ########################################################################################## 

x = seq(from = -5, to = 5, length.out = 100)
pdf_test = spec$f.pdf(x = x, theta = thetaSim, y = y, log = FALSE)
rnd_test = spec$f.rnd(n = nobs, theta = thetaSim, y = y,outputState = TRUE)

f.g = function(x){
  out = x ^ 2 * spec$f.pdf(x = x, theta = thetaSim, y = y, log = FALSE)
  return(out)
}

integrate(f = f.g, lower = -10, upper = 10)
mean(rnd_test ^ 2)

hist(rnd_test, nclass = round(10 * log(nobs)), freq = FALSE)
lines(x, pdf_test)

############################################################################################################ 

  