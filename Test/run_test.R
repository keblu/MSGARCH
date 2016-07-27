###MSGARCH TEST
require(MSGARCH)
model.avail = c("sGARCH","gjrGARCH","eGARCH","tGARCH","GAS")
dist.avail = c("norm","std","ged")

for(i in 1:length(model.avail)){
  
  for(j in 1:length(dist.avail)){
    
    spec = MSGARCH::create.spec(model = model.avail[i], distribution = dist.avail[j])  
    
  }
}

for(i in 1:length(model.avail)){
  
  for(j in 1:length(dist.avail)){
    
    spec = MSGARCH::create.spec(model = model.avail[i], distribution = dist.avail[j], do.mix = TRUE)  
    
  }
}

for(i in 1:length(model.avail)){
  
  for(j in 1:length(dist.avail)){
    
    spec = MSGARCH::create.spec(model = model.avail[i], distribution = dist.avail[j], do.shape.ind = TRUE)  
    
  }
}



###MSGARCH TEST create

##MSGARCH
spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
##SINGLE REGIME
spec = MSGARCH::create.spec(model = c("sGARCH"), distribution = c("norm"),
                              do.skew = c(FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
##MIXTURE
spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(TRUE,TRUE), do.mix = TRUE, do.shape.ind = FALSE) 
#MSGARCH SHAPE IND
spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(TRUE,TRUE), do.mix = FALSE, do.shape.ind = TRUE) 
#MIXTURE SHAPE IND
spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(TRUE,TRUE), do.mix = TRUE, do.shape.ind = TRUE) 

###MSGARCH TEST sim

theta = spec$theta0
y = MSGARCH::sim(spec = spec, n = 1000, theta = theta, burnin = 500, do.state = TRUE)

###MSGARCH TEST bayes
data(sp500ret) 
fit.bay = MSGARCH::fit.bayes(spec = spec, y = sp500ret,
                          ctr = list(N.burn = 100,N.mcmc = 1000, N.thin = 1))
MSGARCH::BIC(fit.bay)
MSGARCH::AIC(fit.bay)
MSGARCH::DIC(fit.bay)
###MSGARCH TEST mle

fit.mle = MSGARCH::fit.mle(spec = spec, y = sp500ret)
theta = fit.mle$theta

MSGARCH::BIC(fit.mle)
MSGARCH::AIC(fit.mle)

fit.mle = MSGARCH::fit.mle(spec = spec, y = sp500ret, ctr = list(do.init = TRUE))
theta = fit.mle$theta


MSGARCH::BIC(fit.mle)
MSGARCH::AIC(fit.mle)
###MSGARCH TEST ht


ht = MSGARCH::ht(spec = spec,theta = theta, y = sp500ret)

###MSGARCH TEST kernel

kernel = MSGARCH::kernel(spec = spec, theta = theta, y = sp500ret, log = TRUE)

###MSGARCH TEST Pstate

Pstate = MSGARCH::Pstate(spec = spec, theta = theta, y = sp500ret)

###MSGARCH TEST Plast
Plast = MSGARCH::Plast(spec = spec, theta = theta, y = sp500ret)

###MSGARCH TEST cdf

set.seed(123)

x = rnorm(100)

cdf = MSGARCH::cdf(spec = spec, x = x, theta = theta, y = sp500ret, log = FALSE)

###MSGARCH TEST pdf
set.seed(123)

x = rnorm(100)

pdf = MSGARCH::pdf(spec = spec, x = x, theta = theta, y = sp500ret, log = FALSE)

###MSGARCH TEST rnd

rnd = MSGARCH::rnd(spec = spec, n = 1000, theta = theta, y = sp500ret, do.state = TRUE)

###MSGARCH TEST risk

risk = MSGARCH::risk(spec = spec, theta = spec$theta0, y = sp500ret,level = c(0.95,0.99), ES = TRUE)

###MSGARCH TEST pred
set.seed(123)

x = rnorm(100)
pred = MSGARCH::pred(spec = spec, x = x, theta = theta, y = sp500ret, log = TRUE)

###MSGARCH TEST unc.vol

unc.vol = MSGARCH::unc.vol(spec = spec, theta = theta)
