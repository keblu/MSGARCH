###MSGARCH TEST
require(MSGARCH)
model.avail = c("sGARCH","gjrGARCH","eGARCH","tGARCH","GAS")
dist.avail = c("norm","std","ged")

for(i in 1:length(model.avail)){
  
  for(j in 1:length(dist.avail)){
    
    spec = MSGARCH::f.create.spec(model = model.avail[i], distribution = dist.avail[j])  
    
  }
}

for(i in 1:length(model.avail)){
  
  for(j in 1:length(dist.avail)){
    
    spec = MSGARCH::f.create.spec(model = model.avail[i], distribution = dist.avail[j], do.mix = TRUE)  
    
  }
}

for(i in 1:length(model.avail)){
  
  for(j in 1:length(dist.avail)){
    
    spec = MSGARCH::f.create.spec(model = model.avail[i], distribution = dist.avail[j], do.shape.ind = TRUE)  
    
  }
}



###MSGARCH TEST create

##MSGARCH
spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
##SINGLE REGIME
spec = MSGARCH::f.create.spec(model = c("sGARCH"), distribution = c("norm"),
                              do.skew = c(FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
##MIXTURE
spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(TRUE,TRUE), do.mix = TRUE, do.shape.ind = FALSE) 
#MSGARCH SHAPE IND
spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(TRUE,TRUE), do.mix = FALSE, do.shape.ind = TRUE) 
#MIXTURE SHAPE IND
spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
                              do.skew = c(TRUE,TRUE), do.mix = TRUE, do.shape.ind = TRUE) 

###MSGARCH TEST sim

theta = spec$theta0
y = spec$f.sim(n = 1000, theta = theta, burnin = 500, do.state = TRUE)

###MSGARCH TEST bayes
data(sp500ret) 
theta = MSGARCH::f.estim.bayes(y = sp500ret, spec = spec, 
                          ctr = list(N.burn = 100,N.mcmc = 1000, N.thin = 1))

###MSGARCH TEST mle

theta = MSGARCH::f.estim.mle(y = sp500ret, spec = spec)
theta = theta$theta

theta = MSGARCH::f.estim.mle(y = sp500ret, spec = spec, ctr = (do.init = TRUE))
theta = theta$theta
###MSGARCH TEST ht


ht = spec$f.ht(theta = theta, y = sp500ret)

###MSGARCH TEST kernel


kernel = spec$f.kernel(theta = theta, y = sp500ret, log = TRUE)

###MSGARCH TEST Pstate


Pstate  = spec$f.Pstate(theta = theta, y = sp500ret)

###MSGARCH TEST Plast


Plast  = spec$f.Plast(theta = theta, y = sp500ret)

###MSGARCH TEST cdf


cdf = spec$f.cdf(x = rnorm(100), theta = theta, y = sp500ret, log = FALSE)
  
###MSGARCH TEST pdf

pdf = spec$f.cdf(x = rnorm(100), theta = theta, y = sp500ret, log = FALSE)

###MSGARCH TEST rnd

rnd = spec$f.rnd(n = 1000, theta = theta, y = sp500ret, do.state = TRUE)
 
###MSGARCH TEST risk


risk = spec$f.risk(theta = spec$theta0, y = sp500ret,level = c(0.95,0.99), ES = TRUE)

###MSGARCH TEST pred


pred = spec$f.pred(x = rnorm(100), theta = theta, y = sp500ret, log = TRUE)

###MSGARCH TEST unc.vol

unc.vol = spec$f.unc.vol(theta = theta)
