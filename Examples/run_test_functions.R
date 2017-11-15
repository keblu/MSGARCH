library("MSGARCH")

####################### FITML #########################
#Loop over all model and all distribution and fit by ML the 
#SMI up to K.max regime with similar conditional variance process and distribution
#store the BIC of each estimation
data("SMI", package = "MSGARCH")
K.max = 2
variance.model = c("sARCH","sGARCH","gjrGARCH","tGARCH","eGARCH")
dist.spec = c("norm","snorm","std","sstd","ged","sged")
out = array(NA,dim = c(5,6,K.max),dimnames = list(variance.model,dist.spec, paste0("K=",1:K.max)))
do.mix = FALSE
for(dist in dist.spec){
  for(model in variance.model){
    for(k in 1:K.max){
      spec <- CreateSpec(variance.spec = list(model = model),
                         distribution.spec = list(distribution = dist),
                         switch.spec = list(do.mix = do.mix, K = k))
      fit <- FitML(spec = spec, data = SMI)
      out[model, dist, k] <- summary(fit)$BIC
    }
  }
}

#fail 4 regime GJR garch sstd 
####################### FITMCMC ########################
#Loop over all model and all distribution and fit by MCMC the 
#SMI up to K.max regime with similar conditional variance process and distribution
#store the DIC of each estimation
data("SMI", package = "MSGARCH")
K.max = 2
variance.model = c("sARCH","sGARCH","gjrGARCH","tGARCH","eGARCH")
dist.spec = c("norm","snorm","std","sstd","ged","sged")
out = array(NA,dim = c(5,6,K.max),dimnames = list(variance.model,dist.spec, paste0("K=",1:K.max)))
do.mix = FALSE
for(dist in dist.spec){
  for(model in variance.model){
    for(k in 1:K.max){
      spec <- CreateSpec(variance.spec = list(model = model),
                         distribution.spec = list(distribution = dist),
                         switch.spec = list(do.mix = do.mix, K = k))
      fit <- FitMCMC(spec = spec, data = SMI)
      out[model, dist, k] <- summary(fit)$DIC
    }
  }
}

####################### PIT #########################

# load data
data("SMI", package = "MSGARCH")

# create model specification
# MS(2)-GARCH(1,1)-Normal (default)
spec <- CreateSpec()

# fit the model on the data by ML
fit <- FitML(spec = spec, data = SMI)

# run PIT method on data in-sample
pit.data.its <- PIT(object = fit, do.norm = TRUE, do.its = TRUE)
length(pit.data.its) == 2500
#see if it fit well with qqplot
qqnorm(pit.data.its)
qqline(pit.data.its)
#PIT test with rugarch
rugarch::HLTest(pit.data.its)
# create a mesh
x <- seq(-3,3,0.01)

# run PIT method on mesh in-sample
pit.x.its <- PIT(object = fit, x = x, do.norm = TRUE, do.its = TRUE)
dim(pit.x.its) == c(2500, 601)
plot(zoo::zoo(t(pit.x.its), order.by = x),plot.type = "single")

# run PIT method on mesh at T + nahead
set.seed(123)
pit.x.ots <- PIT(object = fit, x = x, nahead = 5, do.norm = FALSE, do.its = FALSE, ctr = list(nsim = 50000L))
dim(pit.x.ots)  == c(5, 601)
plot(zoo::zoo(t(pit.x.ots), order.by = x))

#Simulate from the fitted parameter
set.seed(123)
sim.series <- simulate(object = spec, par = fit$par, nahead = 1000L, nsim = 1L)
dim(sim.series$draw) == c(1000,1)
dim(sim.series$state) == c(1000,1)
dim(sim.series$CondVol) == c(1000,1,2)
sim.series <- as.vector(sim.series$draw)

# run PIT method on the simualed serie with the true pararameter 
pit.sim.its <- PIT(object = spec, par = fit$par, data = sim.series, do.norm = TRUE, do.its = TRUE)
length(pit.sim.its) == 1000L
# look if the fit is OK with qqplot
qqnorm(pit.sim.its)
qqline(pit.sim.its)

####################### PRED #########################

# load data
data("SMI", package = "MSGARCH")

# create model specification
# MS(2)-GARCH(1,1)-Normal (default)
spec <- CreateSpec()

# fit the model on the data by ML
fit <- FitML(spec = spec, data = SMI)

# run PredPDF method on data in-sample (log-likelihood)
pred.data.its <- PredPdf(object = fit, log = TRUE, do.its = TRUE)
length(pred.data.its) == 2500
round(sum(pred.data.its),4) == -3391.5866

# create a mesh
x <- seq(-3,3,0.01)

# run PredPDF method on x mesh in-sample
pred.x.its <- PredPdf(object = fit, x = x, log = FALSE, do.its = TRUE)
dim(pred.x.its) == c(2500,601)
#plot all predictive
plot(zoo::zoo(t(pred.x.its), order.by = x), plot.type = "single")

# run PredPDF method on x mesh in-sample
pred.x.ots <- PredPdf(object = fit, x = x, log = FALSE, do.its = FALSE, nahead = 5, ctr = list(nsim = 50000L))
dim(pred.x.ots) == c(5,601)
#plot 1 step ahead predictive
plot(zoo::zoo(t(pred.x.ots), order.by = x), plot.type = "single")


####################### Volatility ########################

# load data
data("SMI", package = "MSGARCH")

# create model specification
# MS(2)-GARCH(1,1)-Normal (default)
spec <- CreateSpec()

# fit the model on the data by ML
fit <- FitML(spec = spec, data = SMI)

# compute the In-sample conditional volatility from the fitted model
cond.vol.its <- Volatility(object = fit)
length(cond.vol.its) == 2500
plot(cond.vol.its)

#standardize residual variance shoyld be near 1
var(SMI/cond.vol.its) - 1 < 0.01


####################### Forecast ########################

# compute the forcast of the conditional volatility from the fitted model up to nahead = 1000
cond.vol.ots <- predict(object = fit, nahead = 5, do.return.draw = TRUE, ctr = list(nsim = 50000L))
length(cond.vol.ots$vol) == 5
dim(cond.vol.ots$draw) == c(5, 50000)
par(mfrow = c(1,1))
plot(cond.vol.ots)
####################### UncVol ########################
par(mfrow = c(1,1))
# load data
data("SMI", package = "MSGARCH")

# create model specification
# MS(2)-GARCH(1,1)-Normal (default)
spec <- CreateSpec()

# fit the model on the data by ML
fit <- FitML(spec = spec, data = SMI)

#evaluate error due to simulation for overall uncvol process
out = NULL
for(i in 1:100){
  out[i] <- UncVol(fit)
}

boxplot(out)
out <- matrix(NA, nrow = 100, ncol = 2)
#evaluate error due to simulation for overall each regime uncvol
for(i in 1:100){
  out[i,] <- sapply(ExtractStateFit(fit), UncVol)
}

boxplot(out)

####################### Sim ########################

# load data
data("SMI", package = "MSGARCH")

# create model specification
# MS(2)-GARCH(1,1)-Normal (default)
spec <- CreateSpec()

# fit the model on the data by ML
fit <- FitML(spec = spec, data = SMI)


simul <- simulate(object = spec, nahead = 100,nsim = 1, par = fit$par)
dim(simul$draw) == c(100,1)
dim(simul$state) == c(100,1)
dim(simul$CondVol) == c(100,1,2)
plot(simul$draw)
plot(simul$state)
plot(zoo::zoo(simul$CondVol[,,]),plot.type = "single",col = c("blue","red"))
####################### State ########################

# load data
data("SMI", package = "MSGARCH")

# create model specification
# MS(2)-GARCH(1,1)-Normal (default)
spec <- CreateSpec()

# fit the model on the data by ML
fit <- FitML(spec = spec, data = SMI)

state = State(fit)
dim(state$FiltProb) == c(2500,1,2)
dim(state$PredProb) == c(2501,1,2)
dim(state$SmoothProb) == c(2501,1,2)
dim(state$Viterbi) == c(2500,1)

plot(state,type.prob = "smoothed")
plot(state,type.prob = "filtered")
plot(state,type.prob = "predictive")
plot(state,type.prob = "viterbi")

####################### RISK ########################

# load data
data("SMI", package = "MSGARCH")

# create model specification
# MS(2)-GARCH(1,1)-Normal (default)
spec <- CreateSpec()

# fit the model on the data by ML
fit <- FitML(spec = spec, data = SMI)

risk <- Risk(fit, do.its = TRUE)
dim(risk$VaR) == c(2500,2)
dim(risk$ES) == c(2500,2)
#Risk test with rugarch
library(rugarch)
rugarch::VaRDurTest(actual = SMI, VaR = risk$VaR[,1],alpha = as.numeric(colnames(risk$VaR)[1]))
rugarch::VaRTest(actual = SMI, VaR = risk$VaR[,1],alpha = as.numeric(colnames(risk$VaR)[1]))
rugarch::ESTest(actual = SMI, VaR = risk$VaR[,1], ES = risk$ES[,1], alpha = as.numeric(colnames(risk$VaR)[1]))
plot(risk)
#step ahead Risk measure calculation
risk <- Risk(fit, nahead = 1000, do.its = FALSE, ctr = list(nsim = 10000L))
dim(risk$VaR) == c(1000,2)
dim(risk$ES) == c(1000,2)
plot(risk)


####################### ExtractStateFit ########################

# load data
data("SMI", package = "MSGARCH")

# create model specification
# MS(2)-GARCH(1,1)-Normal (default)
spec <- CreateSpec()

# fit the model on the data by ML
fit <- FitML(spec = spec, data = SMI)

SR.fit = ExtractStateFit(fit)

print(SR.fit)

