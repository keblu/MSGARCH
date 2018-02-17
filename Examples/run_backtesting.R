#################################################################################
### DESCRIPTION

### This code is used in the illustrations of backtesting with MSGARCH models

rm(list = ls())
library("MSGARCH")
data("SMI", package = "MSGARCH")

## Model specification
models <- list(GARCH   = MSGARCH::CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                                             distribution.spec = list(c("sstd")),
                                             switch.spec = list(K = 1)),
               MSGARCH = MSGARCH::CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                                             distribution.spec = list(c("sstd")),
                                             switch.spec = list(K = 2)))

n_backtest <- 1000 # number of out-of-sample evaluation 
n_its      <- 1500 # fit sample size
alpha      <- 0.05 # risk Level 
k_update   <- 100   # frequence at which we update the model (reestimation)

## Initialization 
VaR   <- matrix(NA, nrow = n_backtest, ncol = length(models))
y_ots <- matrix(NA, nrow = n_backtest, ncol = 1)
model_fit <- vector(mode = "list", length = length(models))

## Backtest 

# iterate over out-of-sample time
for (i in 1:n_backtest) { 
  cat("Backtest - Iteration: ", i, "\n")
  y_its    <- SMI[i:(n_its + i - 1)] # in-sample data 
  y_ots[i] <- SMI[n_its + i]         # out-of-sample data
  
  # iterate over models
  for (j in 1:length(models)) { 
    
    # do we update the model estimation
    if (k_update == 1 || i %% k_update == 1) {
      cat("Model", j, "is reestimated\n")
      model_fit[[j]] <- MSGARCH::FitML(spec = models[[j]], data = y_its) 
    }
    
    # calculate VaR 1-step ahead
    VaR[i,j] <- MSGARCH::Risk(model_fit[[j]]$spec, par = model_fit[[j]]$par,
                              data = y_its,
                              n.ahead = 1,
                              alpha   = alpha,
                              do.es   = FALSE,
                              do.its  = FALSE)$VaR
  }                                
}

## Test the VaR
library("GAS")
CC_pval <- DQ_pval <- NULL
for (j in 1:length(models)) { #iterate over model
  test <- GAS::BacktestVaR(data  = y_ots,
                           VaR   = VaR[,j],
                           alpha = alpha)
  
  CC_pval[j] <- test$LRcc[2]      
  DQ_pval[j] <- test$DQ$pvalue  
}
names(CC_pval) <- names(DQ_pval) <- c("GARCH", "MSGARCH")

print(CC_pval)
print(DQ_pval)

par(mfrow = c(1,1))
plot(y_ots, pch = 20, las = 1, xlab = "time index", ylab = "return [%]")
lines(VaR[,1], type = 'l', col = "red", lwd = 2)
lines(VaR[,2], type = 'l', col = "blue", lwd = 2)
legend("topleft", legend = c("GARCH", "MSGARCH"), col = c("red", "blue"), lwd = 2, cex = 2)
grid()


