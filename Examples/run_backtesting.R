rm(list = ls())
library("MSGARCH")
data("SMI", package = "MSGARCH")

## Model specification
models <- list(GARCH   = CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                                    distribution.spec = list(c("norm")),
                                    switch.spec = list(K = 1)),
               MSGARCH = CreateSpec(variance.spec = list(model = c("gjrGARCH")),
                                    distribution.spec = list(c("norm")),
                                    switch.spec = list(K = 2)))

N_out_of_sample <- 10  # number of out-of-sample evaluation 
N_sample        <- 2000 # fit sample size
alpha           <- 0.05 # risk Level 
k_update        <- 1    # frequence at which we update the model (reestimation)

## Initialization 
VaR_OOS   <- matrix(NA, nrow = N_out_of_sample, ncol = length(models))
y_OOS     <- matrix(NA, nrow = N_out_of_sample, ncol = 1)
model_fit <- vector(mode = "list", length = length(models))

## Backtest 
for (i in 1:N_out_of_sample) { # iterate over out-of-sample time
  cat("Backtest. Iteration: ", i, "\n")
  y_fit    <- SMI[i:(N_sample + i - 1)] # training sample
  y_OOS[i] <- SMI[N_sample + i] # store out-of-sample data
  for (j in 1:length(models)) { # iterate over model
    # do we update the model ?
    if (k_update == 1 || i %% k_update == 1) {
      model_fit[[j]] <- FitML(spec = models[[j]], data = y_fit) # fit the model
    }
    
    # calculate VaR 1-step ahead
    VaR_OOS[i,j] <- Risk(model_fit[[j]],
                         n.ahead = 1,
                         alpha   = alpha,
                         do.es   = FALSE,
                         do.its  = FALSE)$VaR
  }
}

## Test the VaR
library("GAS")
CCtest_pvalue <- NULL
for (j in 1:length(models)) { #iterate over model
  test <- GAS::BacktestVaR(data  = y_OOS,
                           VaR   = VaR_OOS[,j],
                           alpha = alpha)
  
  CCtest_pvalue[j] <- test$LRcc[j] #P-Value 
}
