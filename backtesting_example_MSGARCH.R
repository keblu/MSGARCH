require(MSGARCH)
require(GAS)
data(AMZN)
y = AMZN

N = length(AMZN)

# Model specification
model = list(GARCH  = MSGARCH::create.spec(model = "sGARCH",distribution = "norm"),
             MSGARCH  = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm")))


N_out_of_sample = 100 # number of out-of-sample evaluation 
N_sample = 2000 #training sample size
level_risk = 0.95 #Risk Level 
k_update = 1 # frequence at which we update the model (reestimation)

# Matrices initialization 
VaR_OOS = matrix(NA,nrow = N_out_of_sample,ncol = length(model))
y_OOS = matrix(NA,nrow = N_out_of_sample,ncol = 1)
model_fit = vector(mode = "list", length = length(model))

#Calculation 
for(i in 1:N_out_of_sample){ #iterate over out-of-sample time
  y_train = y[i:(N_sample+i-1)] #training sample
  y_OOS[i] = y[N_sample+i] # Store out-of-sample data
  for(j in 1:length(model)){ #iterate over model
    #Do we update the model ?
    if(k_update == 1 || i %% k_update == 1){
      model_fit[[j]] = MSGARCH::fit.mle(spec = model[[j]], y = y_train) #Fit the model
    }
    #Calculate VaR 1-step ahead (for multistep ahead simulations MSGARCH::simahead must be use)
    VaR_OOS[i,j] = MSGARCH::risk(object = model_fit[[j]]$spec,
                                 theta = model_fit[[j]]$theta,
                                 y = y_train,
                                 level = level_risk,
                                 ES = FALSE,
                                 do.its = FALSE)$VaR
  }
}

#test the VaR
CCtest_pvalue = NULL
for(j in 1:length(model)){ #iterate over model
  test = GAS::BacktestVaR(data = y_OOS,
                                VaR = VaR_OOS[,j],
                                alpha = 1 - level_risk, 
                                Lags = 4)
  
  CCtest_pvalue[j] = test$LRcc[2] #P-Value 
}