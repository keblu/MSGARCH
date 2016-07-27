#' @export
print.MSGARCH_SPEC = function(spec){
  if(spec$K == 1){
    type = "Single-Regime"
    print(paste0("Specification Type: ", type))
    print(paste0("Specification Name: ",paste(spec$name,collapse=" ")))
    print(paste0("Number of parameters in variance model: ", paste(spec$n.params.vol,collapse=" ")))
    print(paste0("Number of parameters in distribution: ",  paste(spec$n.params - spec$n.params.vol,collapse=" ")))
    spec$theta0 = matrix(spec$theta0, ncol= length(spec$theta0))
    colnames(spec$theta0) = spec$label
    print(paste0("Default parameters:"))
    print(spec$theta0)
  } else{
    if(isTRUE(spec$is.mix)){
      type = "Mixture"
    } else {
      type = "Markov-Switching"
    } 
    if(isTRUE(spec$is.shape.ind)){
      type2 = " with Regime-Independent distribution"
    } else {
      type2 = ""
    }
    print(paste0("Specification Type: ", type, type2))
    print(paste0("Specification Name: ",paste(spec$name,collapse=" ")))
    print(paste0("Number of parameters in each variance model: ", paste(spec$n.params.vol,collapse=" ")))
    if(isTRUE(spec$is.shape.ind)){
      print(paste0("Number of parameters in distribution: ",  paste(spec$n.params[1] - spec$n.params.vol[1],collapse=" ")))
    } else {
      print(paste0("Number of parameters in each distribution: ",  paste(spec$n.params - spec$n.params.vol,collapse=" ")))
    }
    spec$theta0 = matrix(spec$theta0, ncol= length(spec$theta0))
    colnames(spec$theta0) = spec$label
    print(paste0("Default parameters:"))
    print(spec$theta0)
  }
}

#' @export
summary.MSGARCH_SPEC = function(spec){
  if(spec$K == 1){
    type = "Single-Regime"
    print(paste0("Specification Type: ", type))
    print(paste0("Specification Name: ",paste(spec$name,collapse=" ")))
    print(paste0("Number of parameters in variance model: ", paste(spec$n.params.vol,collapse=" ")))
    print(paste0("Number of parameters in distribution: ",  paste(spec$n.params - spec$n.params.vol,collapse=" ")))
    spec$theta0 = matrix(spec$theta0, ncol= length(spec$theta0))
    colnames(spec$theta0) = spec$label
    print(paste0("Default parameters:"))
    print(spec$theta0)
  } else{
    if(isTRUE(spec$is.mix)){
      type = "Mixture "
    } else {
      type = "Markov-Switching "
    } 
    if(isTRUE(spec$is.shape.ind)){
      type2 = " with Regime-Independent distribution"
    } else {
      type2 = ""
    }
    print(paste0("Specification Type: ", type, type2))
    print(paste0("Specification Name: ",paste(spec$name,collapse=" ")))
    print(paste0("Number of parameters in each variance model: ", paste(spec$n.params.vol,collapse=" ")))
    if(isTRUE(spec$is.shape.ind)){
      print(paste0("Number of parameters in distribution: ",  paste(spec$n.params[1] - spec$n.params.vol[1],collapse=" ")))
    } else {
      print(paste0("Number of parameters in each distribution: ",  paste(spec$n.params - spec$n.params.vol,collapse=" ")))
    }
    
    spec$theta0 = matrix(spec$theta0, ncol= length(spec$theta0))
    colnames(spec$theta0) = spec$label
    print(paste0("Default parameters:"))
    print(spec$theta0)
  }
}

#'@import ggplot2 reshape2
#'@export
plot.MSGARCH_RND = function(rnd){
  if(nrow(rnd$draws == 1)){
  df.m = data.frame(t(rnd$draws))
  }else{
  df.m = data.frame(rnd$draws)
  }
  df.m$label = 1:nrow(df.m)
  df.m <- reshape2::melt(df.m, id.vars = "label")
  ggplot(df.m,aes(x=value, group=variable)) + geom_density(alpha = 0.05, fill="gray10") + theme(legend.position="none") + ggtitle("Density of the simulated draws")
}

#'@export
summary.MSGARCH_RND= function(rnd){
  if(length(rnd$state) == 0){
    print(paste0("Total number of draws for each parameter estimate :", ncol(rnd$draws)))
    print(paste0("For more information set do.state to TRUE"))
  } else {
  print(paste0("Total number of draws for each parameter estimate :", ncol(rnd$state)))
  if(nrow(rnd$state) == 1){
    n.state = table(rnd$state)
    print("Percentage of draws made in state:")
    print(n.state)
  } else {
    n.state = list()
    for(i in 1:nrow(rnd$state)) {
      n.state[[i]] = table(rnd$state[i,])
    }
    n.state.2 = matrix(NA,nrow = length(n.state),ncol = length(n.state[[1]]))
    for(i in 1:length(n.state[[1]])){
      n.state.2[,i] = sapply(n.state, function (x) x[i])
    }
    n.state.mean = matrix(colMeans(n.state.2))
    rownames(n.state.mean) <- 1:nrow(n.state.mean)
    colnames(n.state.mean) = "State"
    n.state.sd = matrix(sqrt(diag(var(n.state.2/nrow(rnd$state)))))
    rownames(n.state.sd) <- 1:nrow(n.state.sd)
    colnames(n.state.sd) = "State"
    print("Average percentage of draws from state:")
    print(n.state.mean/ncol(rnd$state))
    print("Standard deviation of the percentage of draws from state:")
    print(n.state.sd)
  }
  }
}

#'@import zoo
#'@export
plot.MSGARCH_SIM = function(sim){
  sim$draws = zoo::zoo(t(sim$draws))
  plot(cumsum(sim$draws), plot.type = "single",ylab = "Cummulative draws",main = "Simulated draws")
  if(length(sim$state != 0)){
  sim$state = zoo::zoo(t(sim$state))
  plot(sim$state, plot.type = "single",ylab = "State",main = "Simulated state")
  }
}

#'@import zoo
#'@export
plot.MSGARCH_HT = function(ht){
  for(i in 1:dim(ht)[3]){
  tmp = zoo::zoo(ht[,,i])
  plot(tmp, plot.type = "single",ylab = "Volatility",main = paste0("Volatility of State ",i))
  }
}

plot.MSGARCH_PSTATE = function(Pstate){
  for(i in 1:dim(Pstate)[3]){
    tmp = zoo::zoo(Pstate[,,i])
    plot(tmp, plot.type = "single",ylab = "Probability",main = paste0("Probability to be in State ",i))
  }
}