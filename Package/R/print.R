#' @export
print.MSGARCH_SPEC <- function(x, ...) {
  spec <- x
  if (spec$K == 1) {
    type <- "Single-Regime"
    print(paste0("Specification Type: ", type))
    print(paste0("Specification Name: ", paste(spec$name, collapse = " ")))
    print(paste0("Number of parameters in variance model: ", paste(spec$n.params.vol,
                 collapse = " ")))
    print(paste0("Number of parameters in distribution: ", paste(spec$n.params -
                  spec$n.params.vol, collapse = " ")))
    spec$theta0 <- matrix(spec$theta0, ncol = length(spec$theta0))
    colnames(spec$theta0) <- spec$label
    print(paste0("Default parameters:"))
    print(spec$theta0)
  } else {
    if (isTRUE(spec$is.mix)) {
      type <- "Mixture"
    } else {
      type <- "Markov-Switching"
    }
    if (isTRUE(spec$is.shape.ind)) {
      type2 <- " with Regime-Independent distribution"
    } else {
      type2 <- ""
    }
    print(paste0("Specification Type: ", type, type2))
    print(paste0("Specification Name: ", paste(spec$name, collapse = " ")))
    print(paste0("Number of parameters in each variance model: ", paste(spec$n.params.vol,
                 collapse = " ")))
    if (isTRUE(spec$is.shape.ind)) {
      print(paste0("Number of parameters in distribution: ", paste(spec$n.params[1] -
                   spec$n.params.vol[1], collapse = " ")))
    } else {
      print(paste0("Number of parameters in each distribution: ", paste(spec$n.params -
                    spec$n.params.vol, collapse = " ")))
    }
    spec$theta0 <- matrix(spec$theta0, ncol = length(spec$theta0))
    colnames(spec$theta0) <- spec$label
    print(paste0("Default parameters:"))
    print(spec$theta0)
  }
}
#' @export
summary.MSGARCH_MLE_FIT <- function(object, ...) {
  print(object$spec)
  print(paste0("DEoptim initialization: ", object$is.init))
  print(paste0("Fitted Parameters:"))
  print(object$theta)
  if(object$spec$K > 1){
    if(!object$spec$is.mix){
      print(paste0("Transition matrix:"))
      print(transmat(object))
      stable_prob = transmat(object, n = 100) %*% matrix(rep(1/object$spec$K,object$spec$K), ncol = 1)
    } else {
      stable_prob = t(transmat(object))
    }
    print(paste0("Stable probabilities:"))
    rownames(stable_prob) = paste0("State ", 1:object$spec$K)
    colnames(stable_prob) = "Stable probabilities"
    print(stable_prob)
  }
  print("Unconditional volatility:")
  print(MSGARCH::unc.vol(object = object))
  print(paste0("Log-kernel: ", object$log_kernel))
  print(paste0("AIC: ",MSGARCH::AIC(object)))
  print(paste0("BIC: ",MSGARCH::BIC(object)))
}
#' @export
summary.MSGARCH_BAY_FIT <- function(object, ...) {
  print(object$spec)
  print(paste0("Bayesian posterior mean:"))
  theta_mean = colMeans(object$theta)
  print(theta_mean)
  print(paste0("Posterior variance-covariance matrix"))
  print(var(object$theta))
  if(object$spec$K > 1){
    if(!object$spec$is.mix){
      print(paste0("Posterior mean transition matrix:"))
      print(transmat(object = object$spec,theta = theta_mean))
      stable_prob = transmat(object = object$spec, theta = theta_mean, n = 100) %*% matrix(rep(1/object$spec$K,object$spec$K), ncol = 1)
    } else {
      stable_prob = t(transmat(object = object$spec, theta = theta_mean))
    }
    print(paste0("Posterior mean stable probabilities:"))
    rownames(stable_prob) = paste0("State ", 1:object$spec$K)
    colnames(stable_prob) = "Stable probabilities"
    print(stable_prob)
  }
  print("Posterior mean unconditional volatility:")
  print(MSGARCH::unc.vol(object = object$spec, theta = theta_mean))
  print(paste0("Acceptance rate: ",object$accept))
  print(paste0("AIC: ",MSGARCH::AIC(object)))
  print(paste0("BIC: ",MSGARCH::BIC(object)))
  print(paste0("DIC: ",MSGARCH::DIC(object)$DIC))
}

#' @export
summary.MSGARCH_SPEC <- function(object, ...) {
  spec <- object
  if (spec$K == 1) {
    type <- "Single-Regime"
    print(paste0("Specification Type: ", type))
    print(paste0("Specification Name: ", paste(spec$name, collapse = " ")))
    print(paste0("Number of parameters in variance model: ", paste(spec$n.params.vol,
                 collapse = " ")))
    print(paste0("Number of parameters in distribution: ", paste(spec$n.params -
                 spec$n.params.vol, collapse = " ")))
    spec$theta0 <- matrix(spec$theta0, ncol = length(spec$theta0))
    colnames(spec$theta0) <- spec$label
    print(paste0("Default parameters:"))
    print(spec$theta0)
  } else {
    if (isTRUE(spec$is.mix)) {
      type <- "Mixture "
    } else {
      type <- "Markov-Switching "
    }
    if (isTRUE(spec$is.shape.ind)) {
      type2 <- " with Regime-Independent distribution"
    } else {
      type2 <- ""
    }
    print(paste0("Specification Type: ", type, type2))
    print(paste0("Specification Name: ", paste(spec$name, collapse = " ")))
    print(paste0("Number of parameters in each variance model: ", paste(spec$n.params.vol,
                  collapse = " ")))
    if (isTRUE(spec$is.shape.ind)) {
      print(paste0("Number of parameters in distribution: ", paste(spec$n.params[1] -
                   spec$n.params.vol[1], collapse = " ")))
    } else {
      print(paste0("Number of parameters in each distribution: ", paste(spec$n.params -
                   spec$n.params.vol, collapse = " ")))
    }
    spec$theta0 <- matrix(spec$theta0, ncol = length(spec$theta0))
    colnames(spec$theta0) <- spec$label
    print(paste0("Default parameters:"))
    print(spec$theta0)
  }
}

#'@import zoo
#'@importFrom graphics matplot
#'@export
plot.MSGARCH_SIM <- function(x, ...) {
  sim <- x
  matplot(t(sim$draws), type = "l", ylab = "draws", main = "Simulated draws")
  states = as.numeric(rownames(table(sim$state)))
  for(i in 1:length(states)){
    v = readline(prompt = "Pause. Press <Enter> to continue...")
     percent_sim <- zoo::zoo(colSums(sim$state == states[i])/nrow(sim$state))
     plot(percent_sim, plot.type = "single", ylab = "%", main = paste0("Percentage of simulation done in state ", i))
  }
}

#'@import zoo
#'@importFrom graphics plot
#'@export
plot.MSGARCH_HT <- function(x, ...) {
  ht <- x
  for (i in 1:dim(ht)[3]) {
    tmp <- zoo::zoo(ht[, , i])
    plot(tmp, plot.type = "single", ylab = "Volatility",
         main = paste0("Volatility of State ", i))
    v = readline(prompt = "Pause. Press <Enter> to continue...")
  }
}

#'@import zoo
#'@importFrom graphics plot
#'@export
plot.MSGARCH_PSTATE <- function(x, ...) {
  Pstate <- x
  for (i in 1:dim(Pstate)[3]) {
    tmp <- zoo::zoo(Pstate[, , i])
    plot(tmp, plot.type = "single", ylab = "Probability",
         main = paste0("Probability to be in State ", i))
    v = readline(prompt = "Pause. Press <Enter> to continue...")
  }
}

#'@import zoo
#'@importFrom graphics plot
#'@export
plot.MSGARCH_HT <- function(x, ...) {
  ht <- x
  if (is.matrix(ht)) {
    tmp <- zoo::zoo(t(ht))
    plot(tmp, plot.type = "single", ylab = "Volatility", main = paste0("Conditional volatility"))
  } else {
    for (i in 1:dim(ht)[3]) {
      tmp <- zoo::zoo(ht[, , i])
      plot(tmp, plot.type = "single", ylab = "Volatility",
           main = paste0("Conditional volatility of state ", i))
      v = readline(prompt = "Pause. Press <Enter> to continue...")
    }
  }
}

#'@import ggplot2 reshape2
#'@export
plot.MSGARCH_CDF <- function(x, ...) {
  cdf <- x
  variable <- x <- value <- Var2 <- NULL
  if (isTRUE(cdf$do.its)) {
    stop("no plot method for do.its option")
  }
  ind <- sort(cdf$x, index.return = TRUE)
  if (nrow(cdf$cdf) > 1) {
    df.m <- data.frame(cdf$cdf[, ind$ix])
    df.m <- t(df.m)
    df.m <- reshape2::melt(df.m, id.vars = variable)
    df.m$x <- rep(ind$x, nrow(cdf$cdf))
    ggplot(df.m, aes(x = x, y = value, group = Var2)) + geom_line(size = 1.5) +
      theme(legend.position = "none") + ggtitle("Cummulative of x") + ylab("Cummulative")
  } else {
    df.m <- data.frame(cdf$cdf[, ind$ix])
    df.m$x <- rep(ind$x, nrow(cdf$cdf))
    ggplot(df.m, aes(x = x, y = cdf$cdf[, ind$ix])) + geom_line(size = 1.5) +
      theme(legend.position = "none") + ggtitle("Cummulative of x") + ylab("Cummulative")
  }
}

#'@import ggplot2 reshape2
#'@export
plot.MSGARCH_PDF <- function(x, ...) {
  pdf <- x
  if (isTRUE(pdf$do.its)) {
    stop("no plot method for do.its option")
  }
  variable <- x <- value <- Var2 <- NULL
  ind <- sort(pdf$x, index.return = TRUE)
  if (nrow(pdf$pdf) > 1) {
    df.m <- data.frame(pdf$pdf[, ind$ix])
    df.m <- t(df.m)
    df.m <- reshape2::melt(df.m, id.vars = variable)
    df.m$x <- rep(ind$x, nrow(pdf$pdf))
    ggplot(df.m, aes(x = x, y = value, group = Var2)) + geom_line(size = 1.5) +
      theme(legend.position = "none") + ggtitle("Density") + ylab("Density")
  } else {
    df.m <- data.frame(pdf$pdf[, ind$ix])
    df.m$x <- rep(ind$x, nrow(pdf$pdf))
    ggplot(df.m, aes(x = x, y = pdf$pdf[, ind$ix])) + geom_line(size = 1.5) +
      theme(legend.position = "none") + ggtitle("Density") + ylab("Density")
  }
}

#'@import ggplot2
#'@export
plot.MSGARCH_PIT <- function(x, ...) {
  pit <- x
  pit = pit$pit[2:length(pit$pit)]
  df.m <- data.frame(pit)
  ggplot(data = df.m, aes(x = pit)) + geom_histogram(bins = 100, binwidth = 0.01) +
    theme(legend.position = "none") + ggtitle("Probability integral transform") +
    xlab("PIT")
}

#'@import ggplot2
#'@export
plot.MSGARCH_PRED <- function(x, ...) {
  pred <- x
  if (isTRUE(pred$do.its)) {
    stop("no plot method for do.its option")
  }
  x <- NULL
  ind <- sort(pred$x, index.return = TRUE)
  df.m <- data.frame(pred$pred[ind$ix])
  df.m$x <- ind$x
  ggplot(df.m, aes(x = x, y = pred$pred[ind$ix])) + geom_line(size = 1.5) + theme(legend.position = "none") +
    ggtitle("Predictive") + ylab("Density")
}

#'@import ggplot2
#'@importFrom graphics legend plot
#'@importFrom grDevices rainbow
#'@export
plot.MSGARCH_RISK <- function(x, ...) {
  risk <- x
  ts_rainbow <- rainbow(ncol(risk$VaR))
  plot(zoo::zoo(risk$VaR), plot.type = "single", col = ts_rainbow, ylab = "Return",
    xlab = "T", main = paste0("Value-At-Risk"))
  legend("bottomright", legend = colnames(risk$VaR), col = ts_rainbow, lty = 1)
  if (!is.null(risk$ES)) {
    v = readline(prompt = "Pause. Press <Enter> to continue...")
    plot(zoo::zoo(risk$ES), plot.type = "single", col = ts_rainbow, ylab = "Return",
      xlab = "T", main = paste0("Expected-shortfall"))
    legend("bottomright", legend = colnames(risk$ES), lty = 1, col = ts_rainbow)
  }
}
