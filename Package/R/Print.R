#' @export
print.MSGARCH_SPEC <- function(x, ...) {
  spec <- x
  spec <- f_check_spec(spec)
  if (spec$K == 1L) {
    type <- "Single-regime"
    cat(paste0("Specification type: ", type,"\n"))
    cat(paste0("Specification name: ", paste(spec$name, collapse = " "),"\n"))
    cat(paste0("Number of parameters in variance model: ", paste(spec$n.params.vol, collapse = " "),"\n"))
    cat(paste0("Number of parameters in distribution: ", paste(spec$n.params - spec$n.params.vol, collapse = " "),"\n"))
    cat(paste("------------------------------------------\n"))
    out = list(type = type, spec.name = paste(spec$name, collapse = " "),
               n.params.var = paste(spec$n.params.vol, collapse = " "),
               n.params.dist = paste(spec$n.params - spec$n.params.vol, collapse = " "),
               prior.mean = spec$prior.mean, prior.sd = spec$prior.sd,
               fixed.pars = spec$fixed.pars,
               regime.const.pars = spec$regime.const.pars)
  } else {
    if (isTRUE(spec$is.mix)) {
      type <- "Mixture"
    } else {
      type <- "Markov-switching"
    }
    cat(paste0("Specification type: ", type,"\n"))
    cat(paste0("Specification name: ", paste(spec$name, collapse = " "),"\n"))
    cat(paste0("Number of parameters in each variance model: ", paste(spec$n.params.vol, collapse = " "),"\n"))
    cat(paste0("Number of parameters in each distribution: ", paste(spec$n.params - spec$n.params.vol, collapse = " "),"\n"))
    
    if (isTRUE(spec$fixed.pars.bool)) {
      cat(paste("------------------------------------------\n"))
      cat(paste0("Fixed parameters:","\n"))
      cat(unlist(spec$fixed.pars),"\n")
    } else {
      cat(paste("------------------------------------------\n"))
      cat(paste0("Fixed parameters:","\n"))
      cat("None\n")
    }
    
    if (isTRUE(spec$regime.const.pars.bool)) {
      cat(paste("------------------------------------------\n"))
      cat(paste0("Across regime constrained parameters:","\n"))
      cat(unlist(spec$regime.const.pars),"\n")
    } else {
      cat(paste("------------------------------------------\n"))
      cat(paste0("Across regime constrained parameters:","\n"))
      cat("None\n")
    }
    
    cat(paste("------------------------------------------\n"))
    out = list(type = paste0("Specification type: ", type),
               spec.name = paste0("Specification name: ",
                                  paste(spec$name, collapse = " ")),
               n.params.var = paste(spec$n.params.vol, collapse = " "),
               n.params.dist = paste(spec$n.params[1] - spec$n.params.vol[1], collapse = " "),
               prior.mean = spec$prior.mean, prior.sd = spec$prior.sd)
  }
  return(invisible(out))
}

#' @export
summary.MSGARCH_SPEC <- function(object, ...) {
  out <- print(object)
  return(invisible(out))
}

#' @export
summary.MSGARCH_ML_FIT <- function(object, ...) {
  object$spec <- f_check_spec(object$spec)
  spec <- summary(object$spec)
  cat(paste0("Fitted parameters:","\n"))
  if (is.null(object$Inference)) {
    estimate <- as.matrix(object$par)
    colnames(estimate) <- "Estimate"
    estimate.print <- estimate
  } else {
    estimate <- object$Inference$MatCoef
    estimate.print <- data.frame(estimate)
    colnames(estimate.print) <- colnames(estimate)
    estimate.print[, 1:3] <- round(estimate.print[1:3], 4)
    p.val.below <- estimate.print[, 4] < 1e-16
    p.val.below[is.na(p.val.below)] <- FALSE
    estimate.print[, 4] <- format(estimate.print[, 4], scientific = TRUE, digits = 4)
    estimate.print[p.val.below, 4] <- "<1e-16"
  }
  print(estimate.print)
  cat(paste("------------------------------------------\n"))
  if (object$spec$K > 1L) {
    if (!object$spec$is.mix) {
      cat(paste0("Transition matrix:","\n"))
      trans.mat <- TransMat(object)
      print(round(trans.mat, 4))
      cat(paste("------------------------------------------\n"))
      stable.prob <- t(TransMat(object, nahead = 10000)[1, , drop = FALSE])
    } else {
      stable.prob <- t(TransMat(object))
    }
    cat(paste0("Stable probabilities:","\n"))
    stable.prob <-  as.vector(stable.prob)
    names(stable.prob) <- paste0("State ", 1:object$spec$K)
    print(round(stable.prob, 4))
    cat(paste("------------------------------------------\n"))
    #cat("Unconditional volatility:\n")
    #unc.vol = UncVol(object = object)
    #cat("In each regime:\n")
    #print(round(unc.vol$SR, 4))
    #cat("Overall process:\n")
    #cat(paste0(" ", round(unc.vol$MS, 4),"\n"))
  } else {
    #cat("Unconditional volatility:\n")
    #unc.vol = UncVol(object = object)
    #print(round(unc.vol$SR, 4))
  }
  cat(paste0("LL: ", round(object$loglik, 4), "\n"))
  aic <- stats::AIC(object)
  cat(paste0("AIC: ", round(aic, 4), "\n"))
  bic <- stats::BIC(object)
  cat(paste0("BIC: ", round(bic, 4)))
  cat(paste("\n------------------------------------------\n"))
  if (isTRUE(object$spec$is.mix)) {
    out <- list(spec = spec, estimate = estimate,
                stable.prob = stable.prob,
                LL = object$loglik, AIC = aic, BIC = bic)
  } else {
    if (object$spec$K > 1) {
      out <- list(spec = spec, estimate = estimate,
                  trans.mat = trans.mat,
                  stable.prob = stable.prob,
                  LL = object$loglik, AIC = aic, BIC = bic)
    } else {
      out <- list(spec = spec, estimate = estimate,
                  LL = object$loglik, AIC = aic, BIC = bic)
    }
  }
  return(invisible(out))
}

print.MSGARCH_ML_FIT <- function(x, ...) {
  out <- summary(x, ...)
  return(invisible(out))
}

#' @export
summary.MSGARCH_MCMC_FIT <- function(object, ...) {
  object$spec <- f_check_spec(object$spec)
  spec <- summary(object$spec)
  cat(paste0("Posterior sample (size: ", nrow(object$par), ")","\n"))
  summ <- summary(object$par)$stat[colnames(object$ctr$par0),]
  summ <- cbind(summ,summ[,3]^2 / summ[,4]^2)
  colnames(summ)[c(3,4,5)] <- c("SE","TSSE","RNE")
  print(round(summ, digits = 4))
  cat(paste("------------------------------------------\n"))
  if (object$spec$K > 1) {
    if (!object$spec$is.mix) {
      cat(paste0("Posterior mean transition matrix:","\n"))
      post.trans.mat = TransMat(object = object$spec, par = colMeans(object$par))
      print(round(post.trans.mat, 4))
      cat(paste("------------------------------------------\n"))
      post.stable.prob <- t(TransMat(object = object$spec,
                                     par = colMeans(object$par), nahead = 10000)[1, , drop = FALSE])
    } else {
      post.stable.prob <- t(TransMat(object = object$spec, par =  colMeans(object$par)))
    }
    cat(paste0("Posterior mean stable probabilities:","\n"))
    post.stable.prob <- as.vector(post.stable.prob)
    names(post.stable.prob) <- paste0("State ", 1:object$spec$K)
    print(round(post.stable.prob, 4))
    # cat(paste("------------------------------------------\n"))
    #  cat("Posterior mean unconditional volatility:\n")
    #  post.unc.vol = UncVol(object = object$spec, par = post.mean)
    #  cat("In each regime:\n")
    #  print(round(post.unc.vol$SR, 4))
    #  cat("Overall process:\n")
    #  cat(paste0(" ", round(post.unc.vol$MS, 4),"\n"))
  } else {
    # cat("Posterior mean unconditional volatility:\n")
    #  post.unc.vol = UncVol(object = object, par = post.mean)
    #  print(round(post.unc.vol$SR, 4))
  }
  cat(paste("------------------------------------------\n"))
  cat(paste0("Acceptance rate MCMC sampler: ", round(100 * object$accept, 1), "%\n"))
  cat(paste0("nmcmc: ", object$ctr$nmcmc, "\n"))
  cat(paste0("nburn: ", object$ctr$nburn, "\n"))
  cat(paste0("nthin: ", object$ctr$nthin, "\n"))
  cat(paste("------------------------------------------\n"))
  dic <- DIC(object)$DIC
  cat(paste0("DIC: ", round(dic, 4)))
  cat(paste("\n------------------------------------------\n"))
  if (isTRUE(object$spec$is.mix)) {
    out <- list(spec = spec, summary = summ, post.stable.prob = post.stable.prob,
                accept.rate = object$accept,
                DIC = dic)
  } else {
    if (object$spec$K > 1) {
      out <- list(spec = spec, summary = summ, post.trans.mat = post.trans.mat,
                  post.stable.prob = post.stable.prob,
                  accept.rate = object$accept, DIC = dic)
    } else {
      out <- list(spec = spec, summary = summ, accept.rate = object$accept, DIC = dic)
    }
  }
  return(invisible(out))
}

print.MSGARCH_MCMC_FIT <- function(x, ...) {
  out <- summary(x, ...)
  return(invisible(out))
}

#' @export
print.MSGARCH_FORECAST <- function(x, ...){
  x$draw = unclass(x$draw)
  x$vol = unclass(x$vol)
  print(unclass(x))
}
#' @export
print.MSGARCH_CONDVOL <- function(x, ...){
  if(any(class(x) == "zoo")){
    class(x) = "zoo"
    print(x)
  } else  if(any(class(x) == "ts")){
    class(x) = "ts"
    print(x)
  } else {
    print(unclass(x))
  }
}
#' @export
print.MSGARCH_RISK <- function(x, ...){
  if(any(class(x) == "zoo")){
    class(x) = "zoo"
    print(x)
  } else  if(any(class(x) == "ts")){
    class(x) = "ts"
    print(x)
  } else {
    print(unclass(x))
  }
}
#' @export
print.MSGARCH_PSTATE <- function(x, ...){
  print(unclass(x))
}
#' @export
print.MSGARCH_SIM <- function(x, ...){
  if(any(class(x) == "zoo")){
    class(x) = "zoo"
    print(x)
  } else  if(any(class(x) == "ts")){
    class(x) = "ts"
    print(x)
  } else {
    print(unclass(x))
  }
}
