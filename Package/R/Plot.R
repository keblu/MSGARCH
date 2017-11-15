#'@import zoo
#'@importFrom graphics matplot
#'@export
plot.MSGARCH_SIM <- function(x, main = NULL, type = NULL, xlab = NULL, ylab = NULL, ...) {
  draw <- x$draw
  if (is.null(ylab)) {
    ylab <- "Draws"
  }
  
  if (is.null(xlab)) {
    xlab <- "Index"
  }
  
  if (is.null(main)) {
    main <- "Simulated draws"
  }
  
  if (is.null(type)) {
    type <- "l"
  }
  #graphics::matplot(draw, type = type, ylab = ylab, xlab = xlab, main = main)
  fanplot::fan0(t(draw),type = "interval", ylab = ylab, xlab = xlab, main = main, xlim = c(1, nrow(draw)), ylim = range(draw))
}

#'@import zoo fanplot
#'@importFrom graphics plot axis.Date
#'@export
plot.MSGARCH_PSTATE <- function(x, type.prob = c("Smoothed"), date = NULL, main = NULL,
                                xlab = NULL, ylab = NULL, med.col = NULL, med.ln = NULL, ...) {
  
  if (f_match(type.prob, "SMO")) {
    state    <- x$SmoothProb[1:dim(x$SmoothProb)[[1]] - 1, , , drop = FALSE]
    main.def <- "Smoothed"
  } else if (f_match(type.prob, "FIL")) {
    state    <- x$FiltProb
    main.def <- "Filtered"
  } else if (f_match(type.prob, "PRE")) {
    state    <- x$PredProb[1:dim(x$PredProb)[[1]] - 1, , , drop = FALSE]
    main.def <- "Predictive"
  } else if (f_match(type.prob, "VIT")) {
    state <- x$Viterbi
  } else {
    stop("Please choose a valid type.prob. Valid type.prob are
         c('smoothed','filtered','predictive','viterbi')")
  }
  
  if (is.null(main)) {
    if (f_match(type.prob, "VIT")) {
      generate_main = function(x) {
        return("Most probable path")
      }
    } else {
      generate_main = function(x) {
        paste0(main.def, " for State ", x)
      }
    }
  } else {
    generate_main = function(x) {
      return(main)
    }
  }
  if (is.null(date)) {
    date = dimnames(state)[[1]]
    date = tryCatch({as.Date(date)}, error = function(e) {
      return(as.numeric(date))})
  }
  
  if (is.null(ylab)) {
    if (f_match(type.prob, "VIT")) {
      ylab = "State"
    } else {
      ylab = "Probability"
    }
  }
  
  if (is.null(xlab)) {
    xlab = "Date"
  }
  
  if (is.null(med.ln)) {
    med.ln = TRUE
  }
  
  if (is.null(med.col)) {
    med.col = "blue"
  }
  if (dim(state)[2L] > 1L) {
    plot.func = function(x, generated.main, ...) {
      graphics::axis.Date(side = 1L, x = zoo::index(x), ...)
      fanplot::fan0(x, type = "interval", med.col = med.col, med.ln = med.ln,
                    xlab = xlab, ylab = ylab, xlim = range(zoo::index(x)),
                    ylim = range(x), main = generated.main,
                    ...)
    }
  } else {
    plot.func <- function(x, generated.main, ...) {
    plot(x, plot.type = "single", ylab = ylab, xlab = xlab,
           main = generated.main, ...)
    }
  }
  if (!f_match(type.prob, "VIT")) {
    for (i in 1:dim(state)[3L]) {
      if (ncol(state) == 1L) {
        tmp <- zoo::zoo(state[, , i], order.by = date)
      } else {
        tmp <- zoo::zoo(t(state[, , i]), order.by = date)
      }
      plot.func(tmp, generated.main = generate_main(i), ...)
    }
  } else {
    if (ncol(state) == 1L) {
      tmp <- zoo::zoo(state, order.by = date)
    } else {
      tmp <- zoo::zoo(t(state), order.by = date)
    }
    plot.func(tmp, generated.main = generate_main(i), ...)
  }
  }

plot.MSGARCH_FORECAST <- function(x, ...){
  plot.MSGARCH_CONDVOL(x$vol)
  if(!is.null(x$draw)){
     plot.MSGARCH_SIM(x)
  }
}
#'@import zoo fanplot
#'@importFrom graphics plot axis.Date
#'@export
plot.MSGARCH_CONDVOL<- function(x, date = NULL, main = NULL, xlab = NULL,
                            ylab = NULL, med.col = NULL, med.ln = NULL, ...) {
  
  ht <- x
  
  if (is.null(main)) {
    generate_main = function(x) {
      paste0("Conditional volatility")
    }
  } else {
    generate_main = function(x) {
      return(main)
    }
  }
  if (is.null(date)) {
    if(any(class(ht) == "zoo") || any(class(ht) == "ts")){
      date = zoo::index(ht)
    } else {
    date = 1:length(ht)
    }
  }
  
  if (is.null(ylab)) {
    ylab = "Volatility"
  }
  
  if (is.null(xlab)) {
    xlab = "Date"
  }
  
  if (is.null(med.ln)) {
    med.ln = TRUE
  }
  
  if (is.null(med.col)) {
    med.col = "blue"
  }
  
  plot.func <- function(x, generated.main, ...) {
    graphics::plot(x, plot.type = "single", ylab = ylab, xlab = xlab, main = generated.main, ...)
  }
  tmp <- zoo::zoo(ht, order.by = date)
  plot.func(tmp, generated.main = generate_main(x), ...)
}

#'@importFrom graphics legend plot
#'@importFrom grDevices rainbow
#'@export
plot.MSGARCH_RISK <- function(x, date = NULL, main.VaR = NULL, main.ES = NULL,
                              xlab = NULL, ylab = NULL, col = NULL, plot.type = NULL, ...) {
  risk <- x
  
  if (is.null(date)) {
    if(any(class(risk$VaR) == "zoo") || any(class(risk$VaR) == "ts")){
      date = zoo::index(risk$VaR)
    } else {
      date = 1:length(risk$VaR)
    }
  }
  if (is.null(main.VaR)) {
    main.VaR = "Value-At-Risk"
  }
  
  if (is.null(main.ES)) {
    main.ES = "Expected-shortfall"
  }
  
  if (is.null(col)) {
    col <- rainbow(ncol(risk$VaR))
  }
  
  if (is.null(plot.type)) {
    plot.type = "single"
  }
  
  if (is.null(ylab)) {
    ylab = "Return"
  }
  
  if (is.null(xlab)) {
    xlab = "Index"
  }
  
  tmp.VaR <- zoo::zoo(risk$VaR, order.by = date)
  graphics::plot(tmp.VaR, plot.type = plot.type, col = col, ylab = ylab, xlab = xlab, main = main.VaR, ...)
  graphics::legend("bottomright", legend = colnames(risk$VaR), col = col, lty = 1, ...)
  if (!is.null(risk$ES)) {
    tmp.ES <- zoo::zoo(risk$ES, order.by = date)
    graphics::plot(tmp.ES, plot.type = plot.type, col = col, ylab = ylab, xlab = xlab, main = main.ES, ...)
    graphics::legend("bottomright", legend = colnames(risk$ES), col = col, lty = 1, ...)
  }
}
