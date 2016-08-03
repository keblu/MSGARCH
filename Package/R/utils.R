f.error = function(message) {
  cat(paste0("try/catch ", message, "\n"))
  return(FALSE)
}


f.process.ctr = function(ctr = list()) {
  con = list(theta0 = NULL, do.init = TRUE, N.mcmc = 5000, N.burn = 1000, N.thin = 10, NP = 200, itermax = 200)
  con[names(ctr)] = ctr
  return(con)
}

f.check.y = function(y){
  if(is.null(y)){
    stop("y is NULL")
  }
  if(!is.numeric(y)){
    stop("y must be numeric")
  }
  if(all(is.nan(y))){
    stop("nan dectected in y")
  }
  if(!is.null(dim(y))){
    if(any(dim(y) == 1)){
      y = as.vector(as.matrix(y))
    } else {
      stop("y is not a vector")
    }
  }
  y = as.vector(y)
  y = as.matrix(y)
  return(y)
}

f.check.theta = function(spec, theta){
  
  
  if(is.null(theta)){
    stop("theta is NULL")
  }
  
  if(!is.numeric(theta)){
    stop("theta must be a numeric")
  }
  
  if(all(is.nan(theta))){
    stop("nan dectected in theta")
  }
  
  len.theta = length(spec$theta0)
  if (is.vector(theta)) {
    theta = matrix(theta, nrow = 1)
  }
  if(is.data.frame(theta)){
   theta =  data.matrix(theta)
  }
  
  if(dim(theta)[2] != len.theta){
    stop(paste0("Each parameter estimate in theta must be of length ",len.theta))
  }
  
  if (isTRUE(spec$is.shape.ind)) {
    theta = spec$func$f.do.shape.ind(theta)
  }
  
  if (isTRUE(spec$is.mix)) {
    theta = spec$func$f.do.mix(theta)
  }
  return(theta)
}
