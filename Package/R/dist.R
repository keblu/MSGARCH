#' @export
ddist = function(y, dist = "norm", shape = 100, skew = 1, log = FALSE){
  valid.distribution <- c("norm", "std", "ged")
  if (!any(dist == valid.distribution)) {
    stop(paste0("error: The distribution ", dist, "
                does not appear to be a valid choice."))
  }
  if(dist == "norm"){
    ind = 1
  } else{
    ind = 0
  }
    if (skew == 1) {
      skew_lab <- "sym"
    } else {
      skew_lab <- "skew"
    }
  theta = c(shape, skew)
  dist_func = new(Class = getFromNamespace(paste0(dist,"_",skew_lab),pos = "package:MSGARCH"))
  dist_func$load_theta(theta,ind)
  out = dist_func$f_pdf(y)
  if(isTRUE(log)){
   out =  log(out)
  }
  return(out)
}

#' @export
rdist = function(n, dist = "norm", shape = 100, skew = 1){
  valid.distribution <- c("norm", "std", "ged")
  if (!any(dist == valid.distribution)) {
    stop(paste0("error: The distribution ", dist, "
                does not appear to be a valid choice."))
  }
  if(dist == "norm"){
    ind = 1
  } else{
    ind = 0
  }
  if (skew == 1) {
    skew_lab <- "sym"
  } else {
    skew_lab <- "skew"
  }
  theta = c(shape, skew)
  dist_func = new(Class = getFromNamespace(paste0(dist,"_",skew_lab),pos = "package:MSGARCH"))
  dist_func$load_theta(theta,ind)
  out = dist_func$f_rnd(n)
  return(out)
}

#' @export
qdist = function(y, dist = "norm", shape = 100, skew = 1){
  valid.distribution <- c("norm", "std", "ged")
  if (!any(dist == valid.distribution)) {
    stop(paste0("error: The distribution ", dist, "
                does not appear to be a valid choice."))
  }
  if(dist == "norm"){
    ind = 1
  } else{
    ind = 0
  }
  if (skew == 1) {
    skew_lab <- "sym"
  } else {
    skew_lab <- "skew"
  }
  theta = c(shape, skew)
  dist_func = new(Class = getFromNamespace(paste0(dist,"_",skew_lab),pos = "package:MSGARCH"))
  dist_func$load_theta(theta,ind)
  out =  dist_func$f_invsample(y)
  return(out)
}

#' @export
pdist = function(y, dist = "norm", shape = 100, skew = 1){
  valid.distribution <- c("norm", "std", "ged")
  if (!any(dist == valid.distribution)) {
    stop(paste0("error: The distribution ", dist, "
                does not appear to be a valid choice."))
  }
  if(dist == "norm"){
    ind = 1
  } else{
    ind = 0
  }
  if (skew == 1) {
    skew_lab <- "sym"
  } else {
    skew_lab <- "skew"
  }
  theta = c(shape, skew)
  dist_func = new(Class = getFromNamespace(paste0(dist,"_",skew_lab),pos = "package:MSGARCH"))
  dist_func$load_theta(theta, ind)
  out = dist_func$f_cdf(y)
  return(out)
}
