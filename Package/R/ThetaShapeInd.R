
## DA improve the writing
f.theta.RegIndDist <- function(K, nb_total_params, nb_params_model, theta) {
  if (is.vector(theta))
    theta <- matrix(theta, nrow = 1)
  n_theta <- nrow(theta)
  nb_dist_params <- nb_total_params - nb_params_model
  if (nb_dist_params[1] == 0) {
    return(theta)
  }
  nb_params_transition <- ncol(theta) - sum(nb_params_model) - nb_dist_params[1]
  new_length <- sum(nb_params_model) + sum(nb_dist_params) + nb_params_transition
  new_theta <- matrix(data = NA, nrow = n_theta, ncol = new_length)
  for (i in 1:n_theta) {
    dist_params <- theta[i, (sum(nb_params_model) + 1):(sum(nb_params_model) + nb_dist_params[1])]
    p_params <- theta[i, (sum(nb_params_model) + nb_dist_params[1] + 1):length(theta[i, ])]
    ind1 <- 1
    ind2 <- 1
    for (j in 1:K) {
      new_theta[i, ind2:(ind2 + nb_params_model[j] - 1)] <- theta[i, ind1:(ind1 + nb_params_model[j] - 1)]
      new_theta[i, (ind2 + nb_params_model[j]):(ind2 + nb_params_model[j] + nb_dist_params[1] - 1)] <- dist_params
      ind1 <- ind1 + nb_params_model[j]
      ind2 <- ind2 + nb_params_model[j] + length(dist_params)
    }
    new_theta[i, ind2:new_length] <- p_params
  }
  return(new_theta)
}

f.theta.RegIndDist.reverse <- function(K, nb_total_params, nb_params_model, theta) {
  if (is.vector(theta))
    theta <- matrix(theta, nrow = 1)
  n_theta <- nrow(theta)
  nb_params <- length(theta[1, ])
  nb_dist_params <- nb_total_params - nb_params_model
  if (nb_dist_params[1] == 0) {
    return(theta)
  }
  new_length <- sum(nb_params_model) + nb_dist_params[1] + nb_params - sum(nb_total_params)
  new_theta <- matrix(data = NA, nrow = n_theta, ncol = new_length)
  for (i in 1:n_theta) {
    p_params <- theta[i, (sum(nb_total_params) + 1):nb_params]
    dist_params <- theta[i, (nb_params_model[1] + 1):(nb_params_model[1] + nb_dist_params[1])]
    ind1 <- 1
    ind2 <- 1
    for (j in 1:K) {
      new_theta[i, ind1:(ind1 + nb_params_model[j] - 1)] <- theta[i, ind2:(ind2 + nb_params_model[j] - 1)]
      ind1 <- ind1 + nb_params_model[j]
      ind2 <- ind2 + nb_params_model[j] + length(dist_params)
    }
    new_theta[i, ind1:(ind1 + length(dist_params) - 1)] <- dist_params
    ind1 <- ind1 + length(dist_params)
    new_theta[i, ind1:new_length] <- p_params
  }
  return(new_theta)
}
