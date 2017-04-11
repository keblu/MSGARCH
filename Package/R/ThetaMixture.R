f.theta.mixture <- function(K, nb_params, theta) {
  if (is.vector(theta))
    theta <- matrix(theta, nrow = 1)
  n_total_params = ncol(theta)
  new_theta = matrix(nrow(theta),ncol = nb_params+K*(K-1))
  for (i in 1:nrow(theta)) {
    p_vector <- theta[i, (nb_params + 1):n_total_params]
    new_p = rep(p_vector,K)
    new_theta[i,] = c(theta[i, 1:nb_params],new_p)
  }
  return(new_theta)
}

f.theta.mixture.reverse <- function(K, nb_params, theta) {
  if (is.vector(theta))
    theta <- matrix(theta, nrow = 1)
  n_total_params <- ncol(theta)
  n_theta <- nrow(theta)
  new_theta <- matrix(nrow = n_theta, ncol = nb_params + K - 1)
  for (i in 1:n_theta) {
    p_vector <- theta[i, (nb_params + 1):n_total_params]
    new_theta[i, 1:nb_params] <- theta[i, 1:nb_params]
    idx <- 1
    for (j in 1:(K - 1)) {
      new_theta[i, (nb_params + j)] <- p_vector[idx]
      idx <- idx + K
    }
  }
  return(new_theta)
}