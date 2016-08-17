f.theta.mixture <- function(K, nb_params, theta) {
  if (is.vector(theta))
    theta <- matrix(theta, nrow = 1)
  n_total_params <- ncol(theta)
  n_theta <- nrow(theta)
  p_vector <- theta[(nb_params + 1):n_total_params]
  tmp <- theta[1, ]
  idx <- 1
  for (i in 1:(K - 1)) {
    tmp[(nb_params + idx):(nb_params + idx + K - 1)] <- p_vector[i]
    idx <- idx + K
  }
  new_theta <- matrix(nrow = n_theta, ncol = length(tmp))
  for (i in 1:n_theta) {
    p_vector <- theta[i, (nb_params + 1):n_total_params]
    idx <- 1
    for (j in 1:(K - 1)) {
      new_theta[i, (nb_params + idx):(nb_params + idx + K - 1)] <- p_vector[j]
      idx <- idx + K
    }
  }
  new_theta[, 1:nb_params] <- theta[, 1:nb_params]
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