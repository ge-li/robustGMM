beta_est_func <- function(x, lambda, mu, sigma, beta=0.2, rel.tol=1e-12) {

  valid_normix_input(lambda, mu, sigma)
  N <- length(x) #
  K <- length(lambda)

  gamma <- beta_responsibility(x, lambda, mu, sigma, beta)

  b_d_lambda <- b_beta_d_lambda(lambda, mu, sigma, beta)
  b_d_mu <- b_beta_d_mu(lambda, mu, sigma, beta)
  b_d_var <- b_beta_d_var(lambda, mu, sigma, beta)

  # tau <- sum(gamma) - sum(lambda * b_d_lambda)

  l_d_mu <- matrix(nrow=N, ncol=K)
  l_d_var <- matrix(nrow=N, ncol=K)
  l_d_lambda <- matrix(nrow=N, ncol=K)
  for (k in 1:K) {
    l_d_mu[, k] <- gamma[, k] * (x - mu[k]) / (sigma[k] ^ 2) - b_d_mu[k]
    l_d_lambda[, k] <- gamma[, k] / lambda[k] - b_d_lambda[k]
    l_d_var[, k] <- gamma[, k] * ((x - mu[k])^2 - sigma[k]^2) / (sigma[k]^4) / 2 - b_d_var[k]
  }
  print(colMeans(l_d_mu))
  print(colMeans(l_d_lambda))
  print(colMeans(l_d_var))
}
