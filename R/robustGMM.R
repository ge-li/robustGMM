#' Fitting robust Gaussian Mixture Model
#'
#' @inheritParams beta_likelihood
#' @param tol
#' @param maxiter Max iteration
#' @param verbose
#' @param keep.data
#'
#' @return
#' @export
#'
#' @examples
robustGMM <- function(x,
                      lambda = NULL,
                      mu = NULL,
                      sigma = NULL,
                      beta = 0.2,
                      tol = 1e-6,
                      maxiter = 1000,
                      verbose = FALSE,
                      keep.data = FALSE)
{
  valid_normix_input(lambda, mu, sigma)
  N <- length(x) #
  K <- length(lambda)

  logger <- function(iter, l_beta) {
    s_iter <- format(iter, width=4)
    s_l_beta <- format(round(l_beta, 6), width=12, nsmall=6)
    glue::glue("Iter: {s_iter}, Beta-likelihood: {s_l_beta}")
  }
  l_beta <- beta_likelihood(x, lambda, mu, sigma, beta)
  if (verbose) {
    cli::cli_h1(glue::glue("Fitting Robust Gaussian Mixture Model"))
    cli::cli_alert_info(logger(0, l_beta))
  }

  for (iter in 1:maxiter) {
    # E-step: calculate the responsibilities
    gamma <- beta_responsibility(x, lambda, mu, sigma, beta)
    N_beta <- colSums(gamma)

    # M-step: update lambda, mu, sigma
    b_d_lambda <- b_beta_d_lambda(lambda, mu, sigma, beta)
    b_d_mu <- b_beta_d_mu(lambda, mu, sigma, beta)
    b_d_var <- b_beta_d_var(lambda, mu, sigma, beta)

    mu_new <- (colSums(gamma * x) - N * sigma^2 * b_d_mu) / N_beta
    sigma_new <- numeric(K)
    for (k in 1:K) {
      sigma_new[k] <- sqrt((sum(gamma[, k] * (x - mu[k]) ^ 2) -
                              2 * N * sigma[k] ^ 4 * b_d_var[k]) / N_beta[k])
    }
    lambda_new <- (N_beta -
                     N * lambda * b_d_lambda +
                     N * lambda * sum(lambda * b_d_lambda)) / sum(gamma)

    l_beta_new <- beta_likelihood(x, lambda_new, mu_new, sigma_new, beta)

    if (abs(l_beta_new - l_beta) < tol) {
      res <- list()
      if (keep.data) {res$x <- x}
      res$lambda <- lambda_new
      res$mu <- mu_new
      res$sigma <- sigma_new
      res$beta <- beta
      res$l_beta <- l_beta_new
      if (verbose) {cli::cli_alert_success(logger(iter, l_beta_new))}
      return(res)
    } else {
      if (verbose) {cli::cli_alert_info(logger(iter, l_beta_new))}
      lambda <- lambda_new
      mu <- mu_new
      sigma <- sigma_new
      l_beta <- l_beta_new
    }
  }
}
