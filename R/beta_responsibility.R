#' Generalized (Beta) Responsibility
#'
#' The generalized beta responsibility of a Gaussian Mixture Model (GMM) is \eqn{\gamma_{i, k}(\beta) = \dfrac{\lambda_k N(x_i \mid \mu_k, \Sigma_k)}{f(x_i; \theta)^{1 - \beta}}}. When `beta = 0`, \eqn{\gamma_{i, k}(0) = \Pr(Z_{k} = 1 \mid x_i)}.
#'
#' @param x `numeric(n)` A vector with n observations from a Gaussian Mixture Model.
#' @inheritParams valid_normix_input
#' @param beta `numeric(1)` - parameter corresponds to the beta-divergence induced likelihood.
#'
#' @return `matrix` A matrix with n rows and k columns representing the generalized responsibilities.
#' @export
#'
#' @examples
#' set.seed(42)
#' lambda <- c(0.25, 0.75)
#' mu <- c(0, 4)
#' sigma <- c(1, 1)
#' x <- rnormix(n=4, lambda, mu, sigma)
#' beta_responsibility(x, lambda, mu, sigma, beta=0)
#' beta_responsibility(x, lambda, mu, sigma, beta=0.2)
beta_responsibility <- function(x, lambda, mu, sigma, beta) {
  N <- length(x)
  K <- length(lambda)
  gamma <- matrix(nrow=N, ncol=K)
  f <- dnormix(x, lambda, mu, sigma) # empirical likelihood
  for (k in 1:K) {
    gamma[, k] <- lambda[k] * stats::dnorm(x, mu[k], sigma[k])
  }
  gamma <- gamma / (f ^ (1 - beta))
  gamma
}
