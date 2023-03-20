#' Beta-Likelihood of Gaussian Mixture Model
#'
#' \eqn{l_\beta(\theta) = \dfrac{1}{n \beta} \sum_{i=1}^n f(x_i; \beta)^{\beta} - b_\beta(\theta)}, where
#' \eqn{b_\beta(\theta) = \dfrac{1}{1 + \beta} \int f(x; \theta) ^ {1 + \beta} d x}
#'
#' @param x `numeric(n)` - A numeric vector of observations
#' @inheritParams b_beta
#'
#' @return `numeric(1)` - calculated beta-likelihood
#' @export
#'
#' @examples
#' set.seed(42)
#' lambda <- c(0.25, 0.75)
#' mu <- c(0, 4)
#' sigma <- c(1, 1)
#' x <- rnormix(n=100, lambda, mu, sigma)
#' beta_likelihood(x, lambda, mu, sigma, beta=0.2)
#' beta_likelihood(x, lambda, mu, sigma, beta=1e-3)
beta_likelihood <- function(x, lambda, mu, sigma, beta=1e-3, rel.tol=1e-12) {
  valid_normix_input(lambda, mu, sigma)
  f <- function(x) {
    dnormix(x, lambda, mu, sigma)
  }
  mean(f(x)^beta) / beta - b_beta(lambda, mu, sigma, beta, rel.tol)
}
