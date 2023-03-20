#' B-beta function of GMM model
#'
#' \eqn{b_\beta(\theta) = \dfrac{1}{1 + \beta} \int f(x; \theta) ^ {1 + \beta} d x}
#'
#' @inheritParams valid_normix_input
#' @inheritParams stats::integrate
#' @param beta `numeric(1)` - parameter corresponds to the beta-divergence induced likelihood
#'
#' @return `numeric(1)`
#' @export
#'
#' @examples
#' set.seed(42)
#' lambda <- c(0.25, 0.75)
#' mu <- c(0, 4)
#' sigma <- c(1, 1)
#' b_beta(lambda, mu, sigma, beta=0.2)
#' b_beta(lambda, mu, sigma, beta=1e-3)
b_beta <- function(lambda, mu, sigma, beta=1e-3, rel.tol=1e-12) {
  valid_normix_input(lambda, mu, sigma)
  f <- function(x) {
    dnormix(x, lambda, mu, sigma)
  }
  integral <- stats::integrate(\(x) {f(x)^(1 + beta)}, -Inf, Inf,
                               rel.tol=rel.tol)
  res <- integral$value / (1 + beta)
  attr(res, "abs.error") <- integral$abs.error / (1 + beta)
  res
}

b_beta_d_lambda <- function(lambda, mu, sigma, beta=1e-3, rel.tol=1e-12) {
  integrand_d_b_lambda <- function(x, k) {
    dnormix(x, lambda, mu, sigma)^beta * stats::dnorm(x, mu[k], sigma[k])
  }
  valid_normix_input(lambda, mu, sigma)
  K <- length(lambda)
  res <- numeric(K)
  for (k in 1:K) {
    res[k] <- stats::integrate(\(x){integrand_d_b_lambda(x, k)}, -Inf, Inf, rel.tol=rel.tol)$value
  }
  res
}

b_beta_d_mu <- function(lambda, mu, sigma, beta=1e-3, rel.tol=1e-12) {
  integrand_d_b_mu <- function(x, k) {
    dnormix(x, lambda, mu, sigma)^beta *
      stats::dnorm(x, mu[k], sigma[k]) * (x - mu[k]) *
      lambda[k] / (sigma[k] ^ 2)
  }
  valid_normix_input(lambda, mu, sigma)
  K <- length(lambda)
  res <- numeric(K)
  for (k in 1:K) {
    res[k] <- stats::integrate(\(x){integrand_d_b_mu(x, k)}, -Inf, Inf, rel.tol=rel.tol)$value
  }
  res
}

b_beta_d_var <- function(lambda, mu, sigma, beta=1e-3, rel.tol=1e-12) {
  integrand_d_b_var <- function(x, k) {
    dnormix(x, lambda, mu, sigma)^beta *
      stats::dnorm(x, mu[k], sigma[k]) *
      ((x - mu[k])^2 - sigma[k]^2) *
      lambda[k] / (sigma[k] ^ 4) / 2
  }
  valid_normix_input(lambda, mu, sigma)
  K <- length(lambda)
  res <- numeric(K)
  res <- numeric(K)
  for (k in 1:K) {
    res[k] <- stats::integrate(\(x){integrand_d_b_var(x, k)}, -Inf, Inf, rel.tol=rel.tol)$value
  }
  res
}

# b_beta_d_lambda_fdm <- function(lambda, mu, sigma, beta=1e-3, eps=1e-5) {
#   # finite difference method, do not recommend
#   valid_normix_input(lambda, mu, sigma)
#   K <- length(lambda)
#   res <- numeric(K)
#   for (k in 1:K) {
#     e <- numeric(K)
#     e[k] <- 1
#     res[k] <- (b_beta(lambda + eps * e, mu, sigma, beta) -
#                b_beta(lambda - eps * e, mu, sigma, beta)) / (2 * eps)
#   }
#   res
# }

# b_beta_d_mu_fdm <- function(lambda, mu, sigma, beta=1e-3, eps=1e-5) {
#   valid_normix_input(lambda, mu, sigma)
#   K <- length(lambda)
#   res <- numeric(K)
#   for (k in 1:K) {
#     e <- numeric(K)
#     e[k] <- 1
#     res[k] <- (b_beta(lambda, mu + eps * e, sigma, beta) -
#                b_beta(lambda, mu - eps * e, sigma, beta)) / (2 * eps)
#   }
#   res
# }


