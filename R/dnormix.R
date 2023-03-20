#' Density of normal mixture
#'
#' Calculate probability density function of a mixture of univariate normal distributions.
#'
#' @param x `numeric` - Vector of values at which to evaluate the density
#' @inheritParams valid_normix_input
#' @param log `logical` - Logical value indicating whether to return log density (TRUE) or density (FALSE)
#'
#' @return `numeric`; Vector of densities or log densities
#' @export
#'
#' @examples
#' dnormix(-5:5, c(0.7, 0.3), c(0, 3), c(1, 1))
#'
#' # visualize
#' p <- function(x) {dnormix(x, c(0.7, 0.3), c(0, 3), c(1, 1))}
#' lp <- function(x) {dnormix(x, c(0.7, 0.3), c(0, 3), c(1, 1), log=TRUE)}
#' par(mfrow=c(1, 2))
#' curve(p, -5, 5)
#' curve(lp, -5, 5)
#' par(mfrow=c(1, 1))
#'
#' # supports fast numerical integration without underflow/overflow
#' integrate(p, -Inf, Inf) # equals 1
#' integrate(function(x){- p(x) * lp(x)}, -Inf, Inf) # entropy
#'
dnormix <- function(x, lambda, mu, sigma, log=FALSE) {
  valid_normix_input(lambda, mu, sigma)
  k <- length(lambda)
  if (k == 1) {return(stats::dnorm(x, mu, sigma, log))}
  if (log) {# returns the log density
    ld <- sapply(1:k, \(i) {log(lambda[i]) + stats::dnorm(x, mu[i], sigma[i], log=TRUE)})
    if (is.vector(ld)) {
      return(lse(ld))
    } else{
      return(apply(ld, MARGIN=1, lse))
    }
  } else {
    d <- 0
    for (i in 1:k) { # for loop is actually faster than *apply()
      d <- d + lambda[i] * stats::dnorm(x, mu[i], sigma[i])
    }
    d
  }
}
