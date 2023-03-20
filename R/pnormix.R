#' Cumulative distribution function of normal mixture
#'
#' @param x `numeric` Vector of values at which to evaluate the (log, complementary) cdf.
#' @inheritParams dnormix
#' @param lower.tail `logical` - if TRUE (default), probabilities are `P[X ≤ x]` otherwise, `P[X > x]`.
#' @param log.p `logical`; if TRUE, probabilities p are given as `log(p)`.
#'
#' @return `numeric` Vector of (log) probabilities.
#' @export
#'
#' @examples
#' pnormix(-5:5, c(0.7, 0.3), c(0, 3), c(1, 1))
#'
pnormix <- function(x, lambda, mu, sigma, lower.tail=TRUE, log.p=FALSE) {
  valid_normix_input(lambda, mu, sigma)
  k <- length(lambda)
  if (k == 1) {return (stats::pnorm(x, mu, sigma, lower.tail, log.p))}
  lp <- sapply(1:k, \(i) {log(lambda[i]) +
      stats::pnorm(x, mu[i], sigma[i], lower.tail, log.p=TRUE)})
  if (is.vector(lp)) {
    lp <- lse(lp)
  } else{
    lp <- apply(lp, MARGIN=1, lse)
  }
  if (log.p) {
    lp
  } else{
    exp(lp)
  }
}
