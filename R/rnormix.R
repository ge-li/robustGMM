#' Generate random samples from normal mixture
#'
#' @param n `numeric(1)` - number of samples
#' @inheritParams valid_normix_input
#' @param method "standard" behaves like [rnorm]. "montecarlo" is more efficient
#' and recommended for Monte Carlo computation based on mixture samples; the
#' samples are ordered based on each component. "fixed" will use expected number
#' of samples in each component and produce a `data.frame`.
#'
#' @return `numeric` if "standard" or "montecarlo"; `data.frame` if "fixed"
#' @export
#'
#' @examples
#' set.seed(42)
#' rnormix(10, c(0.5, 0.5), c(-10, 10), c(1, 1), "standard")
#' rnormix(10, c(0.5, 0.5), c(-10, 10), c(1, 1), "montecarlo")
#' rnormix(10, c(0.5, 0.5), c(-10, 10), c(1, 1), "fix")
#'
rnormix <- function(n, lambda, mu, sigma, method=c("standard", "montecarlo", "fix")) {
  valid_normix_input(lambda, mu, sigma)
  method = match.arg(method)
  k <- length(lambda)
  if (k == 1) {return(stats::rnorm(n, mu, sigma))}
  m <- stats::rmultinom(1, n, lambda)
  if (method == "standard") {
    purrr::map(1:k, ~rnorm(m[.], mu[.], sigma[.])) |> unlist() |> sample()
  } else if (method == "montecarlo") {
    purrr::map(1:k, ~rnorm(m[.], mu[.], sigma[.])) |> unlist()
  } else if (method == "fix") {
    m <- as.integer(n * lambda)
    m[-1] <- n - sum(m[1:(k-1)])
    x <- purrr::map(1:k, ~rnorm(m[.], mu[.], sigma[.])) |> unlist()
    y <- rep(0:(k-1), m) |> as.factor()
    data.frame(x, y)
  }
}
