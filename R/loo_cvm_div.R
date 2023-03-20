#' Leave-One-Out Cramer-von Mises Divergence
#'
#' @inheritParams robustGMM
#'
#' @return
#' @export
#'
#' @examples
loo_cvm_div <- function(x, lambda, mu, sigma, beta=0.2,
                        tol=1e-5, maxiter=1000) {
  x <- sort(x)
  N <- length(x)
  div <- numeric(N)
  for (i in 1:N) {
    mod <- robustGMM(x[-i], lambda, mu, sigma, beta, tol)
    div[i] <-  ((i - 0.5) / N - pnormix(x[i], mod$lambda, mod$mu, mod$sigma))^2
  }
  mean(div)
}


