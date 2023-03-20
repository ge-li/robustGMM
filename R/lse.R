#' Log-Sum-Exp
#'
#' \eqn{LSE(x_1, \ldots, x_n) = \log(\exp(x_1) + \dots + \exp(x_n))}
#'
#' @param x numeric
#'
#' @return LSE
#' @export
#'
#' @examples
#' lse(c(-1001, -1002, -1003))
lse <- function(x) {
  # log-sum-exp
  # https://stat.ethz.ch/pipermail/r-help/2011-February/269205.html
  # Or use matrixStats::logSumExp
  xmax <- which.max(x)
  log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
}
