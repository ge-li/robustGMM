#' Check if inputs are valid for *normix()
#'
#' @param lambda `numeric` - Vector of mixing proportions of each normal component
#' @param mu `numeric` - Vector of means of each normal component
#' @param sigma `numeric` - Vector of standard deviations of each normal component
#'
#' @return NULL
#'
valid_normix_input <- function(lambda, mu, sigma) {
  stopifnot("Input length not equal" =
              length(unique(c(length(lambda), length(mu), length(sigma)))) == 1)
  stopifnot(all.equal(sum(lambda), 1))
  stopifnot(min(sigma) > 0)
}
