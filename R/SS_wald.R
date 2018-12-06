#' Sample size calculation based on wald statistics
#'
#' @param k1 estimated ratio of parameters
#' @param dx observation probability in first group
#' @param dy observation probability in second group
#' @param alpha confidence level (type 1 error)
#' @param beta 1 - Power (type 2 error)
#' @param k0 Effect under Null Hypothesis (usually equal to 1)
#' @param c ratio ny/nx [where n* is sample size of first/second group]
#'
#' @return List (3 elements) containing estimated required sample size (observations + censorings)
#'         \itemize{
#'         \item in first group (Nx)
#'         \item in second group (Ny)
#'         \item altogether (N)
#'         }
#' @export
#'
#' @examples
#' SampleSize_wald(2, .9, .9)
#'
SampleSize_wald <- function(k1, dx, dy, alpha = .05, beta = .2, k0 = 1, c = 1){
  checkmate::assert_number(k1, lower = 0, finite = TRUE)
  checkmate::assert_number(dx, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_number(dy, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_number(alpha, lower = 0, upper = .4, finite = TRUE)
  checkmate::assert_number(beta, lower = 0, upper = .4, finite = TRUE)
  checkmate::assert_number(k0, lower = 0, finite = TRUE)
  checkmate::assert_number(c, lower = 0, finite = TRUE)

  c1 <- (dx + c * dy) / (c * dx * dy)
  z1 <- qnorm(1 - alpha / 2)
  z2 <- qnorm(1 - beta)

  nx <- c1 * (z1 * k0 + z2 * k1)^2 / (k1 - k0)^2
  ny <- c * nx
  nx <- ceiling(nx)
  ny <- ceiling(ny)
  N  <- nx + ny

  return(list(N  = N,
              Nx = nx,
              Ny = ny))
}

