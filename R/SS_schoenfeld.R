#' Sample size calculation based on Schoenfeld formula
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
#' SampleSize_schoenfeld(2, .9, .9)
#'
SampleSize_schoenfeld <- function(k1, dx, dy, alpha = .05, beta = .2, k0 = 1, c = 1){
  checkmate::assert_number(k1, lower = 0, finite = TRUE)
  checkmate::assert_number(dx, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_number(dy, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_number(alpha, lower = 0, upper = .4, finite = TRUE)
  checkmate::assert_number(beta, lower = 0, upper = .4, finite = TRUE)
  checkmate::assert_number(k0, lower = 1, upper = 1, finite = TRUE)
  checkmate::assert_number(c, lower = 0, finite = TRUE)

  z1 <- qnorm(1 - alpha / 2)
  z2 <- qnorm(1 - beta)
  c1 <- (1 + c)^3 / (c * dx + c^2 * dy)
  c2 <- log(k1)^2

  N  <- (z1 + z2)^2 / c2 * c1
  nx <- N / (1 + c)
  ny <- N * c / (1 + c)
  nx <- ceiling(nx)
  ny <- ceiling(ny)
  N  <- nx + ny

  return(list(N  = N,
              Nx = nx,
              Ny = ny))
}

