#' Wald Confidence Interval for Ratio of Parameters for Two Exponentially Distributed Samples
#'
#' @import stats
#'
#' @param x_mean Mean of first group
#' @param y_mean Mean of second group
#' @param nx sample size of first group
#' @param ny sample size of second group
#' @param dx observation probability (=1-censoring prob.) of first group
#' @param dy observation probability (=1-censoring prob.) of second group
#' @param alpha confidence level
#'
#' @return a data.frame with columns W1 (lower CI) and W2 (upper CI)
#' @export
#'
#' @examples
#' res <- KI_wald(x_mean = c(3,4),
#'                y_mean = c(4,6),
#'                nx     = c(100, 200),
#'                ny     = c(200, 100),
#'                dx     = c(.4, .5),
#'                dy     = c(.5, .7),
#'                alpha  = .05)
#'
KI_wald <- function(x_mean, y_mean, nx, ny, dx, dy, alpha = .05){
  checkmate::assert_numeric(x_mean)
  checkmate::assert_numeric(y_mean)
  checkmate::assert_numeric(nx,   lower = 0)
  checkmate::assert_numeric(ny,   lower = 0)
  checkmate::assert_numeric(dx,   lower = 0, upper = 1)
  checkmate::assert_numeric(dy,   lower = 0, upper = 1)
  checkmate::assert_number(alpha, lower = 0, upper = .3)
  checkmate::assert_set_equal(c(length(x_mean), length(nx), length(dx)), length(x_mean))
  checkmate::assert_set_equal(c(length(y_mean), length(ny), length(dy)), length(y_mean))

  chi <- qchisq(p = 1 - alpha,
                df = 1)
  if(dx == 0 || dy == 0){
    k_exp <- NA_integer_
  } else{
    k_exp <- dy * x_mean / (dx * y_mean)
  }

  res <- solve_quadratic_equation(
    a = nx * ny * dx * dy - chi * (nx * dx + ny * dy),
    b = -2 * k_exp * nx * ny * dx * dy,
    c = nx * ny * k_exp^2 * dx * dy,
    current_task = "Wald-Delta"
  )
  res <- as.data.frame(res, col.names = paste0("W", c("1","2")))
  return(res)
}
