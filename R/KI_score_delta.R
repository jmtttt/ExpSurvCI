#' Scord-Delta Confidence Interval for Ratio of Parameters for Two Exponentially Distributed Samples
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
#' @return a data.frame with columns Sd1 (lower CI) and Sd2 (upper CI)
#' @export
#'
#' @examples
#' res <- KI_score_delta(x_mean = c(3,4),
#'                       y_mean = c(4,6),
#'                       nx     = c(100, 200),
#'                       ny     = c(200, 100),
#'                       dx     = c(.4, .5),
#'                       dy     = c(.5, .7),
#'                       alpha  = .05)
#'
KI_score_delta <- function(x_mean, y_mean, nx, ny, dx, dy, alpha = .05){
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

  c_sd <- nx * ny * (nx * dx + ny * dy) / (chi * dx * dy)
  res <- solve_quadratic_equation(
    a = c_sd * dx^2 * y_mean^2 - ny^2 * y_mean^2,
    b = -c_sd * 2 * dx * dy * x_mean * y_mean - 2 * nx * ny * x_mean * y_mean,
    c = c_sd * dy^2 * x_mean^2 - nx^2 * x_mean^2,
    current_task = "Score-Delta"
  )
  res <- as.data.frame(res, col.names = paste0("Sd", c("1","2")))
  return(res)
}
