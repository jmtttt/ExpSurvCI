#' helper function: solve quadratic equations (suports vectorized form)
#'
#' @import checkmate
#'
#' @param a numeric (or numeric vector) for parameter
#' @param b numeric (or numeric vector) for parameter
#' @param c numeric (or numeric vector) for parameter
#' @param current_task character-identifier to improve traceback, if equation has no zeros
#'
#' @details This function solves quadratic equations of the form \deqn{a*x^2 + b*x + c = 0}. It also solves multiple equations at once, when the parameters are passed as vectors.
#'
#' @return a two-element list with both zeros (or a vector of such)
#' @export
#'
#' @examples
#' a <- c(1,2,3)
#' b <- c(0,4,6)
#' c <- c(-1,-6,-15)
#' res <- solve_quadratic_equation(a, b, c)
#'
solve_quadratic_equation <- function(a,b,c,
                                     current_task = "<undefined>"){
  checkmate::assert_numeric(a, any.missing = TRUE)
  checkmate::assert_numeric(b, any.missing = TRUE)
  checkmate::assert_numeric(c, any.missing = TRUE)
  n <- length(a)
  checkmate::assertSetEqual(y = c(n, n, n),
                            x = c(length(a), length(b), length(c)))

  root_term <-  b^2 - 4 * a * c
  real_res  <- (b^2 - 4 * a * c) >= 0
  real_res[is.na(real_res)] <- FALSE
  if(any(!real_res, na.rm = TRUE)){
    warning(sprintf(paste("When calculating %s Confidence Interval,",
                          "there were %i out of %i equations",
                          "without solution."),
                    current_task, sum(!real_res), n))
  }


  KI_1 <- rep.int(NA_real_, times = n)
  KI_2 <- rep.int(NA_real_, times = n)

  KI_1[real_res] <- (-b[real_res] - sqrt(root_term[real_res])) / (2 * a[real_res])
  KI_2[real_res] <- (-b[real_res] + sqrt(root_term[real_res])) / (2 * a[real_res])
  return(list(s1 = KI_1,
              s2 = KI_2))
}
