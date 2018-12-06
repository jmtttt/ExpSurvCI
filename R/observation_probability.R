#' Calculate Observation Probability given Distributions for Event and Censoring Times
#'
#' @import checkmate
#'
#' @param dist_E distribution identifier for event distribution, see \code{\link{random_from_dist}}
#' @param param_E parameter set for event distribution, see \code{\link{random_from_dist}}
#' @param dist_C distribution identifier for censoring distribution, see \code{\link{random_from_dist}}
#' @param param_C parameter set for censoring distribution, see \code{\link{random_from_dist}}
#' @param N number of replicates
#'
#' @details For a given distribution of E and C, N samples each are drawn to
#'          approximate the probability for E_i < C_i, say the event happens before censoring
#'
#' @return  a numeric value 0<p<1, denoting observation probability
#' @export
#'
#' @examples
#'   d_o <- observation_probability(dist_E  = "exp",
#'                                  param_E = c(1),
#'                                  dist_C  = "accrual-followup",
#'                                  param_C = c(.05, 5, 2),
#'                                  N       = 1e6)
#'
#'
observation_probability <- function(dist_E, param_E=NULL,
                                    dist_C, param_C=NULL,
                                    N = 1e7){
  checkmate::assert_character(x = dist_E)
  checkmate::assert_choice(x = dist_E, choices = c("exp", "weibull", "norm", "unif", "dirac", "triangle", "followup", "accrual-followup"))
  checkmate::assert_numeric(x = param_E, max.len = 3)
  checkmate::assert_character(x = dist_C)
  checkmate::assert_choice(x = dist_C, choices = c("exp", "weibull", "norm", "unif", "dirac", "triangle", "followup", "accrual-followup"))
  checkmate::assert_numeric(x = param_C, max.len = 3)
  checkmate::assert_number(x = N, lower = 0, upper = 1e10)
  # simulate numbers
  X <- random_from_dist(dist  = dist_E,
                        param = param_E,
                        N     = N)
  Y <- random_from_dist(dist  = dist_C,
                        param = param_C,
                        N     = N)
  # evaluate X<Y
  res <- mean(X < Y)
  return(res)
}
