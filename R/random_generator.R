#' Generate Random Numbers from a specified distribution
#'
#' @import stats
#' @import checkmate
#'
#' @details This function is a wrapper for various random number generators
#'          (for different distributions). Its main purpose is to deliver
#'          an interface for mixed distributions.
#'
#'          Supported are
#'          \itemize{
#'            \item exp (exponential distribution), param[1] = parameter
#'            \item weibull (weibull distribution), param[1] = shape, param[2] = scale
#'            \item norm (normal distribution), param[1] = mu, param[2] = sd, if negative set to 0
#'            \item unif (uniform distribution from 0 to param[1])
#'            \item gamma (gamma distribution with shape = param[1] and scale = param[2])
#'            \item dirac not really random, but N replicates of param[1]
#'            \item triangle (random numbers with parabolic distribution), param[1] = length (= max - min), param[2] = max
#'          }
#'
#' @param dist a character specifying the desired distribution
#' @param param A Vector giving the parameters of the distributions respectively, see details
#' @param N Number of Randoms
#'
#' @return N random numbers
#' @export
#'
random_from_dist <- function(dist, param=NULL, N = 1e7){
  checkmate::assert_character(x = dist)
  checkmate::assert_choice(x = dist, choices = c("exp", "weibull", "norm", "unif", "gamma", "dirac", "triangle", "followup", "accrual-followup"))
  checkmate::assert_numeric(x = param, max.len = 3)
  checkmate::assert_number(x = N, lower = 0, upper = 1e10)
  # if dist is elementary, just simulate
  rnd <- switch (dist,
          "exp" = rexp(n    = N,
                       rate = param[1]),
          "weibull" = rweibull(n = N,
                               shape = param[1],
                               scale = param[2]),
          "norm" = mapply(max, 0, rnorm(n = N,
                                        mean = param[1],
                                        sd = param[2])),
          "unif" = runif(n = N,
                         min = 0,
                         max = param[1]),
          "gamma" = rgamma(n = N,
                           shape = param[1],
                           scale = param[2]),
          "dirac" = rep.int(x = param[1],
                            times = N),
          "triangle" = param[2] - param[1] * abs(runif(n = N,
                                                       min = -.5,
                                                       max = .5) +
                                                   runif(n = N,
                                                         min = -.5,
                                                         max = .5)),
          c()
  )
  # if dist is a compound, split it up
  if(length(rnd) == 0){
    if(dist == "followup"){
      nc <- ceiling(param[2] * param[1] * N)
      ne <- N - nc
      rnd <- c(random_from_dist(dist = "unif", param = param[2], N = nc),
               random_from_dist(dist = "dirac", param = param[2], N = ne))
    }
    if(dist == "accrual-followup"){
      nc <- ceiling(param[2] * param[1] * N)
      ne <- N - nc
      rnd <- c(random_from_dist(dist = "unif", param = param[2], N = nc),
               random_from_dist(dist = "unif", param = c(param[3]), N = ne) + (param[2] - param[3]))
    }
  }
  # return anyway
  return(rnd)
}
