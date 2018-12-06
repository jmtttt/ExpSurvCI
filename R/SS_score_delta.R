#' Sample size calculation based on score statistics with delta approximation
#'
#' @param k1 estimated ratio of parameters
#' @param dx observation probability in first group
#' @param dy observation probability in second group
#' @param alpha confidence level (type 1 error)
#' @param beta 1 - Power (type 2 error)
#' @param k0 Effect under Null Hypothesis (usually equal to 1)
#' @param c ratio ny/nx [where n* is sample size of first/second group]
#'
#' @details Warning: Simulations show, that this formula might be unreliable!
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
#' SampleSize_score_delta(2, .9, .9)
#'
SampleSize_score_delta <- function(k1, dx, dy, alpha = .05, beta = .2, k0 = 1, c = 1){
  checkmate::assert_number(k1, lower = 0, finite = TRUE)
  checkmate::assert_number(dx, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_number(dy, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_number(alpha, lower = 0, upper = .4, finite = TRUE)
  checkmate::assert_number(beta, lower = 0, upper = .4, finite = TRUE)
  checkmate::assert_number(k0, lower = 0, finite = TRUE)
  checkmate::assert_number(c, lower = 0, finite = TRUE)

  z1 <- qnorm(1 - alpha / 2)
  z2 <- qnorm(1 - beta)

  # output from wolfram (Mathematica) for call
  # FortranForm[Solve[{z1 == n * c * dx * dy * (dx + c * dy) * (k - k0)^2 / (dx * k + c * dy * k0)^2, z2 == n * c * dx * dy * (dx + c * dy) * (k - k1)^2 / (dx * k + c * dy * k1)^2},{k, n}]]
  # using the largest (4th.) solution
  # slightly casted by hand from fortran to R
  #  - replaced starting "   -   " by nothing
  #  - replaced "Sqrt" by "sqrt"
  #  - removed linebreaks "\<\n>"

  nx <- (c**2*dy**2*k0**3*z1 - 3*c**2*dy**2*k0**2*k1*z1 + 3*c**2*dy**2*k0*k1**2*z2 - c**2*dy**2*k1**3*z2 +
           2*c*dx*dy*k0**2*z1*(((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
                                 sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                                        (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                                        (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                                         (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                                         (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                                        (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.) + 2*c**2*dy**2*k0**2*z1*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.) - 6*c*dx*dy*k0*k1*z1*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.) + 6*c*dx*dy*k0*k1*z2*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.) - 2*c*dx*dy*k1**2*z2*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.) - 2*c**2*dy**2*k1**2*z2*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.) + dx**2*k0*z1*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.)**2 + 4*c*dx*dy*k0*z1*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.)**2 - 3*dx**2*k1*z1*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.)**2 + 3*dx**2*k0*z2*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.)**2 - dx**2*k1*z2*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.)**2 - 4*c*dx*dy*k1*z2*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.)**2 + 2*dx**2*z1*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.)**3 - 2*dx**2*z2*
           (((dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2))/(2.*dx*(z1 - z2)) + (-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)/(2.*dx*(z1 - z2)) +
              sqrt((2*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**2)/(dx**2*(z1 - z2)**2) - (-2*c*dx*dy*k0*k1*z1 + 2*c*dx*dy*k0*k1*z2)/(dx**2*z1 - dx**2*z2) -
                     (c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2)/(dx**2*(z1 - z2)) +
                     (dx*(z1 - z2)*((8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)**3)/(dx**3*(z1 - z2)**3) -
                                      (8*(-(c*dy*k0*z1) + dx*k1*z1 - dx*k0*z2 + c*dy*k1*z2)*(c**2*dy**2*k0**2*z1 - 4*c*dx*dy*k0*k1*z1 + dx**2*k1**2*z1 - dx**2*k0**2*z2 + 4*c*dx*dy*k0*k1*z2 - c**2*dy**2*k1**2*z2))/
                                      (dx**3*(z1 - z2)**2) - (16*(-(c**2*dy**2*k0**2*k1*z1) + c*dx*dy*k0*k1**2*z1 - c*dx*dy*k0**2*k1*z2 + c**2*dy**2*k0*k1**2*z2))/(dx**2*(z1 - z2))))/
                     (4.*(dx + c*dy)*(k0 - k1)*sqrt(z1)*sqrt(z2)))/2.)**3)/
    (c*dx**2*dy*k0**3 + c**2*dx*dy**2*k0**3 - 3*c*dx**2*dy*k0**2*k1 - 3*c**2*dx*dy**2*k0**2*k1 + 3*c*dx**2*dy*k0*k1**2 + 3*c**2*dx*dy**2*k0*k1**2 - c*dx**2*dy*k1**3 - c**2*dx*dy**2*k1**3)

  ny <- c * nx
  nx <- ceiling(nx)
  ny <- ceiling(ny)
  N  <- nx + ny

  return(list(N  = N,
              Nx = nx,
              Ny = ny))
}

