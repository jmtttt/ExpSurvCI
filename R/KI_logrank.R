#' Pseudo (!) Confidence Interval for Logrank test for the Ratio of Parameters for Two Exponentially Distributed Samples
#'
#' @import survival
#' @import stats
#'
#' @details this function just performs a logrank test and returns
#'          \itemize{
#'          \item (0, 1), if the test is significant and mean(Y) > mean(X)
#'          \item (1, Inf), if the test is significant and mean(X) > mean(Y)
#'          \item (0, Inf), if the test is not significant.
#'          }
#'          Warning! It shall be noted, that theese pseudo confidence intervalls do not hold any confidence levels as it may be intended.
#'          Rather, this procedure is implemented to fit a power simulation for the logrank test into the framwork of the other confidence intervals.
#'
#' @param X numeric matrix -- Observed survival times in first group
#' @param Y numeric matrix -- Observed survival times in second group
#' @param Dx numeric matrix -- Observation indicator (1 for events, 0 for censorings) for first group
#' @param Dy numeric matrix -- Observation indicator (1 for events, 0 for censorings) for second group
#' @param alpha confidence level
#'
#' @return a data.frame with columns L1 (lower CI) and L2 (upper CI), giving pseudo confidence levels.
#' @export
#'
KI_logrank <- function(X, Y, Dx, Dy,
                       alpha = .05){
  checkmate::assert_numeric(X,  lower = 0)
  checkmate::assert_numeric(Y,  lower = 0)
  checkmate::assert_numeric(Dx, lower = 0, upper = 1)
  checkmate::assert_numeric(Dy, lower = 0, upper = 1)
  checkmate::assert_number(alpha, lower = 0, upper = .3)
  checkmate::assert_set_equal(dim(X), dim(Dx))
  checkmate::assert_set_equal(dim(Y), dim(Dy))
  checkmate::assert_set_equal(c(nrow(X), nrow(Y), nrow(Dx)), nrow(Dy))

  nx <- ncol(X)
  ny <- ncol(Y)

  if(is.vector(X)){
    X <- matrix(data = X, nrow = 1)
  }
  if(is.vector(Dx)){
    Dx <- matrix(data = Dx, nrow = 1)
  }
  if(is.vector(Y)){
    Y <- matrix(data = Y, nrow = 1)
  }
  if(is.vector(Dy)){
    Dy <- matrix(data = Dy, nrow = 1)
  }

  surv_data <- lapply(X   = 1:nrow(X),
                      FUN = function(t){
                        res <- data.frame(time  = survival::Surv(time  = c(X[t, ], Y[t, ]),
                                                                 event = c(Dx[t, ], Dy[t, ])),
                                          group = c(rep.int(x = 0, times = nx),
                                                    rep.int(x = 1, times = ny)))
                        return(res)
                      })

  res <- sapply(X = surv_data,
                FUN = function(x){
                  fm   <- survival::survdiff(formula = time ~ group,
                                            data = x)
                  isSignificant <- fm$chisq > qchisq(p = 1 - alpha,
                                                     df = 1)
                  if(isSignificant){
                    if(fm$exp[1] > fm$obs[1]){
                      res <- c(lower = 1+1e-5,
                               upper = 1e10)
                    } else{
                      res <- c(lower = 0,
                               upper = 1-1e-5)
                    }
                  } else{
                    res <- c(lower = 0,
                             upper = 1e10)
                  }
                  return(res)
                },
                simplify = TRUE)
  res <- as.data.frame(x = t(res))
  colnames(res) <- paste0("L", c("1","2"))
  return(res)
}
