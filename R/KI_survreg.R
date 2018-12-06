#' Confidence Interval for Regression Parameter for the Ratio of Parameters for Two Exponentially Distributed Samples
#'
#' @import survival
#'
#' @param X numeric matrix -- Observed survival times in first group
#' @param Y numeric matrix -- Observed survival times in second group
#' @param Dx numeric matrix -- Observation indicator (1 for events, 0 for censorings) for first group
#' @param Dy numeric matrix -- Observation indicator (1 for events, 0 for censorings) for second group
#' @param alpha confidence level
#'
#' @return a data.frame with columns R1 (lower CI) and R2 (upper CI)
#' @export
#'
KI_survreg <- function(X, Y, Dx, Dy,
                       alpha = .05){
  checkmate::assert_numeric(X,  lower = 0)
  checkmate::assert_numeric(Y,  lower = 0)
  checkmate::assert_numeric(Dx, lower = 0, upper = 1)
  checkmate::assert_numeric(Dy, lower = 0, upper = 1)
  checkmate::assert_number(alpha, lower = 0, upper = .3)
  checkmate::assert_set_equal(dim(X), dim(Dx))
  checkmate::assert_set_equal(dim(Y), dim(Dy))
  checkmate::assert_set_equal(c(nrow(X), nrow(Y), nrow(Dx)), nrow(Dy))

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

  nx <- ncol(X)
  ny <- ncol(Y)

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
                  fm   <- survival::survreg(formula = time ~ group,
                                          data = x,
                                          dist = "exponential")
                  est  <- exp(-fm$coefficients[[2]])
                  varD <- exp(-2*fm$coefficients[[2]]) * fm$var[2,2]
                  res <- c(lower = est - sqrt(varD) * qnorm(p = 1 - alpha / 2),
                           upper = est + sqrt(varD) * qnorm(p = 1 - alpha / 2))
                  return(res)
                },
                simplify = TRUE)
  res <- as.data.frame(x = t(res))
  colnames(res) <- paste0("R", c("1","2"))
  return(res)
}
