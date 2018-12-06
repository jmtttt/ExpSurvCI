#' Confidence Interval for Hazard Ratio as approximation for the Ratio of Parameters for Two Exponentially Distributed Samples
#'
#' @import stats
#'
#' @param X numeric matrix -- Observed survival times in first group
#' @param Y numeric matrix -- Observed survival times in second group
#' @param Dx numeric matrix -- Observation indicator (1 for events, 0 for censorings) for first group
#' @param Dy numeric matrix -- Observation indicator (1 for events, 0 for censorings) for second group
#' @param alpha confidence level
#'
#' @return a data.frame with columns C1 (lower CI) and C2 (upper CI)
#' @export
#'
KI_coxph <- function(X, Y, Dx, Dy,
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
                            fm <- survival::coxph(formula = time ~ group,
                                        data = x)
                            return(unlist(exp(confint(object = fm,
                                                      level = 1 - alpha))))
                          }, simplify = TRUE)
  res <- as.data.frame(x = t(res))
  colnames(res) <- paste0("C", c("1","2"))
  return(res)
}
