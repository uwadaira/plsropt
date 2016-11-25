#' @title Calculates selectivity ratio for a given regression model
#'
#' @param mvr an \code{mvr} object. The fitted model built by \code{\link{plsr}} function (see \code{\link{plsr}}).
#' @param ncomp a positive integer. The number of components included in the model.
#'
#' @return  selectivity ratio for each variable included in the model
#'
#' @author Yasuhiro Uwadaira (supported by Vipavee Trivittayasil)
#'
#' @examples
#' data(peach)
#' peach.pls <- plsr(Brix ~ NIR, ncomp = 5, data = peach, validation = "CV")
#' selectivityRatio(peach.pls)
#'
#' @export

selectivityRatio <- function(mvr, ncomp = mvr$ncomp){

  # X-variable
  X <- mvr$model[[2]]
  X <- apply(X, 2, function(x) x - mean(x)) # mean-centering
  if(!is.null(mvr$scale)){
    X <- apply(X, 2, function(x) x/sd(x)) # scaling
  }

  # projection of the rows of X onto the normalized regression coefficients vector
  # t_TP is proportional on the predicted values, yhat
  b <- as.matrix(mvr$coefficients[, 1, ncomp])
  t_TP <- X %*% (b/norm(b))

  # loading, p_TP, are obtained by projecting the columns of X onto the obtained score vector, t_TP
  p_TP <- t(X) %*% (t_TP/c(t(t_TP) %*% t_TP))

  SR <- c()
  for (i in 1:ncol(X)){
    E_TPi <- as.matrix(X[, i]) - as.matrix(p_TP[i] * t_TP)
    V_ex <- var(t_TP %*% p_TP[i]) # explained variance
    V_re <- var(E_TPi) # residual variance

    if(is.nan(V_ex/V_re)){
      SR <- c(SR, 0)
    }else SR <- c(SR, V_ex/V_re)
  }
  return(SR)
}
