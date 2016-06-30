#' @title Bland-Altman plot
#'
#' @param y actual values
#' @param yhat predicted values
#' @param sample Name of samples
#' @param nsd a numeric value multiplied to the standard deviations. The product is used for determining the limits of agreement.
#' @param ... additional arguments passed to the \code{\link{plot}}.
#'
#' @return Bland-Altman plot
#'
#' @examples
#' data(peach)
#' peach.pls <- plsr(Brix ~ NIR, ncomp = 5, data = peach, validation = "CV")
#' yhat <- peach.pls$fitted.values[, , 5]
#' baplot(peach$Brix, yhat, sample = rownames(peach), nsd = 3)
#'
#' @name plsrauto
#' @docType package
#' @export
#'
baplot <- function(y, yhat, sample=seq(1,length(y),by=1), nsd=3, ...){

  residual <- y - yhat
  ylim <- c( -1.5*max(abs(residual)), 1.5*max(abs(residual)) )
  plot((y+yhat)/2, residual, ylim=ylim, xlab="(y + y.pred)/2", ylab="residual(y - y.pred)", ...)
  title(paste("Bland-Altman plot (border: Mean+-", nsd, "SD)", sep=""))
  abline(h = mean(residual)+nsd*sd(residual))
  abline(h = mean(residual)-nsd*sd(residual))
  text((y+yhat)/2, residual, labels = sample, adj = c(0.7, 1.5), cex = 0.7)
}
