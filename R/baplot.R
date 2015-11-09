#' @title Bland-Altman plot
#'
#' @param y actual values
#' @param yhat predicted values
#' @param sample Name of samples
#' @param nsd Number of standard deviations for determining the limits of agreement.
#'
#' @return Bland-Altman plot
#'
#' @examples
#' baplot(y, yhat, nsd = 3)
#'
#' @export
#'
baplot <- function(y, yhat, sample=seq(1,length(y),by=1), nsd){

  residual <- y - yhat
  ylim <- c( -1.5*max(abs(residual)), 1.5*max(abs(residual)) )
  plot((y+yhat)/2, residual, ylim=ylim, xlab="(y + y.pred)/2", ylab="residual(y - y.pred)")
  title(paste("Bland-Altman plot (border: Mean+-", nsd, "SD)", sep=""))
  abline(h = mean(residual)+nsd*sd(residual))
  abline(h = mean(residual)-nsd*sd(residual))
  text((y+yhat)/2, residual, labels = sample, adj = c(0.7, 1.5), cex = 0.7)
}
