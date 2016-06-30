#' @title Apply Savitzky-Golay smoothing or 1st/2nd derivative to the matrix
#'
#' @description The function \code{matsgolay} applys Savitzky-Golay smoothing, 1st derivative or 2nd derivative to
#' an object of class \code{matrix} by using \code{\link{sgolayfilt}} in \code{signal} package.
#'
#' @param x x m-by-n data of class \code{data.frame} or \code{matrix} to be filtered. m is the number of samples (observations) and n is the number of variables.
#' @param ... additional arguments passed to the \code{\link{sgolayfilt}} in 'signal' package.
#'
#' @return an fitted matrix
#'
#' @seealso \code{\link{sgolayfilt}}
#'
#' @examples
#' data(peach)
#' nir.1d <- matsgolay(peach$NIR, p = 2, n = 11, m = 1)
#'
#' @name plsrauto
#' @docType package
#' @import signal
#' @export

matsgolay <- function(x, ...){

  x.sg <- x

  for(i in 1:nrow(x)){
    x.sg[i,] <- sgolayfilt(as.numeric(x[i,]), ...)
  }
  colnames(x.sg) <- colnames(x)

  return(x.sg)
}
