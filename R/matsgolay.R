#' @title Apply Savitzky-Golay smoothing or 1st/2nd derivative to the matrix or dataframe
#'
#' @description The function \code{matsgolay} applys Savitzky-Golay smoothing, 1st derivative or 2nd derivative to
#' an object of class \code{matrix} or \code{data.frame} by using \link{\code{sgolayfilt}} in 'signal' package.
#'
#' @param x x m-by-n data of class \code{data.frame} or \code{matrix} to be filtered. m is the number of samples (observations) and n is the number of variables.
#' @param ... additional arguments passed to the \code{\link{sgolayfilt}} in 'signal' package.
#'
#' @return an object of class \code{list} containing the smoothed, 1st derivative and 2nd derivative data
#'
#' @seealso \code{\link{sgolayfilt}}
#'
#' @examples
#' # Apply a Savitzky-Golay smoothing & 1st derivative to a matirx \code{dat}
#' dat.sg1d <- mysgolay(dat, p = 2, n = 11, m = 1)
#'
#' @export

matsgolay <- function(x, ...){

  x.sg <- x

  # Smoothing
  for(i in 1:nrow(x)){
    x.sg[i,] <- sgolayfilt(as.numeric(x[i,]), ...)
  }
  colnames(x.sg) <- colnames(x)

  return(x.sg)
}
