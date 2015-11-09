#' @title Savitzky-Golay smoothing, 1st/2nd derivative
#'
#' @description The function \code{matsgolay} applys Savitzky-Golay smoothing, 1st/2nd derivative to the data by using \link{\code{sgolayfilt}} in 'signal' package.
#'
#' @param x x m-by-n data of class \code{data.frame} or \code{matrix} to be filtered. m is the number of samples (observations) and n is the number of variables.
#' @param p filter order for Savitzky-Golay smoothing
#' @param n filter length (window size) for Savitzky-Golay smoothing (must be odd. default value is 11)
#'
#' @return an object of class \code{list} containing the smoothed, 1st derivative and 2nd derivative data
#'
#' @seealso \code{\link{sgolayfilt}}
#'
#' @examples
#' x.sg <- mysgolay(x, p = 2, n = 11)
#'
#' @export

matsgolay <- function(x, p = 2, n = 11){

  X.sgsm <- x
  X.sg1d <- x
  X.sg2d <- x

  # Smoothing
  for(i in 1:nrow(x)){
    X.sgsm[i,] <- sgolayfilt(as.numeric(x[i,]), p=p, n=n, m=0)
  }
  colnames(X.sgsm) <- colnames(x)

  # 1st derivative
  for(i in 1:nrow(x)){
    X.sg1d[i,] <- sgolayfilt(as.numeric(x[i,]), p=p, n=n, m=1)
  }
  colnames(X.sg1d) <- colnames(x)

  # 2nd derivative
  for(i in 1:nrow(x)){
    X.sg2d[i,] <- sgolayfilt(as.numeric(x[i,]), p=p, n=n, m=2)
  }
  colnames(X.sg2d) <- colnames(x)

  return(list(window.size=n, sm=X.sgsm, first.d=X.sg1d, second.d=X.sg2d))
}
