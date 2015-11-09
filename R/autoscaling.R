#' @title Auto-scaling
#'
#' @param x m-by-n data of class \code{data.frame} or \code{matrix}. m is the number of samples (observations) and n is the number of variables.
#'
#' @return Auto-scaled data
#'
#' @examples
#' dat.as <- autoscaling(dat)
#'
#' @export

autoscaling <- function(x){

  x.as <- apply(x, 2, function(x) (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))
  colnames(x.as) <- colnames(x)

  return(x.as)
}
