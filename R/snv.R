#' @title Standard normal variate (SNV)
#'
#' @description The function \code{snv} applys the standard normal variate to the data.
#'
#' @param x x m-by-n data of class \code{data.frame} or \code{matrix} to be filtered. m is the number of samples (observations) and n is the number of variables.
#'
#' @return data applied the standard normal variate
#'
#' @examples
#' dat.snv <- snv(dat)
#'
#' @export

snv <- function(x){

  x.snv <-t(apply(x, 1, function(x) (x - mean(x))/sd(x)))
  colnames(x.snv) <- colnames(x)

  return(x.snv)
}
