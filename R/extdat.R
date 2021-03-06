#' @title Extract the specific range of column from a data set.
#'
#' @param x an object of class \code{data.frame} or \code{matrix}.
#' @param start name of the 1st column of new data set to be extracted.
#' @param end name of the last column of new data set to be extracted.
#'
#' @return a new data set of class \code{data.frame} or \code{matrix}
#'
#' @examples
#' data(peach)
#' x <- extdat(peach$NIR, start = 700, end = 1050)
#'
#' @export

extdat <- function(x, start, end){

  x2 <- x[, which(colnames(x)==start) : which(colnames(x)==end)]

  return(x2)
}
