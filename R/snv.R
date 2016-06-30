#' @title Standard normal variate (SNV)
#'
#' @description The function \code{snv} applys the standard normal variate to the spectral data.
#'
#' @param x \code{data.frame} or \code{matrix} with spectra in columns.
#' @param range range of \code{x} used for calculating mean value and sd value. Specify as the column numbers.
#' @return spectral data applied the standard normal variate
#'
#' @examples
#' data(peach)
#' nir.snv <- snv(peach$NIR)
#'
#' @name plsrauto
#' @docType package
#' @export

snv <- function(x, range = c(1, ncol(x))){

  x.snv.all <- c()
  for(i in 1:nrow(x)){
    x.use <- as.numeric(x[i, c(range[1]:range[2])])
    x.snv <- (x[i,] - mean(x.use))/sd(x.use)
    x.snv.all <- rbind(x.snv.all, x.snv)
  }
  colnames(x.snv.all) <- colnames(x)
  rownames(x.snv.all) <- rownames(x)

  return(x.snv.all)
}
