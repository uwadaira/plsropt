#' @title Determine the optimal number of components to include in the model
#'
#' @description The function \code{ncompopt} determines the optimal number of components to include the model.
#'
#' @param x RMSECV values calculated by the partial least squares regression
#'
#' @details If \code{ncomp = "auto"}, the number of components with the minimum RMSECV is selected.
#' If the selected n-th component with the minimum RMSECV does not improve the RMSECV by at least 2 percent, then (n-1)th compoment is selected and so forth.
#'
#' @return the optimal number of components
#'
#' @seealso \code{\link{rowpeak}}
#'
#' @examples
#' data(peach)
#' peach.pls <- plsr(Brix ~ NIR, ncomp = 15, data = peach, validation = "CV")
#' ncompopt(RMSEP(peach.pls)$val[2,,])
#'
#' @name plsropt
#' @docType package
#' @export

ncompopt <- function(x){
  row.peak <- rowpeak(x)
  if(length(row.peak$low) == 0){
    ncomp.opt <- length(x)
  }else if(length(row.peak$low) == 1){
    ncomp.opt  <- row.peak$low[1]
  }else{
    ncomp.opt <- row.peak$low[1]
    for(i in 1:(length(row.peak$low) - 1)){
      diff.x <- x[ncomp.opt] - x[row.peak$low[i+1]]
      if(diff.x > 0 && diff.x >= x[row.peak$low[i]]*0.02){
        ncomp.opt <- row.peak$low[i+1]
      }
    }
  }
  if(ncomp.opt > 1){
    ncomp.opt2 <- ncomp.opt - 1
    diff.x2 <- x[ncomp.opt2] - x[ncomp.opt]
    while(diff.x2 < x[ncomp.opt]*0.02){
      ncomp.opt <- ncomp.opt2
      if(ncomp.opt == 1){ break }
      ncomp.opt2 <- ncomp.opt - 1
      diff.x2 <- x[ncomp.opt2] - x[ncomp.opt]
    }
  }
  if(ncomp.opt >= 2) ncomp.opt <- ncomp.opt - 1

  return(ncomp.opt)
}
