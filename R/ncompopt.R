#' @title Determine the optimun number of components to include in the model
#'
#' @description The function \code{ncompopt} determines the optimum number of components to include the model.
#'
#' @param x RMSECV values calculated by the partial least squares regression
#'
#' @details If \code{ncomp = "auto"}, the number of components with the minimum RMSECV is selected.
#' If there is no component with the minimum RMSECV, the maximum number of components is selected.
#' Then, if the last component does not improve the RMSECV by at least 2 percent, the last one excluded and so forth.
#'
#' @return the optimun number of components
#'
#' @seealso \code{\link{rowpeak}}
#'
#' @examples
#' ncomp <- ncompopt(x)
#'
#' @export

ncompopt <- function(x){
    row.peak <- rowpeak(x)
    if(length(row.peak$low) == 0){
      ncomp.opt <- length(x)
    }else if(length(row.peak$low) == 1){
      ncomp.opt  <- row.peak$low[1]
    }else{
      for(i in 1:(length(row.peak$low) - 1)){
        diff.x <- x[row.peak$low[i]] - x[row.peak$low[i+1]]
        if(diff.x > 0 && diff.x >= x[row.peak$low[i]]*0.02){
          ncomp.opt <- row.peak$low[i+1]
        }else{
          ncomp.opt <- row.peak$low[i]
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
    return(ncomp.opt)
}
