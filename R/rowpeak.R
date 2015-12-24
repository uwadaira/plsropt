#' @title Number of rows with high and low peak value
#'
#' @param x an object of class \code{vector}
#'
#' @return an object of class \code{list} containing the number of rows which have high and low peak values
#'
#' @examples
#' rowpeak(x)
#'
#' @export

rowpeak <- function(x){

  row.peak.high <- c()
  row.peak.low  <- c()
  i <- 1
  if(x[i+1]-x[i] >= 0){
    while(i <= (length(x)-1)){
      while(x[i] <= x[i+1]){
        i <- i + 1
        if(i == (length(x))){break}
      }
      if(i == (length(x))){break
      }else{row.peak.high <- c(row.peak.high, i)}

      while(x[i] >= x[i+1]){
        i <- i + 1
        if(i == (length(x))){break}
      }
      if(i == (length(x))){break
      }else{row.peak.low <- c(row.peak.low, i)}
    }
  }else{
    while(i <= (length(x)-1)){
      while(x[i] >= x[i+1]){
        i <- i + 1
        if(i == (length(x))){break}
      }
      if(i == (length(x))){break
      }else{row.peak.low <- c(row.peak.low, i)}

      while(x[i] <= x[i+1]){
        i <- i + 1
        if(i == (length(x))){break}
      }
      if(i == (length(x))){break
      }else{row.peak.high <- c(row.peak.high, i)}
    }
  }
  return(list(high=row.peak.high, low=row.peak.low))
}
