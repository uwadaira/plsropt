#' @title PLS regressions under different combinations of X-variable range and preprocessing method
#'
#' @description The function \code{plsrauto} performs the partial least squares (PLS) regressions under different combinations of X-variable range and preprocessing method automatically.
#'
#' @param formula	a model formula like \code{y ~ x}. See \code{\link{plsr}} and \code{\link{lm}} for a detailed description.
#' @param testdata data set for prediction.
#' @param xrange an object of class \code{list} which contains ranges of X-variables (see below).
#' @param p filter order for Savitzky-Golay smoothing (default value is 2).
#' @param n filter length (window size) for Savitzky-Golay smoothing (must be odd. default value is 11).
#' @param ... additional arguments passed to \code{\link{plsrPlot}} and \code{\link{plsr}} in 'pls' package.
#'
#' @details Three steps of preprocessing are automatically applied to the X-variable data set.
#' First, standard normal variate (SNV) is applied or not.
#' Second, Savitzky-Golay smoothing, 1st derivative or 2nd derivative is applied or not, respectively.
#' Finally, auto-scaling is applied or not.
#' Total 16 (2*4*2) kinds of preprocessing methods are applied and the results of PLS regressions are returned as an object of class \code{data.frame} and output as a CSV file.
#' If \code{output == TRUE}, PLS regression result per combination of a X-variable range and a preprocessing method is output as PDF and CSV files in one directory.
#'
#' @return an object of class \code{data.frame} containing the statistics of PLS regressions under the different combinations of X-variable range and preprocessing method is returned.
#'
#' @seealso \code{\link{plsrPlot}}, \code{\link{plsr}}, \code{\link{lm}}
#'
#' @examples
#' # import CSV file
#' datAll <- read.csv(file = "./peachNIR.csv", row.names = 1, check.names = F)
#'
#' datAll[1:5, 1:7]
#' Brix firmness       700       702       704       706       708
#' ak01 12.4     1.37 0.1264722 0.1206582 0.1162688 0.1129761 0.1105212
#' ak02 13.1     0.83 0.0894380 0.0878140 0.0866534 0.0858628 0.0853562
#' ak03 12.7     1.28 0.0887944 0.0860848 0.0840838 0.0826256 0.0815916
#' ak04 10.2     1.50 0.1637433 0.1554227 0.1490126 0.1440929 0.1403139
#' ak05 11.4     1.00 0.0983719 0.0939621 0.0906616 0.0882200 0.0864351
#'
#' # extract X-variables
#' x <- extdat(datAll, start = 700, end = 2498)
#'
#' # merge Y-variables and X-variables
#' dat <- data.frame(datAll[, 1:2], NIR = I(x))
#'
#' # divide the data set into training and test data set
#' datTrain <- dat[1:50, ]
#' datTest  <- dat[51:74, ]
#'
#' result.all <- plsrauto(Brix ~ NIR, data = datTrain, testdata = datTest,
#'                        xrange = list(c(700, 1098), c(1100, 2498)))
#'
#' @export

plsrauto <- function(formula, data, testdata=NULL, xrange=NULL, p=2, n=11, ...){

  dir1 <- paste("./PLSR_auto", terms(formula)[[2]], sep = "/")
  dir.create(dir1, showWarnings = FALSE, recursive = TRUE)

  dTrain <- model.frame(formula, data = data)
  yTrain <- dTrain[[1]]
  xTrain <- dTrain[[2]]

  if(!is.null(testdata)){
    dTest <- model.frame(formula, data = testdata)
    yTest <- dTest[[1]]
    xTest <- dTest[[2]]
  }

  if(is.null(testdata)){
    xAll <- xTrain
  }else{
    set <- c(rep("train", nrow(xTrain)), rep("test", nrow(xTest)))
    xAll <- rbind(xTrain, xTest)
    colnames(xAll) <- colnames(xTrain)
  }

  ### Run PLS regressions automatically
  result.all <- c()

  if(is.null(xrange)) xrange <- list(c(colnames(x)[1], colnames(x)[ncol(x)]))

  # Preprocessing methods
  sgpara <- paste("(wsize-", n, "pt_forder-", p, ")", sep = "")
  preproc1 <- c("Raw", "SNV")
  preproc2 <- c("Raw", paste(c("SGSM", "SG1D","SG2D"), sgpara, sep = ""))
  preproc3 <- c("Raw", "Auto-scaling")

  for(h in 1:length(xrange)){
    xAllSub <- extdat(xAll, start = xrange[[h]][1], end = xrange[[h]][2])
    rname <- paste(xrange[[h]][1], xrange[[h]][2], sep = "-")

    for(i in preproc1){
      if(i == "SNV") xAllSub <- snv(xAllSub)
      xAllSub.sgsm <- matsgolay(x = xAllSub, p = p, n = n, m = 0)
      xAllSub.sg1d <- matsgolay(x = xAllSub, p = p, n = n, m = 1)
      xAllSub.sg2d <- matsgolay(x = xAllSub, p = p, n = n, m = 2)
      prename1 <- i

      for(j in preproc2){
        if(j == preproc2[1]){
          xUse <- xAllSub
        }else if(j == preproc2[2]){
          xUse <- xAllSub.sgsm
        }else if(j == preproc2[3]){
          xUse <- xAllSub.sg1d
        }else if(j == preproc2[4]){
          xUse <- xAllSub.sg2d
        }
        if(j != preproc2[1]){
          prename2 <- paste(prename1, j, sep = " + ")
        }else{
          prename2 <- prename1
        }

        for(k in preproc3){
          if(k == "Auto-scaling"){
            xUse <- autoscaling(xUse)
            prename3 <- paste(prename2, k, sep = " + ")
          }else{
            prename3 <- prename2
          }

          if(is.null(testdata)){
            xTrain   <- xUse
            datTrain <- data.frame(y = yTrain, x = I(as.matrix(xTrain)))
            datTest  <- NULL
          }else{
            xTrain  <- subset(xUse, set == "train")
            xTest   <- subset(xUse, set == "test")
            datTrain <- data.frame(y = yTrain, x = I(as.matrix(xTrain)))
            datTest  <- data.frame(y = yTest,  x = I(as.matrix(xTest)))
          }

          dir <- paste(dir1, rname, prename3, sep="/")
          result <- plsrPlot(y ~ x, data = datTrain, testdata = datTest, return.stats=TRUE, dir = dir, ...)
          result.all <- rbind.data.frame(result.all, data.frame(Xrange=rname, Preprocessing=prename3, result))
        }
      }
    }
  }

  # Sort by the correlation coefficient for validation set or test set
  if(is.null(testdata)){
    sortlist <- result.all$R.val
    fname <- "result_train.csv"
  }else{
    sortlist <- result.all$R.test
    fname <- "result_test.csv"
  }
  result.all <- result.all[sort.list(sortlist, decreasing = TRUE), ]
  rownames(result.all) <- seq(1, nrow(result.all), by = 1)
  write.csv(result.all, paste(dir1, fname, sep = "/"))

  return(result.all)
}
