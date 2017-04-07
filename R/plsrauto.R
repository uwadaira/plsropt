#' @title PLS regressions under different combinations of X-variable range and preprocessing method
#'
#' @description The function \code{plsrauto} performs the partial least squares (PLS) regressions under different combinations of X-variable range and preprocessing method automatically.
#'
#' @param formula	a model formula like \code{y ~ x}. See \code{\link{plsr}} for a detailed description.
#' @param data an optional data frame with the data to fit the model from.
#' @param testdata data set for prediction.
#' @param xrange an object of class \code{list} which contains ranges of X-variables (see below).
#' @param p filter order for Savitzky-Golay smoothing (default value is 2).
#' @param n filter length (window size) for Savitzky-Golay smoothing (must be odd. default value is 11).
#' @param output if \code{TRUE}, the results are output as PDF and CSV files in the 'PLSR_auto' directory.
#' @param ... additional arguments passed to \code{\link{plsr}} in 'pls' package.
#'
#' @details Three steps of preprocessing are automatically applied to the X-variable data set.
#' First, standard normal variate (SNV) is applied or not.
#' Second, Savitzky-Golay smoothing, 1st derivative or 2nd derivative is applied or not, respectively.
#' Finally, auto-scaling is applied or not.
#' Total 16 (2*4*2) kinds of preprocessing methods are applied and the results of PLS regression are returned as an object of class \code{data.frame}.
#' If \code{output == TRUE}, each PLS regression result is output as PDF and CSV files in one directory.
#'
#' @return an object of class \code{data.frame} containing the statistics of PLS regressions under the different combinations of X-variable range and preprocessing method is returned.
#'
#' @seealso \code{\link{plsrPlot}}, \code{\link{plsr}}
#'
#' @references
#' Vignette \url{https://www.gitbook.com/book/uwadaira/plsropt_vignette_ver1-2-0}
#'
#' @examples
#' data(peach)
#' datTrain <- peach[1:50, ]
#' datTest  <- peach[51:74, ]
#' result.all <- plsrauto(Brix ~ NIR, data = datTrain, testdata = datTest, xrange = list(c(700, 1098), c(1100, 2498)))
#'
#' @import signal
#' @export

plsrauto <- function(formula = NULL, data = NULL, testdata = NULL,
                     yTrain = NULL, xTrain = NULL, xTest = NULL, yTest = NULL, yname = NULL,
                     xrange = NULL, p = 2, n = 11, maxcomp = 10, plot = FALSE, output = FALSE, ...){

  if(!is.null(formula)){
    if(is.null(data)) stop("'data' is not specified.")
    if(is.null(rownames(data))) stop("Set the sample name as the row name of 'data'.")
    dTrain <- model.frame(formula, data = data, na.action = na.omit)
    yTrain <- dTrain[[1]]
    xTrain <- dTrain[[2]]
    sampleTrain <- rownames(xTrain)
  }else{
    if(is.null(yTrain)) stop("'yTrain' is not specified.")
    if(is.null(xTrain)) stop("'xTrain' is not specified.")
    if(is.null(rownames(xTrain))) stop("Set the sample name as the row name of 'xTrain'.")
    yTrain <- yTrain
    xTrain <- as.matrix(xTrain)
    sampleTrain <- rownames(xTrain)

    # Remove the observations of [yTrain = NA]
    sampleTrain <- sampleTrain[!is.na(yTrain)]
    xTrain <- xTrain[!is.na(yTrain), ]
    yTrain <- yTrain[!is.na(yTrain)]
  }

  # Create a directory for storing the resulting files
  if(output == TRUE){
    if(is.null(yname)){
      if(!is.null(formula)){
        dir1 <- paste("./PLSR_auto", terms(formula)[[2]], sep = "/")
      }else{
        stop("'yname' must be specified.")
      }
    }else{
      dir1 <- paste("./PLSR_auto", yname, sep = "/")
    }
    dir.create(dir1, showWarnings = FALSE, recursive = TRUE)
  }

  if(!is.null(testdata) || !is.null(xTest)){ # In the case of using a test set
    if(!is.null(testdata)){
      if(is.null(rownames(testdata))) stop("Set the sample name as the row name of 'testdata'.")
      dTest <- model.frame(formula, data = testdata, na.action = na.omit)
      yTest <- dTest[[1]]
      xTest <- dTest[[2]]
      sampleTest <- rownames(xTest)
    }else if(!is.null(xTest)){
      if(is.null(yTest)) stop("'yTest' is not specified.")
      if(is.null(rownames(xTest))) stop("Set the sample name as the row name of 'xTest'.")
      sampleTest <- rownames(xTest)
      yTest <- yTest
      xTest <- as.matrix(xTest)

      # Remove the observations of [yTest = NA]
      sampleTest <- sampleTest[!is.na(yTest)]
      xTest <- xTest[!is.na(yTest), ]
      yTest <- yTest[!is.na(yTest)]
    }

    set <- c(rep("train", nrow(xTrain)), rep("test", nrow(xTest)))
    xAll <- rbind(xTrain, xTest)
    colnames(xAll) <- colnames(xTrain)
  }else{
    xAll <- xTrain
  }

  ### Run PLS regressions automatically
  result.all <- c()

  if(is.null(xrange)) xrange <- list(range(as.numeric(colnames(xTrain))))

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
            xUse <- scale(xUse, center = TRUE, scale = TRUE)
            prename3 <- paste(prename2, k, sep = " + ")
          }else{
            prename3 <- prename2
          }

          if(!is.null(testdata) || !is.null(xTest)){
            xTrain  <- subset(xUse, set == "train")
            xTest   <- subset(xUse, set == "test")
            datTrain <- data.frame(y = yTrain, x = I(as.matrix(xTrain)))
            datTest  <- data.frame(y = yTest,  x = I(as.matrix(xTest)))
            rownames(datTest)  <- sampleTest
          }else{
            xTrain   <- xUse
            datTrain <- data.frame(y = yTrain, x = I(as.matrix(xTrain)))
            datTest  <- NULL
          }
          rownames(datTrain) <- sampleTrain

          if(output == TRUE){
            dir <- paste(dir1, rname, prename3, sep="/")
            plot <- TRUE
          }else dir <- NULL

          result <- plsrPlot(y ~ x, data = datTrain, testdata = datTest, maxcomp = maxcomp, plot = plot, return.stats=TRUE, dir = dir, output = output, ...)
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

  if(output == TRUE){
    dir.create(dir1, showWarnings = F, recursive = T)
    write.csv(result.all, paste(dir1, fname, sep = "/"))
  }

  return(result.all)
}
