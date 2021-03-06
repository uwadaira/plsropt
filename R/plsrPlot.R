#' @title Perform PLS regression and draw important plots at once
#'
#' @description The function \code{plsrPlot} performs the partial least squares (PLS) regression by using the function \code{\link{plsr}} in 'pls' package and draws four kinds of plots
#' (cross-validation, regression vector, variable importance [selectivity ratio and variable importance in projection], actual vs. predicted value) at once.
#'
#' @param formula	a model formula (see below).
#' @param data an optional data frame with the data to fit the model form.
#' @param testdata data set for prediction.
#' @param ncomp the number of components (latent variables) to include in the model (see below).
#' @param maxcomp Maximum number of components (latent variables) to be used for cross-validation.
#' @param plot if \code{TRUE}, the four plots are drawn.
#' @param output if \code{TRUE}, the results are output as PDF and CSV files.
#' @param dir path to the directory where the results are output (default value is \code{"PLSR_result"}).
#' @param return.stats if \code{TRUE}, an object of class \code{data.frame} containing the statistics of PLS regression is returned.
#' @param ... additional arguments passed to the \code{\link{plsr}} in 'pls' package.
#'
#' @details The \code{formula} argument should be a symbolic formula of the form \code{response ~ terms},
#' where \code{response} is the name of the response vector
#' and \code{terms} is the name of the predictor matrix, e.g., water ~ FTIR.
#' See \code{\link{plsr}} for a detailed description.
#'
#' If \code{ncomp = "auto"}, the optimum number of components is automatically selected (see \code{\link{ncompopt}}).
#'
#' @return If \code{return.stats = TRUE}, an object of class \code{data.frame} containing the statistics of PLS regression is returned.
#' Otherwise, an object of class \code{mvr} is returned (see \code{\link{plsr}}).
#'
#' @references
#' Vignette \url{https://www.gitbook.com/book/uwadaira/plsropt_vignette_ver1-2-0}
#'
#' @seealso \code{\link{ncompopt}}, \code{\link{plsr}}
#'
#' @examples
#'
#' # alternative way
#' result <- plsrPlot(yTrain = datTrain$Brix, xTrain = datTrain$NIR, yTest=datTest$Brix, xTest=datTest$NIR)
#'
#' @import pls
#' @export

plsrPlot <- function(formula = NULL, data = NULL, testdata = NULL,
                     yTrain = NULL, xTrain = NULL, xTest = NULL, yTest = NULL, yname = NULL,
                     ncomp = "auto", maxcomp = 10, plot = TRUE, validation = "CV", segment.type = "interleaved",
                     output = FALSE, dir = NULL, return.stats = FALSE, ...){

  if(!is.null(formula)){
    if(is.null(data)) stop("data is not specified")
    d <- model.frame(formula, data = data)
    y <- d[[1]]
    x <- d[[2]]
    if(is.null(rownames(data))) stop("Set the sample name as the row name of data.")
    sampleTrain <- rownames(d)
  }else{
    if(is.null(yTrain)) stop("yTrain is not specified.")
    y <- yTrain
    if(is.null(xTrain)) stop("xTrain is not specified.")
    x <- xTrain
    if(is.null(rownames(xTrain))) stop("Set the sample name as the row name of xTrain.")
    sampleTrain <- rownames(xTrain)

    # Remove the observations of [y = NA]
    sampleTrain <- sampleTrain[!is.na(y)]
    x <- x[!is.na(y), ]
    y <- y[!is.na(y)]
  }
  x <- as.matrix(x)

  xvar <- as.numeric(colnames(x))
  if(is.na(range(xvar)[1])) xvar <- 1:ncol(x)

  # Number of components considered
  if(maxcomp > length(y) - 2){maxcomp <- length(y) - 2}

  # Boolean for executing prediction
  if(!is.null(testdata) || !is.null(xTest)){
    pred <- TRUE
  }else{
    pred <- FALSE
  }

### Perform PLS regression

  result <- plsr(y ~ x, ncomp=maxcomp, method="oscorespls", validation = validation, segment.type = segment.type, ...)
  vip <- VIP(result)

  # Number of components used in the model
  if(ncomp == "auto"){
    ncomp <- ncompopt(RMSEP(result)$val[2,,])
  }

  result.ncomp <- plsr(y ~ x, ncomp=ncomp, method="oscorespls", validation = validation, segment.type = segment.type, ...)

  y <- result$model$y
  yhat.cal <- result$fitted.values[, , ncomp]
  result.cal.lm <- lm(yhat.cal ~ y)
  Slope.cal <- result.cal.lm$coefficients[2]

  yhat.val <- result$validation$pred[, , ncomp]
  result.val.lm <- lm(yhat.val ~ y)
  Slope.val <- result.val.lm$coefficients[2]

  Bias.val <- mean(y - yhat.val)
  SEC      <- sqrt(sum((y - yhat.cal)^2)/length(y))
  SECV     <- sqrt(sum((y - yhat.val - Bias.val)^2)/(length(y)-1))
  RPD.val  <- sd(y)/SECV

  stats <- data.frame(N=length(y), ncomp, R.cal=cor(y, yhat.cal), Slope.cal, SEC, R.val=cor(y, yhat.val), Slope.val, SECV, Bias.val, RPD.val)

### Prediction with test set

  if(pred == TRUE){
    if(!is.null(testdata)){
      d.test <- model.frame(formula, data = testdata)
      y.test <- d.test[[1]]
      x.test <- d.test[[2]]
      x.test <- as.matrix(x.test)
      if(is.null(rownames(testdata))) stop("Set the sample name as the row name of testdata.")
      sampleTest <- rownames(testdata)
    }else{
      if(is.null(yTest)) stop("yTest is not specified")
      y.test <- yTest
      x.test <- as.matrix(xTest)
      if(is.null(rownames(xTest))) stop("Set the sample name as the row name of xTest.")
      sampleTest <- rownames(xTest)
    }

    # Remove the observation with [y.test = NA]
    x.test <- x.test[!is.na(y.test), ]
    y.test <- y.test[!is.na(y.test)]

    yhat.test <- predict(result, ncomp = ncomp, newdata = as.matrix(x.test))

    result.test.lm <- lm(yhat.test ~ y.test)
    Slope.test <- result.test.lm$coefficients[2]

    Bias.test <- mean(y.test - yhat.test)
    SEP       <- sqrt(sum((y.test - yhat.test-Bias.test)^2)/(length(y.test)-1))
    RPD.test  <- sd(y.test)/SEP

    stats <- data.frame(stats, N.test=length(y.test), R.test=cor(y.test, yhat.test), Slope.test, SEP, Bias.test, RPD.test)
  }

### Plot graphics

  if(plot == TRUE){
    if(output == TRUE){
      if(is.null(dir)){
        if(!is.null(formula)){
          yname <- terms(formula)[[2]]
        }else{
          if(is.null(yname)) stop("yname must be specified")
          yname <- yname
        }
        dir <- paste("./PLSR_result", yname, sep = "/")
      }
      if(is.null(testdata)){
        subdir <- paste(ncomp, "comps(n=", nrow(x), ")_train", sep = "")
      }else{
        subdir <- paste(ncomp, "comps(n=", nrow(x), ")_test", sep = "")
      }
      dir <- paste(dir, subdir, sep = "/")
      dir.create(dir, showWarnings = FALSE, recursive = TRUE)
      pdf(paste(dir, "PLS_plots.pdf", sep="/"), width=9, height=9)
    }

    par(mfrow=c(2,2))

    # Cross-validation plot
    plot(result, "validation", main="Cross-validation", legendpos = "topright", type = "b", lty = 1, pch = 16, col = c("royal blue", "red"), cex = 1.2)

    # Regression vector
    plot(result, "coef", ncomp = ncomp, labels = "numbers", main=paste0("Regression vector (", ncomp, " comps)"))
    abline(h=0, lty=2, col="red")
    axis(1, tck = 1, col = "lightgrey", lty = "dotted")
    box()

    # Predicted values vs. actual values
    if(pred == TRUE){
      x.range <- range(y, y.test)
      y.range <- range(yhat.cal, yhat.test)
      plot(result$model$y, yhat.cal, xlim = x.range, ylim = y.range, asp = 1, pch = 21, bg = 3, cex = 1.2, axes = F, xlab = "", ylab = "")
      par(new=T)
      plot(y.test, yhat.test, xlim = x.range, ylim = y.range, asp = 1, pch = 21, bg = 2, cex = 1.2)
      abline(a = 0, b = 1)
      legend("topleft", legend = c("Calibration", "Prediction"), pch = 16, col = c(3, 2))
      legend("bottomright", legend = as.expression(bquote(italic(paste({R^2}, "test") == .(round(cor(y.test, yhat.test)^2, 2))))), cex = 1.2, bty = "n")
    }else{
      x.range <- range(y)
      y.range <- range(yhat.cal, yhat.val)
      plot(y, yhat.cal, xlim = x.range, ylim = y.range, asp = 1, pch = 21, bg = 3, cex = 1.2, axes = F, xlab = "", ylab = "")
      par(new=T)
      plot(y, yhat.val, xlim = x.range, ylim = y.range, asp = 1, pch = 21, bg = 2, cex = 1.2)
      abline(a = 0, b = 1)
      abline(a = 0, b = 1)
      legend("topleft", legend = c("Calibration", "Cross-valudation"), pch = 16, col = c(3, 2))
      legend("bottomright", legend = as.expression(bquote(italic(paste({R^2}, "cal") == .(round(cor(y, yhat.cal)^2, 2))))), cex = 1.2, bty = "n")
    }

    # VIP & Selectivity ratio
    if(xvar[1] > xvar[length(xvar)]){
      x.lim <- rev(range(xvar))
    }else{
      x.lim <- range(xvar)
    }
    sratio <- selectivityRatio(result, ncomp = ncomp)
    plot(xvar, sratio, type="l", xlim=x.lim, col="royal blue", xlab="variable", ylab="", main=paste0("Variable importance (", ncomp, " comps)"))
    axis(2, col = "royal blue")
    par(new=T)
    plot(xvar, vip[, ncomp], type="l", xlim=x.lim, col="red", axes=FALSE, xlab="", ylab="", main="")
    axis(4, col = "red")
    abline(h=1, lty=2, col="red")
    axis(1, tck = 1, col = "lightgrey", lty = "dotted")
    box()
    legend("topleft", legend=c("Selectivity ratio", "VIP (right axis)"), lty=1, col=c("royal blue", "red"), bty="n")

    par(mar=c(5,4,4,2) + 0.1)
    par(xpd=F)
    par(mfrow=c(1,1)) # reset panel number

    if(output == TRUE) dev.off()
  }

  ### Save the result
  if(output == TRUE){
    write.csv(stats,
              paste(dir, "stats.csv", sep="/"), row.names=FALSE)
    write.csv(data.frame(sample=sampleTrain, y, yhat.cal, yhat.val),
              paste(dir, "fittedvalue.csv", sep="/"), row.names=FALSE)
    write.csv(coef(result, ncomp=ncomp, intercept=TRUE),
              paste(dir, "regcoef.csv", sep="/"))
    write.csv(vip,
              paste(dir, "vip.csv", sep="/"))
    write.csv(sratio,
              paste(dir, "selectivityRatio.csv", sep="/"))
    write.csv(loading.weights(result)[,1:ncomp],
              paste(dir,"loading.csv", sep="/"))

    if(pred == TRUE){
      write.csv(data.frame(sample=sampleTest, y.test, yhat.test),
                paste(dir, "fittedvalue_test.csv", sep="/"), row.names=FALSE)
    }

    # Save Bland-Altman plot for outlier detection
    pdf(paste(dir, "BA_plot.pdf", sep="/"), width=6, height=5)
    baplot(y, yhat.val, sample=seq(1,length(y),by=1), nsd=3)
    dev.off()
  }

  if(return.stats == TRUE){return(stats)
  }else return(result.ncomp)
}
