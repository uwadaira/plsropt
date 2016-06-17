#' @title Draw important plots for the partial least squares regression at once
#'
#' @description The function \code{plsrPlot} performs the partial least squares (PLS) regression by using the function \code{\link{plsr}} in 'pls' package and draws four kinds of plots
#' (cross-validation, variable importance [regression coefficient and variable importance in projection], loading, actual vs. predicted value) at once.
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
#' See \code{\link{plsr}} and \code{\link{lm}} for a detailed description.
#'
#' If \code{ncomp = "auto"}, the optimum number of components is automatically selected (see \code{\link{ncompopt}}).
#'
#' @return If \code{return.stats = TRUE}, an object of class \code{data.frame} containing the statistics of PLS regression is returned.
#' Otherwise, an object of class \code{mvr} is returned (see \code{\link{plsr}}).
#'
#' @seealso \code{\link{ncompopt}}, \code{\link{plsr}}, \code{\link{lm}}
#'
#' @examples
#' result <- plsrPlot(y = dat$Brix, x = x, ncomp = "auto", output = TRUE)
#'
#' result <- plsrPlot(y = dat.train$Brix, x = x.train,
#'                    y.test=dat.test$Brix, x.test=x.test, ncomp = 5)
#'
#' @export

plsrPlot <- function(formula = NULL, data = NULL, testdata = NULL,
                     yTrain = NULL, xTrain = NULL, xTest = NULL, yTest = NULL, yname = NULL,
                     ncomp = "auto", maxcomp = 10, plot = TRUE, validation = "CV", segment.type ="interleaved",
                     output = FALSE, dir = NULL, return.stats = FALSE, ...){

  if(!is.null(formula)){
    if(is.null(data)) stop("data is not specified")
    d <- model.frame(formula, data = data)
    y <- d[[1]]
    x <- d[[2]]
    x <- as.matrix(x)
  }else{
    if(is.null(yTrain)) stop("yTrain is not specified")
    y <- yTrain
    if(is.null(xTrain)) stop("xTrain is not specified")
    x <- as.matrix(xTrain)
  }

  # Remove the observation of [y = NA]
  x <- x[!is.na(y), ]
  y <- y[!is.na(y)]

  xvar <- as.numeric(colnames(x))

  # Maximum number of components
  if(maxcomp > length(y) - 2){maxcomp <- length(y) - 2}

### Perform PLS regression

  result <- plsr(y ~ x, ncomp=maxcomp, method="oscorespls", ...)
  vip <- VIP(result)

  # Number of components to be included in the model
  if(ncomp == "auto"){
    RMSECV <- sqrt(c(result$validation$PRESS0, result$validation$PRESS)/(length(y)-1))
    ncomp.opt <- ncompopt(RMSECV)
    if(ncomp.opt == 1){
      ncomp <- ncomp.opt
    }else{
      ncomp <- ncomp.opt - 1
    }
  }

  result.ncomp <- plsr(y ~ x, ncomp=ncomp, method="oscorespls", ...)

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

  if(!is.null(testdata) || !is.null(xTest)){

    if(!is.null(testdata)){
      d.test <- model.frame(formula, data = testdata)
      y.test <- d.test[[1]]
      x.test <- d.test[[2]]
      x.test <- as.matrix(x.test)
    }else{
      if(is.null(yTest)) stop("yTest is not specified")
      y.test <- yTest
      x.test <- as.matrix(xTest)
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
    plot(result, "validation", main="Cross-validation", legendpos = "topright")

    # Plot predicted values vs. actual values
    if(is.null(testdata)){  # Training set
      xylim <- c(min(y, yhat.cal, yhat.val), max(y, yhat.cal, yhat.val))
      plot(y, yhat.cal, type="p", pch=16, col="blue", asp=1, xlim=xylim, ylim=xylim, xlab="actual", ylab="predicted", main=paste("Actual vs. Predicted (", ncomp, " comps, Training set)", sep=""))
      abline(a=0,b=1)
      par(new = TRUE)
      plot(y, yhat.val, type="p", pch=16, col="red", asp=1, xlim=xylim, ylim=xylim, axes=F, xlab="", ylab="", main="")
      legend("topleft", legend=c("Calibration", "Validation"), pch=16, col=c("blue", "red"), bg="white")
    }else{  # Test set
      xylim <- c(min(y, yhat.cal, y.test, yhat.test), max(y, yhat.cal, y.test, yhat.test))
      plot(y, yhat.cal, pch=16, col="blue", asp=1, xlim=xylim, ylim=xylim, xlab="actual value", ylab="predicted value", main=paste("Actual vs. Predicted (", ncomp, " comps, Test set)", sep=""))
      abline(a=0,b=1)
      par(new = TRUE)
      plot(y.test, yhat.test, pch=16, col="red", asp=1, xlim=xylim, ylim=xylim, axes=F, main="", xlab="", ylab="")
      legend("topleft", legend=c("Calibration", "Prediction"), pch=16, col=c("blue", "red"), bg="white")
      box()
    }

    # Regression coefficient & VIP plot
    plot(xvar, result$coefficients[, , ncomp], type="l", col="blue", xlab="variable", ylab="", main=paste("Variable importance (", ncomp, " comps)", sep=""))
    abline(h=0, lty=2, col="blue")
    par(new=T)
    plot(xvar, vip[, ncomp], type="l", col="red", axes=FALSE, xlab="", ylab="", main="")
    axis(4)
    abline(h=1, lty=2, col="red")
    legend("topleft", legend=c("Regression coefficient", "VIP (right axis)"), lty=1, col=c("blue", "red"), bty="n")

    # Loading plot
    colsld <- seq(1, ncomp, by=1)
    par(mar=c(5,4,4,5) + 0.1)
    matplot(xvar, loading.weights(result)[,1:ncomp], type="l", lty=1, col=colsld, xlab="variable", ylab="", main="Loading")
    abline(h=0, lty=2)
    par(xpd=T)
    label <- seq(1, ncomp, by=1)
    legend(par()$usr[2], par()$usr[4], legend=label, col=colsld, lty=1, ncol=1)
    par(mar=c(5,4,4,2) + 0.1)
    par(xpd=F)
    par(mfrow=c(1,1)) # reset panel number

    if(output == TRUE) dev.off()
  }

  ### Save the result
  if(output == TRUE){
    write.csv(stats,
              paste(dir, "stats.csv", sep="/"), row.names=FALSE)
    write.csv(data.frame(y, yhat.cal, yhat.val),
              paste(dir, "fittedvalue.csv", sep="/"), row.names=FALSE)
    write.csv(coef(result, ncomp=ncomp, intercept=TRUE),
              paste(dir, "regcoef.csv", sep="/"))
    write.csv(vip,
              paste(dir, "vip.csv", sep="/"))
    write.csv(loading.weights(result)[,1:ncomp],
              paste(dir,"loading.csv", sep="/"))

    if(!is.null(testdata)){
      write.csv(data.frame(y.test, yhat.test),
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
