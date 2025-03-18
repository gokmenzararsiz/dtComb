#' @title Standardization according to the chosen method.
#'
#' @description The \code{std.train} Standardization (range, zScore etc.) can be
#' estimated from the training data and applied to any dataset with the same
#' variables.
#'
#' @param data a \code{numeric} data frame of biomarkers
#'
#' @param standardize a \code{character} string indicating the name of the
#' standardization method. The default option is no standardization applied.
#' Available options are:
#' \itemize{
#' \item \bold{Z-score} \code{(zScore)}: This method scales the data to have a mean
#' of 0 and a standard deviation of 1. It subtracts the mean and divides by the standard
#'  deviation for each feature. Mathematically,
#'  \deqn{ Z-score = \frac{x - (\overline x)}{sd(x)}}
#'
#'   where \eqn{x} is the value of a marker, \eqn{\overline{x}} is the mean of the marker and \eqn{sd(x)} is the standard deviation of the marker.
#' \item \bold{T-score} \code{(tScore)}: T-score is commonly used
#' in data analysis to transform raw scores into a standardized form.
#'  The standard formula for converting a raw score \eqn{x} into a T-score is:
#'  \deqn{T-score = \Biggl(\frac{x - (\overline x)}{sd(x)}\times 10 \Biggl) +50}
#'   where \eqn{x} is the value of a marker, \eqn{\overline{x}} is the mean of the marker
#'    and \eqn{sd(x)} is the standard deviation of the marker.
#'
#' \item \bold{Range (a.k.a. min-max scaling)} \code{(range)}: This method transforms data to
#' a specific range, between 0 and 1. The formula for this method is:
#' \deqn{Range = \frac{x - min(x)}{max(x) - min(x)}}
#'
#' \item \bold{Mean} \code{(mean)}: This method, which helps
#' to understand the relative size of a single observation concerning
#' the mean of dataset, calculates the ratio of each data point to the mean value
#' of the dataset.
#' \deqn{Mean =  \frac{x}{\overline{x}}}
#' where \eqn{x} is the value of a marker and \eqn{\overline{x}} is the mean of the marker.
#'
#' \item \bold{Deviance} \code{(deviance)}: This method, which allows for
#' comparison of individual data points in relation to the overall spread of
#' the data, calculates the ratio of each data point to the standard deviation
#' of the dataset.
#' \deqn{Deviance = \frac{x}{sd(x)}}
#' where \eqn{x} is the value of a marker and \eqn{sd(x)} is the standard deviation of the marker.
#' }
#'
#' @return A \code{numeric} data.frame of standardized biomarkers
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' # call data
#' data(laparoscopy)
#'
#' # define the function parameters
#' markers <- laparoscopy[, -1]
#' markers2 <- std.train(markers, "deviance")
#'
#' @export


std.train <- function(data, standardize = NULL) {
  std <- matrix(, 2, 4)
  colnames(std) <- c("mean", "sd", "min", "max")

  for (j in 1:2) {
    std[, j]
    for (i in 1:ncol(data)) {
      std[i, ] <- cbind(
        mean(data[, i]), sd(data[, i]),
        min(data[, i]), max(data[, i])
      )
    }
  }

  if (any(standardize == "range")) {
    for (i in 1:ncol(data)) {
      data[, i] <- (data[, i] - min(data[, i])) /
        (max(data[, i]) - min(data[, i]))
    }
  } else if (any(standardize == "zScore")) {
    for (i in 1:ncol(data)) {
      data[, i] <- (data[, i] - mean(data[, i])) / sd(data[, i])
    }
  } else if (any(standardize == "tScore")) {
    for (i in 1:ncol(data)) {
      data[, i] <- (10 * ((data[, i] - mean(data[, i]))
      / sd(data[, i]))) + 50
    }
  } else if (any(standardize == "mean")) {
    for (i in 1:ncol(data)) {
      data[, i] <- data[, i] / mean(data[, i])
    }
  } else if (any(standardize == "deviance")) {
    for (i in 1:ncol(data)) {
      data[, i] <- data[, i] / sd(data[, i])
    }
  }
  std.model <- list(
    data = data,
    std = std
  )

  return(std.model)
}



#' @title Standardization according to the training model parameters.
#'
#' @description The \code{std.test} Standardization parameters will be taken
#' from the fitted training model and applied to the new data set.
#'
#' @param newdata a \code{numeric} data frame of biomarkers
#'
#' @param model a \code{list} of parameters from the output of linComb,
#' nonlinComb, mlComb or mathComb functions.
#
#' @return A \code{numeric} dataframe of standardized biomarkers
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz



std.test <- function(newdata, model) {
  if (any(model$fit$Standardize == "range")) {
    for (i in 1:ncol(newdata)) {
      newdata[, i] <- ((newdata[, i] - model$fit$Std.model[i, 3]) /
        (model$fit$Std.model[i, 4] - model$fit$Std.model[i, 3]))
    }
  } else if (any(model$fit$Standardize == "zScore")) {
    for (i in 1:ncol(newdata)) {
      newdata[, i] <- (newdata[, i] - model$fit$Std.model[i, 1]) /
        model$fit$Std.model[i, 2]
    }
  } else if (any(model$fit$Standardize == "tScore")) {
    for (i in 1:ncol(newdata)) {
      newdata[, i] <- (10 * ((newdata[, i] - model$fit$Std.model[i, 1]) /
        model$fit$Std.model[i, 2])) + 50
    }
  } else if (any(model$fit$Standardize == "mean")) {
    for (i in 1:ncol(newdata)) {
      newdata[, i] <- newdata[, i] / model$fit$Std.model[i, 1]
    }
  } else if (any(model$fit$Standardize == "deviance")) {
    for (i in 1:ncol(newdata)) {
      newdata[, i] <- newdata[, i] / model$fit$Std.model[i, 2]
    }
  }

  return(newdata)
}
