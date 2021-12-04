# TODO: Add comment
#
# Author: serra
###############################################################################
#' @title Standardization with respect to range.
#'
#' @description The \code{std.range} standardizes to a range between 0 and 1
#'
#' @param markers a \code{numeric} data frame of biomarkers
#'
#' @return A \code{numeric} data frame of standardized biomarkers
#'
#' @author Serra Bersan Gengec, Ilayda Serra Yerlitas
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' markers2 <- std.range(markers)
#'
#' @export

std.range <- function(markers){

  for (i in 1:ncol(markers)){

   markers[ , i] <- (markers[ , i] - min(markers[ , i])) /
     (max(markers[ , i]) - min(markers[ , i]))

  }
  return(markers)}

#' @title Standardization with respect to z score.
#'
#' @description The \code{std.zscore} standardizes using z scores with
#' mean = 0 and standard deviation = 1
#'
#' @param markers a \code{numeric} data frame of biomarkers
#'
#' @return A \code{numeric} data frame of standardized biomarkers
#'
#' @author Serra Bersan Gengec, Ilayda Serra Yerlitas
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' markers2 <- std.zscore(markers)
#'
#' @export

std.zscore <- function(markers){

  for (i in 1:ncol(markers)){

    markers[ , i] <- (markers[ , i] - mean(markers[ , i])) / sd(markers[ , i])

  }
  return(markers)}


#' @title Standardization with respect to the sample mean.
#'
#' @description The \code{std.mean} standardizes with sample mean = 1
#'
#' @param markers a \code{numeric} data frame of biomarkers
#'
#' @return A \code{numeric} data frame of standardized biomarkers
#'
#' @author Serra Bersan Gengec, Ilayda Serra Yerlitas
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' markers2 <- std.mean(markers)
#'
#' @export

std.mean <- function(markers){

  for (i in 1:ncol(markers)){

    markers[ , i] <- markers[ , i] / mean(markers[ , i])

  }
  return(markers)}

#' @title Standardization with respect to the sample standard deviation.
#'
#' @description The \code{std.deviance} standardizes with sample standard
#'   deviation = 1
#'
#' @param markers a \code{numeric} data frame of biomarkers
#'
#' @return A \code{numeric} data frame of standardized biomarkers
#'
#' @author Serra Bersan Gengec, Ilayda Serra Yerlitas
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' markers2 <- std.deviance(markers)
#'
#' @export

std.deviance <- function(markers){

  for (i in 1:ncol(markers)){

    markers[ , i] <- markers[ , i] / sd(markers[ , i])

  }
  return(markers)}


#' @title Standardization with respect to the t score.
#'
#' @description The \code{std.tscore} standardizes using T scores. The range
#' varies between usually 20 and 80
#'
#' @param markers a \code{numeric} dataframe of biomarkers
#'
#' @return A \code{numeric} dataframe of standardized biomarkers
#'
#' @author Serra Bersan Gengec, Ilayda Serra Yerlitas
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' markers2 <- std.tscore(markers)
#'
#' @export

std.tscore <- function(markers){

  for (i in 1:ncol(markers)){

    markers[ , i] <- 10 * ((markers[ , i] - mean(markers[ , i]))
                           / sd(markers[ , i])) + 50

  }
  return(markers)}