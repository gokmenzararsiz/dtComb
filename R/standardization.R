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
#' @author Serra Ilayda  Yerlitas, Serra Bersan Gengec 
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

std.range <- function(testMark, trainMark){

  for (i in 1:ncol(testMark)){

    testMark[ , i] <- (testMark[ , i] - min(trainMark[ , i])) /
     (max(trainMark[ , i]) - min(trainMark[ , i]))

  }
  return(testMark)}

#' @title Standardization with respect to z score.
#'
#' @description The \code{std.zscore} standardizes using z scores with
#' mean = 0 and standard deviation = 1
#'
#' @param markers a \code{numeric} data frame of biomarkers
#'
#' @return A \code{numeric} data frame of standardized biomarkers
#'
#' @author Serra Ilayda  Yerlitas, Serra Bersan Gengec 
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

std.zscore <- function(testMark, trainMark){

  for (i in 1:ncol(testMark)){

    testMark[ , i] <- (testMark[ , i] - mean(trainMark[ , i])) / sd(trainMark[ , i])

  }
  return(testMark)}


#' @title Standardization with respect to the sample mean.
#'
#' @description The \code{std.mean} standardizes with sample mean = 1
#'
#' @param markers a \code{numeric} data frame of biomarkers
#'
#' @return A \code{numeric} data frame of standardized biomarkers
#'
#' @author Serra Ilayda  Yerlitas, Serra Bersan Gengec 
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

std.mean <- function(testMark, trainMark){

  for (i in 1:ncol(testMark)){

    testMark[ , i] <- testMark[ , i] / mean(trainMark[ , i])

  }
  return(testMark)}

#' @title Standardization with respect to the sample standard deviation.
#'
#' @description The \code{std.deviance} standardizes with sample standard
#'   deviation = 1
#'
#' @param markers a \code{numeric} data frame of biomarkers
#'
#' @return A \code{numeric} data frame of standardized biomarkers
#'
#' @author Serra Ilayda  Yerlitas, Serra Bersan Gengec 
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

std.deviance <- function(testMark, trainMark){

  for (i in 1:ncol(testMark)){

    testMark[ , i] <- testMark[ , i] / sd(trainMark[ , i])

  }
  return(testMark)}


#' @title Standardization with respect to the t score.
#'
#' @description The \code{std.tscore} standardizes using T scores. The range
#' varies between usually 20 and 80
#'
#' @param markers a \code{numeric} dataframe of biomarkers
#'
#' @return A \code{numeric} dataframe of standardized biomarkers
#'
#' @author Serra Ilayda  Yerlitas, Serra Bersan Gengec 
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

std.tscore <- function(testMark, trainMark){

  for (i in 1:ncol(testMark)){

    testMark[ , i] <- 10 * ((testMark[ , i] - mean(trainMark[ , i]))
                           / sd(trainMark[ , i])) + 50

  }
  return(testMark)}



std <- function(data1, data2, standardize) {
  
 if (any(standardize == "range")){
  
    data1 <- std.range(data1, data2)
  
  }
  else if (any(standardize == "zScore")){
  
    data1 <- std.zscore(data1, data2)
  
  }
  else if (any(standardize == "tScore")){
  
    data1 <- std.tscore(data1, data2)
  
  }
  else if (any(standardize == "mean")){
  
    data1 <- std.mean(data1, data2)
  
  }
  else if (any(standardize == "deviance")){
  
    data1 <- std.deviance(data1, data2)
  
  }
  return(data1)
}  
