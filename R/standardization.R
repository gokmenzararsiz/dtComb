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

std.range <- function(newdata, model, type = TRUE){
  if(type == TRUE){
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- (newdata[ , i] - min(model[ , i])) /
        (max(model[ , i]) - min(model[ , i]))
      
    }
    return(newdata)
    
} else {
  
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- (newdata[ , i] - model$fit$Std[i, 3] /
                          model$fit$Std[i, 4] - model$fit$Std[i, 3])
      
    }
    return(newdata)}
  }


 

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

std.zscore <- function(newdata, model, type = TRUE){

 if(type == TRUE){ 
   
   for (i in 1:ncol(newdata)){

     newdata[ , i] <- (newdata[ , i] - mean(model[ , i])) / sd(model[ , i])

  }
  return(newdata)
   
 }
  else {
    
      
      for (i in 1:ncol(newdata)){
        
        newdata[ , i] <- (newdata[ , i] - model$fit$Std[i, 1] / model$fit$Std[i, 2])
        
      }
      return(newdata)
  }
}


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

std.mean <- function(newdata, model, type = TRUE){
  if (type == TRUE){
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- newdata[ , i] / mean(model[ , i])
      
    }
    return(newdata)
    }
  else {
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- newdata[ , i] / model$fit$Std[i, 1]
      
    }
    return(newdata)}
  }
  

  

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

std.deviance <- function(newdata, model, type = TRUE){
  if (type == TRUE){
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- newdata[ , i] / sd(model[ , i])
      
    }
    return(newdata)
    
  } 
  else {
    for (i in 1:ncol(newdata)){
      
      newdata[, i] <- newdata[ , i] / model$fit$Std[i, 2]
      
    }
    return(newdata)
  }

  }


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

std.tscore <- function(newdata, model,type = TRUE){

  if(type == TRUE){

    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- 10 * ((newdata[ , i] - mean(model[ , i]))
                              / sd(model[ , i])) + 50
      
    }
    return(newdata)
    
  }
  else{

    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- 10 * ((newdata[ , i] - model$fit$Std[i, 1])
                             / model$fit$Std[i, 2]) + 50
      
    }
    return(newdata)
  }
 }



std <- function(newdata, model, standardize = NULL, type = TRUE) {

  if (any(standardize == "range")){
  
   newdata <- std.range(newdata, model,type)
  
  }
  else if (any(standardize == "zScore")){
  
    newdata <- std.zscore(newdata, model,type)
  
  }
  else if (any(standardize == "tScore")){
    
    newdata <- std.tscore(newdata, model,type)
  
  }
  else if (any(standardize == "mean")){
  
    newdata <- std.mean(newdata, model,type)
  
  }
  else if (any(standardize == "deviance")){
  
    newdata <- std.deviance(newdata, model,type)
  
  }
  
  return(newdata)
  }


