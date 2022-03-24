# TODO: Add comment
#
# Author: serra
###############################################################################
#' @title Standardization with respect to the t score.
#'
#' @description The \code{std.train} standardizes using T scores. The range
#' varies between usually 20 and 80
#'
#' @param data a \code{numeric} dataframe of biomarkers
#' 
#' @param standardize a a \code{string} cdscfvf
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
#' markers <- exampleData1[, -1]
#' markers2 <- std.train(markers, tScore)
#'
#' @export


std.train <- function(data, standardize = NULL) {
  
  std = matrix(,2,4)
  colnames(std) <- c("mean", "sd", "min", "max")
  
  for (j in 1:2) {
    
    std[, j]
    for (i in 1:ncol(data)) {
      
      std[i, ] = cbind(mean(data[, i]),sd(data[, i]), 
                       min(data[, i]),max(data[, i]))
      
    }
  }

  if (any(standardize == "range")){
  
    for (i in 1:ncol(data)){
      
      data[ , i] <- (data[ , i] - min(data[ , i])) /
        (max(data[ , i]) - min(data[ , i]))
      
    }
  
  }
  else if (any(standardize == "zScore")){
  
    for (i in 1:ncol(data)){
      
      data[ , i] <- (data[ , i] - mean(data[ , i])) / sd(data[ , i])
      
    }
  
  }
  else if (any(standardize == "tScore")){
    
    for (i in 1:ncol(data)){
      
      data[ , i] <- (10 * ((data[ , i] - mean(data[ , i]))
                              / sd(data[ , i]))) + 50
      
    }
  
  }
  else if (any(standardize == "mean")){
  
    for (i in 1:ncol(data)){
      
      data[ , i] <- data[ , i] / mean(data[ , i])
      
    }
  
  }
  else if (any(standardize == "deviance")){
  
    for (i in 1:ncol(data)){
      
      data[ , i] <- data[ , i] / sd(data[ , i])
      
    }
  
  }
  return_model <- list(data = data,
                       std = std)
  
  return(return_model)
  }



#' @title Standardization with respect to the t score.
#'
#' @description The \code{std.test} standardizes using T scores. The range
#' varies between usually 20 and 80
#'
#' @param newdata a \code{numeric} dataframe of biomarkers
#' 
#' @param model a \code{string} dcdcdcdz
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
#' markers <- exampleData1[, -1]
#' markers2 <- std.test(markers, score1)
#'
#' @export

std.test <- function(newdata, model) {
  
  if (any(model$fit$Standardize == "range")){
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- ((newdata[ , i] - std.model$std[i, 3]) /
                          (std.model$std[i, 4] - std.model$std[i, 3]))
      
    }
    
  }
  else if (any(model$fit$Standardize == "zScore")){
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- (newdata[ , i] - std.model$std[i, 1]) /  std.model$std[i, 2]
      
    }
    
  }
  else if (any(model$fit$Standardize == "tScore")){
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- (10 * ((newdata[ , i] - std.model$std[i, 1])
                              / std.model$std[i, 2])) + 50
      
    }
    
  }
  else if (any(model$fit$Standardize == "mean")){
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- newdata[ , i] / std.model$std[i, 1]
      
    }
    
  }
  else if (any(model$fit$Standardize == "deviance")){
    
    for (i in 1:ncol(newdata)){
      
      newdata[, i] <- newdata[ , i] / std.model$std[i, 2]
      
    }
    
  }
  
  return(newdata)
}
