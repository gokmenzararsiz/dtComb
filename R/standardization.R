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
  
  return(data)
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
      
      newdata[ , i] <- ((newdata[ , i] - model$fit$Std[i, 3]) /
                          (model$fit$Std[i, 4] - model$fit$Std[i, 3]))
      
    }
    
  }
  else if (any(model$fit$Standardize == "zScore")){
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- (newdata[ , i] - model$fit$Std[i, 1]) /  model$fit$Std[i, 2]
      
    }
    
  }
  else if (any(model$fit$Standardize == "tScore")){
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- (10 * ((newdata[ , i] - model$fit$Std[i, 1])
                              / model$fit$Std[i, 2])) + 50
      
    }
    
  }
  else if (any(model$fit$Standardize == "mean")){
    
    for (i in 1:ncol(newdata)){
      
      newdata[ , i] <- newdata[ , i] / model$fit$Std[i, 1]
      
    }
    
  }
  else if (any(model$fit$Standardize == "deviance")){
    
    for (i in 1:ncol(newdata)){
      
      newdata[, i] <- newdata[ , i] / model$fit$Std[i, 2]
      
    }
    
  }
  
  return(newdata)
}
