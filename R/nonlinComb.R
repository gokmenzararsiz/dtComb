# TODO: Add comment
#
# Author: serra ilayda yerlitas
###############################################################################
#' @title Combine two diagnostic tests with several non-linear combination methods.
#'
#' @description The \code{nonlinComb} function calculates the combination
#' scores of two diagnostic tests selected among several non-linear combination
#' methods and standardization options
#'
#' @param markers a \code{numeric} data frame that includes two diagnostic tests
#'  results
#'
#' @param status a \code{factor} vector that includes the actual disease
#'  status of the patients
#'
#' @param event a \code{character} string that indicates the event in the status
#'  to be considered as positive event
#'
#' @param method a \code{character} string specifying the method used for
#'  combining the markers. The available methods are:
#'  \itemize{
#'  \item \code{polyreg}: The method builds a logistic regression model with the
#'   feature space created and returns the probability of a positive event for 
#'   each observation.
#'  \item \code{ridgereg}: Ridge regression is a shrinkage method used to 
#'  estimate the coefficients of highly correlated variables and in this case 
#'  the polynomial feature space created from two biomarkers. For the 
#'  implementation of the method, glmnet library [8] is used with two functions:
#'  cv.glmnet() to run a cross validation  model to determine the tuning 
#'  parameter λ and glmnet() to fit the model with the selected tuning parameter. 
#'  \item \code{lassoreg}: Lasso regression is also a shrinkage method with one
#'  difference is that at the end this method returns the coefficients of some 
#'  features as 0, makes this method useful for feature elimination as well. 
#'  The implementation is similar to Ridge regression, cross validation for 
#'  parameter selection and model fit are implemented with glmnet library.
#'  \item \code{elasticreg}: Elastic Net regression is obtained by combining the 
#'  penalties of Ridge and Lasso regression to get the best of both models. The 
#'  model again includes a tuning parameter λ as well as a mixing parameter α 
#'  taken form the user which takes a value between 0 (ridge) and 1 (lasso) to 
#'  determine the weights of the loss functions of Ridge and Lasso regressions.
#'  \item \code{splines}: With the applications of regression models in a 
#'  polynomial feature space the second non-linear approach to combining 
#'  biomarkers comes from applying several regression models to the dataset 
#'  using a function derived from piecewise polynomials. Splines are implemented
#'  with degrees of freedom and degrees of the fitted polynomials taken from the 
#'  user. For the implementation splines library is used to build piecewise 
#'  logistic regression models with base splines.
#'  \item \code{sgam}: In addition to the basic spline structure, Generalized 
#'  Additive Models are applied with smoothing splines using the gam 
#'  library R.
#'  \item \code{nsgam}: In addition to the basic spline structure, Generalized 
#'  Additive Models are applied with natural cubic splines using the gam 
#'  library R.
#'  }
#' 
#' @param degree1 a \code{numeric} value for polynomial based methods indicates 
#'  the degree of the feature space created for marker 1, for spline based 
#'  methods the degree of the fitted polynomial between each node for marker 1. 
#'  (3, default)
#' 
#' @param degree2 a \code{numeric} value for polynomial based methods indicates 
#'  the degree of the feature space created for marker 2, for spline based 
#'  methods the degree of the fitted polynomial between each node for marker 2
#'  (3, default)
#' 
#'  @param df1 a \code{numeric} value that indicates the number of knots as the
#'   degrees of freedom in spline based methods for marker 1 (4, default)
#'  
#'  @param df2 a \code{numeric} value that indicates the number of knots as the
#'   degrees of freedom in spline based methods for marker 2 (4, default)
#' 
#' 
#' @param standardize a \code{character} string indicating the name of the
#'  standardization method. The default option is no standardization applied.
#'  Available options are:
#'  \itemize{
#'  \item \code{range}: Standardization to a range between 0 and 1
#'  \item \code{zScore}: Standardization using z scores with mean = 0
#'  and standard deviation = 1
#'  \item \code{tScore}: Standardization using T scores. The range varies between
#'  usually 20 and 80
#'  \item \code{mean}: Standardization with sample mean = 1
#'  \item \code{deviance}: Standardization with sample standard deviation = 1
#' }
#' 
#' @param interact a \code{logical} indicator that specifies whether to include
#'  the interaction between the markers to the feature space created for 
#'  polynomial based methods (FALSE, default)
#' 
#' @param alpha a \code{numeric} value as the mixing parameter in Elastic Net 
#' Regression method (0.5, default)
#' 
#' @param direction a \code{character} string determines in which direction the 
#'  comparison will be made.  “>”: if the predictor values for the control group 
#'  are higher than the values of the case group (controls > cases). 
#'  “<”: if the predictor values for the control group are lower or equal than 
#'  the values of the case group (controls < cases). 
#'
#' @param conf.level a \code{numeric} values determines the confidens interval
#'  for the roc curve(0.95, default).
#' 
#' @param cutoff.method  a \code{character} string determines the cutoff method
#'  for the roc curve.
#' 
#' @return A list of \code{numeric} nonlinear combination scores calculated
#'  according to the given method and standardization option
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec
#'
#' @examples
#' data("exampleData1")
#' data <- exampleData1
#'
#' markers <- data[, -1]
#' status <- factor(data$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#' direction <- "<"
#' cutoff.method <- "youden"
#' 
#' score1 <- nonlinComb(markers = markers, status = status, event = event,
#' method = "ridgereg", resample = "boot", include.interact = TRUE, 
#' cutoff.method = "youden", standardize = "zScore")
#'  
#'  
#' 
#' score3 <- nonlinComb(markers = markers, status = status, event = event, 
#' method = "polyreg", resample = "cv", include.interact = FALSE, direction = "<", 
#' cutoff.method = "youden", standardize = "zScore")

nonlinComb <- function(markers = NULL, status = NULL, event = NULL,
                    method = c("polyreg", "ridgereg", "lassoreg", "elasticreg",
                               "splines", "sgam", "nsgam"),
                    degree1 = 3, degree2 = 3, df1 = 4, df2 = 4,
                    resample = c("none", "cv", "repeatedcv", "boot"),
                    nfolds = 5, nrepeats = 3, niters = 10,
                    standardize = c("none", "range","zScore", "tScore", "mean", 
                                    "deviance"),
                    include.interact = FALSE, alpha = 0.5,  
                    direction = c("auto" ,"<", ">"), conf.level = 0.95, 
                    cutoff.method = c("youden", "roc01")){
  match.arg(method)
  match.arg(direction)
  match.arg(cutoff.method)
  
  if (!is.data.frame(markers)) {
    markers <- as.data.frame(markers)
  }
  
  for(i in 1:ncol(markers)) if(!is.numeric(markers[, i]))
    stop("at least one variable is not numeric")
  
  if(!ncol(markers) == 2)
    stop("the number of markers should be 2")
  
  if(!is.factor(status)) status <- as.factor(status)
  
  if(!length(levels(status)) == 2)
    stop("the number of status levels should be 2")
  
  stopifnot(event %in% status)
  levels(status)[levels(status) == "NA"] <- NA
  stopifnot(nrow(markers) == length(status))
  
  comp <- complete.cases(markers)
  markers <- markers[comp, ]
  status <- status[comp]
  
  if (length(method) != 1){
    stop("No method provided")
  }
  
  if(any(resample == "none")){
    
    resample <- "none"
  }
  
  if(any(resample == "cv")){
    
    nrepeats = 1
  }
  
  if(any(standardize == "none")){
    
    standardize <- "none"
  }

  markersData <- markers
  colnames(markersData) <- c("m1", "m2")
  data <- cbind(status,markersData) 
  
  resample_results <- vector(mode = "list", length = 2)
  names(resample_results) <- c("parameters", "AUC")
  repeated_results <- vector(mode = "list", length = 2)
  names(repeated_results) <- c("parameters", "AUC")
  
  if(method == "polyreg"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = data[iters[[i]], ]
        trainMarkBaseStatus = trainMarkBase[, 1]
        trainMarkBase = trainMarkBase[, -1]
        train = std(trainMarkBase, trainMarkBase, standardize)
        
        train = cbind(trainMarkBaseStatus, train)
        colnames(train) <- c("status","m1", "m2")
        
        test = std(data[, -1], trainMarkBase, standardize)
        
        test = cbind(data[, 1], test)    
        colnames(test) <- c("status","m1", "m2")
        
        if(include.interact == TRUE){
          
          interact <- train$m1 * train$m2
          
          model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2) 
                   + interact, data = train, family = binomial(link = "logit"))
        } else {
          
          model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2), 
                       data = train, family = binomial(link = "logit"))
        }
        
        interact <- test$m1 * test$m2
        comb.score <- predict(model, newdata = test, type = "response")
        interact <- data$m1 * data$m2
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newdata = data, 
                                      type = "response")
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = data[-folds[[i]], ]
          trainMarkBaseStatus = trainMarkBase[, 1]
          trainMarkBase = trainMarkBase[, -1]
          train = std(trainMarkBase, trainMarkBase, standardize)
          test = data[folds[[i]], ]
          testStatus = test[, 1]
          test = std(test[, -1], trainMarkBase, standardize)     
          train = cbind(trainMarkBaseStatus, train)
          colnames(train) <- c("status","m1", "m2")
          
          test = cbind(testStatus, test)  
          colnames(test) <- c("status","m1", "m2")
          
          if( include.interact == TRUE){
            
            interact <- train$m1 * train$m2
            
            model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2) 
                    + interact, data = train, family = binomial(link = "logit"))
          } else {
            
            model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2), 
                         data = train, family = binomial(link = "logit"))
          }
          
          interact <- test$m1 * test$m2
          comb.score <- predict(model, newdata = test, type = "response")
          interact <- data$m1 * data$m2
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(test$status, as.numeric(comb.score))))
          
          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- 
          resample_results$parameters[[max_AUC]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newdata = data, 
                                      type = "response")
      
    }
    
    else {
      
      data2 <- std(data[, -1], data[, -1], standardize)
      
      if( include.interact == TRUE){
        
        interact <- data2$m1 * data2$m2
        
        model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2) 
                     + interact, data = data2, family = binomial(link = "logit"))
      } else {
        
        model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2), 
                     data = data2, family = binomial(link = "logit"))
      }
      
      comb.score <- predict(model, newdata = data2, type = "response")
      
      parameters <- model
      
    }
  
  } 
  
  else if (method == "ridgereg"){
    
    if(any(resample== "boot")){
      
      markersData <- markers
      colnames(markersData) <- c("m1", "m2")
      data <- cbind(status,markersData)
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
   
        trainMarkBase = data[iters[[i]], ]
        trainMarkBaseStatus = trainMarkBase[, 1]
        trainMarkBase = trainMarkBase[, -1]
        train = std(trainMarkBase, trainMarkBase, standardize)
        
        train = cbind(trainMarkBaseStatus, train)
        colnames(train) <- c("status","m1", "m2")
        
        test = std(data[, -1], trainMarkBase, standardize)
        
        test = cbind(data[, 1], test)    
        colnames(test) <- c("status","m1", "m2")
        
        if( include.interact == TRUE){
          
          interact <- train$m1 * train$m2
          
          space <- cbind.data.frame(poly(train$m1, degree1),
                                    poly(train$m2, degree2), interact)
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = 0, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                            family = "binomial", lambda = cv.model$lambda.min)
          interact <- test$m1 * test$m2
          
          testspace <- cbind.data.frame(poly(test$m1, degree1),
                                        poly(test$m2, degree2), interact)
          interact <- data$m1 * data$m2
          
          dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                        poly(data$m2, degree2), interact)
        } else {
          
          space <- cbind.data.frame(poly(train$m1, degree1),
                                    poly(train$m2, degree2))
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = 0, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                            family = "binomial", lambda = cv.model$lambda.min) 
          testspace <- cbind.data.frame(poly(test$m1, degree1),
                                        poly(test$m2, degree2))
          dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                        poly(data$m2, degree2))
        } 
        
        comb.score<-predict(model, newx = as.matrix(testspace), type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newx = as.matrix(dataspace), 
                                      type = "response")
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(data$status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = data[-folds[[i]], ]
          trainMarkBaseStatus = trainMarkBase[, 1]
          trainMarkBase = trainMarkBase[, -1]
          train = std(trainMarkBase, trainMarkBase, standardize)
          test = data[folds[[i]], ]
          testStatus = test[, 1]
          test = std(test[, -1], trainMarkBase, standardize)     
          train = cbind(trainMarkBaseStatus, train)
          colnames(train) <- c("status","m1", "m2")
          
          test = cbind(testStatus, test)  
          colnames(test) <- c("status","m1", "m2")
            
          if( include.interact == TRUE){
            
            interact <- train$m1 * train$m2
            
            space <- cbind.data.frame(poly(train$m1, degree1),
                                      poly(train$m2, degree2), interact)
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = 0, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                            family = "binomial", lambda = cv.model$lambda.min)
            
            interact <- test$m1 * test$m2
            
            testspace <- cbind.data.frame(poly(test$m1, degree1),
                                          poly(test$m2, degree2), interact)
            interact <- data$m1 * data$m2
            
            dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                          poly(data$m2, degree2), interact)
          } else {
            
            space <- cbind.data.frame(poly(train$m1, degree1),
                                      poly(train$m2, degree2))
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = 0, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                             family = "binomial", lambda = cv.model$lambda.min)
            
            testspace <- cbind.data.frame(poly(test$m1, degree1),
                                          poly(test$m2, degree2))
            
            dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                          poly(data$m2, degree2))
          } 
          
          status <- test$status
          comb.score<-predict(model, newx = as.matrix(testspace), type="response")
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(status, as.numeric(comb.score))))
          
          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- 
          resample_results$parameters[[max_AUC]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newx = as.matrix(dataspace), 
                                      type = "response")
      
    }
    
    else {
      
      data2 <- std(data[, -1], data[, -1], standardize)
      
      if( include.interact == TRUE){
        
        interact <- data2$m1 * data2$m2
        
        space <- cbind.data.frame(poly(data2$m1, degree1),
                                  poly(data2$m2, degree2), interact)
        cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                      alpha = 0, family = "binomial")
        model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                           family = "binomial", lambda = cv.model$lambda.min)
      } else {
        
        space <- cbind.data.frame(poly(data2$m1, degree1),
                                  poly(data2$m2, degree2))
        cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                      alpha = 0, family = "binomial")
        model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                           family = "binomial", lambda = cv.model$lambda.min) 
      } 
      
      comb.score<-predict(model, newx = as.matrix(space), type="response")
      
      parameters <- model
      
    }
    
  }
  
  else if (method == "lassoreg"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(data$status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = data[iters[[i]], ]
        trainMarkBaseStatus = trainMarkBase[, 1]
        trainMarkBase = trainMarkBase[, -1]
        train = std(trainMarkBase, trainMarkBase, standardize)
        
        train = cbind(trainMarkBaseStatus, train)
        colnames(train) <- c("status","m1", "m2")
        
        test = std(data[, -1], trainMarkBase, standardize)
        
        test = cbind(data[, 1], test)    
        colnames(test) <- c("status","m1", "m2")
        
        if( include.interact == TRUE){
          
          interact <- train$m1 * train$m2
          
          space <- cbind.data.frame(poly(train$m1, degree1),
                                    poly(train$m2, degree2), interact)
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = 1, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                                  family = "binomial", lambda = cv.model$lambda.min)
          interact <- test$m1 * test$m2
          
          testspace <- cbind.data.frame(poly(test$m1, degree1),
                                        poly(test$m2, degree2), interact)
          interact <- data$m1 * data$m2
          
          dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                        poly(data$m2, degree2), interact)
        } else {
          
          space <- cbind.data.frame(poly(train$m1, degree1),
                                    poly(train$m2, degree2))
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = 1, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                                  family = "binomial", lambda = cv.model$lambda.min) 
          testspace <- cbind.data.frame(poly(test$m1, degree1),
                                        poly(test$m2, degree2))
          dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                        poly(data$m2, degree2))
        } 
        
        status <- test$status
        comb.score<-predict(model, newx = as.matrix(testspace), type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- (predict(parameters, newx = as.matrix(dataspace), 
                                      type = "response"))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(data$status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = data[-folds[[i]], ]
          trainMarkBaseStatus = trainMarkBase[, 1]
          trainMarkBase = trainMarkBase[, -1]
          train = std(trainMarkBase, trainMarkBase, standardize)
          test = data[folds[[i]], ]
          testStatus = test[, 1]
          test = std(test[, -1], trainMarkBase, standardize)     
          train = cbind(trainMarkBaseStatus, train)
          colnames(train) <- c("status","m1", "m2")
          
          test = cbind(testStatus, test)  
          colnames(test) <- c("status","m1", "m2")
          
          if( include.interact == TRUE){
            
            interact <- train$m1 * train$m2
            
            space <- cbind.data.frame(poly(train$m1, degree1),
                                      poly(train$m2, degree2), interact)
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = 1, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                                    family = "binomial", lambda = cv.model$lambda.min)
            
            interact <- test$m1 * test$m2
            
            testspace <- cbind.data.frame(poly(test$m1, degree1),
                                          poly(test$m2, degree2), interact)
            interact <- data$m1 * data$m2
            
            dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                          poly(data$m2, degree2), interact)
          } else {
            
            space <- cbind.data.frame(poly(train$m1, degree1),
                                      poly(train$m2, degree2))
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = 1, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                                    family = "binomial", lambda = cv.model$lambda.min)
            
            testspace <- cbind.data.frame(poly(test$m1, degree1),
                                          poly(test$m2, degree2))
            dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                          poly(data$m2, degree2))
          } 
          
          status <- test$status
          comb.score<-predict(model, newx = as.matrix(testspace), type ="response")
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(test$status, as.numeric(comb.score))))
          
          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- 
          resample_results$parameters[[max_AUC]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newx = as.matrix(dataspace), 
                                      type = "response")
      
    }
    
    else {
      
      data2 <- std(data[, -1], data[, -1], standardize)
      
      if( include.interact == TRUE){
        
        interact <- data2$m1 * data2$m2
        
        space <- cbind.data.frame(poly(data2$m1, degree1),
                                  poly(data2$m2, degree2), interact)
        cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                      alpha = 1, family = "binomial")
        model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                           family = "binomial", lambda = cv.model$lambda.min)
      } else {
        
        space <- cbind.data.frame(poly(data2$m1, degree1),
                                  poly(data2$m2, degree2))
        cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                      alpha = 1, family = "binomial")
        model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                           family = "binomial", lambda = cv.model$lambda.min) 
      }
      
      comb.score <- predict(model, newx = as.matrix(space), type="response")
      
      parameters <- model
      
    }
    
  } 
  
  else if (method == "elasticreg"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(data$status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = data[iters[[i]], ]
        trainMarkBaseStatus = trainMarkBase[, 1]
        trainMarkBase = trainMarkBase[, -1]
        train = std(trainMarkBase, trainMarkBase, standardize)
        
        train = cbind(trainMarkBaseStatus, train)
        colnames(train) <- c("status","m1", "m2")
        
        test = std(data[, -1], trainMarkBase, standardize)
        
        test = cbind(data[, 1], test)    
        colnames(test) <- c("status","m1", "m2")
        
        if( include.interact == TRUE){
          
          interact <- data$m1 * data$m2
          
          space <- cbind.data.frame(poly(data$m1, degree1),
                                    poly(data$m2, degree2), interact)
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = alpha, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                             family = "binomial", lambda = cv.model$lambda.min)
          interact <- test$m1 * test$m2
          
          testspace <- cbind.data.frame(poly(test$m1, degree1),
                                        poly(test$m2, degree2), interact)
          interact <- data$m1 * data$m2
          
          dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                        poly(data$m2, degree2), interact)
        } else {
          
          space <- cbind.data.frame(poly(data$m1, degree1),
                                    poly(data$m2, degree2))
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = alpha, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                            family = "binomial", lambda = cv.model$lambda.min) 
          testspace <- cbind.data.frame(poly(test$m1, degree1),
                                        poly(test$m2, degree2))
          dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                        poly(data$m2, degree2))
        }
        
        status <- test$status
        comb.score <- predict(model, newx = as.matrix(testspace), type="response")
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newx = as.matrix(dataspace), 
                                      type = "response")
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(data$status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = data[-folds[[i]], ]
          trainMarkBaseStatus = trainMarkBase[, 1]
          trainMarkBase = trainMarkBase[, -1]
          train = std(trainMarkBase, trainMarkBase, standardize)
          test = data[folds[[i]], ]
          testStatus = test[, 1]
          test = std(test[, -1], trainMarkBase, standardize)     
          train = cbind(trainMarkBaseStatus, train)
          colnames(train) <- c("status","m1", "m2")
          
          test = cbind(testStatus, test)  
          colnames(test) <- c("status","m1", "m2")
          
          if( include.interact == TRUE){
            
            interact <- train$m1 * train$m2
            
            space <- cbind.data.frame(poly(train$m1, degree1),
                                      poly(train$m2, degree2), interact)
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = alpha, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                                    family = "binomial", lambda = cv.model$lambda.min)
            
            interact <- test$m1 * test$m2
            
            testspace <- cbind.data.frame(poly(test$m1, degree1),
                                          poly(test$m2, degree2), interact)
            interact <- data$m1 * data$m2
            
            dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                          poly(data$m2, degree2), interact)
          } else {
            
            space <- cbind.data.frame(poly(train$m1, degree1),
                                      poly(train$m2, degree2))
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = alpha, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                                    family = "binomial", lambda = cv.model$lambda.min)
            
            testspace <- cbind.data.frame(poly(test$m1, degree1),
                                          poly(test$m2, degree2))
            
            dataspace <- cbind.data.frame(poly(data$m1, degree1),
                                          poly(data$m2, degree2))
          } 
          
          status <- test$status
          comb.score<-predict(model, newx = as.matrix(testspace), type="response")
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(test$status, as.numeric(comb.score))))
          
          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- 
          resample_results$parameters[[max_AUC]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC]]
      comb.score <- as.matrix(predict(parameters, newx = as.matrix(dataspace), 
                                      type = "response"))
      
    }
    
    else {
      
      data2 <- std(data[, -1], data[, -1], standardize)
      
      if( include.interact == TRUE){
        
        interact <- data2$m1 * data2$m2
        
        space <- cbind.data.frame(poly(data2$m1, degree1),
                                  poly(data2$m2, degree2), interact)
        cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                      alpha = alpha, family = "binomial")
        model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                          family = "binomial", lambda = cv.model$lambda.min)
      } else {
        
        space <- cbind.data.frame(poly(data2$m1, degree1),
                                  poly(data2$m2, degree2))
        cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                      alpha = alpha, family = "binomial")
        model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                          family = "binomial", lambda = cv.model$lambda.min) 
      }
      
      comb.score <- predict(model, newx = as.matrix(space), type="response")
      
      parameters <- model
      
    }
    
  }
  
  else if (method == "splines"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = data[iters[[i]], ]
        trainMarkBaseStatus = trainMarkBase[, 1]
        trainMarkBase = trainMarkBase[, -1]
        train = std(trainMarkBase, trainMarkBase, standardize)
        
        train = cbind(trainMarkBaseStatus, train)
        colnames(train) <- c("status","m1", "m2")
        
        test = std(data[, -1], trainMarkBase, standardize)
        
        test = cbind(data[, 1], test)    
        colnames(test) <- c("status","m1", "m2")
        
        model <- glm(status ~ splines::bs(m1,degree = degree1, df = df1) + 
                       splines::bs(m2,degree = degree2, df = df2), 
                     data = train, family = binomial)
        
        comb.score <- predict(model,newdata = test,type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newdata = data, 
                                      type = "response")
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = data[-folds[[i]], ]
          trainMarkBaseStatus = trainMarkBase[, 1]
          trainMarkBase = trainMarkBase[, -1]
          train = std(trainMarkBase, trainMarkBase, standardize)
          test = data[folds[[i]], ]
          testStatus = test[, 1]
          test = std(test[, -1], trainMarkBase, standardize)     
          train = cbind(trainMarkBaseStatus, train)
          colnames(train) <- c("status","m1", "m2")
          
          test = cbind(testStatus, test)  
          colnames(test) <- c("status","m1", "m2")
          
          model <- glm(status ~ splines::bs(m1,degree = degree1, df = df1) + 
                         splines::bs(m2,degree = degree2, df = df2), 
                       data = train, family = binomial)
          
          comb.score <- predict(model,newdata = test,type="response")
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(test$status, as.numeric(comb.score))))
          
          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- 
          resample_results$parameters[[max_AUC]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newdata = data,
                                      type = "response")
      
    }
    
    else {
      
      data2 <- std(data[, -1], data[, -1], standardize)
      
      model <- glm(status ~ splines::bs(m1,degree = degree1, df = df1) + 
                     splines::bs(m2,degree = degree2, df = df2), 
                   data = data2, family = binomial)
      
      comb.score <- predict(model,newdata = markers, type="response")
      
      parameters <- model
      
    }
    
  }
  
  else if (method == "sgam"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = data[iters[[i]], ]
        trainMarkBaseStatus = trainMarkBase[, 1]
        trainMarkBase = trainMarkBase[, -1]
        train = std(trainMarkBase, trainMarkBase, standardize)
        
        train = cbind(trainMarkBaseStatus, train)
        colnames(train) <- c("status","m1", "m2")
        
        test = std(data[, -1], trainMarkBase, standardize)
        
        test = cbind(data[, 1], test)    
        colnames(test) <- c("status","m1", "m2")
        
        model <- gam::gam(status ~ gam::s(m1, df = df1) + 
                       gam::s(m2, df = df2), 
                     data = train, family = binomial)
        
        comb.score <- predict(model,newdata = test,type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newdata = data, 
                                      type = "response")
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = data[-folds[[i]], ]
          trainMarkBaseStatus = trainMarkBase[, 1]
          trainMarkBase = trainMarkBase[, -1]
          train = std(trainMarkBase, trainMarkBase, standardize)
          test = data[folds[[i]], ]
          testStatus = test[, 1]
          test = std(test[, -1], trainMarkBase, standardize)     
          train = cbind(trainMarkBaseStatus, train)
          colnames(train) <- c("status","m1", "m2")
          
          test = cbind(testStatus, test)  
          colnames(test) <- c("status","m1", "m2")
          
          model <- gam::gam(status ~ gam::s(m1, df = df1) + 
                         gam::s(m2, df = df2), 
                       data = train, family = binomial)
          
          comb.score <- predict(model,newdata = test,type="response")
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(test$status, as.numeric(comb.score))))
          
          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- 
          resample_results$parameters[[max_AUC]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newdata = data, 
                                      type = "response")
      
    }
    
    else {
      
      data2 <- std(data[, -1], data[, -1], standardize)
      
      model <- gam::gam(status ~ gam::s(m1, df = df1) + 
                     gam::s(m2, df = df2), 
                   data = data2, family = binomial)
      
      comb.score <- predict(model, newdata = markers, type="response")
      
      parameters <- model
      
    }
    
  }
  
  else if (method == "nsgam"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = data[iters[[i]], ]
        trainMarkBaseStatus = trainMarkBase[, 1]
        trainMarkBase = trainMarkBase[, -1]
        train = std(trainMarkBase, trainMarkBase, standardize)
        
        train = cbind(trainMarkBaseStatus, train)
        colnames(train) <- c("status","m1", "m2")
        
        test = std(data[, -1], trainMarkBase, standardize)
        
        test = cbind(data[, 1], test)    
        colnames(test) <- c("status","m1", "m2")
        
        model <- gam::gam(status ~ splines::ns(m1, df = df1) + 
                       splines::ns(m2, df = df2), 
                     data = train, family = binomial)
        
        comb.score <- predict(model, newdata = test, type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newdata = data, 
                                      type = "response")
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = data[-folds[[i]], ]
          trainMarkBaseStatus = trainMarkBase[, 1]
          trainMarkBase = trainMarkBase[, -1]
          train = std(trainMarkBase, trainMarkBase, standardize)
          test = data[folds[[i]], ]
          testStatus = test[, 1]
          test = std(test[, -1], trainMarkBase, standardize)     
          train = cbind(trainMarkBaseStatus, train)
          colnames(train) <- c("status","m1", "m2")
          
          test = cbind(testStatus, test)  
          colnames(test) <- c("status","m1", "m2")
          
          model <- gam::gam(status ~ splines::ns(m1, df = df1) + 
                         splines::ns(m2, df = df2), 
                       data = train, family = binomial)
          
          comb.score <- predict(model, newdata = test, type="response")
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(test$status, as.numeric(comb.score))))
          
          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- 
          resample_results$parameters[[max_AUC]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC]]
      comb.score <- predict(parameters, newdata = data, 
                                      type = "response")
      
    }
    
    else {
      
      data2 <- std(data[, -1], data[, -1], standardize)
      
      model <- gam::gam(status ~ splines::ns(m1, df = df1) + 
                     splines::ns(m2, df = df2), 
                   data = data2, family = binomial)
      
      comb.score <- predict(model, newdata = markers, type="response")
      
      parameters <- model
      
    }
    
  }
  
  comb.score <- as.matrix(comb.score)
  status <- data$status
  
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status,
                 event = event, direction = direction, conf.level = conf.level,
                 cutoff.method = cutoff.method)
  
   model_fit <- list(CombType = "nonlinComb",
                  Method = method,
                  Standardize = standardize,
                  Parameters = parameters
                 )
  
  allres$fit <- model_fit

  return(allres)


}


