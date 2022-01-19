# data("exampleData1")
# data <- exampleData1
# 
# markers <- data[, -1]
# status <- factor(data$group, levels = c("not_needed", "needed"))
# event <- "needed"
# direction <- "<"
# cutoff.method <- "youden"
# 
#  score3 <- nonlinComb(markers = markers, status = status, event = event, 
#  method = "ridgereg", interact = FALSE, resample = "cv", nfolds = 5,
#  direction = "<", cutoff.method = "youden")


nonlinComb <- function(markers = NULL, status = NULL, event = NULL,
                    method = c("polyreg", "ridgereg", "lassoreg", "elasticreg",
                               "splines", "sgam", "nsgam"),
                    degree1 = 3, degree2 = 3, df1 = 4, df2 = 4,
                    resample = c("none", "cv", "repeatedcv", "boot"),
                    nfolds = 5, nrepeats = 3,
                    standardize = c("none", "range","zScore", "tScore", "mean", 
                                    "deviance"),
                    interact = FALSE, alpha = 0,  
                    direction = c("<", ">"), conf.level = 0.95, 
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
  
  if (is.null(resample)){
    resample <- "none"
  }
  
  if(any(resample == "cv")){
    nrepeats = 1
  }
  
  
  if (is.null(standardize)){
    standardize <- "none"
  }
  
  if(any(standardize == "none")){
    
    markers <- markers
    
  }
  else if (any(standardize == "range")){
    
    markers <- std.range(markers)
    
  }
  else if (any(standardize == "zScore")){
    
    markers <- std.zscore(markers)
    
  }
  else if (any(standardize == "tScore")){
    
    markers <- std.tscore(markers)
    
  }
  else if (any(standardize == "mean")){
    
    markers <- std.mean(markers)
    
  }
  else if (any(standardize == "deviance")){
    
    markers <- std.deviance(markers)
    
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
      
      folds = caret::createDataPartition(status, nfolds)
      
      for(i in (1:nfolds)){
        
        train = data[folds[[i]], ]
        test = data[-folds[[i]], ]
        
        if( interact == TRUE){
          
          interact <- data$m1 * data$m2
          
          model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2) 
                   + interact, data = data, family = binomial(link = "logit"))
        } else {
          
          model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2), 
                       data = data, family = binomial(link = "logit"))
        }
        
        comb.score <- predict(model, newdata = data, type = "response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          train = data[folds[[i]], ]
          test = data[-folds[[i]], ]
          
          if( interact == TRUE){
            
            interact <- data$m1 * data$m2
            
            model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2) 
                    + interact, data = data, family = binomial(link = "logit"))
          } else {
            
            model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2), 
                         data = data, family = binomial(link = "logit"))
          }
          
          comb.score <- predict(model, newdata = data, type = "response")
          
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
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else{
      
      if( interact == TRUE){
        
        interact <- data$m1 * data$m2
        
        model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2) 
                     + interact, data = data, family = binomial(link = "logit"))
      } else {
        
        model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2), 
                     data = data, family = binomial(link = "logit"))
      }
      
      comb.score <- predict(model, newdata = data, type = "response")
      
      parameters <- model
      
    }
  
  } 
  
  else if (method == "ridgereg"){
    
    if(any(resample== "boot")){
      
      folds = caret::createDataPartition(status, nfolds)
      
      for(i in (1:nfolds)){
        
        train = data[folds[[i]], ]
        test = data[-folds[[i]], ]
        
        if( interact == TRUE){
          
          interact <- data$m1 * data$m2
          
          space <- cbind.data.frame(poly(data$m1, degree1),
                                    poly(data$m2, degree2), interact)
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = 0, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                            family = "binomial", lambda = cv.model$lambda.min)
        } else {
          
          space <- cbind.data.frame(poly(data$m1, degree1),
                                    poly(data$m2, degree2))
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = 0, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                            family = "binomial", lambda = cv.model$lambda.min) 
        } 
        
        comb.score<-predict(model, newx = as.matrix(space), type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          train = data[folds[[i]], ]
          test = data[-folds[[i]], ]
          
          if( interact == TRUE){
            
            interact <- data$m1 * data$m2
            
            space <- cbind.data.frame(poly(data$m1, degree1),
                                      poly(data$m2, degree2), interact)
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = 0, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                            family = "binomial", lambda = cv.model$lambda.min)
          } else {
            
            space <- cbind.data.frame(poly(data$m1, degree1),
                                      poly(data$m2, degree2))
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = 0, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                             family = "binomial", lambda = cv.model$lambda.min) 
          } 
          
          comb.score<-predict(model, newx = as.matrix(space), type="response")
          
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
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else{
      
      if( interact == TRUE){
        
        interact <- data$m1 * data$m2
        
        space <- cbind.data.frame(poly(data$m1, degree1),
                                  poly(data$m2, degree2), interact)
        cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                      alpha = 0, family = "binomial")
        model <- glmnet::glmnet(x = space, y = status, alpha = 0, 
                           family = "binomial", lambda = cv.model$lambda.min)
      } else {
        
        space <- cbind.data.frame(poly(data$m1, degree1),
                                  poly(data$m2, degree2))
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
      
      folds = caret::createDataPartition(status, nfolds)
      
      for(i in (1:nfolds)){
        
        train = data[folds[[i]], ]
        test = data[-folds[[i]], ]
        
        if( interact == TRUE){
          
          interact <- data$m1 * data$m2
          
          space <- cbind.data.frame(poly(data$m1, degree1),
                                    poly(data$m2, degree2), interact)
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = 1, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                            family = "binomial", lambda = cv.model$lambda.min)
        } else {
          
          space <- cbind.data.frame(poly(data$m1, degree1),
                                    poly(data$m2, degree2))
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = 1, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                             family = "binomial", lambda = cv.model$lambda.min) 
        }
        
        comb.score <- predict(model, newx = as.matrix(space), type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          train = data[folds[[i]], ]
          test = data[-folds[[i]], ]
          
          if( interact == TRUE){
            
            interact <- data$m1 * data$m2
            
            space <- cbind.data.frame(poly(data$m1, degree1),
                                      poly(data$m2, degree2), interact)
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = 1, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                             family = "binomial", lambda = cv.model$lambda.min)
          } else {
            
            space <- cbind.data.frame(poly(data$m1, degree1),
                                      poly(data$m2, degree2))
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = 1, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                             family = "binomial", lambda = cv.model$lambda.min) 
          }
          
          comb.score <- predict(model, newx = as.matrix(space), type="response")
          
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
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else{
      
      if( interact == TRUE){
        
        interact <- data$m1 * data$m2
        
        space <- cbind.data.frame(poly(data$m1, degree1),
                                  poly(data$m2, degree2), interact)
        cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                      alpha = 1, family = "binomial")
        model <- glmnet::glmnet(x = space, y = status, alpha = 1, 
                           family = "binomial", lambda = cv.model$lambda.min)
      } else {
        
        space <- cbind.data.frame(poly(data$m1, degree1),
                                  poly(data$m2, degree2))
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
      
      folds = caret::createDataPartition(status, nfolds)
      
      for(i in (1:nfolds)){
        
        train = data[folds[[i]], ]
        test = data[-folds[[i]], ]
        
        if( interact == TRUE){
          
          interact <- data$m1 * data$m2
          
          space <- cbind.data.frame(poly(data$m1, degree1),
                                    poly(data$m2, degree2), interact)
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = alpha, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                             family = "binomial", lambda = cv.model$lambda.min)
        } else {
          
          space <- cbind.data.frame(poly(data$m1, degree1),
                                    poly(data$m2, degree2))
          cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                        alpha = alpha, family = "binomial")
          model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                            family = "binomial", lambda = cv.model$lambda.min) 
        }
        
        comb.score <- predict(model, newx = as.matrix(space), type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          train = data[folds[[i]], ]
          test = data[-folds[[i]], ]
          
          if( interact == TRUE){
            
            interact <- data$m1 * data$m2
            
            space <- cbind.data.frame(poly(data$m1, degree1),
                                      poly(data$m2, degree2), interact)
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = alpha, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                            family = "binomial", lambda = cv.model$lambda.min)
          } else {
            
            space <- cbind.data.frame(poly(data$m1, degree1),
                                      poly(data$m2, degree2))
            cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                          alpha = alpha, family = "binomial")
            model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                            family = "binomial", lambda = cv.model$lambda.min) 
          }
          
          comb.score <- predict(model, newx = as.matrix(space), type="response")
          
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
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else{
      
      if( interact == TRUE){
        
        interact <- data$m1 * data$m2
        
        space <- cbind.data.frame(poly(data$m1, degree1),
                                  poly(data$m2, degree2), interact)
        cv.model <- glmnet::cv.glmnet(x = as.matrix(space), y = status, 
                                      alpha = alpha, family = "binomial")
        model <- glmnet::glmnet(x = space, y = status, alpha = alpha, 
                          family = "binomial", lambda = cv.model$lambda.min)
      } else {
        
        space <- cbind.data.frame(poly(data$m1, degree1),
                                  poly(data$m2, degree2))
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
      
      folds = caret::createDataPartition(status, nfolds)
      
      for(i in (1:nfolds)){
        
        train = data[folds[[i]], ]
        test = data[-folds[[i]], ]
        
        model <- glm(status ~ splines::bs(m1,degree = degree1, df = df1) + 
                       splines::bs(m2,degree = degree2, df = df2), 
                     data = data, family = binomial)
        
        comb.score <- predict(model,newdata = markers,type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          train = data[folds[[i]], ]
          test = data[-folds[[i]], ]
          
          model <- glm(status ~ splines::bs(m1,degree = degree1, df = df1) + 
                         splines::bs(m2,degree = degree2, df = df2), 
                       data = data, family = binomial)
          
          comb.score <- predict(model,newdata = markers,type="response")
          
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
      comb.score <- as.matrix(predict(parameters, newdata = data,
                                      type = "response"))
      
    }
    
    else{
      
      model <- glm(status ~ splines::bs(m1,degree = degree1, df = df1) + 
                     splines::bs(m2,degree = degree2, df = df2), 
                   data = data, family = binomial)
      
      comb.score <- predict(model,newdata = markers,type="response")
      
      parameters <- model
      
    }
    
  }
  
  else if (method == "sgam"){
    
    if(any(resample== "boot")){
      
      folds = caret::createDataPartition(status, nfolds)
      
      for(i in (1:nfolds)){
        
        train = data[folds[[i]], ]
        test = data[-folds[[i]], ]
        
        model <- glm(status ~ gam::s(m1, df = df1) + 
                       gam::s(m2, df = df2), 
                     data = data, family = binomial)
        
        comb.score <- predict(model,newdata = markers,type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          train = data[folds[[i]], ]
          test = data[-folds[[i]], ]
          
          model <- glm(status ~ gam::s(m1, df = df1) + 
                         gam::s(m2, df = df2), 
                       data = data, family = binomial)
          
          comb.score <- predict(model,newdata = markers,type="response")
          
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
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else{
      
      model <- glm(status ~ gam::s(m1, df = df1) + 
                     gam::s(m2, df = df2), 
                   data = data, family = binomial)
      
      comb.score <- predict(model,newdata = markers,type="response")
      
      parameters <- model
      
    }
    
  }
  
  else if (method == "nsgam"){
    
    if(any(resample== "boot")){
      
      folds = caret::createDataPartition(status, nfolds)
      
      for(i in (1:nfolds)){
        
        train = data[folds[[i]], ]
        test = data[-folds[[i]], ]
        
        model <- glm(status ~ splines::ns(m1, df = df1) + 
                       splines::ns(m2, df = df2), 
                     data = data, family = binomial)
        
        comb.score <- predict(model,newdata = markers,type="response")
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- model
        resample_results$AUC[[i]] <- auc_value
        
      }
      
      max_AUC <- which(resample_results$AUC == 
                         max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC]]
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          train = data[folds[[i]], ]
          test = data[-folds[[i]], ]
          
          model <- glm(status ~ splines::ns(m1, df = df1) + 
                         splines::ns(m2, df = df2), 
                       data = data, family = binomial)
          
          comb.score <- predict(model,newdata = markers,type="response")
          
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
      comb.score <- as.matrix(predict(parameters, newdata = data, 
                                      type = "response"))
      
    }
    
    else{
      
      model <- glm(status ~ splines::ns(m1, df = df1) + 
                     splines::ns(m2, df = df2), 
                   data = data, family = binomial)
      
      comb.score <- predict(model,newdata = markers,type="response")
      
      parameters <- model
      
    }
    
  }
  
  
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status,
                 event = event, direction = direction, conf.level = conf.level,
                 cutoff.method = cutoff.method)
  
  allres$fit <- parameters

  return(allres)


}


