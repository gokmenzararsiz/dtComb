# TODO: conf.level, ci.auc.method, nBoot, cutoff.method
#
# Author: serra ilayda yerlitas
###############################################################################
#' @title Combine two diagnostic tests with several linear combination methods.
#'
#' @description The \code{linComb} function calculates the combination
#' scores of two diagnostic tests selected among several linear combination
#' methods and standardization options
#'
#' @param markers a \code{numeric} data frame that includes two diagnostic tests
#' results
#'
#' @param status a \code{factor} vector that includes the actual disease
#' status of the patients
#'
#' @param event a \code{character} string that indicates the event in the status
#' to be considered as positive event
#'
#' @param method a \code{character} string specifying the method used for
#' combining the markers. The available methods are:
#' \itemize{
#' \item \code{scoring}: Combination score obtained using the slope values of
#' the relevant logistic regression model
#' \item \code{SL}: Su and Liu combination score obtained by using Fisher's
#' discriminant function under the assumption of multivariate normal
#' distribution model and proportionate covariance matrices
#' \item \code{logistic}: Combination score obtained by fitting a logistic
#' regression model
#' \item \code{minmax}: Linearly combines the minimum and maximum values of
#' the markers
#' \item \code{PT}: Pepe and Thompson combination score obtained by
#' proportioning the slope values the relevant logistic regression model
#' \item \code{PCT}: Pepe, Cai and Langton combination score obtained by using
#' AUC as the parameter of a logistic regression model
#' \item \code{minimax}: Combination score obtained with the Minimax procedure
#' \item \code{TS}: Combination score obtained by using the trigonometric
#' functions of the theta value that optimizes the AUC
#' }
#' \bold{IMPORTANT}: See Details for further information.
#'
#' @param standardize a \code{character} string indicating the name of the
#' standardization method. The default option is no standardization applied.
#' Available options are:
#' \itemize{
#' \item \code{range}: Standardization to a range between 0 and 1
#' \item \code{zScore}: Standardization using z scores with mean = 0
#' and standard deviation = 1
#' \item \code{tScore}: Standardization using T scores. The range varies between
#'  usually 20 and 80
#' \item \code{mean}: Standardization with sample mean = 1
#' \item \code{deviance}: Standardization with sample standard deviation = 1
#' }
#' 
#' @param ndigits a \code{integer} value to indicate the number of decimal places
#' to be used for rounding in Scoring method
#'
#' @param init.param a \code{numeric} initial value to be used for optimization
#' in minmax, PCL, minimax and TS methods
#' 
#' @return A list of \code{numeric} linear combination scores calculated
#' according to the given method and standardization option
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- exampleData1[, -1]
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#'
#' score1 <- linComb(markers = markers, status = status, event = event,
#' method = "TS", resample= "cv",
#' standardize = "tScore", direction = "<", cutoff.method = "youden")
#'
#' score2 <- linComb(markers = markers, status = status, event = event,
#' method = "TS", resample = "boot", standardize = "zScore", direction = "<", 
#' cutoff.method = "youden")
#'
#' score3 <- linComb(markers = markers, status = status, event = event,
#' method = "logistic", resample = "cv", direction = "<", cutoff.method = "youden")
#' 
#' @export

linComb <- function(markers = NULL, status = NULL, event = NULL,
                    method = c("scoring", "SL", "logistic", "minmax", 
                               "PT", "PCL", "minimax", "TS"), 
                    resample = c("none", "cv", "repeatedcv", "boot"),
                    nfolds = 5, nrepeats = 3, niters = 10,
                    standardize = c("none", "range", 
                                    "zScore", "tScore", "mean", "deviance"),
                    ndigits = 0,
                    direction = c("auto","<", ">"), conf.level = 0.95, 
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
  
  status <- factor(ifelse(status == event, 1, 0))
  
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
  
  if (method %in% c("minmax", "PT", "PCL") && (!standardize == "range")){
    
    warning("The used combination method requires range standardization. 
            All biomarker values are standardized to a range between 0 and 1.")
    standardize <- "range"
    
  }
  
  markers <- std(markers, markers, standardize) 
  
  neg.markers <- markers[status != 1, ]
  pos.markers <- markers[status == 1, ]
  
  resample_results <- vector(mode = "list", length = 3)
  names(resample_results) <- c("parameters", "AUC","trainMarks")
  repeated_results <- vector(mode = "list", length = 3)
  names(repeated_results) <- c("parameters", "AUC","trainMarks")
  
  
  if (method == "scoring"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = markers[iters[[i]], ] 
        trainMark = std(trainMarkBase, trainMarkBase, standardize)
        testMark = std(markers, trainMarkBase, standardize)
        
        trainStat = status[iters[[i]] ]
        testStat = status
        
        res <- glm(trainStat ~ trainMark[ , 1] + trainMark[ , 2],
                   family = binomial((link = "logit")))
        
        round.coef <- round(res$coefficients, digits = ndigits)
        comb.score <- as.matrix(testMark) %*% as.matrix(round.coef[-1])
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(testStat, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- round.coef
        resample_results$AUC[[i]] <- auc_value
        resample_results$trainMarks[[i]] <- trainMarkBase
        
      }
      
      max_AUC <- which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% as.matrix(parameters[-1])
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = markers[-folds[[i]], ]
          trainMark = std(trainMarkBase, trainMarkBase, standardize)
          testMark = markers[folds[[i]], ]
          testMark = std(testMark, trainMarkBase, standardize)
          
          trainStat = status[-folds[[i]] ]
          testStat = status[folds[[i]] ]
          
          res <- glm(trainStat ~ trainMark[ , 1] + trainMark[ , 2],
                     family = binomial((link = "logit")))
          
          round.coef <- round(res$coefficients, digits = ndigits)
          comb.score <- as.matrix(testMark) %*% as.matrix(round.coef[-1])
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(testStat, as.numeric(comb.score))))
          
          resample_results$parameters[[i]] <- round.coef
          resample_results$AUC[[i]] <- auc_value
          resample_results$trainMarks[[i]] <- trainMarkBase
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC[1]]]
        repeated_results$trainMarks[[r]] <- resample_results$trainMarks[[max_AUC[1]]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% as.matrix(parameters[-1])
      
    }
    
    else {
  
      res <- glm(status ~ markers[ , 1] + markers[ , 2],
                 family = binomial((link = "logit")))
      
      round.coef <- round(res$coefficients, digits = ndigits)
      parameters <- round.coef
      comb.score <- as.matrix(markers) %*% as.matrix(round.coef[-1])
      
    }
    
  } 
  
  else if (method == "SL" ) {
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = markers[iters[[i]], ]
        trainMark = std(trainMarkBase, trainMarkBase, standardize)
        testMark = std(markers, trainMarkBase, standardize)
        
        trainStat = status[iters[[i]] ]
        testStat = status
        
        neg.markers <- trainMark[trainStat != 1, ]
        pos.markers <- trainMark[trainStat == 1, ]
        
        sum.var <- var(pos.markers) + var(neg.markers)
        subs_mean <- colMeans(pos.markers) - colMeans(neg.markers)
        est.coef <- as.numeric(abs(solve(sum.var) %*% subs_mean))
        comb.score <- as.matrix(testMark) %*% est.coef
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(testStat, as.numeric(comb.score))))
        
        names(est.coef) <- c("alpha", "beta")
        resample_results$parameters[[i]] <- est.coef
        resample_results$AUC[[i]] <- auc_value
        resample_results$trainMarks[[i]] <- trainMarkBase
        
      }
      
      max_AUC <- which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% parameters
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = markers[-folds[[i]], ]
          trainMark = std(trainMarkBase, trainMarkBase, standardize)
          testMark = markers[folds[[i]], ]
          testMark = std(testMark, trainMarkBase, standardize)
          
          trainStat = status[-folds[[i]] ]
          testStat = status[folds[[i]] ]
          
          neg.markers <- trainMark[trainStat != 1, ]
          pos.markers <- trainMark[trainStat == 1, ]
          
          sum.var <- var(pos.markers) + var(neg.markers)
          subs_mean <- colMeans(pos.markers) - colMeans(neg.markers)
          est.coef <- as.numeric(abs(solve(sum.var) %*% subs_mean))
          comb.score <- as.matrix(testMark) %*% est.coef
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(testStat, as.numeric(comb.score))))
          
          names(est.coef) <- c("alpha", "beta")
          resample_results$parameters[[i]] <- est.coef
          resample_results$AUC[[i]] <- auc_value
          resample_results$trainMarks[[i]] <- trainMarkBase
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC[1]]]
        repeated_results$trainMarks[[r]] <- resample_results$trainMarks[[max_AUC[1]]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% parameters
      
    }
    
    else {

      sum.var <- var(pos.markers) + var(neg.markers)
      subs_mean <- colMeans(pos.markers) - colMeans(neg.markers)
      est.coef <- as.numeric(abs(solve(sum.var) %*% subs_mean))
      names(est.coef) <- c("alpha", "beta")
      parameters <- est.coef
      comb.score <- as.matrix(markers) %*% est.coef
      
    }
  } 
  
  else if (method == "logistic"){
    
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
        test = std(data[, -1], trainMarkBase, standardize)
        train = cbind(trainMarkBaseStatus, train)
        colnames(train) <- c("status","m1", "m2")
        
        test = cbind(data[, 1], test)    
        colnames(test) <- c("status","m1", "m2")
        
        res <- glm(status ~ m1 + m2,
                   family = binomial((link = "logit")), data = train)
        
        comb.score <- as.matrix(predict(res, newdata = test, type = "response"))
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(test$status, as.numeric(comb.score))))
        
        resample_results$parameters[[i]] <- res
        resample_results$AUC[[i]] <- auc_value
        resample_results$trainMarks[[i]] <- trainMarkBase
        
      }
      
      max_AUC <- which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(predict(parameters, newdata = data, type = "response"))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      markersData <- markers
      colnames(markersData) <- c("m1", "m2")
      data <- cbind(status,markersData)
      
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
          
          res <- glm(status ~ m1 + m2,
                     family = binomial((link = "logit")), data = train)
          
          comb.score <- as.matrix(predict(res, newdata = test, type = "response"))
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(test$status, as.numeric(comb.score))))
          
          resample_results$parameters[[i]] <- res
          resample_results$AUC[[i]] <- auc_value
          resample_results$trainMarks[[i]] <- trainMarkBase
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC[1]]]
        repeated_results$trainMarks[[r]] <- resample_results$trainMarks[[max_AUC[1]]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(predict(parameters, newdata = data, type = "response"))
      
    }
    
    else {

      res <- glm(status ~ markers[ , 1] + markers[ , 2],
                 family = binomial((link = "logit")))
      parameters <- res
      comb.score <- as.matrix(predict(res, newdata = markers, type = "response"))
      
    }
    
  } 
  
  else if (method == "minmax"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = markers[iters[[i]], ]
        trainMark = std(trainMarkBase, trainMarkBase, standardize)
        testMark = std(markers, trainMarkBase, standardize)
        
        trainStat = status[iters[[i]] ]
        testStat = status
        
        neg.markers <- trainMark[trainStat != 1, ]
        pos.markers <- trainMark[trainStat == 1, ]
        
        init.param <- runif(1, 0, 1)
        
        opt.func <- optim(par = init.param, fn = helper_minmax, neg.set = neg.markers,
                          pos.set = pos.markers, method = "Brent",
                          lower = 0, upper = 1)
        lambda <- as.numeric(opt.func$par)
        
        comb.score <- as.matrix(apply(testMark, 1, max) 
                                + lambda * apply(testMark, 1, min))
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(testStat, as.numeric(comb.score))))
        
        names(lambda) <- "lambda"
        resample_results$parameters[[i]] <- lambda
        resample_results$AUC[[i]] <- auc_value
        resample_results$trainMarks[[i]] <- trainMarkBase
        
      }
      
      max_AUC <- which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(apply(markers, 1, max) 
                              + parameters * apply(markers, 1, min))
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = markers[-folds[[i]], ]
          trainMark = std(trainMarkBase, trainMarkBase, standardize)
          testMark = markers[folds[[i]], ]
          testMark = std(testMark, trainMarkBase, standardize)
          
          trainStat = status[-folds[[i]] ]
          testStat = status[folds[[i]] ]
          
          neg.markers <- trainMark[trainStat != 1, ]
          pos.markers <- trainMark[trainStat == 1, ]
          
          init.param <- runif(1, 0, 1)
          
          opt.func <- optim(par = init.param, fn = helper_minmax, neg.set = neg.markers,
                            pos.set = pos.markers, method = "Brent",
                            lower = 0, upper = 1)
          lambda <- as.numeric(opt.func$par)
          
          comb.score <- as.matrix(apply(testMark, 1, max) 
                                  + lambda * apply(testMark, 1, min))
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(testStat, as.numeric(comb.score))))
          
          names(lambda) <- "lambda"
          resample_results$parameters[[i]] <- lambda
          resample_results$AUC[[i]] <- auc_value
          resample_results$trainMarks[[i]] <- trainMarkBase
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC[1]]]
        repeated_results$trainMarks[[r]] <- resample_results$trainMarks[[max_AUC[1]]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(apply(markers, 1, max) 
                              + parameters * apply(markers, 1, min))
      
    }
    
    else {

      init.param <- runif(1, 0, 1)
      
      opt.func <- optim(par = init.param, fn = helper_minmax, neg.set = neg.markers,
                        pos.set = pos.markers, method = "Brent",
                        lower = 0, upper = 1)
      lambda <- as.numeric(opt.func$par)
      
      names(lambda) <- "lambda"
      parameters <- lambda
      
      comb.score <- as.matrix(apply(markers, 1, max) 
                              + lambda * apply(markers, 1, min))
      
    }
    
  } 
  
  else if(method == "PT"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = as.matrix(markers[iters[[i]], ])
        trainMark = std(trainMarkBase, trainMarkBase, standardize)
        testMark = std(as.matrix(markers), trainMarkBase, standardize)
        
        trainStat = status[iters[[i]] ]
        testStat = status
        
        model <- glm(trainStat ~ trainMark[ , 1] + trainMark[ , 2],
                     family = binomial(link = "logit"))
        lambda <- as.numeric(model$coefficients[3] / model$coefficients[2])
        
        comb.score <- as.matrix(testMark[, 1] + lambda * testMark[, 2])
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(testStat, as.numeric(comb.score))))
        
        names(lambda) <- "lambda"
        resample_results$parameters[[i]] <- lambda
        resample_results$AUC[[i]] <- auc_value
        resample_results$trainMarks[[i]] <- trainMarkBase
        
      }
      
      max_AUC <- which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers[, 1] + parameters * markers[, 2])
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = as.matrix(markers[-folds[[i]], ])
          trainMark = std(trainMarkBase, trainMarkBase, standardize)
          testMark = as.matrix(markers[folds[[i]], ])
          testMark = std(testMark, trainMarkBase, standardize)
          
          trainStat = status[-folds[[i]] ]
          testStat = status[folds[[i]] ]
          
          model <- glm(trainStat ~ trainMark[ , 1] + trainMark[ , 2],
                       family = binomial(link = "logit"))
          lambda <- as.numeric(model$coefficients[3] / model$coefficients[2])
          
          comb.score <- as.matrix(testMark[, 1] + lambda * testMark[, 2])
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(testStat, as.numeric(comb.score))))
          
          names(lambda) <- "lambda"
          resample_results$parameters[[i]] <- lambda
          resample_results$AUC[[i]] <- auc_value
          resample_results$trainMarks[[i]] <- trainMarkBase
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC[1]]]
        repeated_results$trainMarks[[r]] <- resample_results$trainMarks[[max_AUC[1]]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers[, 1] + parameters * markers[, 2])
      
    }
    
    else {
      
      markers <- as.matrix(markers)
      model <- glm(status ~ markers, family = binomial(link = "logit"))
      lambda <- model$coefficients[3] / model$coefficients[2]
      
      names(lambda) <- "lambda"
      parameters <- lambda
      
      comb.score <- as.matrix(markers[, 1] + lambda * markers[, 2])
      
    }
    
  } 
  
  else if(method == "PCL"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = markers[iters[[i]], ]
        trainMark = std(trainMarkBase, trainMarkBase, standardize)
        testMark = std(markers, trainMarkBase, standardize)
        
        trainStat = status[iters[[i]] ]
        testStat = status
        
        neg.markers <- trainMark[trainStat != 1, ]
        pos.markers <- trainMark[trainStat == 1, ]
        
        init.param <- runif(1, 0, 1)
        
        opt.func <- optim(par = init.param, fn = helper_PCL,
                          neg.set = neg.markers , pos.set = pos.markers,
                          method = "Brent", lower = 0, upper = 1)
        
        lambda <- as.numeric(opt.func$par)
        markers <- as.matrix(markers)
        
        comb.score <-  as.matrix(testMark[ ,1] + testMark[ ,2] * lambda)
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(testStat, as.numeric(comb.score))))
        
        names(lambda) <- "lambda"
        resample_results$parameters[[i]] <- lambda
        resample_results$AUC[[i]] <- auc_value
        resample_results$trainMarks[[i]] <- trainMarkBase
        
      }
      
      max_AUC <- which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers[ ,1] + markers[ ,2] * parameters)
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = markers[-folds[[i]], ]
          trainMark = std(trainMarkBase, trainMarkBase, standardize)
          testMark = markers[folds[[i]], ]
          testMark = std(testMark, trainMarkBase, standardize)
          
          trainStat = status[-folds[[i]] ]
          testStat = status[folds[[i]] ]
          
          neg.markers <- trainMark[trainStat != 1, ]
          pos.markers <- trainMark[trainStat == 1, ]
          
          init.param <- runif(1, 0, 1)
          
          opt.func <- optim(par = init.param, fn = helper_PCL,
                            neg.set = neg.markers , pos.set = pos.markers,
                            method = "Brent", lower = 0, upper = 1)
          
          lambda <- as.numeric(opt.func$par)
          markers <- as.matrix(markers)
          
          comb.score <-  as.matrix(testMark[ ,1] + testMark[ ,2] * lambda)
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(testStat, as.numeric(comb.score))))
          
          names(lambda) <- "lambda"
          resample_results$parameters[[i]] <- lambda
          resample_results$AUC[[i]] <- auc_value
          resample_results$trainMarks[[i]] <- trainMarkBase
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC[1]]]
        repeated_results$trainMarks[[r]] <- resample_results$trainMarks[[max_AUC[1]]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers[ ,1] + markers[ ,2] * parameters)
      
    }
    
    else {

      init.param <- runif(1, 0, 1)
      
      opt.func <- optim(par = init.param, fn = helper_PCL,
                        neg.set = neg.markers , pos.set = pos.markers,
                        method = "Brent", lower = 0, upper = 1)
      
      lambda <- as.numeric(opt.func$par)
      markers <- as.matrix(markers)
      
      names(lambda) <- "lambda"
      parameters <- lambda
      
      comb.score <-  as.matrix(markers[ ,1] + markers[ ,2] * lambda)
      
    }
    
  } 
  
  else if(method == "minimax"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = markers[iters[[i]], ]
        trainMark = std(trainMarkBase, trainMarkBase, standardize)
        testMark = std(markers, trainMarkBase, standardize)
        
        trainStat = status[iters[[i]] ]
        testStat = status
        
        neg.markers <- trainMark[trainStat != 1, ]
        pos.markers <- trainMark[trainStat == 1, ]
        
        init.param <- runif(1, 0, 1)
        
        opt.func <- optim(par = init.param, fn = helper_minimax,
                          neg.set = neg.markers , pos.set = pos.markers,
                          markers = trainMark, status = trainStat,
                          method = "Brent", lower = 0, upper = 1)
        t <- as.numeric(opt.func$par)
        
        b.coef <- as.numeric((solve(t * var(pos.markers)) + (1 - t) * var(neg.markers)) %*%
          (colMeans(pos.markers) - colMeans(neg.markers)))
        comb.score <- as.matrix(testMark) %*% b.coef
        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(testStat, as.numeric(comb.score))))
        
        names(b.coef) <- c("b1", "b2")
        resample_results$parameters[[i]] <- b.coef
        resample_results$AUC[[i]] <- auc_value
        resample_results$trainMarks[[i]] <- trainMarkBase
        
      }
      
      max_AUC <- which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% parameters
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = markers[-folds[[i]], ]
          trainMark = std(trainMarkBase, trainMarkBase, standardize)
          testMark = markers[folds[[i]], ]
          testMark = std(testMark, trainMarkBase, standardize)
          
          trainStat = status[-folds[[i]] ]
          testStat = status[folds[[i]] ]
          
          neg.markers <- trainMark[trainStat != 1, ]
          pos.markers <- trainMark[trainStat == 1, ]
          
          init.param <- runif(1, 0, 1)
          
          opt.func <- optim(par = init.param, fn = helper_minimax,
                            neg.set = neg.markers , pos.set = pos.markers,
                            markers = trainMark, status = trainStat,
                            method = "Brent", lower = 0, upper = 1)
          t <- as.numeric(opt.func$par)
          
          b.coef <- as.numeric((solve(t * var(pos.markers)) + (1 - t) * var(neg.markers)) %*%
                                 (colMeans(pos.markers) - colMeans(neg.markers)))
          comb.score <- as.matrix(testMark) %*% b.coef
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(testStat, as.numeric(comb.score))))
          
          names(b.coef) <- c("b1", "b2")
          resample_results$parameters[[i]] <- b.coef
          resample_results$AUC[[i]] <- auc_value
          resample_results$trainMarks[[i]] <- trainMarkBase
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC[1]]]
        repeated_results$trainMarks[[r]] <- resample_results$trainMarks[[max_AUC[1]]]
        
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% parameters
      
    }
    
    else {
      
      init.param <- runif(1, 0, 1)
      
      opt.func <- optim(par = init.param, fn = helper_minimax,
                        neg.set = neg.markers , pos.set = pos.markers,
                        markers = markers, status = status,
                        method = "Brent", lower = 0, upper = 1)
      t <- as.numeric(opt.func$par)
      
      b.coef <- (solve(t * var(pos.markers)) + (1 - t) * var(neg.markers)) %*%
        (colMeans(pos.markers) - colMeans(neg.markers))
      
      names(b.coef) <- c("b1", "b2")
      parameters <- b.coef
      
      comb.score <- as.matrix(markers) %*% b.coef
      }
    
  } 
  
  else if(method == "TS"){
    
    if(any(resample== "boot")){
      
      iters = caret::createResample(status, niters)
      
      for(i in (1:niters)){
        
        trainMarkBase = markers[iters[[i]], ]
        trainMark = std(trainMarkBase, trainMarkBase, standardize)
        testMark = std(markers, trainMarkBase, standardize)
        
        trainStat = status[iters[[i]] ]
        testStat = status
        
        init.param <- runif(1, -1.57079633, 1.57079633)
        
        opt.func <- optim(par = init.param, fn = helper_TS, markers = trainMark,
                          status = trainStat, method = "Brent",
                          lower = -1.57079633, upper = 1.57079633)
        theta <- as.numeric(opt.func$par)
        
        a1 <- sin(theta)
        a2 <- cos(theta)
        
        comb.score <- as.matrix(a1 * testMark[, 1] + a2 * testMark[, 2])

        
        auc_value <- suppressMessages(as.numeric(
          pROC::auc(testStat, as.numeric(comb.score))))
        
        trigFuncs <- as.numeric(c(a1, a2))
        names(trigFuncs) <- c("sin(theta)", "cos(theta)")
        resample_results$parameters[[i]] <- trigFuncs
        resample_results$AUC[[i]] <- auc_value
        resample_results$trainMarks[[i]] <- trainMarkBase
      }
      
      max_AUC <- which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(parameters[1] * markers[, 1] + parameters[2] * markers[, 2])
      
    }
    
    else if(any(resample == "cv") || any(resample == "repeatedcv")){
      
      for(r in (1:nrepeats)){
        
        folds = caret::createFolds(status, nfolds)
        
        for(i in (1:nfolds)){
          
          trainMarkBase = markers[-folds[[i]], ]
          trainMark = std(trainMarkBase, trainMarkBase, standardize, TRUE)
          testMark = markers[folds[[i]], ]
          testMark = std(testMark, trainMarkBase, standardize, TRUE)
          
          trainStat = status[-folds[[i]] ]
          testStat = status[folds[[i]] ]
          
          init.param <- runif(1, -1.57079633, 1.57079633)
          
          opt.func <- optim(par = init.param, fn = helper_TS, markers = trainMark,
                            status = trainStat, method = "Brent",
                            lower = -1.57079633, upper = 1.57079633)
          theta <- as.numeric(opt.func$par)
          
          a1 <- sin(theta)
          a2 <- cos(theta)
          
          comb.score <- as.matrix(a1 * testMark[, 1] + a2 * testMark[, 2])
          
          
          auc_value <- suppressMessages(as.numeric(
            pROC::auc(testStat, as.numeric(comb.score))))
          
          trigFuncs <- as.numeric(c(a1, a2))
          names(trigFuncs) <- c("sin(theta)", "cos(theta)")
          resample_results$parameters[[i]] <- trigFuncs
          resample_results$AUC[[i]] <- auc_value
          resample_results$trainMarks[[i]] <- trainMarkBase
          
        }
        
        max_AUC <- which(resample_results$AUC ==
                           max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <- resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC[1]]]
        repeated_results$trainMarks[[r]] <- resample_results$trainMarks[[max_AUC[1]]]
      }
      
      max_AUC <- which(repeated_results$AUC ==
                         max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      
      comb.score <- as.matrix(parameters[1] * markers[, 1] + parameters[2] * markers[, 2])
    }
    
    else {
      
      init.param <- runif(1, -1.57079633, 1.57079633)
      
      opt.func <- optim(par = init.param, fn = helper_TS, markers = markers,
                        status = status, method = "Brent",
                        lower = -1.57079633, upper = 1.57079633)
      theta <- as.numeric(opt.func$par)
      
      a1 <- sin(theta)
      a2 <- cos(theta)
      
      trigFuncs <- as.numeric(c(a1, a2))
      names(trigFuncs) <- c("sin(theta)", "cos(theta)")
      parameters <- trigFuncs
      
      comb.score <- as.matrix(a1 * markers[, 1] + a2 * markers[, 2])
      
    }
    
  }
  
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status,
                   event = event, direction = direction, conf.level = conf.level,
                   cutoff.method = cutoff.method)
  std = matrix(,2,4)
  colnames(std) <- c("mean", "sd", "min", "max")

  max_AUC <- c()
  trainMarkBase <- data.frame()
  if(resample == "boot"){
    
    max_AUC <- which(resample_results$AUC ==
                       max(unlist(resample_results$AUC)))
    trainMarkBase <- resample_results$trainMarks[[max_AUC[1]]]
    
  } else {
    
    max_AUC <- which(repeated_results$AUC ==
                       max(unlist(repeated_results$AUC)))
    trainMarkBase <- repeated_results$trainMarks[[max_AUC[1]]]
    
  }
  
  
  for (j in 1:2) {
    
    std[, j]
    for (i in 1:ncol(trainMarkBase)) {
      
      std[i, ] = cbind(mean(trainMarkBase[, i]),sd(trainMarkBase[, i]), 
                       min(trainMarkBase[, i]),max(trainMarkBase[, i]))
      
    }
  }

   model_fit <- list(CombType = "linComb",
                  Method = method,
                  Standardize = standardize,
                  Parameters = parameters,
                  Std = std
                  
                  )
   
  
  allres$fit <- model_fit
  
  return(allres)
  
}

#' @title Mann Whitney AUC estimator for minmax method.
#'
#' @description The \code{helper_minmax} function estimates non-parametric AUC
#' for the given biomarkers
#'
#' @param lambda a \code{numeric} parameter that will be estimated in minmax
#' method for the combination score
#'
#' @param neg.set a \code{numeric} data frame that contains the observations
#' with negative status
#'
#' @param pos.set a \code{numeric} data frame that contains the observations
#' with positive status
#'
#' @return A \code{numeric} value for the estimated AUC
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' neg.set <- markers[status == levels(status)[1], ]
#' pos.set <- markers[status == levels(status)[2], ]
#'
#' lambda = 0.5
#'
#' stat <- helper_minmax(lambda, neg.set = neg.set, pos.set = pos.set)
#'
#' @export

helper_minmax <- function(lambda, neg.set, pos.set){
  
  Xmax <- as.matrix(apply(neg.set, 1, max))
  Xmin <- as.matrix(apply(neg.set, 1, min))
  
  Ymax <- as.matrix(apply(pos.set, 1, max))
  Ymin <- as.matrix(apply(pos.set, 1, min))
  
  W.lambda <- 0
  n <- dim(neg.set)[1]
  m <- dim(pos.set)[1]
  
  for (i in 1:n){
    for(j in 1:m){
      W.lambda <- W.lambda + as.numeric(Ymax[j, ] + lambda * Ymin[j, ] >
                                          Xmax[i, ] + lambda * Xmin[i, ])
    }
  }
  
  return(-W.lambda / (n * m))
}

#' @title Mann Whitney AUC estimator for PCL method.
#'
#' @description The \code{helper_PCL} function estimates non-parametric
#' AUC for the given biomarkers
#'
#' @param lambda a \code{numeric} parameter that will be estimated in minmax
#' method for the combination score
#'
#' @param neg.set a \code{numeric} data frame that contains the observation with
#' negative status
#'
#' @param pos.set a \code{numeric} data frame that contains the observation with
#' positive status
#'
#' @return A \code{numeric} value for the estimated AUC
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' neg.set <- markers[status == levels(status)[1], ]
#' pos.set <- markers[status == levels(status)[2], ]
#'
#' lambda = 0.5
#'
#' stat <- helper_PCL(lambda, neg.set = neg.set, pos.set = pos.set)
#'
#' @export

helper_PCL <- function(lambda, neg.set, pos.set){
  
  YD1 <- as.matrix(pos.set[ , 1])
  YD2 <- as.matrix(pos.set[ , 2])
  
  YDN1 <- as.matrix(neg.set[ , 1])
  YDN2 <- as.matrix(neg.set[ , 2])
  
  W.lambda <- 0
  
  n <- dim(pos.set)[1]
  m <- dim(neg.set)[1]
  
  for (i in 1:n){
    for(j in 1:m){
      W.lambda <- W.lambda +
        as.numeric(YD1[i, ] + lambda * YD2[i, ] >
                     YDN1[j, ] + lambda * YDN2[j, ])+
        as.numeric(YD1[i, ] + lambda * YD2[i, ] == YDN1[j, ] +
                     lambda * YDN2[j, ])/2
    }
  }
  
  return(-W.lambda / (n * m))
}

#' @title AUC calculator for minimax method.
#'
#' @description The \code{helper_minimax} function calculates the combination
#' coefficient and AUC value of given biomarkers for minimax method
#'
#' @param t a \code{numeric} parameter that will be estimated in minimax
#' method for the combination score
#'
#' @param neg.set a \code{numeric} data frame that contains the observation with
#' negative status
#'
#' @param pos.set a \code{numeric} data frame that contains the observation with
#' positive status
#'
#' @param marker.set a \code{numeric} data frame that contains the biomarkers
#'
#' @param status a \code{factor} data frame that includes the actual disease
#' status of the patients
#'
#' @return A \code{numeric} AUC value calculated with combination scores using t
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' neg.set <- markers[status == levels(status)[1], ]
#' pos.set <- markers[status == levels(status)[2], ]
#'
#' t <- 0.5
#'
#' stat <- helper_minimax(t, neg.set = neg.set, pos.set = pos.set,
#' marker.set = markers, status)
#'
#' @importFrom pROC auc
#'
#' @export

helper_minimax <- function(t, neg.set, pos.set, markers, status){
  
  b.coef <- (solve(t * var(pos.set)) + (1 - t) * var(neg.set)) %*%
    (colMeans(pos.set) - colMeans(neg.set))
  comb.score <- as.matrix(markers) %*% b.coef
  comb.score <- as.numeric(comb.score)
  
  auc_value <- suppressMessages(as.numeric(pROC::auc(status, comb.score)))
  
  return(-(auc_value))
}

#' @title AUC calculator for minimax method.
#'
#' @description The \code{helper_minimax} function calculates the combination
#' coefficient and AUC value of given biomarkers for minimax method
#'
#' @param t a \code{numeric} parameter that will be estimated in minimax
#' method for the combination score
#'
#' @param neg.set a \code{numeric} data frame that contains the observation with
#' negative status
#'
#' @param pos.set a \code{numeric} data frame that contains the observation with
#' positive status
#'
#' @param marker.set a \code{numeric} data frame that contains the biomarkers
#'
#' @param status a \code{factor} data frame that includes the actual disease
#' status of the patients
#'
#' @return A \code{numeric} AUC value calculated with combination scores using t
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' neg.set <- markers[status == levels(status)[1], ]
#' pos.set <- markers[status == levels(status)[2], ]
#'
#' t <- 0.5
#'
#' stat <- helper_minimax(t, neg.set = neg.set, pos.set = pos.set,
#' marker.set = markers, status)
#'
#' @importFrom pROC auc
#'
#' @export

helper_TS <- function(theta, markers, status){
  
  a1 <- sin(theta)
  a2 <- cos(theta)
  z <- a1 * markers[, 1] + a2 * markers[, 2]
  
  roc_obj <- suppressMessages(pROC::roc(status, z))
  auc_value <- as.numeric(pROC::auc(roc_obj))
  
  return(-(auc_value))
}
