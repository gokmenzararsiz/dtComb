# TODO: predict function
#
# Author: serra berþan gengeç
###############################################################################
#' @title Combine two diagnostic tests with Machine Learning Algorithms.
#'
#' @description The \code{mlComb} function calculates the combination
#' scores of two diagnostic tests selected among several Machine Learning 
#' Algorithms
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
#' combining the markers. For the available methods see availableMethods()
#' 
#' \bold{IMPORTANT}: See https://topepo.github.io/caret/available-models.html 
#' for further information about the methods used in this function.
#' 
#' @param resample a \code{character} string that indicates the resampling 
#' method used while training the model. The available methods are "boot", 
#' "boot632", "optimism_boot", "boot_all", "cv", "repeatedcv", "LOOCV", "LGOCV",
#' "none", "oob", "adaptive_cv", "adaptive_boot" and "adaptive_LGOCV". for 
#' details of these resampling methods see ?caret::trainControl
#' 
#' @param nfolds a \code{numeric} value that indicates the number of folds or 
#' the number of resampling iterations (5, default)
#' 
#' @param nrepeats a \code{numeric} value that indicates the number of repeats 
#' for "repeatedcv" option of resampling methods (3, default)
#' 
#' @param preProcess a \code{character} string that indicates the pre-processing
#' options to be applied in the data before training the model. Available 
#' pre-processing methods are: "BoxCox", "YeoJohnson", "expoTrans", "center", 
#' "scale", "range", "knnImpute", "bagImpute", "medianImpute", "pca", "ica", 
#' "spatialSign", "corr", "zv", "nzv", and "conditionalX". For detailed 
#' information about the methods see ?caret::preProcess
#' 
#' @param B a \code{numeric} value that is the number of bootstrap samples for 
#' bagging classifiers, "bagFDA", "bagFDAGCV", "bagEarth" and "bagEarthGCV". 
#' (25, default)
#' 
#' @param direction a \code{character} string determines in which direction the 
#' comparison will be made.  “>”: if the predictor values for the control group 
#' are higher than the values of the case group (controls > cases). 
#' “<”: if the predictor values for the control group are lower or equal than 
#' the values of the case group (controls < cases). 
#'
#' @param conf.level a \code{numeric} value to  determine the confidence interval
#' for the ROC curve(0.95, default).
#' 
#' @param cutoff.method  a \code{character} string determines the cutoff method
#' for the ROC curve.
#' 
#' @param \dots optional arguments passed to selected classifiers.
#' 
#' @return A \code{list} of AUC values, diagnostic statistics,
#' coordinates of the ROC curve for the combination score obtained using 
#' Machine Learning Algorithms as well as the given biomarkers individually, a 
#' comparison table for the AUC values of individual biomarkers and combination 
#' score obtained and the fitted model.
#'
#' @author Serra Ilayda  Yerlitas, Serra Bersan Gengec
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
#' score1 <- mlComb(markers = markers, status = status, event = event,
#' method = "knn", resample = "boot632",  niters = 15,
#' preProcess = "center", direction = "<", cutoff.method ="youden")
#' 
#' @export


mlComb <- function(markers = NULL, status = NULL, event = NULL,
                            method = NULL, resample = NULL, 
                            niters = 5, nfolds = 5, nrepeats = 3,
                            preProcess = NULL, B = 25,
                            direction = c("auto", "<", ">"), conf.level = 0.95, 
                            cutoff.method = c("youden", "roc01"), ...){
  
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
  
  data <- cbind(status,markers)
  
  if (is.null(resample)){
    resample <- "none"
  }
  
  
  BMethods <- c("bagFDA", "bagFDAGCV", "bagEarth", "bagEarthGCV")
  
  verboseMethods <- c("gbm", "mlpKerasDecay",	"mlpKerasDecayCost",	
                      "mlpKerasDropout", "mlpKerasDropoutCost",	"deepboost",
                      "hda",	"mxnet", "plsRglm",	"sda",	"ORFlog",	"ORFpls",
                      "ORFridge",	"ORFsvm",	"bartMachine")
  
  if(!(method %in% allMethods[,1])){
    
    stop("The response given method is not available for mlComb function. 
         See availableMethods() for the list of methods available.")
  }
 
  if(method %in% BMethods){
    
    if(any(resample == "repeatedcv")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                   trControl = caret::trainControl(method = resample, 
                                                   number = nfolds, 
                                                   repeats = nrepeats,
                                                   classProbs =  TRUE), 
                               preProc = preProcess, B = B, ...)
      
    } else if (resample %in%  c("boot", "boot632", "optimism_boot", "boot_all")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                               trControl = caret::trainControl(method = resample, 
                                                               number = niters,
                                                               classProbs =  TRUE), 
                               preProc = preProcess, B = B, ...)
      
    } else {
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                    trControl = caret::trainControl(method = resample, 
                                                    number = nfolds,
                                                    classProbs =  TRUE), 
                               preProc = preProcess, B = B, ...)
      
    }
    
    score <- predict(modelFit, newdata = markers, type = "prob")
  }
  
  else if(method %in% verboseMethods){
    
    if(any(resample == "repeatedcv")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                     trControl = caret::trainControl(method = resample, 
                                                     number = nfolds, 
                                                     repeats = nrepeats,
                                                     classProbs =  TRUE), 
                     preProc = preProcess, verbose = FALSE, ...)
      
    } else if (resample %in%  c("boot", "boot632", "optimism_boot", "boot_all")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                               trControl = caret::trainControl(method = resample, 
                                                               number = niters,
                                                               classProbs =  TRUE), 
                               preProc = preProcess, verbose = FALSE, ...)
      
    } else {
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                     trControl = caret::trainControl(method = resample, 
                                                     number = nfolds,
                                                     classProbs =  TRUE), 
                     preProc = preProcess, verbose = FALSE, ...)
      
    }
    
    score <- predict(modelFit, newdata = markers, type = "prob")
    
  }
  
  else{
    
    if(any(resample == "repeatedcv")){
     
       modelFit <- caret::train(status ~ ., data = data, method = method,
                     trControl = caret::trainControl(method = resample, 
                                                     number = nfolds, 
                                                     repeats = nrepeats,
                                                     classProbs =  TRUE), 
                               preProc = preProcess, ...)
    } else if (resample %in%  c("boot", "boot632", "optimism_boot", "boot_all")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                               trControl = caret::trainControl(method = resample, 
                                                               number = niters,
                                                               classProbs =  TRUE), 
                               preProc = preProcess, ...)
      
    } else {
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                    trControl = caret::trainControl(method = resample, 
                                                    number = nfolds, 
                                                    classProbs =  TRUE), 
                               preProc = preProcess, ...)
      
    }
    
    score <- predict(modelFit, newdata = markers, type = "prob")
    
  }
  
  comb.score <- as.numeric(score[,levels(status)==event])
  
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status,
                   event = event, direction = direction, conf.level = conf.level,
                   cutoff.method = cutoff.method)
  
   model_fit <- list(CombType = "mlComb",
                     Model = modelFit)
   
  allres$fit <- model_fit
  
  return(allres)
}

availableMethods <- function(){
  
  message("The available methods are listed below. For more information 
about the methods see https://topepo.github.io/caret/available-models.html") 
  print(allMethods)           
  
}
