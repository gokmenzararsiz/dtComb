##implementation

# data(exampleData1)
# markers <- exampleData1[, -1]
# status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
# method <- "C5.0Tree"
# resample <- "repeatedcv"
# preProcess = NULL#c("center","scale")
# nrepeats <- 3
# event <- "needed"
# nfolds <- 10
# direction <- "<"
# cutoff.method <- "youden"

# model <- mlComb(markers = markers, status = status, event = event,
#                 method = method,
#                 resample = resample, nrepeats = nrepeats,
#                 preProcess = preProcess,
#                 direction = direction, cutoff.method = cutoff.method)


mlComb <- function(markers = NULL, status = NULL, event = NULL,
                            method = NULL, resample = NULL, 
                            nfolds = 5, nrepeats = 3,
                            preProcess = NULL, B = 25,
                            direction = c("<", ">", "auto"), conf.level = 0.95, 
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
  
  excludedMethods <- c("glmStepAIC", "ordinalRF", "pcaNNet", "polr", 
                       "stepLDA", "stepQDA", "vglmAdjCat", "vglmContRatio", 
                       "vglmCumulative", "amdai", "awnb", "awtan", "binda", 
                       "chaid", "elm", "logicBag", "logreg", "manb", "mxnet",
                       "mxnetAdam", "nbDiscrete", "nbSearch", "randomGLM", 
                       "rbf", "rrlda", "SLAVE", "svmBoundrangeString", 
                       "svmExpoString", "svmSpectrumString", "tan", "tanSearch",
                       "vbmpRadial")
  
  BMethods <- c("bagFDA", "bagFDAGCV", "bagEarth", "bagEarthGCV")
  
  verboseMethods <- c("gbm", "mlpKerasDecay",	"mlpKerasDecayCost",	
                      "mlpKerasDropout", "mlpKerasDropoutCost",	"deepboost",
                      "hda",	"mxnet", "plsRglm",	"sda",	"ORFlog",	"ORFpls",
                      "ORFridge",	"ORFsvm",	"bartMachine")
  
  probMethods <- c("BstLm", "bstSm", "bstTree", "C5.0Cost", "CSimca", 
                   "FRBCS.CHI", "FH.GBML", "FRBCS.W", "lssvmLinear", 
                   "lssvmPoly", "lssvmRadial", "lvq", "ownn", "partDSA", 
                   "protoclass", "rFerns", "RFlda", "rfRules", "rocc", 
                   "rpartCost", "rpartScore","RSimca", "smda", "snn",
                   "svmLinear3", "svmLinearWeights2")
  
  if(method %in% excludedMethods){
    
    stop("The response variable or the given markers are not 
         compatible with the given method.")
  }
 
  if(method %in% bagMethods){
    
    if(any(resample == "repeatedcv")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                   trControl = caret::trainControl(method = resample, 
                                                   number = nfolds, 
                                                   repeats = nrepeats,
                                                   classProbs =  TRUE), 
                               preProc = preProcess, B = B, ...)
      
    } else{
      
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
      
    } else{
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                     trControl = caret::trainControl(method = resample, 
                                                     number = nfolds,
                                                     classProbs =  TRUE), 
                     preProc = preProcess, verbose = FALSE, ...)
      
    }
    
    score <- predict(modelFit, newdata = markers, type = "prob")
    
  }
  if(method %in% probMethods){
    
    if(any(resample == "repeatedcv")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                    trControl = caret::trainControl(method = resample, 
                                                    number = nfolds, 
                                                    repeats = nrepeats), 
                               preProc = preProcess, ...)
      
    } else{
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                    trControl = caret::trainControl(method = resample, 
                                                    number = nfolds), 
                               preProc = preProcess, ...)
      
    }
    
    score <- predict(modelFit, newdata = markers)
    
  }
  if(method == "deepboost"){
    
    if(any(resample == "repeatedcv")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                    trControl = caret::trainControl(method = resample, 
                                                    number = nfolds, 
                                                    repeats = nrepeats), 
                               preProc = preProcess, verbose = FALSE, ...)
      
    } else{
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                    trControl = caret::trainControl(method = resample, 
                                                    number = nfolds), 
                               preProc = preProcess, verbose = FALSE, ...)
      
    }
    
    score <- predict(modelFit, newdata = markers)
    
  }
  
  else{
    
    if(any(resample == "repeatedcv")){
     
       modelFit <- caret::train(status ~ ., data = data, method = method,
                     trControl = caret::trainControl(method = resample, 
                                                     number = nfolds, 
                                                     repeats = nrepeats,
                                                     classProbs =  TRUE), 
                               preProc = preProcess, ...)
    } else{
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                    trControl = caret::trainControl(method = resample, 
                                                    number = nfolds, 
                                                    classProbs =  TRUE), 
                               preProc = preProcess, ...)
      
    }
    
    score <- predict(modelFit, newdata = markers, type = "prob")
    
  }
  
  
  print(modelFit)
  
  comb.score <- as.numeric(score[,levels(status)==event])
  
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status,
                   event = event, direction = direction, conf.level = conf.level,
                   cutoff.method = cutoff.method)

  allres$fit <- modelFit
  
  return(allres)
}


