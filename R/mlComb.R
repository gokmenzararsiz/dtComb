##implementation

# data(exampleData1)
# markers <- exampleData1[, -1]
# status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
# method <- "C5.0Tree"
# train_method <- "repeatedcv"
# preProcess = NULL#c("center","scale")
# train_repeats <- 3
# event <- "needed"
# train_number <- 10
# direction <- "<"
# cutoff.method <- "youden"

# model <- mlComb(markers = markers, status = status, event = event,
#                 method = method,
#                 train_method = train_method, train_repeats = train_repeats,
#                 preProcess = preProccess,
#                 direction = direction, cutoff.method = cutoff.method)


mlComb <- function(markers = NULL, status = NULL, event = NULL,
                            method = NULL, train_method = NULL, 
                            train_number = 10, train_repeats = 1,
                            preProcess = NULL, B = 10,
                            direction = c("<", ">"), conf.level = 0.95, 
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
  
  if (is.null(train_method)){
    train_method <- "none"
  }
  
  ordinalMethods <- c("rpartScore", "ordinalRF", "polr", "vglmAdjCat",	
                      "vglmContRatio",	"vglmCumulative")
  
  bagMethods <- c("bagFDA", "bagFDAGCV", "bagEarth", "bagEarthGCV", "logicBag")
  
  verboseMethods <- c("mlpKerasDecay",	"mlpKerasDecayCost",	"mlpKerasDropout",	
                      "mlpKerasDropoutCost",	"deepboost",	"hda",	"mxnet",
                      "plsRglm",	"sda",	"binda",	"gbm",	"ORFlog",	"ORFpls",
                      "ORFridge",	"ORFsvm",	"bartMachine")
  
  if(method %in% ordinalMethods){
    
    stop("The response variable is not compatible with the given method.")
  }
 
  if(method %in% bagMethods){
    
    if(any(train_method == "repeatedcv")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                               trControl = caret::trainControl(method = train_method, 
                                                               number = train_number, 
                                                               repeats = train_repeats,
                                                               classProbs =  TRUE), 
                               preProc = preProcess, metric = "ROC", B = B, ...)
      
    } else{
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                               trControl = caret::trainControl(method = train_method, 
                                                               number = train_number, 
                                                               classProbs =  TRUE), 
                               preProc = preProcess, metric = "ROC", B = B, ...)
      
    }
    
  }
  
  else if(method %in% verboseMethods){
    
    if(any(train_method == "repeatedcv")){
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                     trControl = caret::trainControl(method = train_method, 
                                                     number = train_number, 
                                                     repeats = train_repeats,
                                                     classProbs =  TRUE), 
                     preProc = preProcess, metric = "ROC", verbose = FALSE, ...)
      
    } else{
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                     trControl = caret::trainControl(method = train_method, 
                                                     number = train_number, 
                                                     classProbs =  TRUE), 
                     preProc = preProcess, metric = "ROC", verbose = FALSE, ...)
      
    }
    
  }
  
  else{
    
    if(any(train_method == "repeatedcv")){
     
       modelFit <- caret::train(status ~ ., data = data, method = method,
                     trControl = caret::trainControl(method = train_method, 
                                                     number = train_number, 
                                                     repeats = train_repeats,
                                                     classProbs =  TRUE), 
                               preProc = preProcess, metric = "ROC", ...)
    } else{
      
      modelFit <- caret::train(status ~ ., data = data, method = method,
                    trControl = caret::trainControl(method = train_method, 
                                                    number = train_number, 
                                                    classProbs =  TRUE), 
                               preProc = preProcess, metric = "ROC", ...)
      
    }
  }
  
  
  print(modelFit)
  
  score <- predict(modelFit, newdata = markers, type = "prob")
  
  comb.score <- as.numeric(score[,levels(status)==event])
  
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status,
                   event = event, direction = direction, conf.level = conf.level,
                   cutoff.method = cutoff.method)

  return(allres)
}


