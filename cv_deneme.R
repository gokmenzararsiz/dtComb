
cv.linComb <- function(markers = NULL, status = NULL,  
                       train_method = "cv", 
                       train_number = 10){
  
  
  set.seed(100)
  folds = caret::createFolds(status, train_number)
  
  for(i in (1:train_number)){
    
    trainMark = markers[-folds[[i]], ]
    testMark = markers[folds[[i]], ]
    
    trainStat = status[-folds[[i]] ]
    testStat = status[folds[[i]] ]
    
  }
  
  
}