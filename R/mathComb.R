##implementation

# data(exampleData1)
# markers <- exampleData1[, -1]
# status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
# event <- "needed"
# direction <- "<"
# cutoff.method <- "youden"

# score1 <- mathComb(markers = markers, status = status, event = event,
# method = "distance", distance ="euclidean", direction = direction, 
# cutoff.method = cutoff.method)

# score2 <- mathComb(markers = markers, status = status, event = event,
# method = "sec^first", transform = "log", direction = direction, cutoff.method = cutoff.method)

# score3 <- mathComb(markers = markers, status = status, event = event,
# method = "subtract", power.transform = TRUE, direction = direction, 
# cutoff.method = cutoff.method)

mathComb <- function(markers = NULL, status = NULL, event = NULL,

                     method = c("add", "multiply", "divide", "subtract",
                                  "distance", "first^sec", "sec^first"),
                     distance = c("euclidean", "manhattan", "chebyshev",
                                    "kulczynski_d", "lorentzian", "taneja",
                                      "kumar-johnson", "avg"),
                     standardize = c("none", "range", 
                                     "zScore", "tScore", "mean", "deviance"),
                     transform = c("none", "log", "exp", 
                                    "sin", "cos"), 
                     power.transform = FALSE, direction = c("<", ">"), 
                     conf.level = 0.95, cutoff.method = c("youden", "roc01")){


  match.arg(method)
  match.arg(distance)
  match.arg(transform)
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
  if (is.null(transform)){
    
    transform <- "none"
  }
  
  if(any(transform == "none")){
    
    markers <- markers
    
  }
   else if(any(transform == "log")){
     
    markers <- log(markers)
  }
    else if(any(transform == "exp")){
    
    markers <- exp(markers)
    
  }
    else if(any(transform == "sin")){
    
    markers <- sin(markers)
    
  }
    else {
    
    markers <- cos(markers)
    
  }
  
  x <- as.matrix(seq(-3, 3, 0.1))
  n <- seq(1:nrow(x))
  get_roc <- function(n){
    values <- suppressMessages(pROC::roc(status , power[,n], direction = direction))
    auc <- values$auc
    return(auc)
  }
  
  if (method == "add"){
    
   if(power.transform == TRUE){
      
      power <- apply(x, 1, power.add)
      
      auc_list <- sapply(n, get_roc)
      max_index <- which(auc_list == max(auc_list))
      
        if(length(max_index)>1){
          max_index <- max_index[1]
          }
      
        comb.score <- power[,max_index]
        } 
    
        else {comb.score <- markers[ ,1] + markers[ ,2]}
  
  } else if (method == "multiply") {
    
   comb.score <- markers[ ,1] * markers[ ,2]
    
  } else if (method == "divide"){
    
    comb.score <- markers[ ,1] / markers[ ,2]
    
  } else if (method == "subtract"){
    
     if(power.transform == TRUE){
    
        power <- apply(x, 1, power.subt)
    
      auc_list <- sapply(n, get_roc)
      max_index <- which(auc_list == max(auc_list))
      
      if(length(max_index)>1){
        max_index <- max_index[1]
        }
      comb.score <- power[,max_index]
      }
      else{comb.score <- (markers[ ,1] - markers[ ,2])}
    
  }  else if(method == "distance"){
    
      distMethod <- function(params){
       origin <-c(0,0)
       suppressMessages(philentropy::distance(rbind(origin, params), 
                                method = distance, 
                                  use.row.names = TRUE))
    }
    
    comb.score <- as.matrix(unlist(apply(markers, 1, distMethod)))
    rownames(comb.score) <- NULL
    
  } else if (method == "first^sec") {
    
    comb.score <- markers[ ,1] ^ markers[ ,2]
    
  } else if (method == "sec^first") {
    
    comb.score <- markers[ ,2] ^ markers[ ,1]
    
  }
  
  comb.score <- as.matrix(comb.score)
  
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status, 
                   event = event, direction = direction, conf.level = conf.level,
                   cutoff.method = cutoff.method)
  
  return(allres)
}


##### Helper Functions#####

power.add <- function(x){
  power1 <- markers[,1] ^ x
  power2 <- markers[,2] ^ x
  
  return (power1 + power2)
}

power.subt <- function(x){
  power1 <- markers[,1] ^ x
  power2 <- markers[,2] ^ x
  
  return (power1 - power2)
}
