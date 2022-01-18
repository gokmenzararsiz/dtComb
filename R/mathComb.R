# TODO: Add comment
#
# Author: serra ilayda yerlitas
###############################################################################
#' @title Combine two diagnostic tests with several mathematical operators with 
#'  methods.
#'
#' @description The code{mathCom} function returns the combination results of 
#'  two diagnostic tests with different mathematical operators, standardization,
#'  and transform options.
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
#' combining the markers. The available methods are:
#' \itemize{
#' \item \code{add}: 
#' \item \code{multiply}: 
#' \item \code{divide}: 
#' \item \code{subtract}: 
#' \item \code{distance}: 
#' \item \code{baseinexp}: 
#' \item \code{expinbase}: 
#' \item \code{TS}: Combination score obtained by using the trigonometric
#' functions of the theta value that optimizes the AUC
#' }
#' 
#' @param distance a \code{character} string specifying the method used for
#' combining the markers. The available methods are:
#' \itemize{
#' \item \code{euclidean}: 
#' \item \code{manhattan}: 
#' \item \code{chebyshev}: 
#' \item \code{kulczynski_d}: 
#' \item \code{lorentzian}: 
#' \item \code{avg}: 
#' }
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
#' @param transform a \code{character} string indicating the name of the
#' standardization method. The default option is no standardization applied.
#' Available options are:
#' \itemize{
#' \item \code{log}: 
#' \item \code{exp}: 
#' \item \code{sin}: 
#' \item \code{cos}: 
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


##implementation

#' data(exampleData1)
#' markers <- exampleData1[, -1]
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#' direction <- "<"
#' cutoff.method <- "youden"

# 
# score1 <- mathComb(markers = markers, status = status, event = event,
# method = "distance", distance ="avg", direction = direction,
# cutoff.method = cutoff.method)

score2 <- mathComb(markers = markers, status = status, event = event,
method = "baseinexp", transform = "exp", direction = direction,
cutoff.method = cutoff.method)
# 
# score3 <- mathComb(markers = markers, status = status, event = event,
# method = "add", power.transform = FALSE, direction = direction,
# cutoff.method = cutoff.method)




mathComb <- function(markers = NULL, status = NULL, event = NULL,

                     method = c("add", "multiply", "divide", "subtract",
                                  "distance", "baseinexp", "expinbase"),
                     distance = c("euclidean", "manhattan", "chebyshev",
                                    "kulczynski_d", "lorentzian", "avg"),
                     standardize = c("none", "range", 
                                     "zScore", "tScore", "mean", "deviance"),
                     transform = c("none", "log", "exp", "sin", "cos"), 
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
  
  if (method %in% c("baseinexp", "expinbase") && transform == "exp" 
     
       && standardize == "none"){
    
      var <- readline(prompt = "Please enter a standardize for this method (range, zScore, tScore, mean, deviance): ")
   
      standardize <- var
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
    
    values <- suppressMessages(pROC::roc(status , power[,n], 
                                         direction = direction))
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

        if(distance == "euclidean"){

          comb.score <- sqrt(markers[, 1] ^ 2 + markers[, 2] ^ 2)

        } else if (distance == "manhattan"){
    
          comb.score <- markers[, 1] + markers[, 2]
    
        } else if (distance == "chebyshev"){
          
          comb.score <- apply(markers, 1, max)
          
        } else if (distance == "kulczynski_d"){
          
          a <- abs(markers[, 1] + markers[, 2])
          b <- min(markers)
          comb.score <- a / b
          
        } else if (distance == "lorentzian"){
          
          comb.score <- log(abs(markers[, 1]) + 1) + log(abs(markers[, 2]) + 1)
          
        } else {
          
          comb.score <- (markers[, 1] + markers[, 2] + apply(markers, 1, max)) / 2
          
        }     
  } else if (method == "baseinexp") {
    
    comb.score <- markers[ ,1] ^ markers[ ,2]
    
  } else if (method == "expinbase") {
    
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
