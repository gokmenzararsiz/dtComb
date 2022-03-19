# TODO: Add comment
#
# Author: serra ilayda yerlitas
###############################################################################
#' @title Combine two diagnostic tests with several mathematical operators and 
#'  distance measures.
#'
#' @description The code{mathComb} function returns the combination results of 
#'  two diagnostic tests with different mathematical operators, distance 
#'  measures, standardization, and transform options.
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
#'  \item \code{add}: Combination score obtained by adding markers
#'  \item \code{multiply}: Combination score obtained by multiplying markers
#'  \item \code{divide}: Combination score obtained by dividing markers
#'  \item \code{subtract}: Combination score obtained by subtracting markers
#'  \item \code{distance}: Combination score obtained with the help of 
#'  distance measures.
#'  \item \code{baseinexp}: Combination score obtained by marker1 power marker2.
#'  \item \code{expinbase}: Combination score obtained by marker2 power marker1.
#' 
#' @param distance a \code{character} string specifying the method used for
#'  combining the markers. The available methods are:
#'  \itemize{
#'  \item \code{euclidean}: The euclidean distance is (named after Euclid) a 
#'  straight line distance between two points. Euclid argued that the shortest 
#'  distance between two points is always a line. 
#'  Euclidean: d = sqrt( ∑ | P_i - Q_i |^2)
#'  \item \code{manhattan}: The Manhattan distance, also called the Taxicab 
#'  distance or the City Block distance, calculates the distance between two 
#'  real-valued vectors.
#'  Manhattan : d = ∑ | P_i - Q_i |
#'  \item \code{chebyshev}: Chebyshev distance (or Tchebychev distance), maximum 
#'  metric, or L∞ metric is a metric defined on a vector space where the 
#'  distance between two vectors is the greatest of their differences along any 
#'  coordinate dimension.
#'  Chebyshev : $d = max | P_i - Q_i |$
#'  \item \code{kulczynski_d}: 
#'  Kulczynski d : d = ∑ | P_i - Q_i | / ∑ min(P_i , Q_i)
#'  \item \code{lorentzian}: 
#'  Lorentzian : d = ∑ ln(1 + | P_i - Q_i |)
#'  \item \code{avg}: 
#'  Avg(L_1, L_n) : d = ∑ | P_i - Q_i| + max{ | P_i - Q_i |} / 2
#'  \item \code{taneja}:
#'  Taneja : d = ∑ ( P_i + Q_i / 2) * log( P_i + Q_i / ( 2 * sqrt( P_i * Q_i)) )
#'  \item \code{kumar-johnson}: 
#'  Kumar-Johnson : d = ∑ (P_i^2 - Q_i^2)^2 / 2 * (P_i * Q_i)^1.5
#' }
#' 
#' @param standardize a \code{character} string indicating the name of the
#'  standardization method. The default option is no standardization applied.
#'  Available options are:
#'  \itemize{
#'  \item \code{range}: Standardization to a range between 0 and 1
#'  \item \code{zScore}: Standardization using z scores with mean = 0
#'  and standard deviation = 1
#'  \item \code{tScore}: Standardization using T scores. The range varies 
#'  between usually 20 and 80
#'  \item \code{mean}: Standardization with sample mean = 1
#'  \item \code{deviance}: Standardization with sample standard deviation = 1
#' }
#' 
#' @param transform a \code{character} string indicating the name of the
#'  standardization method. The default option is no standardization applied.
#'  Available options are:
#'  \itemize{
#'  \item \code{log}: Applies logarithm transform to markers before calculating 
#'  combination score
#'  \item \code{exp}: Applies exponential transform to markers before 
#'  calculating combination score
#'  \item \code{sin}: Applies sinus trigonometric transform to markers before 
#'  calculatin combination score
#'  \item \code{cos}: Applies cosinus trigonometric transform to markers before 
#'  calculating combination score
#' }
#' 
#' @param power.transform a \code{logical} variable that determines whether to 
#'  apply power to the markers giving the optimum AUC value in the [-3, 3] 
#'  range, before calculating the combination score (FALSE, default).
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
#' @return A list of \code{numeric} linear combination scores calculated
#'  according to the given method and standardization option
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec
#'
#' @examples
#'
#' data(exampleData1)
#' markers <- exampleData1[, -1]
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#' direction <- "<"
#' cutoff.method <- "youden"
#'
#' score1 <- mathComb(markers = markers, status = status, event = event,
#' method = "distance", distance ="taneja", direction = "auto", standardize = "range",
#' cutoff.method = cutoff.method)
#'
#' score2 <- mathComb(markers = markers, status = status, event = event,
#' method = "expinbase", transform = "exp", direction = direction,
#' cutoff.method = cutoff.method)
#'
#' score3 <- mathComb(markers = markers, status = status, event = event,
#' method = "add", transform = "log", direction = direction,
#' cutoff.method = cutoff.method)
#' 
#' @export



mathComb <- function(markers = NULL, status = NULL, event = NULL,
                     
                     method = c("add", "multiply", "divide", "subtract",
                                "distance", "baseinexp", "expinbase"),
                     distance = c("euclidean", "manhattan", "chebyshev",
                                  "kulczynski_d", "lorentzian", "avg", 
                                  "taneja","kumar-johnson"),
                     standardize = c("none", "range", 
                                     "zScore", "tScore", "mean", "deviance"),
                     transform = c("none", "log", "exp", "sin", "cos"), 
                     power.transform = FALSE, direction = c("auto", "<", ">"), 
                     conf.level = 0.95, cutoff.method = c("youden", "roc01")){

  match.arg(method)
  match.arg(distance)
  match.arg(direction)
  match.arg(cutoff.method)
  
  raw.markers <- markers
  raw.status <- status
  
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
  
  status_levels <- levels(status)
  status <- factor(ifelse(status == event, 1, 0))
  
  comp <- complete.cases(markers)
  markers <- markers[comp, ]
  status <- status[comp]
  
  if (length(method) != 1){
    stop("No method provided")
  }
  
  if (method != "distance"){
    
    distance <- NULL
    
  }
  
  if (power.transform != TRUE){
    
    max_power <- NULL
    
  }
  
  if(any(standardize == "none")){
    
    standardize <- "none"
  }
  
  std = matrix(,2,4)
  colnames(std) <- c("mean", "sd", "min", "max")
  
  for (j in 1:2) {
    
    std[, j]
    for (i in 1:ncol(markers)) {
      
      std[i, ] = cbind(mean(markers[, i]),sd(markers[, i]), 
                       min(markers[, i]),max(markers[, i]))
      
    }
  }
  
  
  markers <- std.train(markers, standardize) 
  
  if(any(transform == "none")){
    
    markers <- markers
    transform <- "none"
    
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
  
  n <- as.matrix(seq(-3, 3, 0.1))
  p <- seq(1:nrow(n))
  
  get_roc <- function(p){
    
    values <- suppressMessages(pROC::roc(status , power[,p], 
                                         direction = direction))
    auc <- values$auc
    
    return(auc)
    
  }
  if (method == "add"){
    
    if(power.transform == TRUE){
      
      power <- apply(n, 1, power.add)
      
      auc_list <- sapply(p, get_roc)
      max_index <- which(auc_list == max(auc_list))
      
      if(length(max_index)>1){
        
        max_index <- max_index[1]
        
      }
      max_power <- (n[max_index])
      comb.score <- power[,max_index]
    } 
    else {comb.score <- markers[ ,1] + markers[ ,2]}
    
  } else if (method == "multiply") {
    
    comb.score <- markers[ ,1] * markers[ ,2]
    
  } else if (method == "divide"){
    
    comb.score <- markers[ ,1] / markers[ ,2]
    
  } else if (method == "subtract"){
    
    if(power.transform == TRUE){
      
      power <- apply(n, 1, power.subt)
      
      auc_list <- sapply(p, get_roc)
      max_index <- which(auc_list == max(auc_list))
      
      if(length(max_index)>1){ 
        
        max_index <- max_index[1]
        
      }
      max_power <- (n[max_index])
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
      
    } else if (distance == "taneja"){
      
      epsilon <- 0.00001
      z1 <- (markers[, 1] + 0.00001) / 2
      z2 <- (markers[, 2] + 0.00001) / 2
      comb.score <- (z1 / 2) * log(z1 * sqrt(markers[, 1] * epsilon)) +
        (z2 / 2) * log(z2 * sqrt(markers[, 2] * epsilon))
      
    } else if (distance == "kumar-johnson"){
      
      epsilon <- 0.00001
      z1 <- (((markers[, 1] ^ 2) - (epsilon ^ 2)) ^ 2) /
        2 * ((markers[, 1] * epsilon) ^ 1.5) 
      z2 <- (((markers[, 2] ^ 2) - (epsilon ^ 2)) ^ 2) /
        2 * ((markers[, 2] * epsilon) ^ 1.5) 
      comb.score <- z1 + z2
      
    } else {
      
      comb.score <- (markers[, 1] + markers[, 2] + apply(markers, 1, max)) / 2
      
    }     
  } else if (method == "baseinexp") {
    
    comb.score <- markers[ ,1] ^ markers[ ,2]
    
  } else if (method == "expinbase") {
     
    comb.score <- markers[ ,2] ^ markers[ ,1]
    
  }
  if(length(which(is.infinite(comb.score))) > 0 && standardize != "none"){
    warning("Since inifinity is generated in markers, standardize changed to 'none'.")
    return(mathComb(markers = raw.markers, status = raw.status, event = event, 
                    method = method, distance = distance, direction = direction,
                    standardize = "none", cutoff.method = cutoff.method, 
                    transform = transform, power.transform = power.transform,
                    conf.level = conf.level))
  }
  
  if(length(which(is.infinite(comb.score)) ) > 0 && transform != "none"){
    warning("Since inifinity is generated in markers, transform changed to 'none'.")
    return(mathComb(markers = raw.markers, status = raw.status, event = event,
                    method = method, distance = distance, direction = direction,
                    standardize = standardize, cutoff.method = cutoff.method, 
                    transform = "none", power.transform = power.transform,
                    conf.level = conf.level))
  }
  
  comb.score <- as.matrix(comb.score)
  
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status, 
                   event = event, direction = direction, conf.level = conf.level,
                   cutoff.method = cutoff.method)

  model_fit <- list(CombType = "mathComb",
                    Method = method,
                    Distance = distance,
                    Standardize = standardize,
                    Transform = transform,
                    PowerTransform = power.transform,
                    MaxPower = max_power,
                    Std = std)
  
  allres$fit <- model_fit

  xtab <- as.table(cbind(as.numeric(allres$DiagStatCombined$tab$`   Outcome +`),
                         as.numeric(allres$DiagStatCombined$tab$`   Outcome -`)))
  xtab <- xtab[-3,]
  diagonal.counts <- diag(xtab)
  N <- sum(xtab)
  row.marginal.props <- rowSums(xtab)/N
  col.marginal.props <- colSums(xtab)/N

  Po <- sum(diagonal.counts)/N
  Pe <- sum(row.marginal.props*col.marginal.props)
  k <- (Po - Pe)/(1 - Pe)
  
  accuracy = sum(diagonal.counts) / N
  
  print_model = list(CombType = "mathComb",
                     Method = method,
                     Distance = distance,
                     rowcount = nrow(markers),
                     colcount = ncol(markers),
                     classification = status_levels,
                     Pre_processing = standardize,
                     Accuracy = accuracy,
                     Kappa = k,
                     AUC_table = allres$AUC_table,
                     MultComp_table = allres$MultComp_table,
                     DiagStatCombined = allres$DiagStatCombined
                     
  )
  
  print_allres(print_model)

  return(allres)
}


#' @title The power giving optimum AUC value for the addition operator
#'
#' @description The \code{power.add} function calculates the nth power of the 
#' markers for n selected in the range -3 <= n <= 3 for the addition function.
#'
#' @param n a \code{numeric} parameter to estimate the best power.transform for 
#' the combination score
#' 
#' @param marker.set a \code{numeric} data frame that contains the biomarkers
#'
#' @return A \code{numeric} combination score value for calculated with the 
#' help of exponent
#'
#' @author Ilayda Serra Yerlitas, Serra Bersan Gengec
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' n = 0.5
#'
#' comb.score <- power.add(n, marker.set = markers)
#'
#' @export

power.add <- function(n, marker.set){
  
  power1 <- markers[,1] ^ n
  power2 <- markers[,2] ^ n
  
  return (power1 + power2)
}

#' @title The power giving optimum AUC value for the subtraction operator
#'
#' @description The \code{power.subt} function calculates the nth power of the 
#' markers for n selected in the range -3 <= n <= 3 for the subtraction function.
#'
#' @param n a \code{numeric} parameter to estimate the best power.transform for 
#' the combination score
#'
#' @param marker.set a \code{numeric} data frame that contains the biomarkers
#'
#' @return A \code{numeric} combination score value for calculated with the 
#' help of exponent
#'
#' @author Ilayda Serra Yerlitas, Serra Bersan Gengec
#'
#' @examples
#' #call data
#' data(exampleData1)
#'
#' #define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' n = 2.1
#'
#' comb.score <- power.subt(n, marker.set = markers)
#'
#' @export

power.subt <- function(n, marker.set){
  
  power1 <- markers[,1] ^ n
  power2 <- markers[,2] ^ n
  
  return (power1 - power2)
}
