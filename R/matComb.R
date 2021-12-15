##implementation

# data(exampleData1)
# markers <- exampleData1[, -1]
# status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
# method <- "multiplication"
# event <- "needed"
# direction <- "<"
# cutoff.method <- "youden"

# score <-matComb(markers = markers, status = status, event = event,
# method = "addition", direction = "<", cutoff.method = "youden")

matComb <- function(markers = NULL, status = NULL, event = NULL,
                    method = c("addition", "multiplication", "log.division",
                               "subtraction", "power2", "power1"), 
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
  
  status <- factor(ifelse(status == event, 1, 0))
  
  comp <- complete.cases(markers)
  markers <- markers[comp, ]
  status <- status[comp]
  

  
  if (method == "addition"){
    
    add <- markers[ ,1] + markers[ ,2]
    comb.score <- as.matrix(add)
    
  } else if (method == "multiplication") {
    
    mult <- markers[ ,1] * markers[ ,2]
    comb.score <- as.matrix(mult)
    
  } else if (method == "log.division"){
    
    log.div <- log(markers[ ,1]) / log(markers[ ,2])
    comb.score <- as.matrix(div)
    
  } else if (method == "subtraction"){
    
    sub <- markers[ ,1] - markers[ ,2]
    comb.score <- as.matrix(sub)
    
  } else if(method == "power2"){
    second.power <- round((markers[, 2] ^ markers[,1]), 3)
    comb.score <- as.matrix(second.power)
    
  } else if(method == "power1"){
    
    first.power <- round((markers[, 1] ^ markers[,2]), 3)    
    comb.score <-  as.matrix(first.power)
    
  }
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status, 
                   event = event, direction = direction, conf.level = conf.level,
                   cutoff.method = cutoff.method)
  
  return(allres)
}



markers <- data[, -1]
status <- factor(data$group, levels = c("not_needed", "needed"))
event <- "needed"
score3 <- matComb(markers = markers, status = status, event = event,
 method = "addition", direction = "<", cutoff.method = "youden")
