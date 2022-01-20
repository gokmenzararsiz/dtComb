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
#  method = "ridgereg", interact = FALSE, direction = "<", cutoff.method = "youden")


nonlinComb <- function(markers = NULL, status = NULL, event = NULL,
                    method = c("polyreg", "ridgereg", "lassoreg", "elasticreg",
                               "splines", "sgam", "nsgam"),
                    degree1 = 3, degree2 = 3, df1 = 4, df2 = 4,
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

  colnames(markers) <- c("m1", "m2")
  data <- cbind(status,markers) 
  if(method == "polyreg"){
    
    if( interact == TRUE){
      
      interact <- data$m1 * data$m2
      
      model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2) 
                   + interact, data = data, family = binomial(link = "logit"))
    } else {
      
      model <- glm(status ~ poly(m1, degree1) + poly(m2, degree2), 
                   data = data, family = binomial(link = "logit"))
    }
    
    comb.score<-predict(model,newdata=data,type="response")
  
  } else if (method == "ridgereg"){
    
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
    
  }
  else if (method == "lassoreg"){
    
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
    
  } 
  
  else if (method == "elasticreg"){
    
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
    
  }
  else if (method == "splines"){
    
    model <- glm(status ~ splines::bs(m1,degree = degree1, df = df1) + 
                    splines::bs(m2,degree = degree2, df = df2), 
                  data = data, family = binomial)
    
    comb.score <- predict(model,newdata = markers,type="response")
    
  } 
  else if (method == "sgam"){
    
    model <- glm(status ~ gam::s(m1, df = df1) + 
                   gam::s(m2, df = df2), 
                 data = data, family = binomial)
    
    comb.score <- predict(model,newdata = markers,type="response")
    
  }
  else if (method == "nsgam"){
    
    model <- glm(status ~ splines::ns(m1, df = df1) + 
                   splines::ns(m2, df = df2), 
                 data = data, family = binomial)
    
    comb.score <- predict(model,newdata = markers,type="response")
    
  }
  allres <- rocsum(markers = markers, comb.score = comb.score, status = status,
                 event = event, direction = direction, conf.level = conf.level,
                 cutoff.method = cutoff.method)

  return(allres)


}


