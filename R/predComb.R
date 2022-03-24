# score1 <- nonlinComb(markers = markers, status = status, event = event,
# method = "nsgam", resample = "boot", include.interact = FALSE,
# standardize = "zScore", cutoff.method = "youden")
# newdata <- markers
#'
#'comb.score <-  comb.predict(score1, markers)

comb.predict <- function(model, newdata){

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  combtype <- model$fit$CombType
  
  if(combtype != "mlComb"){
    
    colnames(newdata) <- c("m1", "m2")
    newdata = std.test(newdata, model)
    
  }
  
  if(combtype == "linComb"){
    
    method <- model$fit$Method
    parameters <- model$fit$Parameters
    
    if (method == "scoring"){
      
      comb.score <- as.matrix(newdata) %*% as.matrix(model$fit$Parameters[-1])
    } 
    else if (method == "SL"){
      
      comb.score <- as.matrix(newdata) %*% model$fit$Parameters
    } 
    else if (method == "logistic"){
      
      comb.score <- as.matrix(predict(model$fit$Parameters, newdata = newdata, 
                                      type = "response"))
    }
    else if (method == "minmax"){
      
      comb.score <- as.matrix(apply(newdata, 1, max) 
                              + model$fit$Parameters * apply(newdata, 1, min))
    }
    else if (method == "PT"){
      
      comb.score <- as.matrix(newdata[, 1] + model$fit$Parameters * newdata[, 2])
    }
    else if (method == "PCL"){
      
      comb.score <- as.matrix(newdata[ ,1] + newdata[ ,2] * model$fit$Parameters)
    }
    else if (method == "minimax"){
      
      comb.score <- as.matrix(newdata) %*% model$fit$Parameters
    }
    else {

      comb.score <- as.matrix(model$fit$Parameters[1] * newdata[, 1] + model$fit$Parameters[2] * 
                                newdata[, 2])
    }
  } 
  
  else if(combtype == "nonlinComb"){
    
    method <- model$fit$Method
    parameters <- model$fit$Parameters
    interact <- model$fit$Interact
    
    if(method == "polyreg"){
      
      comb.score <- predict(model$fit$Parameters, newdata = newdata, 
                            type = "response")
    }
    else if(method %in% c("ridgereg", "lassoreg", "elasticreg")){
      
      if(interact == TRUE){
        
        interact <- newdata$m1 * newdata$m2
      
        dataspace <- cbind.data.frame(poly(newdata$m1, model$fit$Degree1),
                                    poly(newdata$m2, model$fit$Degree2), interact)
      
        comb.score <- predict(model$fit$Parameters, newx = as.matrix(dataspace), 
                            type = "response")
      } else {
        
        dataspace <- cbind.data.frame(poly(newdata$m1, model$fit$Degree1),
                                      poly(newdata$m2, model$fit$Degree2))
        
        comb.score <- predict(model$fit$Parameters, newx = as.matrix(dataspace), 
                            type = "response")
      }
    }
    else {
      
      comb.score <- predict(model$fit$Parameters, newx = as.matrix(newdata), 
                            type="response")
      
    }
    
  } 
  else if(combtype == "mlComb"){
    
    model_fit <- model$fit$Model
    
      
    comb.score <- predict(model_fit, newdata = newdata, type = "prob")

  } 
  else {
    
    method = model$fit$Method
    distance = model$fit$Distance
    transform = model$fit$Transform
    power.transform = model$fit$PowerTransform
    max_power = model$fit$MaxPower
    
    if(any(transform == "none")){
      
      newdata <- newdata
      
    }
    else if(any(transform == "log")){
      
      newdata <- log(newdata)
    }
    else if(any(transform == "exp")){
      
      newdata <- exp(newdata)
      
    }
    else if(any(transform == "sin")){
      
      newdata <- sin(newdata)
      
    }
    else {
      
      newdata <- cos(newdata)
      
    }
    
    if (method == "add"){
      
      if(power.transform == TRUE){
        
        comb.score <- apply(max_power, 1, power.add)
        
      } 
      else {comb.score <- newdata[ ,1] + newdata[ ,2]}
      
    } else if (method == "multiply") {
      
      comb.score <- newdata[ ,1] * newdata[ ,2]
      
    } else if (method == "divide"){
      
      comb.score <- newdata[ ,1] / newdata[ ,2]
      
    } else if (method == "subtract"){
      
      if(power.transform == TRUE){
        
        comb.score <- apply(max_power, 1, power.subt)
      }
      else{comb.score <- (newdata[ ,1] - newdata[ ,2])}
      
    }  else if(method == "distance"){
      
      if(distance == "euclidean"){
        
        comb.score <- sqrt(newdata[, 1] ^ 2 + newdata[, 2] ^ 2)
        
      } else if (distance == "manhattan"){
        
        comb.score <- newdata[, 1] + newdata[, 2]
        
      } else if (distance == "chebyshev"){
        
        comb.score <- apply(newdata, 1, max)
        
      } else if (distance == "kulczynski_d"){
        
        a <- abs(newdata[, 1] + newdata[, 2])
        b <- min(newdata)
        comb.score <- a / b
        
      } else if (distance == "lorentzian"){
        
        comb.score <- log(abs(newdata[, 1]) + 1) + log(abs(newdata[, 2]) + 1)
        
      } else if (distance == "taneja"){
        
        epsilon <- 0.00001
        z1 <- (newdata[, 1] + 0.00001) / 2
        z2 <- (newdata[, 2] + 0.00001) / 2
        comb.score <- (z1 / 2) * log(z1 * sqrt(newdata[, 1] * epsilon)) +
          (z2 / 2) * log(z2 * sqrt(newdata[, 2] * epsilon))
        
      } else if (distance == "kumar-johnson"){
        
        epsilon <- 0.00001
        z1 <- (((newdata[, 1] ^ 2) - (epsilon ^ 2)) ^ 2) /
          2 * ((newdata[, 1] * epsilon) ^ 1.5) 
        z2 <- (((newdata[, 2] ^ 2) - (epsilon ^ 2)) ^ 2) /
          2 * ((newdata[, 2] * epsilon) ^ 1.5) 
        comb.score <- z1 + z2
        
      } else {
        
        comb.score <- (newdata[, 1] + newdata[, 2] + apply(newdata, 1, max)) / 2
        
      }     
    } else if (method == "baseinexp") {
      
      comb.score <- newdata[ ,1] ^ newdata[ ,2]
      
    } else if (method == "expinbase") {
      
      comb.score <- newdata[ ,2] ^ newdata[ ,1]
      
    }
    
  }
  
  if(combtype != "mlComb"){

    comb.score <- as.matrix(comb.score)

    labels <- (comb.score > model$ThresholdCombined)
    labels[labels==TRUE] <- "event"
    labels[labels=="FALSE"] <- "not_event"
    comb.score <- data.frame(comb.score, labels)
    colnames(comb.score) <- c("comb.score", "labels")
  }
  
  return(comb.score)
  
}



power.add <- function(n, marker.set){
  
  power1 <- markers[,1] ^ n
  power2 <- markers[,2] ^ n
  
  return (power1 + power2)
}


power.subt <- function(n, marker.set){
  
  power1 <- markers[,1] ^ n
  power2 <- markers[,2] ^ n
  
  return (power1 - power2)
}
