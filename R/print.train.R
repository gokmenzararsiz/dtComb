print.train <- function(print_model) {


 cat("Method:", print_model$Method,"\n")
  
if(print_model$Method == "distance"){
  
  cat(paste("Distance", print_model$Distance, sep = ": "),"\n")
  
}
 cat(paste("Samples", print_model$rowcount, sep = ": "),"\n")
 cat(paste("Markers", print_model$colcount, sep = ": "),"\n")
 cat("Event:", paste(print_model$classification,collapse = ", "),"\n")
 cat(paste("Standardization", print_model$Pre_processing, sep = ": "),"\n")
 
 if(print_model$CombType == "mathComb"){
   
   cat(paste("Transform", print_model$Transform, sep = ": "),"\n")
   if(print_model$PowerTransform == TRUE)
   cat(paste("MaxPower", print_model$MaxPower, sep = ": "),"\n")
   
 }

if(print_model$CombType != "mathComb"){
  if(print_model$Resampling == "boot"){
   
   cat("Resampling: boot (niters:",paste(print_model$niters,")",sep=""))
   
  }else if(print_model$Resampling == "cv"){
   
   cat("Resampling: cv (nfolds:",paste(print_model$nfolds,")",sep=""))
   
   
  }else if(print_model$Resampling == "repeatedcv"){
   
   cat("Resampling: repeatedcv (nfolds:",print_model$nfolds,",","nrepeats:",paste(print_model$nrepeats,")",sep=""))
   
  }
 }
 
 cat("\n")
 cat("\n")
 cat("Area Under the Curves of markers and combined score: ","\n")
 print(print_model$AUC_table)
 cat("\n")
 
 cat("Area Under the Curve comparison of markers and combined score: ","\n")
 print(print_model$MultComp_table)
 cat("\n")
 cat("Confusion matrix: ","\n")
 print(print_model$DiagStatCombined)
 cat("\n")
 
 cat(" Kappa ","   "," Accuracy ","\n",print_model$Kappa," ",print_model$Accuracy)
 cat("\n")
 
 
}





