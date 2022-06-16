# TODO: Add comment
#
# Author: serra ilayda yerlitas
###############################################################################
#' @title Combine two diagnostic tests with several mathematical operators and
#'  distance measures.
#'
#' @description The \code{mathComb} function returns the combination results of
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
#'  }
#' @param distance a \code{character} string specifying the method used for
#'  combining the markers. The available methods are:
#'  \itemize{
#'  \item \code{euclidean}:
#'  Euclidean: d = sqrt( ∑ | P_i - Q_i |^2)
#'  \item \code{manhattan}:
#'  Manhattan : d = ∑ | P_i - Q_i |
#'  \item \code{chebyshev}:
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
#' @param show.plot a \code{logical} a \code{logical}. If TRUE, a ROC curve is 
#' plotted. Default is TRUE
#' 
#' @param direction a \code{character} string determines in which direction the
#'  comparison will be made.  “>”: if the predictor values for the control group
#'  are higher than the values of the case group (controls > cases).
#'  “<”: if the predictor values for the control group are lower or equal than
#'  the values of the case group (controls < cases).
#'
#' @param conf.level a \code{numeric} values determines the confidence interval
#'  for the roc curve(0.95, default).
#'
#' @param cutoff.method  a \code{character} string determines the cutoff method
#'  for the roc curve.
#'  
#' @param \dots further arguments. Currently has no effect on the results.
#'
#' @return A list of \code{numeric} mathematical combination scores calculated
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
#' method = "distance", distance = "avg", direction = direction, show.plot = FALSE,
#' standardize = "none", cutoff.method = cutoff.method)
#'
#' score2 <- mathComb(markers = markers, status = status, event = event,
#' method = "baseinexp", transform = "exp", direction = direction,
#' cutoff.method = cutoff.method)
#'
#' score3 <- mathComb(markers = markers, status = status, event = event,
#' method = "add",  direction = direction, power.transform = TRUE,
#' cutoff.method = cutoff.method)
#'
#' @export

mathComb <- function(markers = NULL,
                     status = NULL,
                     event = NULL,
                     
                     method = c("add",
                                "multiply",
                                "divide",
                                "subtract",
                                "distance",
                                "baseinexp",
                                "expinbase"),
                     distance = c(
                       "euclidean",
                       "manhattan",
                       "chebyshev",
                       "kulczynski_d",
                       "lorentzian",
                       "avg",
                       "taneja",
                       "kumar-johnson"
                     ),
                     standardize = c("none", "range",
                                     "zScore", "tScore", "mean", "deviance"),
                     transform = c("none", "log", "exp", "sin", "cos"),
                     power.transform = FALSE, show.plot = TRUE,
                     direction = c("auto", "<", ">"),
                     conf.level = 0.95,
                     cutoff.method = c("youden", "roc01"), ...) {
  methods <-
    c("add",
      "multiply",
      "divide",
      "subtract",
      "distance",
      "baseinexp",
      "expinbase")
  
  distances <-
    c(
      "euclidean",
      "manhattan",
      "chebyshev",
      "kulczynski_d",
      "lorentzian",
      "avg",
      "taneja",
      "kumar-johnson"
    )
  
  standardizes <-
    c("none", "range", "zScore", "tScore", "mean", "deviance")
  
  transforms <-  c("none", "log", "exp", "sin", "cos")
  
  directions <- c("auto", "<", ">")
  
  cutoff.methods <- c("youden", "roc01")
  
  raw.markers <- markers
  raw.status <- status
  
  if (!is.data.frame(markers)) {
    markers <- as.data.frame(markers)
  }
  
  for (i in 1:ncol(markers))
    if (!is.numeric(markers[, i]))
      stop("at least one variable is not numeric")
  
  if (!ncol(markers) == 2)
    stop("the number of markers should be 2")
  
  if (!is.factor(status))
    status <- as.factor(status)
  
  if (!length(levels(status)) == 2)
    stop("the number of status levels should be 2")
  
  if(!(event %in% status))
    stop("status does not include event")
  
  levels(status)[levels(status) == "NA"] <- NA
  
  if (nrow(markers) != length(status))
    stop(
      paste(
        "the number of rows of markers is not equal to the number of",
        "elements of the status"
      )
    )
  
  status_levels <- levels(status)
  status <- factor(ifelse(status == event, 1, 0))
  
  if (length(which(is.na(markers))) > 0) {
    comp <- complete.cases(markers)
    markers <- markers[comp, ]
    status <- status[comp]
    warning(paste(
      "Rows with NA removed from the dataset since markers",
      "include NA"
    ))
  }
  
  if (length(which(is.na(status))) > 0) {
    comp <- complete.cases(status)
    status <- status[comp]
    markers <- markers[comp, ]
    warning(paste(
      "Rows with NA removed from the dataset since status",
      "include NA"
    ))
  }
  
  if (length(which(methods == method)) == 0 || length(method) != 1)
    stop(
      paste(
        "method should be one of “add”, “multiply”, “divide”, “subtract”"
        ,
        ",“distance”, “baseinexp”, “expinbase”"
      )
    )
  
  if (method != "distance") {
    distance <- NULL
    
  } else {
    if (length(which(distances == distance)) == 0 ||
        length(distance) != 1) {
      stop(
        paste(
          "distance should be one of “euclidean”, “manhattan”,",
          "“chebyshev”, “kulczynski_d”, “lorentzian”, “avg”, “taneja”,",
          "“kumar-johnson”"
        )
      )
    }
    
  }
  
  if (length(which(standardizes == standardize)) == 0)
    stop(
      paste(
        "standardize should be one of “range”, “zScore”, “tScore”,",
        "“mean”, “deviance”"
      )
    )
  
  if (length(which(transforms == transform)) == 0)
    stop("transforms should be one of “none”, “log”, “exp”, “sin”, “cos”")
  
  if (length(which(directions == direction)) == 0 )
    stop("direction should be one of “auto”, “<”, “>”")
  
  if(length(direction) != 1)
    warning("Direction is set to “auto”")
  
  if (length(which(cutoff.methods == cutoff.method)) == 0 ||
      length(cutoff.method) != 1)
    stop("cutoff.method should be one of “youden”, “roc01”")
  
  if (power.transform != TRUE) {
    max_power <- NULL
    
  }
  
  if (length(standardize) != 1) {
    standardize <- "none"
  }
  if (length(transform) != 1) {
    transform <- "none"
  }
  
  std.model <- std.train(markers, standardize)
  markers <- std.model$data
  
  markers <- transform.math(markers, transform)
  
  
  n <- as.matrix(seq(-3, 3, 0.1))
  p <- seq(1:nrow(n))
  
  get_roc <- function(p) {
    values <- suppressMessages(pROC::roc(status , power[, p],
                                         direction = direction))
    auc <- values$auc
    
    return(auc)
    
  }
  if (method == "add") {
    if (power.transform == TRUE) {
      power <- apply(n, 1, power.add)
      
      auc_list <- sapply(p, get_roc)
      max_index <- which(auc_list == max(auc_list))
      
      if (length(max_index) > 1) {
        max_index <- max_index[1]
        
      }
      max_power <- (n[max_index])
      comb.score <- power[, max_index]
    }
    else {
      comb.score <- markers[, 1] + markers[, 2]
    }
    
  } else if (method == "multiply") {
    comb.score <- markers[, 1] * markers[, 2]
    
  } else if (method == "divide") {
    comb.score <- markers[, 1] / markers[, 2]
    
  } else if (method == "subtract") {
    if (power.transform == TRUE) {
      power <- apply(n, 1, power.subt)
      
      auc_list <- sapply(p, get_roc)
      max_index <- which(auc_list == max(auc_list))
      
      if (length(max_index) > 1) {
        max_index <- max_index[1]
        
      }
      max_power <- (n[max_index])
      comb.score <- power[, max_index]
      
    }
    else{
      comb.score <- (markers[, 1] - markers[, 2])
    }
    
  }  else if (method == "distance") {
    if (distance == "euclidean") {
      comb.score <- sqrt((markers[, 1] ^ 2) + (markers[, 2] ^ 2))
      
    } else if (distance == "manhattan") {
      comb.score <- markers[, 1] + markers[, 2]
      
    } else if (distance == "chebyshev") {
      comb.score <- apply(markers, 1, max)
      
    } else if (distance == "kulczynski_d") {
      comb.score <- (abs(markers[, 1]) + abs(markers[, 2])) / 0.00002
      
    } else if (distance == "lorentzian") {
      comb.score <-
        log(abs(markers[, 1]) + 1) + log(abs(markers[, 2]) + 1)
      
    } else if (distance == "taneja") {
      epsilon <- 0.00001
      z1 <- (markers[, 1]) / 2
      z2 <- (markers[, 2]) / 2
      comb.score <-
        (z1 * (log(z1 / sqrt(
          markers[, 1] * epsilon
        )))) +
        (z2 * (log(z2 / sqrt(
          markers[, 2] * epsilon
        ))))
      
    } else if (distance == "kumar-johnson") {
      epsilon <- 0.00001
      z1 <- ((markers[, 1] ^ 2) ^ 2) /
        (2 * ((markers[, 1] * epsilon) ^ 1.5))
      z2 <- ((markers[, 2] ^ 2) ^ 2) /
        (2 * ((markers[, 2] * epsilon) ^ 1.5))
      comb.score <- z1 + z2
      
    } else {
      comb.score <-
        (markers[, 1] + markers[, 2] + apply(markers, 1, max)) / 2
      
    }
  } else if (method == "baseinexp") {
    comb.score <- markers[, 1] ^ markers[, 2]
    
  } else if (method == "expinbase") {
    comb.score <- markers[, 2] ^ markers[, 1]
    
  }
  if ((length(which(is.infinite(comb.score))) ||
       length(which(is.nan(comb.score))))  > 0 &&
      standardize != "none") {
    warning("Infinity or NaNs values generated in markers, standardization changed to 'none'.")
    return(
      mathComb(
        markers = raw.markers,
        status = raw.status,
        event = event,
        method = method,
        distance = distance,
        direction = direction,
        standardize = "none",
        cutoff.method = cutoff.method,
        transform = transform,
        power.transform = power.transform,
        conf.level = conf.level
      )
    )
  }
  
  if ((length(which(is.infinite(comb.score))) ||
       length(which(is.nan(comb.score)))) > 0 &&
      transform != "none") {
    warning("Infinity or NaNs values generated in markers, transformation changed to 'none'.")
    return(
      mathComb(
        markers = raw.markers,
        status = raw.status,
        event = event,
        method = method,
        distance = distance,
        direction = direction,
        standardize = standardize,
        cutoff.method = cutoff.method,
        transform = "none",
        power.transform = power.transform,
        conf.level = conf.level
      )
    )
  }
  
  comb.score <- as.matrix(comb.score)
  
  allres <-
    rocsum(
      markers = markers,
      comb.score = comb.score,
      status = status,
      event = event,
      direction = direction,
      conf.level = conf.level,
      cutoff.method = cutoff.method,
      show.plot = show.plot
    )
  
  model_fit <- list(
    CombType = "mathComb",
    Method = method,
    Distance = distance,
    Transform = transform,
    PowerTransform = power.transform,
    MaxPower = max_power,
    Classification = status_levels,
    Std.model = std.model$std,
    Standardize = standardize
  )
  
  allres$fit <- model_fit
  
  print_model = list(
    CombType = "mathComb",
    Method = method,
    Distance = distance,
    rowcount = nrow(markers),
    colcount = ncol(markers),
    classification = status_levels,
    Pre_processing = standardize,
    Transform = transform,
    PowerTransform = power.transform,
    MaxPower = max_power,
    AUC_table = allres$AUC_table,
    MultComp_table = allres$MultComp_table,
    DiagStatCombined = allres$DiagStatCombined
  )
  
  print_train(print_model)
  
  invisible(allres)
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
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec
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

power.add <- function(n, marker.set) {
  power1 <- markers[, 1] ^ n
  power2 <- markers[, 2] ^ n
  
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
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec
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

power.subt <- function(n, marker.set) {
  power1 <- markers[, 1] ^ n
  power2 <- markers[, 2] ^ n
  
  return (power1 - power2)
}



#'
#'
#'
#'
#'
#'
#'
#' @export
transform.math <- function(markers, transform) {
  if (any(transform == "none")) {
    markers <- markers
    transform <- "none"
    
  }
  else if (any(transform == "log")) {
    markers <- log(markers)
  }
  else if (any(transform == "exp")) {
    markers <- exp(markers)
    
  }
  else if (any(transform == "sin")) {
    markers <- sin(markers)
    
  }
  else {
    markers <- cos(markers)
    
  }
  return(markers)
}