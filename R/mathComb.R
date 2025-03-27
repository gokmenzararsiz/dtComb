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
#'  \item \bold{Euclidean} (\code{euclidean}):
#'  \eqn{c_i = {\sqrt{(x_{i1}-0)^2+(x_{i2}-0)^2}}}
#'  \item \bold{Manhattan}(\code{manhattan}):
#'  \eqn{c_i = |x_{i1}-0|+|x_{i2}-0|}
#'  \item \bold{Chebyshev} (\code{chebyshev}):
#'  \eqn{c_i = max{|x_{i1}-0|,|x_{i2}-0|}}
#'  \item \bold{Kulczynski} (\code{kulczynski_d}):
#'  \eqn{c_i = \frac{|x_{i1}-0|+|x_{i2}-0|}{min(x_{i1},x_{i2})}}
#'  \item \bold{Lorentzian} (\code{lorentzian}):
#'  \eqn{c_i = (ln(1+|x_{i1}-0|))+ (ln(1+|x_{i2}-0|))}
#'  \item \bold{Taneja} (\code{taneja}):
#'  \eqn{c_i = z_1\times\Biggl(log\frac{z_1}{\sqrt{(x_{i1}\times \epsilon )}}\Biggl)+z_2\times\Biggl(log\frac{z_2}{\sqrt{(x_{i2}\times\epsilon)}}\Biggl)}
#'  \item \bold{Kumar-Johnson} (\code{kumar-johnson}):
#'  \eqn{c_i = {\frac{(x_{i1}-0)^2}{2(x_{i1}\times\epsilon)}}+{\frac{(x_{i2}-0)^2}{2(x_{i2}\times\epsilon)}}, \epsilon = 0.00001}
#'  \item \bold{Avg} (\code{avg}):
#'  \deqn{(L_1, L_n) = \frac{|x_{i1}-0|+|x_{i2}-0| + max{(x_{i1}-0),(x_{i2}-0)}}{2}}
#' }
#'
#' @param standardize a \code{character} string indicating the name of the
#'  standardization method. The default option is no standardization applied.
#'  Available options are:
#'  \itemize{
#' \item \bold{Z-score} \code{(zScore)}: This method scales the data to have a mean
#' of 0 and a standard deviation of 1. It subtracts the mean and divides by the standard
#'  deviation for each feature. Mathematically,
#'  \deqn{ Z-score = \frac{x - (\overline x)}{sd(x)}}
#'
#'   where \eqn{x} is the value of a marker, \eqn{\overline{x}} is the mean of the marker and \eqn{sd(x)} is the standard deviation of the marker.
#' \item \bold{T-score} \code{(tScore)}: T-score is commonly used
#' in data analysis to transform raw scores into a standardized form.
#'  The standard formula for converting a raw score \eqn{x} into a T-score is:
#'  \deqn{T-score = \Biggl(\frac{x - (\overline x)}{sd(x)}\times 10 \Biggl) +50}
#'   where \eqn{x} is the value of a marker, \eqn{\overline{x}} is the mean of the marker
#'    and \eqn{sd(x)} is the standard deviation of the marker.
#'
#' \item \bold{min_max_scale} \code{(min_max_scale)}: This method transforms data to
#' a specific scale, between 0 and 1. The formula for this method is:
#' \deqn{min_max_scale = \frac{x - min(x)}{max(x) - min(x)}}
#'
#' \item \bold{scale_mean_to_one} \code{(scale_mean_to_one)}: This method scales
#' the arithmetic mean to 1. The formula for this method is:
#' \deqn{scale_mean_to_one =  \frac{x}{\overline{x}}}
#' where \eqn{x} is the value of a marker and \eqn{\overline{x}} is the mean of the marker.
#'
#' \item \bold{scale_sd_to_one} \code{(scale_sd_to_one)}: This method, which allows for
#' comparison of individual data points in relation to the overall spread of
#' the data, scales the standard deviation to 1. The formula for this method is:
#' \deqn{scale_sd_to_one = \frac{x}{sd(x)}}
#' where \eqn{x} is the value of a marker and \eqn{sd(x)} is the standard deviation of the marker.
#' }
#'
#' @param transform A \code{character} string specifying the mathematical
#' transformation method.
#' Available options: \code{"log"}, \code{"exp"}, \code{"sin"}, \code{"cos"}.
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
#' @param show.plot a \code{logical}. If TRUE, a ROC curve is
#' plotted. Default is TRUE
#'
#' @param direction a \code{character} string determines in which direction the
#'  comparison will be made.  ">": if the predictor values for the control group
#'  are higher than the values of the case group (controls > cases).
#'  "<": if the predictor values for the control group are lower or equal than
#'  the values of the case group (controls < cases).
#'
#' @param conf.level a \code{numeric} values determines the confidence interval
#'  for the roc curve(0.95, default).
#'
#' @param cutoff.method  a \code{character} string determines the cutoff method
#'  for the roc curve.
#'
#' @param show.result A \code{logical} value indicating whether the results
#' should be printed to the console.
#'
#' @param \dots further arguments. Currently has no effect on the results.
#'
#' @return A list containing the computed combination scores and, optionally,
#'  diagnostic performance metrics.
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#'
#' data(laparotomy)
#' markers <- laparotomy[, -1]
#' status <- factor(laparotomy$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#' direction <- "<"
#' cutoff.method <- "Youden"
#'
#' score1 <- mathComb(
#'   markers = markers, status = status, event = event,
#'   method = "distance", distance = "avg", direction = direction, show.plot = FALSE,
#'   standardize = "none", cutoff.method = cutoff.method
#' )
#'
#' score2 <- mathComb(
#'   markers = markers, status = status, event = event,
#'   method = "baseinexp", transform = "exp", direction = direction,
#'   cutoff.method = cutoff.method
#' )
#'
#' score3 <- mathComb(
#'   markers = markers, status = status, event = event,
#'   method = "subtract", direction = "auto", cutoff.method = "MinValueSp", transform = "sin"
#' )
#'
#' @export

mathComb <- function(markers = NULL,
                     status = NULL,
                     event = NULL,
                     method = c(
                       "add",
                       "multiply",
                       "divide",
                       "subtract",
                       "distance",
                       "baseinexp",
                       "expinbase"
                     ),
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
                     standardize = c(
                       "none", "min_max_scale",
                       "zScore", "tScore", "scale_mean_to_one", "scale_sd_to_one"
                     ),
                     transform = c("none", "log", "exp", "sin", "cos"),
                     show.plot = TRUE,
                     direction = c("auto", "<", ">"),
                     conf.level = 0.95,
                     cutoff.method = c(
                       "CB", "MCT", "MinValueSp", "MinValueSe", "ValueSp",
                       "ValueSe", "MinValueSpSe", "MaxSp", "MaxSe",
                       "MaxSpSe", "MaxProdSpSe", "ROC01", "SpEqualSe",
                       "Youden", "MaxEfficiency", "Minimax", "MaxDOR",
                       "MaxKappa", "MinValueNPV", "MinValuePPV", "ValueNPV",
                       "ValuePPV", "MinValueNPVPPV", "PROC01", "NPVEqualPPV",
                       "MaxNPVPPV", "MaxSumNPVPPV", "MaxProdNPVPPV",
                       "ValueDLR.Negative", "ValueDLR.Positive", "MinPvalue",
                       "ObservedPrev", "MeanPrev", "PrevalenceMatching"
                     ), show.result = FALSE, ...) {
  methods <-
    c(
      "add",
      "multiply",
      "divide",
      "subtract",
      "distance",
      "baseinexp",
      "expinbase"
    )

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
    c("none", "min_max_scale", "zScore", "tScore", "scale_mean_to_one", "scale_sd_to_one")

  transforms <- c("none", "log", "exp", "sin", "cos")

  directions <- c("auto", "<", ">")

  cutoff.methods <- c(
    "CB", "MCT", "MinValueSp", "MinValueSe", "ValueSp",
    "ValueSe", "MinValueSpSe", "MaxSp", "MaxSe",
    "MaxSpSe", "MaxProdSpSe", "ROC01", "SpEqualSe",
    "Youden", "MaxEfficiency", "Minimax", "MaxDOR",
    "MaxKappa", "MinValueNPV", "MinValuePPV", "ValueNPV",
    "ValuePPV", "MinValueNPVPPV", "PROC01", "NPVEqualPPV",
    "MaxNPVPPV", "MaxSumNPVPPV", "MaxProdNPVPPV",
    "ValueDLR.Negative", "ValueDLR.Positive", "MinPvalue",
    "ObservedPrev", "MeanPrev", "PrevalenceMatching"
  )

  raw.markers <- markers
  raw.status <- status

  if (!is.data.frame(markers)) {
    markers <- as.data.frame(markers)
  }

  for (i in 1:ncol(markers)) {
    if (!is.numeric(markers[, i])) {
      stop("at least one variable is not numeric")
    }
  }

  if (!ncol(markers) == 2) {
    stop("the number of markers should be 2")
  }

  if (!is.factor(status)) {
    status <- as.factor(status)
  }

  if (!length(levels(status)) == 2) {
    stop("the number of status levels should be 2")
  }

  if (!(event %in% status)) {
    stop("status does not include event")
  }

  levels(status)[levels(status) == "NA"] <- NA

  if (nrow(markers) != length(status)) {
    stop(
      paste(
        "the number of rows of markers is not equal to the number of",
        "elements of the status"
      )
    )
  }

  status_levels <- levels(status)
  if (status_levels[1] == event) {
    firstStatus <- status_levels[1]
    secondStatus <- status_levels[2]
    status_levels[1] <- secondStatus
    status_levels[2] <- firstStatus
  }
  status <- factor(ifelse(status == event, 1, 0), ordered = TRUE)

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

  if (length(which(methods == method)) == 0 || length(method) != 1) {
    stop(
      paste(
        "method should be one of 'add', 'multiply', 'divide', 'subtract'",
        ",'distance', 'baseinexp', 'expinbase'"
      )
    )
  }

  if (method != "distance") {
    distance <- NULL
  } else {
    if (length(which(distances == distance)) == 0 ||
      length(distance) != 1) {
      stop(
        paste(
          "distance should be one of 'euclidean', 'manhattan',",
          "'chebyshev', 'kulczynski_d', 'lorentzian', 'avg', 'taneja',",
          "'kumar-johnson'"
        )
      )
    }
  }

  if (length(which(standardizes == standardize)) == 0) {
    stop(
      paste(
        "standardize should be one of 'min_max_scale', 'zScore', 'tScore',",
        "'scale_mean_to_one', 'scale_sd_to_one'"
      )
    )
  }

  if (length(which(transforms == transform)) == 0) {
    stop("transforms should be one of 'none', 'log', 'exp', 'sin', 'cos'")
  }

  if (length(which(directions == direction)) == 0) {
    stop("direction should be one of 'auto', '<', '>'")
  }

  if (length(direction) != 1) {
    warning("Direction is set to 'auto'")
  }

  if (length(which(cutoff.methods == cutoff.method)) == 0 ||
    length(cutoff.method) != 1) {
    stop("The entered cutoff.method is invalid")
  }

  if (length(standardize) != 1) {
    standardize <- "none"
  }
  if (length(transform) != 1) {
    transform <- "none"
  }

  std.model <- std.train(markers, standardize)
  markers <- std.model$data

  markers <- transform_math(markers, transform)


  n <- as.matrix(seq(-3, 3, 0.1))
  p <- seq(1:nrow(n))

  get_roc <- function(p) {
    tryCatch(
      {
        values <- suppressMessages(pROC::roc(status, power[, p],
          direction = direction
        ))
        auc <- values$auc

        return(auc)
      },
      error = function(cond) {
        return(0)
      }
    )
  }
  power.add <- function(n) {
    power1 <- markers[, 1]^n
    power2 <- markers[, 2]^n

    return(power1 + power2)
  }
  power.subt <- function(n) {
    power1 <- markers[, 1]^n
    power2 <- markers[, 2]^n

    return(power1 - power2)
  }
  max_power <- NA

  if (method == "add") {
    power <- apply(n, 1, power.add)

    auc_list <- sapply(p, get_roc)
    max_index <- which(auc_list == max(auc_list))

    if (length(max_index) > 1) {
      max_index <- max_index[1]
    }
    max_power <- (n[max_index])
    comb.score <- power[, max_index]
  } else if (method == "multiply") {
    comb.score <- markers[, 1] * markers[, 2]
  } else if (method == "divide") {
    comb.score <- markers[, 1] / markers[, 2]
  } else if (method == "subtract") {
    power <- apply(n, 1, power.subt)

    auc_list <- sapply(p, get_roc)
    max_index <- which(auc_list == max(auc_list))

    if (length(max_index) > 1) {
      max_index <- max_index[1]
    }
    max_power <- (n[max_index])
    comb.score <- power[, max_index]
  } else if (method == "distance") {
    if (distance == "euclidean") {
      comb.score <- sqrt((markers[, 1]^2) + (markers[, 2]^2))
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
        (z1 * (log(z1 / sqrt(markers[, 1] * epsilon)))) +
        (z2 * (log(z2 / sqrt(markers[, 2] * epsilon))))
    } else if (distance == "kumar-johnson") {
      epsilon <- 0.00001
      z1 <- ((markers[, 1]^2)^2) /
        (2 * ((markers[, 1] * epsilon)^1.5))
      z2 <- ((markers[, 2]^2)^2) /
        (2 * ((markers[, 2] * epsilon)^1.5))
      comb.score <- z1 + z2
    } else {
      comb.score <-
        (markers[, 1] + markers[, 2] + apply(markers, 1, max)) / 2
    }
  } else if (method == "baseinexp") {
    comb.score <- markers[, 1]^markers[, 2]
  } else if (method == "expinbase") {
    comb.score <- markers[, 2]^markers[, 1]
  }
  if ((length(which(is.infinite(comb.score))) ||
    length(which(is.nan(comb.score)))) > 0 &&
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
    MaxPower = max_power,
    Classification = status_levels,
    Std.model = std.model$std,
    Standardize = standardize
  )

  allres$fit <- model_fit

  if (show.result) {
    print_model <- list(
      CombType = "mathComb",
      Method = method,
      Distance = distance,
      rowcount = nrow(markers),
      colcount = ncol(markers),
      classification = status_levels,
      Pre_processing = standardize,
      Transform = transform,
      MaxPower = max_power,
      AUC_table = allres$AUC_table,
      MultComp_table = allres$MultComp_table,
      DiagStatCombined = allres$DiagStatCombined,
      Cutoff_method = cutoff.method,
      ThresholdCombined = allres$ThresholdCombined,
      Criterion = allres$Criterion.c
    )
    print_train(print_model)
  }
  invisible(allres)
}


#' @title Mathematical transformations for biomarkers
#'
#' @description Applies a selected mathematical transformation (\code{"log"}, \code{"exp"}, \code{"sin"}, or \code{"cos"}) to biomarker data.
#'
#' @param markers A \code{numeric} data frame that contains the biomarkers.
#'
#' @param transform A character string specifying the transformation method to be applied.
#' Supported methods are: \code{"log"}, \code{"exp"}, \code{"sin"}, \code{"cos"}.
#'
#' @return A \code{numeric} data frame containing the transformed biomarkers.
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' data(laparotomy)
#' markers <- laparotomy[, -1]
#' transform_math(markers, transform = "log")
#'
#' @export


transform_math <- function(markers, transform) {
  if (any(transform == "none")) {
    markers <- markers
    transform <- "none"
  } else if (any(transform == "log")) {
    markers <- log(markers)
  } else if (any(transform == "exp")) {
    markers <- exp(markers)
  } else if (any(transform == "sin")) {
    markers <- sin(markers)
  } else {
    markers <- cos(markers)
  }
  return(markers)
}
