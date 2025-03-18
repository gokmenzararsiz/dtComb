#' @title Combine two diagnostic tests with Machine Learning Algorithms.
#'
#' @description The \code{mlComb} function calculates the combination
#' scores of two diagnostic tests selected among several Machine Learning
#' Algorithms
#'
#' @param markers a \code{numeric} data frame that includes two diagnostic tests
#' results
#'
#' @param status a \code{factor} vector that includes the actual disease
#' status of the patients
#'
#' @param event a \code{character} string that indicates the event in the status
#' to be considered as positive event
#'
#' @param method a \code{character} string specifying the method used for
#' combining the markers. For the available methods see availableMethods()
#'
#' \bold{IMPORTANT}: See https://topepo.github.io/caret/available-models.html
#' for further information about the methods used in this function.
#'
#' @param resample a \code{character} string that indicates the resampling
#' method used while training the model. The available methods are "boot",
#' "boot632", "optimism_boot", "boot_all", "cv", "repeatedcv", "LOOCV", "LGOCV",
#' "none", "oob", "adaptive_cv", "adaptive_boot" and "adaptive_LGOCV". for
#' details of these resampling methods see ?caret::trainControl
#'
#' @param niters a \code{numeric} value that indicates the number of
#' bootstrapped resampling iterations (10, default)
#'
#' @param nfolds a \code{numeric} value that indicates the number of folds for
#' cross validation based resampling methods  (5, default)
#'
#' @param nrepeats a \code{numeric} value that indicates the number of repeats
#' for "repeatedcv" option of resampling methods (3, default)
#'
#' @param preProcess a \code{character} string that indicates the pre-processing
#' options to be applied in the data before training the model. Available
#' pre-processing methods are: "BoxCox", "YeoJohnson", "expoTrans", "center",
#' "scale", "range", "knnImpute", "bagImpute", "medianImpute", "pca", "ica",
#' "spatialSign", "corr", "zv", "nzv", and "conditionalX". For detailed
#' information about the methods see ?caret::preProcess
#'
#' @param B a \code{numeric} value that is the number of bootstrap samples for
#' bagging classifiers, "bagFDA", "bagFDAGCV", "bagEarth" and "bagEarthGCV".
#' (25, default)
#'
#' @param show.plot a \code{logical}. If TRUE, a ROC curve is
#' plotted. Default is TRUE
#'
#' @param direction a \code{character} string determines in which direction the
#' comparison will be made.  ">": if the predictor values for the control group
#' are higher than the values of the case group (controls > cases).
#' "<": if the predictor values for the control group are lower or equal than
#' the values of the case group (controls < cases).
#'
#' @param conf.level a \code{numeric} value to  determine the confidence interval
#' for the ROC curve(0.95, default).
#'
#' @param cutoff.method  a \code{character} string determines the cutoff method
#' for the ROC curve.
#'
#' @param show.result a \code{logical} string indicating whether the results
#' should be printed to the console.
#'
#' @param \dots optional arguments passed to selected classifiers.
#'
#' @return A \code{list} of AUC values, diagnostic statistics,
#' coordinates of the ROC curve for the combination score obtained using
#' Machine Learning Algorithms as well as the given biomarkers individually, a
#' comparison table for the AUC values of individual biomarkers and combination
#' score obtained and the fitted model.
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' # call data
#' data(laparoscopy)
#'
#' # define the function parameters
#' markers <- laparoscopy[, -1]
#' status <- factor(laparoscopy$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#'
#' model <- mlComb(
#'   markers = markers, status = status, event = event,
#'   method = "knn", resample = "repeatedcv", nfolds = 10, nrepeats = 5,
#'   preProcess = c("center", "scale"), direction = "<", cutoff.method = "Youden"
#' )
#'
#' @export


mlComb <- function(markers = NULL,
                   status = NULL,
                   event = NULL,
                   method = NULL,
                   resample = NULL,
                   niters = 5,
                   nfolds = 5,
                   nrepeats = 3,
                   preProcess = NULL,
                   show.plot = TRUE,
                   B = 25,
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

  data <- cbind(status, markers)

  if (is.null(resample)) {
    resample <- "none"
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


  BMethods <- c("bagFDA", "bagFDAGCV", "bagEarth", "bagEarthGCV")

  verboseMethods <- c(
    "gbm",
    "mlpKerasDecay",
    "mlpKerasDecayCost",
    "mlpKerasDropout",
    "mlpKerasDropoutCost",
    "deepboost",
    "hda",
    "mxnet",
    "plsRglm",
    "sda",
    "ORFlog",
    "ORFpls",
    "ORFridge",
    "ORFsvm",
    "bartMachine"
  )

  if (is.null(method) || !(method %in% allMethods[, 1])) {
    stop(paste(
      "The response given method is not available for mlComb function.",
      "See availableMethods function for the list of methods available."
    ))
  }

  if (method %in% BMethods) {
    if (any(resample == "repeatedcv")) {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = resample,
          number = nfolds,
          repeats = nrepeats,
          classProbs = TRUE
        ),
        preProc = preProcess,
        B = B,
        ...
      )
    } else if (resample %in% c("boot", "boot632", "optimism_boot", "boot_all")) {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = resample,
          number = niters,
          classProbs =  TRUE
        ),
        preProc = preProcess,
        B = B,
        ...
      )
    } else if (resample == "none") {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        preProc = preProcess,
        B = B,
        ...
      )
    } else {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = resample,
          number = nfolds,
          classProbs =  TRUE
        ),
        preProc = preProcess,
        B = B,
        ...
      )
    }

    score <- predict(modelFit, newdata = markers, type = "prob")
  } else if (method %in% verboseMethods) {
    if (any(resample == "repeatedcv")) {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = resample,
          number = nfolds,
          repeats = nrepeats,
          classProbs = TRUE
        ),
        preProc = preProcess,
        verbose = FALSE,
        ...
      )
    } else if (resample %in% c("boot", "boot632", "optimism_boot", "boot_all")) {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = resample,
          number = niters,
          classProbs =  TRUE
        ),
        preProc = preProcess,
        verbose = FALSE,
        ...
      )
    } else if (resample == "none") {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        preProc = preProcess,
        verbose = FALSE,
        ...
      )
    } else {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = resample,
          number = nfolds,
          classProbs =  TRUE
        ),
        preProc = preProcess,
        verbose = FALSE,
        ...
      )
    }

    score <- predict(modelFit, newdata = markers, type = "prob")
  } else {
    if (any(resample == "repeatedcv")) {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = resample,
          number = nfolds,
          repeats = nrepeats,
          classProbs = TRUE
        ),
        preProc = preProcess,
        ...
      )
    } else if (resample %in% c("boot", "boot632", "optimism_boot", "boot_all")) {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = resample,
          number = niters,
          classProbs =  TRUE
        ),
        preProc = preProcess,
        ...
      )
    } else if (resample == "none") {
      print(method)
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        preProc = preProcess,
        ...
      )
    } else {
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = resample,
          number = nfolds,
          classProbs =  TRUE
        ),
        preProc = preProcess,
        ...
      )
    }

    score <- predict(modelFit, newdata = markers, type = "prob")
  }

  comb.score <- as.numeric(score[, levels(status) == event])
  status <- factor(ifelse(status == event, 1, 0), ordered = TRUE)

  allres <-
    rocsum(
      markers = markers,
      comb.score = as.matrix(comb.score),
      status = status,
      event = event,
      direction = direction,
      conf.level = conf.level,
      cutoff.method = cutoff.method,
      show.plot = show.plot
    )

  model_fit <- list(
    CombType = "mlComb",
    Model = modelFit
  )

  allres$fit <- model_fit

  if (show.result) {
    print_model <- list(
      CombType = "mlComb",
      Model = modelFit,
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

###############################################################################
#' @title Available classification/regression methods in \code{dtComb}
#'
#' @description This function returns a data.frame of available classification
#' methods in \code{dtComb}. These methods are imported from the caret package.
#'
#' @return \code{No return value} contains the method names and explanations of the
#' machine-learning models available for the dtComb package.
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#'
#' availableMethods()
#'
#' @export

availableMethods <- function() {
  message(
    paste(
      "The available methods are listed below. For more information",
      "about the methods see https://topepo.github.io/caret/available-models.html"
    )
  )
  print(allMethods)
}
