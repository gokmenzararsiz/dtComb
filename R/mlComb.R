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
#' @param resample A character string specifying the resampling method used while
#' training the model. Available methods are \code{"boot"}, \code{"boot632"},
#' \code{"optimism_boot"}, \code{"boot_all"}, \code{"cv"},
#' \code{"repeatedcv"}, \code{"LOOCV"}, \code{"LGOCV"}, \code{"none"},
#' \code{"oob"}, \code{"adaptive_cv"}, \code{"adaptive_boot"}, and
#' \code{"adaptive_LGOCV"}. The \code{"timeslice"} method is not supported by
#' \code{mlComb}. If \code{NULL}, \code{"none"} is used.
#'
#' @param niters A positive integer indicating the number of bootstrap
#' resampling iterations. Used for \code{"boot"}, \code{"boot632"},
#' \code{"optimism_boot"}, \code{"boot_all"}, and \code{"adaptive_boot"}.
#' Default is 10.
#'
#' @param nfolds A positive integer. For \code{"cv"}, \code{"repeatedcv"},
#' and \code{"adaptive_cv"}, it indicates the number of folds. For
#' \code{"LGOCV"} and \code{"adaptive_LGOCV"}, it indicates the number of
#' repeated training/test splits. Default is 5.
#'
#' @param nrepeats A positive integer indicating the number of repeats for
#' \code{"repeatedcv"}. Default is 3.
#'
#' @param p A numeric value between 0 and 1 specifying the training proportion
#' used for \code{"LGOCV"} and \code{"adaptive_LGOCV"}. Default is 0.75.
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
#' @param show.plot a \code{logical}. If TRUE, a ROC curve is
#' plotted. Default is TRUE
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
#' @author Serra Ilayda Yerlitas Tastan, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @seealso \code{\link{availableMethods}}, \code{caret::train},
#' \code{caret::trainControl}
#'
#' @examples
#' # call data
#' data(laparotomy)
#'
#' # define the function parameters
#' markers <- laparotomy[, -1]
#' status <- factor(laparotomy$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#'
#' model <- mlComb(
#'   markers = markers, status = status, event = event,
#'   method = "knn", resample = "cv", nfolds = 5,
#'   preProcess = c("center", "scale"), direction = "<", cutoff.method = "Youden"
#' )
#'
#' @export


mlComb <- function(markers = NULL,
                   status = NULL,
                   event = NULL,
                   method = NULL,
                   resample = NULL,
                   niters = 10,
                   nfolds = 5,
                   nrepeats = 3,
                   preProcess = NULL,
                   show.plot = TRUE,
                   B = 25,
                   direction = c("<", ">"),
                   conf.level = 0.95,
                   p = 0.75,
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
  if (is.null(resample)) {
    resample <- "none"
  }

  params <- validateParameters(
    markers = markers,
    status = status,
    event = event,
    method = method,
    direction = direction,
    conf.level = conf.level,
    cutoff.method = cutoff.method,
    show.plot = show.plot,
    show.result = show.result
  )

  markers <- params$markers
  status <- params$status
  event <- params$event
  method <- params$method
  direction <- params$direction
  conf.level <- params$conf.level
  cutoff.method <- params$cutoff.method
  show.plot <- params$show.plot
  show.result <- params$show.result

  validateMlComb(
    niters = niters,
    nfolds = nfolds,
    nrepeats = nrepeats,
    preProcess = preProcess,
    B = B,
    p = p,
    resample = resample,
    method = method
  )

  data <- cbind(status, markers)

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

  if (method %in% BMethods) {
    if (resample %in% c("repeatedcv", "adaptive_cv")) {
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
    } else if (resample %in% c("boot", "boot632", "optimism_boot", "boot_all", "adaptive_boot")) {
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
        trControl = caret::trainControl(
          method = "none",
          classProbs = TRUE
        ),
        preProc = preProcess,
        B = B,
        ...
      )
    } else {
      ctrl <- if (resample %in% c("LGOCV", "adaptive_LGOCV")) {
        caret::trainControl(method = resample, number = nfolds, p = p, classProbs = TRUE)
      } else if (identical(resample, "LOOCV")) {
        caret::trainControl(method = resample, classProbs = TRUE)
      } else if (resample == "oob") {
        caret::trainControl(method = "oob", classProbs = TRUE)
      } else {
        caret::trainControl(method = resample, number = nfolds, classProbs = TRUE)
      }
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = ctrl,
        preProc = preProcess,
        B = B,
        ...
      )
    }

    score <- tryCatch(
      predict(modelFit, newdata = markers, type = "prob"),
      error = function(e) {
        stop(
          paste(
            "The selected method does not provide class probabilities.",
            "Please choose a caret method that supports type = 'prob'."
          )
        )
      }
    )
  } else if (method %in% verboseMethods) {
    if (resample %in% c("repeatedcv", "adaptive_cv")) {
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
    } else if (resample %in% c("boot", "boot632", "optimism_boot", "boot_all", "adaptive_boot")) {
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
        trControl = caret::trainControl(
          method = "none",
          classProbs = TRUE
        ),
        preProc = preProcess,
        verbose = FALSE,
        ...
      )
    } else {
      ctrl <- if (resample %in% c("LGOCV", "adaptive_LGOCV")) {
        caret::trainControl(method = resample, number = nfolds, p = p, classProbs = TRUE)
      } else if (identical(resample, "LOOCV")) {
        caret::trainControl(method = resample, classProbs = TRUE)
      } else if (resample == "oob") {
        caret::trainControl(method = "oob", classProbs = TRUE)
      } else {
        caret::trainControl(method = resample, number = nfolds, classProbs = TRUE)
      }
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = ctrl,
        preProc = preProcess,
        verbose = FALSE,
        ...
      )
    }

    score <- tryCatch(
      predict(modelFit, newdata = markers, type = "prob"),
      error = function(e) {
        stop(
          paste(
            "The selected method does not provide class probabilities.",
            "Please choose a caret method that supports type = 'prob'."
          )
        )
      }
    )
  } else {
    if (resample %in% c("repeatedcv", "adaptive_cv")) {
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
    } else if (resample %in% c("boot", "boot632", "optimism_boot", "boot_all", "adaptive_boot")) {
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
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = caret::trainControl(
          method = "none",
          classProbs = TRUE
        ),
        preProc = preProcess,
        ...
      )
    } else {
      ctrl <- if (resample %in% c("LGOCV", "adaptive_LGOCV")) {
        caret::trainControl(method = resample, number = nfolds, p = p, classProbs = TRUE)
      } else if (identical(resample, "LOOCV")) {
        caret::trainControl(method = resample, classProbs = TRUE)
      } else if (resample == "oob") {
        caret::trainControl(method = "oob", classProbs = TRUE)
      } else {
        caret::trainControl(method = resample, number = nfolds, classProbs = TRUE)
      }
      modelFit <- caret::train(
        status ~ .,
        data = data,
        method = method,
        trControl = ctrl,
        preProc = preProcess,
        ...
      )
    }

    score <- tryCatch(
      predict(modelFit, newdata = markers, type = "prob"),
      error = function(e) {
        stop(
          paste(
            "The selected method does not provide class probabilities.",
            "Please choose a caret method that supports type = 'prob'."
          )
        )
      }
    )
  }

  if (!(event %in% colnames(score))) {
    stop("The predicted probability table does not include the event class")
  }
  comb.score <- as.numeric(score[, event])
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
#' @author Serra Ilayda Yerlitas Tastan, Serra Bersan Gengec, Necla Kochan,
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
