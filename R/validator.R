#' @title Validate and prepare common diagnostic-combination parameters
#'
#' @description
#' Internal helper used to validate and prepare parameters shared by diagnostic
#' combination functions. It checks marker data, disease status, event class,
#' ROC direction, confidence level, cutoff method, and display options. It also
#' converts \code{markers} to a data frame and \code{status} to a factor when
#' needed.
#'
#' @param markers A numeric data frame or matrix with exactly two columns. Each
#' column should contain one diagnostic marker.
#'
#' @param status A factor or vector indicating the true disease status. It must
#' contain exactly two non-missing classes, and both classes must be represented
#' in the data.
#'
#' @param event A single non-missing value indicating which level of
#' \code{status} should be treated as the positive/event class.
#'
#' @param method A character string specifying the model or combination method.
#' This argument is returned unchanged; method-specific validation should be
#' handled by the calling function.
#'
#' @param direction A character string specifying the ROC direction. Possible
#' values are \code{"<"} and \code{">"}.
#'
#' @param conf.level A numeric value between 0 and 1 specifying the confidence
#' level used for ROC-related summaries.
#'
#' @param cutoff.method A character string specifying the cutoff selection
#' method used for ROC analysis.
#'
#' @param show.plot A logical value indicating whether the ROC plot should be
#' displayed.
#'
#' @param show.result A logical value indicating whether results should be
#' printed to the console.
#'
#' @return
#' A list containing the validated and prepared values:
#' \code{markers}, \code{status}, \code{event}, \code{method},
#' \code{direction}, \code{conf.level}, \code{cutoff.method},
#' \code{show.plot}, and \code{show.result}.
#'
#' @keywords internal

validateParameters <- function(markers = NULL,
                               status = NULL,
                               event = NULL,
                               method = NULL,
                               direction = c("<", ">"),
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
                               ),
                               show.plot = TRUE,
                               show.result = FALSE) {
  directions <- c("<", ">")

  if (length(direction) != 1) {
    warning("Direction is set to '<'")
    direction <- "<"
  }
  if (length(which(directions == direction)) == 0) {
    stop("direction should be one of '<', '>'")
  }

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

  if (is.null(markers)) {
    stop("markers should not be NULL")
  }

  if (!is.data.frame(markers)) {
    markers <- as.data.frame(markers)
  }

  if (ncol(markers) != 2) {
    stop("the number of markers should be 2")
  }

  if (is.null(names(markers)) || any(names(markers) == "")) {
    stop("markers should have non-empty column names")
  }

  if (anyDuplicated(names(markers)) > 0) {
    stop("marker names should be unique")
  }

  if ("status" %in% names(markers)) {
    stop("markers should not include a column named 'status'")
  }
  if (length(cutoff.method) != 1 ||
    length(which(cutoff.methods == cutoff.method)) == 0) {
    stop("The entered cutoff.method is invalid")
  }
  for (i in seq_len(ncol(markers))) {
    if (!is.numeric(markers[, i])) {
      stop("at least one variable is not numeric")
    }
  }

  if (anyNA(markers)) {
    stop("markers should not include missing values")
  }

  if (any(!is.finite(as.matrix(markers)))) {
    stop("markers should include only finite numeric values")
  }

  if (is.null(status)) {
    stop("status should not be NULL")
  }

  if (anyNA(status)) {
    stop("status should not include missing values")
  }

  if (nrow(markers) != length(status)) {
    stop("the number of rows of markers is not equal to the number of elements of the status")
  }

  if (any(as.character(status) == "NA")) {
    stop("status should not include 'NA' as a class label")
  }

  if (!is.factor(status)) {
    status <- as.factor(status)
  }

  if (length(levels(status)) != 2) {
    stop("the number of status levels should be 2")
  }

  if (length(unique(status)) != 2) {
    stop("status should include observations from both classes")
  }

  if (is.null(event) || length(event) != 1 || is.na(event)) {
    stop("event should be a single non-missing value")
  }

  event <- as.character(event)

  if (!(event %in% levels(status))) {
    stop("status does not include event")
  }

  if (!is.logical(show.plot) || length(show.plot) != 1 || is.na(show.plot)) {
    stop("show.plot should be TRUE or FALSE")
  }

  if (!is.logical(show.result) || length(show.result) != 1 || is.na(show.result)) {
    stop("show.result should be TRUE or FALSE")
  }

  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
    is.na(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level should be a numeric value between 0 and 1")
  }

  list(
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
}
###############################################################################
#' @title Validate parameters specific to \code{mlComb}
#'
#' @description
#' Internal helper used by \code{mlComb} to validate machine-learning-specific
#' arguments before fitting the \code{caret} model. It checks the selected
#' method, resampling method, resampling control parameters, preprocessing
#' options, bagging parameter, and \code{oob}/adaptive resampling restrictions.
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
#' @param preProcess A character vector specifying preprocessing methods passed
#' to \code{caret::train}. Available options are \code{"BoxCox"},
#' \code{"YeoJohnson"}, \code{"expoTrans"}, \code{"center"}, \code{"scale"},
#' \code{"range"}, \code{"knnImpute"}, \code{"bagImpute"},
#' \code{"medianImpute"}, \code{"pca"}, \code{"ica"}, \code{"spatialSign"},
#' \code{"corr"}, \code{"zv"}, \code{"nzv"}, and \code{"conditionalX"}.
#' Default is \code{NULL}.
#'
#' @param B A positive integer indicating the number of bootstrap samples for
#' bagging classifiers such as \code{"bagFDA"}, \code{"bagFDAGCV"},
#' \code{"bagEarth"}, and \code{"bagEarthGCV"}. Default is 25.
#'
#' @param p A numeric value between 0 and 1 specifying the training proportion
#' used for \code{"LGOCV"} and \code{"adaptive_LGOCV"}. Default is 0.75.
#'
#' @param resample A character string specifying the resampling method. Allowed
#' values are \code{"boot"}, \code{"boot632"}, \code{"optimism_boot"},
#' \code{"boot_all"}, \code{"cv"}, \code{"repeatedcv"}, \code{"LOOCV"},
#' \code{"LGOCV"}, \code{"none"}, \code{"oob"}, \code{"adaptive_cv"},
#' \code{"adaptive_boot"}, and \code{"adaptive_LGOCV"}. The
#' \code{"timeslice"} method is not supported by \code{mlComb}.
#'
#' @param method A character string specifying the \code{caret} model used by
#' \code{mlComb}. It must be available in \code{allMethods[, 1]}.
#'
#' @details
#' The \code{"oob"} resampling method is restricted to suitable random forest,
#' bagged tree, bagged earth, bagged flexible discriminant analysis, and
#' conditional tree forest methods.
#'
#' For \code{"adaptive_cv"} and \code{"adaptive_LGOCV"}, \code{nfolds} must be
#' greater than 5. For \code{"repeatedcv"}, \code{nrepeats} must be at least 2.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass. Otherwise, an error
#' is thrown.
#'
#' @keywords internal

validateMlComb <- function(niters = 10,
                           nfolds = 5,
                           nrepeats = 3,
                           preProcess = NULL,
                           B = 25,
                           p = 0.75,
                           resample = "none",
                           method = NULL) {
  available_resampling_methods <- c(
    "boot", "boot632", "optimism_boot", "boot_all",
    "cv", "repeatedcv", "LOOCV", "LGOCV",
    "none", "oob",
    "adaptive_cv", "adaptive_boot", "adaptive_LGOCV"
  )

  if (identical(resample, "timeslice")) {
    stop("resample = 'timeslice' is not supported by mlComb")
  }

  available_preprocess_methods <- c(
    "BoxCox", "YeoJohnson", "expoTrans", "center", "scale", "range",
    "knnImpute", "bagImpute", "medianImpute", "pca", "ica",
    "spatialSign", "corr", "zv", "nzv", "conditionalX"
  )

  if (is.null(method) || length(method) != 1 || is.na(method) ||
    !(method %in% allMethods[, 1])) {
    stop(
      paste(
        "The given method is not available for mlComb function.",
        "See availableMethods() for the list of available methods."
      )
    )
  }

  if (is.null(resample) || length(resample) != 1 || is.na(resample) ||
    !(resample %in% available_resampling_methods)) {
    stop(
      paste(
        "resample should be one of:",
        paste(available_resampling_methods, collapse = ", ")
      )
    )
  }

  bootstrap_resampling_methods <- c(
    "boot", "boot632", "optimism_boot", "boot_all", "adaptive_boot"
  )

  if (resample %in% bootstrap_resampling_methods &&
    (!is.numeric(niters) || length(niters) != 1 ||
      is.na(niters) || niters < 1 || niters != as.integer(niters))) {
    stop("niters should be a positive integer for bootstrap-based resampling methods")
  }

  if (resample == "repeatedcv" &&
    (!is.numeric(nrepeats) || length(nrepeats) != 1 ||
      is.na(nrepeats) || nrepeats < 1 || nrepeats != as.integer(nrepeats))) {
    stop("nrepeats should be a positive integer when resample = 'repeatedcv'")
  }

  if (resample == "repeatedcv" && nrepeats < 2) {
    stop("nrepeats should be at least 2 when resample = 'repeatedcv'")
  }

  bagging_methods <- c("bagFDA", "bagFDAGCV", "bagEarth", "bagEarthGCV")

  if (method %in% bagging_methods &&
    (!is.numeric(B) || length(B) != 1 ||
      is.na(B) || B < 1 || B != as.integer(B))) {
    stop("B should be a positive integer for bagging methods")
  }

  if (resample %in% c("cv", "repeatedcv", "LGOCV", "adaptive_cv", "adaptive_LGOCV") &&
    (!is.numeric(nfolds) || length(nfolds) != 1 ||
      is.na(nfolds) || nfolds < 2 || nfolds != as.integer(nfolds))) {
    stop("nfolds should be an integer greater than or equal to 2")
  }

  if (!is.null(preProcess)) {
    if (!is.character(preProcess) ||
      anyNA(preProcess) ||
      any(!(preProcess %in% available_preprocess_methods))) {
      stop("preProcess includes at least one invalid preprocessing method")
    }
  }

  oob_methods <- c(
    "rf",
    "ranger",
    "treebag",
    "bagEarth",
    "bagEarthGCV",
    "bagFDA",
    "bagFDAGCV",
    "cforest",
    "parRF",
    "extraTrees",
    "Rborist"
  )

  if (resample == "oob" && !(method %in% oob_methods)) {
    stop(
      paste(
        "resample = 'oob' can only be used with random forest,",
        "bagged tree, bagged earth, bagged FDA, or conditional tree forest methods"
      )
    )
  }

  if (resample %in% c("LGOCV", "adaptive_LGOCV") &&
    (!is.numeric(p) || length(p) != 1 || is.na(p) || p <= 0 || p >= 1)) {
    stop("p should be a numeric value between 0 and 1 for LGOCV-based resampling methods")
  }


  if (resample %in% c("adaptive_cv", "adaptive_LGOCV") && nfolds <= 5) {
    stop("nfolds should be greater than 5 for adaptive CV/LGOCV resampling methods")
  }

  invisible(TRUE)
}


##############################################################################
#' @title Validate parameters specific to \code{linComb}
#'
#' @description
#' Internal helper used by \code{linComb} to validate linear-combination-specific
#' arguments. It checks the selected combination method, resampling method,
#' standardization method, and method-specific standardization requirements.
#' It also adjusts \code{nrepeats} when \code{resample = "cv"} and forces the
#' required standardization for selected methods.
#'
#' @param method A character string specifying the linear combination method.
#' Available methods are \code{"scoring"}, \code{"SL"}, \code{"logistic"},
#' \code{"minmax"}, \code{"PT"}, \code{"PCL"}, \code{"minimax"}, and
#' \code{"TS"}.
#'
#' @param resample A character string specifying the resampling method. Available
#' options are \code{"none"}, \code{"cv"}, \code{"repeatedcv"}, and
#' \code{"boot"}.
#'
#' @param nfolds A positive integer indicating the number of folds for
#' cross-validation-based resampling methods. Default is 5.
#'
#' @param nrepeats A positive integer indicating the number of repeats for
#' \code{"repeatedcv"}. If \code{resample = "cv"}, this is set to 1.
#' Default is 3.
#'
#' @param niters A positive integer indicating the number of bootstrap
#' resampling iterations when \code{resample = "boot"}. Default is 10.
#'
#' @param standardize A character string specifying the standardization method.
#' Available options are \code{"none"}, \code{"min_max_scale"}, \code{"zScore"},
#' \code{"tScore"}, \code{"scale_mean_to_one"}, and
#' \code{"scale_sd_to_one"}.
#'
#' @param ndigits A non-negative integer indicating the number of decimal places
#' used for rounding coefficients in the \code{"scoring"} method. Default is 0.
#'
#' @details
#' The \code{"minmax"} and \code{"PCL"} methods require
#' \code{standardize = "min_max_scale"}. If another standardization method is
#' supplied, it is replaced with \code{"min_max_scale"} and a warning is issued.
#'
#' The \code{"PT"} method requires \code{standardize = "zScore"}. If another
#' standardization method is supplied, it is replaced with \code{"zScore"} and a
#' warning is issued.
#'
#' @return
#' A list containing the validated and prepared values:
#' \code{method}, \code{resample},
#' \code{nfolds}, \code{nrepeats}, \code{niters},
#' \code{standardize}, and \code{ndigits}.
#'
#' @keywords internal


validateLinComb <- function(method = c(
                              "scoring",
                              "SL",
                              "logistic",
                              "minmax",
                              "PT",
                              "PCL",
                              "minimax",
                              "TS"
                            ), resample = c("none", "cv", "repeatedcv", "boot"),
                            nfolds = 5,
                            nrepeats = 3,
                            niters = 10,
                            standardize = c(
                              "none", "min_max_scale",
                              "zScore", "tScore", "scale_mean_to_one", "scale_sd_to_one"
                            ),
                            ndigits = 0) {
  methods <-
    c(
      "scoring",
      "SL",
      "logistic",
      "minmax",
      "PT",
      "PCL",
      "minimax",
      "TS"
    )

  if (length(which(methods == method)) == 0 || length(method) != 1) {
    stop(
      paste(
        "method should be one of 'scoring', 'SL', 'logistic', 'minmax',",
        "'PT', 'PCL', 'minimax', 'TS'"
      )
    )
  }

  resamples <- c("none", "cv", "repeatedcv", "boot")

  if (length(which(resamples == resample)) == 0) {
    stop(paste("resample should be one of 'none', 'cv', 'repeatedcv', 'boot'"))
  }

  if (any(resample == "none")) {
    resample <- "none"
  }

  if (resample %in% c("cv", "repeatedcv") &&
    (!is.numeric(nfolds) || length(nfolds) != 1 ||
      is.na(nfolds) || nfolds < 2 || nfolds != as.integer(nfolds))) {
    stop("nfolds should be an integer greater than or equal to 2")
  }

  if (resample == "repeatedcv" &&
    (!is.numeric(nrepeats) || length(nrepeats) != 1 ||
      is.na(nrepeats) || nrepeats < 2 || nrepeats != as.integer(nrepeats))) {
    stop("nrepeats should be an integer greater than or equal to 2 when resample = 'repeatedcv'")
  }

  if (resample == "cv") {
    nrepeats <- 1
  }

  if (resample == "boot" &&
    (!is.numeric(niters) || length(niters) != 1 ||
      is.na(niters) || niters < 1 || niters != as.integer(niters))) {
    stop("niters should be a positive integer when resample = 'boot'")
  }

  if (!is.numeric(ndigits) || length(ndigits) != 1 ||
    is.na(ndigits) || ndigits < 0 || ndigits != as.integer(ndigits)) {
    stop("ndigits should be a non-negative integer")
  }

  standardizes <-
    c("none", "min_max_scale", "zScore", "tScore", "scale_mean_to_one", "scale_sd_to_one")

  if (length(which(standardizes == standardize)) == 0) {
    stop(
      paste(
        "standardize should be one of 'min_max_scale', 'zScore', 'tScore',",
        "'scale_mean_to_one', 'scale_sd_to_one'"
      )
    )
  }

  if (length(standardize) != 1) {
    standardize <- "none"
  }

  if (method %in% c("minmax", "PCL") &&
    (!standardize == "min_max_scale")) {
    warning(
      paste(
        "The used combination method requires min_max_scale standardization.",
        "All biomarker values are standardized to a scale between 0 and 1."
      )
    )
    standardize <- "min_max_scale"
  }
  if (method %in% "PT" &&
    (!standardize == "zScore")) {
    warning(
      paste("The used combination method requires zScore standardization.")
    )
    standardize <- "zScore"
  }
  list(
    method = method,
    resample = resample,
    nfolds = nfolds,
    nrepeats = nrepeats,
    niters = niters,
    standardize = standardize,
    ndigits = ndigits
  )
}

##############################################################################
#' @title Validate parameters specific to \code{nonlinComb}
#'
#' @description
#' Internal helper used by \code{nonlinComb} to validate
#' nonlinear-combination-specific arguments. It checks the selected nonlinear
#' combination method, polynomial degrees, spline/GAM degrees of freedom,
#' resampling method and related resampling parameters, standardization method,
#' interaction option, and elastic net mixing parameter.
#'
#' @param method A character string specifying the nonlinear combination method.
#' Available methods are \code{"polyreg"}, \code{"ridgereg"},
#' \code{"lassoreg"}, \code{"elasticreg"}, \code{"splines"}, \code{"sgam"},
#' and \code{"nsgam"}.
#'
#' @param degree1 A positive integer specifying the polynomial degree for the
#' first marker in polynomial-based methods. Default is 3.
#'
#' @param degree2 A positive integer specifying the polynomial degree for the
#' second marker in polynomial-based methods. Default is 3.
#'
#' @param df1 A positive integer specifying the degrees of freedom for the first
#' marker in spline/GAM-based methods. Default is 4.
#'
#' @param df2 A positive integer specifying the degrees of freedom for the second
#' marker in spline/GAM-based methods. Default is 4.
#'
#' @param resample A character string specifying the resampling method. Available
#' options are \code{"none"}, \code{"cv"}, \code{"repeatedcv"}, and
#' \code{"boot"}.
#'
#' @param nfolds A positive integer indicating the number of folds for
#' cross-validation-based resampling methods. Default is 5.
#'
#' @param nrepeats A positive integer indicating the number of repeats for
#' \code{"repeatedcv"}. If \code{resample = "cv"}, this is set to 1.
#' Default is 3.
#'
#' @param niters A positive integer indicating the number of bootstrap
#' resampling iterations when \code{resample = "boot"}. Default is 10.
#'
#' @param standardize A character string specifying the standardization method.
#' Available options are \code{"none"}, \code{"min_max_scale"}, \code{"zScore"},
#' \code{"tScore"}, \code{"scale_mean_to_one"}, and
#' \code{"scale_sd_to_one"}.
#'
#' @param include.interact A logical value indicating whether the interaction
#' term between the two markers should be included in polynomial-based methods.
#' Default is \code{FALSE}.
#'
#' @param alpha A numeric value between 0 and 1 specifying the elastic net mixing
#' parameter when \code{method = "elasticreg"}. Default is 0.5.
#'
#' @return
#' A list containing the validated and prepared values:
#' \code{method}, \code{degree1}, \code{degree2}, \code{df1}, \code{df2},
#' \code{resample}, \code{nfolds}, \code{nrepeats}, \code{niters},
#' \code{standardize}, \code{include.interact}, and \code{alpha}.
#'
#' @keywords internal
#'
validateNonlinComb <- function(method = c(
                                 "polyreg",
                                 "ridgereg",
                                 "lassoreg",
                                 "elasticreg",
                                 "splines",
                                 "sgam",
                                 "nsgam"
                               ),
                               degree1 = 3,
                               degree2 = 3,
                               df1 = 4,
                               df2 = 4,
                               resample = c("none", "cv", "repeatedcv", "boot"),
                               nfolds = 5,
                               nrepeats = 3,
                               niters = 10,
                               standardize = c(
                                 "none", "min_max_scale", "zScore", "tScore",
                                 "scale_mean_to_one", "scale_sd_to_one"
                               ),
                               include.interact = FALSE,
                               alpha = 0.5) {
  methods <-
    c(
      "polyreg",
      "ridgereg",
      "lassoreg",
      "elasticreg",
      "splines",
      "sgam",
      "nsgam"
    )

  resamples <- c("none", "cv", "repeatedcv", "boot")

  standardizes <-
    c("none", "min_max_scale", "zScore", "tScore", "scale_mean_to_one", "scale_sd_to_one")

  if (length(which(methods == method)) == 0 || length(method) != 1) {
    stop(
      paste(
        "method should be one of 'polyreg', 'ridgereg', 'lassoreg',",
        "'elasticreg', 'splines', 'sgam', 'nsgam'"
      )
    )
  }

  if (length(which(resamples == resample)) == 0) {
    stop(paste("resample should be one of 'none', 'cv', 'repeatedcv', 'boot'"))
  }

  if (any(resample == "none")) {
    resample <- "none"
  }

  if (resample %in% c("cv", "repeatedcv") &&
    (!is.numeric(nfolds) || length(nfolds) != 1 ||
      is.na(nfolds) || nfolds < 2 || nfolds != as.integer(nfolds))) {
    stop("nfolds should be an integer greater than or equal to 2")
  }

  if (resample == "repeatedcv" &&
    (!is.numeric(nrepeats) || length(nrepeats) != 1 ||
      is.na(nrepeats) || nrepeats < 2 || nrepeats != as.integer(nrepeats))) {
    stop("nrepeats should be an integer greater than or equal to 2 when resample = 'repeatedcv'")
  }

  if (resample == "cv") {
    nrepeats <- 1
  }

  if (resample == "boot" &&
    (!is.numeric(niters) || length(niters) != 1 ||
      is.na(niters) || niters < 1 || niters != as.integer(niters))) {
    stop("niters should be a positive integer when resample = 'boot'")
  }

  if (length(which(standardizes == standardize)) == 0) {
    stop(
      paste(
        "standardize should be one of 'min_max_scale', 'zScore', 'tScore',",
        "'scale_mean_to_one', 'scale_sd_to_one'"
      )
    )
  }

  if (length(standardize) != 1) {
    standardize <- "none"
  }
  if (method %in% c("polyreg", "ridgereg", "lassoreg", "elasticreg") &&
    (!is.numeric(degree1) || length(degree1) != 1 ||
      is.na(degree1) || degree1 < 1 || degree1 != as.integer(degree1))) {
    stop("degree1 should be a positive integer for polynomial-based methods")
  }

  if (method %in% c("polyreg", "ridgereg", "lassoreg", "elasticreg") &&
    (!is.numeric(degree2) || length(degree2) != 1 ||
      is.na(degree2) || degree2 < 1 || degree2 != as.integer(degree2))) {
    stop("degree2 should be a positive integer for polynomial-based methods")
  }

  if (method %in% c("splines", "sgam", "nsgam") &&
    (!is.numeric(df1) || length(df1) != 1 ||
      is.na(df1) || df1 < 1 || df1 != as.integer(df1))) {
    stop("df1 should be a positive integer for spline/GAM-based methods")
  }

  if (method %in% c("splines", "sgam", "nsgam") &&
    (!is.numeric(df2) || length(df2) != 1 ||
      is.na(df2) || df2 < 1 || df2 != as.integer(df2))) {
    stop("df2 should be a positive integer for spline/GAM-based methods")
  }

  if (!is.logical(include.interact) || length(include.interact) != 1 ||
    is.na(include.interact)) {
    stop("include.interact should be TRUE or FALSE")
  }

  if (method == "elasticreg" &&
    (!is.numeric(alpha) || length(alpha) != 1 ||
      is.na(alpha) || alpha < 0 || alpha > 1)) {
    stop("alpha should be a numeric value between 0 and 1 when method = 'elasticreg'")
  }

  list(
    method = method,
    degree1 = degree1,
    degree2 = degree2,
    df1 = df1,
    df2 = df2,
    resample = resample,
    nfolds = nfolds,
    nrepeats = nrepeats,
    niters = niters,
    standardize = standardize,
    include.interact = include.interact,
    alpha = alpha
  )
}

##############################################################################
#' @title Validate parameters specific to \code{mathComb}
#'
#' @description
#' Internal helper used by \code{mathComb} to validate mathematical-combination
#' specific arguments. It checks the selected mathematical operator, distance
#' metric, standardization method, transformation method, and marker-value
#' requirements for selected transformations.
#'
#' @param markers A numeric data frame or matrix containing the marker values.
#' Marker values must be finite. If \code{transform = "log"}, all marker values
#' must be positive. If \code{transform = "exp"}, transformed values must remain
#' finite.
#'
#' @param method A character string specifying the mathematical combination
#' method. Available methods are \code{"add"}, \code{"multiply"},
#' \code{"divide"}, \code{"subtract"}, \code{"distance"},
#' \code{"baseinexp"}, and \code{"expinbase"}.
#'
#' @param distance A character string specifying the distance metric used when
#' \code{method = "distance"}. Available options are \code{"euclidean"},
#' \code{"manhattan"}, \code{"chebyshev"}, \code{"kulczynski_d"},
#' \code{"lorentzian"}, \code{"avg"}, \code{"taneja"}, and
#' \code{"kumar-johnson"}. For non-distance methods, this is set to
#' \code{NULL}.
#'
#' @param standardize A character string specifying the standardization method.
#' Available options are \code{"none"}, \code{"min_max_scale"}, \code{"zScore"},
#' \code{"tScore"}, \code{"scale_mean_to_one"}, and
#' \code{"scale_sd_to_one"}.
#'
#' @param transform A character string specifying the transformation applied to
#' marker values. Available options are \code{"none"}, \code{"log"},
#' \code{"exp"}, \code{"sin"}, and \code{"cos"}.
#'
#' @return
#' A list containing the validated and prepared values:
#' \code{standardize} and \code{transform}.
#'
#' @keywords internal
#'
validateMathComb <- function(
  markers, method = c(
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
  transform = c("none", "log", "exp", "sin", "cos")
) {
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

  if (length(standardize) != 1) {
    standardize <- "none"
  }
  if (length(transform) != 1) {
    transform <- "none"
  }
  if (any(!is.finite(as.matrix(markers)))) {
    stop("markers should include only finite numeric values")
  }

  if (transform == "log") {
    if (any(as.matrix(markers) <= 0)) {
      stop("markers should include only positive values when transform = 'log'")
    }
  }

  if (transform == "exp") {
    transformed <- exp(as.matrix(markers))

    if (any(!is.finite(transformed))) {
      stop("transform = 'exp' produced non-finite values; marker values may be too large")
    }
  }

  list(
    standardize = standardize,
    transform = transform
  )
}
