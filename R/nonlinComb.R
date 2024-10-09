#' @title Combine two diagnostic tests with several non-linear combination methods.
#'
#' @description The \code{nonlinComb} function calculates the combination
#' scores of two diagnostic tests selected among several non-linear combination
#' methods and standardization options
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
#'  \item \bold{Logistic Regression with Polynomial Feature Space} \code{(polyreg)}:  The method
#'  builds a logistic regression model with the polynomial feature space and returns the probability
#'  of a positive event for each observation.
#'  \item \bold{Ridge Regression with Polynomial Feature Space} \code{(ridgereg)}: Ridge regression is a
#'  shrinkage method used to estimate the coefficients of highly correlated variables and in this case
#'   the polynomial feature space created from two markers. For the implementation of the method,
#'   glmnet() library is used with two functions: cv.glmnet() to run a cross
#'   validation model to determine the tuning parameter \eqn{\lambda} and glmnet() to fit the
#'    model with the selected tuning parameter. For Ridge regression,
#'    the glmnet() package is integrated into the dtComb package to facilitate the implementation
#'    of this method.
#'  \item \bold{Lasso Regression with Polynomial Feature Space} \code{(lassoreg)}: Lasso regression,
#'   like Ridge regression, is a type of shrinkage method. However, a notable difference is that
#'   Lasso tends to set some feature coefficients to zero, making it useful for feature elimination.
#'    It also employs cross-validation for parameter selection and model fitting using the glmnet library.
#'  \item \bold{Elastic Net Regression with Polynomial Feature Space} \code{(elasticreg)}: Elastic Net
#'   regression is a hybrid model that merges the penalties from Ridge and Lasso regression, aiming
#'   to leverage the strengths of both approaches. This model involves two parameters: \eqn{\lambda},
#'   similar to Ridge and Lasso, and \eqn{\alpha}, a user-defined mixing parameter ranging between 0 (representing Ridge)
#'    and 1 (representing Lasso). The \eqn{\alpha} parameter determines the balance or weights between the loss functions
#'     of Ridge and Lasso regressions.
#'  \item \bold{Splines} \code{(splines)}: Another non-linear approach to combine markers
#'   involves employing regression models within a polynomial feature space. This approach
#'    applies multiple regression models to the dataset using a function derived from
#'     piecewise polynomials. This implementation uses splines with user-defined degrees
#'      of freedom and degrees for the fitted polynomials. The splines library
#'       is employed to construct piecewise logistic regression models using base splines.
#'  \item \bold{Generalized Additive Models with Smoothing Splines and Generalized Additive Models
#'  with Natural Cubic Splines} \code{(sgam & nsgam)}: In addition to the basic spline structure,
#'   Generalized Additive Models are applied with natural cubic splines and smoothing splines
#'    using the gam library in R.
#'  }
#'
#' @param degree1 a \code{numeric} value for polynomial based methods indicates
#'  the degree of the feature space created for marker 1, for spline based
#'  methods the degree of the fitted polynomial between each node for marker 1.
#'  (3, default)
#'
#' @param degree2 a \code{numeric} value for polynomial based methods indicates
#'  the degree of the feature space created for marker 2, for spline based
#'  methods the degree of the fitted polynomial between each node for marker 2
#'  (3, default)
#'
#' @param df1 a \code{numeric} value that indicates the number of knots as the
#' degrees of freedom in spline based methods for marker 1 (4, default)
#'
#' @param df2 a \code{numeric} value that indicates the number of knots as the
#' degrees of freedom in spline based methods for marker 2 (4, default)
#'
#' @param resample a \code{character} string indicating the name of the
#' resampling options. Bootstrapping Cross-validation and repeated cross-validation
#' are given as the options for resampling, along with the number
#' of folds and number of repeats.
#' \itemize{
#'  \item \code{boot}: Bootstrapping is performed similarly; the dataset
#'   is divided into folds with replacement and models are trained and tested
#'   in these folds to determine the best parameters for the given method and
#'   dataset.
#'  \item \code{cv}: Cross-validation resampling, the dataset is divided into the
#'  number of folds given without replacement; in each iteration, one fold is
#'  selected as the test set, and the model is built using the remaining folds
#'  and tested on the test set. The corresponding AUC values and the parameters
#'  used for the combination are kept in a list. The best-performed model is
#'  selected, and the combination score is returned for the whole dataset.
#'  \item \code{repeatedcv}: Repeated cross-validation the process is repeated,
#'  and the best-performed models selected at each step are stored in another
#'  list; the best performed among these models is selected to be applied to
#'  the entire dataset.
#' }
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
#' @param standardize a \code{character} string indicating the name of the
#' standardization method. The default option is no standardization applied.
#' Available options are:
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
#' \item \bold{Range (a.k.a. min-max scaling)} \code{(range)}: This method transforms data to
#' a specific range, between 0 and 1. The formula for this method is:
#' \deqn{Range = \frac{x - min(x)}{max(x) - min(x)}}
#'
#' \item \bold{Mean} \code{(mean)}: This method, which helps
#' to understand the relative size of a single observation concerning
#' the mean of dataset, calculates the ratio of each data point to the mean value
#' of the dataset.
#' \deqn{Mean =  \frac{x}{\overline{x}}}
#' where \eqn{x} is the value of a marker and \eqn{\overline{x}} is the mean of the marker.
#'
#' \item \bold{Deviance} \code{(deviance)}: This method, which allows for
#' comparison of individual data points in relation to the overall spread of
#' the data, calculates the ratio of each data point to the standard deviation
#' of the dataset.
#' \deqn{Deviance = \frac{x}{sd(x)}}
#' where \eqn{x} is the value of a marker and \eqn{sd(x)} is the standard deviation of the marker.
#' }
#'
#' @param include.interact a \code{logical} indicator that specifies whether to
#' include the interaction between the markers to the feature space created for
#' polynomial based methods (FALSE, default)
#'
#' @param alpha a \code{numeric} value as the mixing parameter in Elastic Net
#' Regression method (0.5, default)
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
#' @param conf.level a \code{numeric} values determines the confidence interval
#' for the ROC curve(0.95, default).
#'
#' @param cutoff.method  a \code{character} string determines the cutoff method
#' for the ROC curve.
#'
#' @param show.result a \code{logical} string indicating whether the results
#' should be printed to the console.
#'
#' @param \dots further arguments. Currently has no effect on the results.
#'
#' @return A list of \code{numeric} nonlinear combination scores calculated
#' according to the given method and standardization option
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' data("exampleData1")
#' data <- exampleData1
#'
#' markers <- data[, -1]
#' status <- factor(data$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#' cutoff.method <- "Youden"
#'
#' score1 <- nonlinComb(
#'   markers = markers, status = status, event = event,
#'   method = "lassoreg", include.interact = FALSE, resample = "boot", niters = 5,
#'   degree1 = 4, degree2 = 4, cutoff.method = cutoff.method,
#'   direction = "<"
#' )
#'
#' score2 <- nonlinComb(
#'   markers = markers, status = status, event = event,
#'   method = "splines", resample = "none", cutoff.method = cutoff.method,
#'   standardize = "tScore", direction = "<"
#' )
#'
#' score3 <- nonlinComb(
#'   markers = markers, status = status, event = event,
#'   method = "lassoreg", resample = "repeatedcv", include.interact = TRUE,
#'   cutoff.method = "ROC01", standardize = "zScore", direction = "auto"
#' )
#'
#' @export

nonlinComb <- function(markers = NULL,
                       status = NULL,
                       event = NULL,
                       method = c(
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
                         "none", "range", "zScore", "tScore",
                         "mean", "deviance"
                       ),
                       include.interact = FALSE,
                       alpha = 0.5, show.plot = TRUE,
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
    c("none", "range", "zScore", "tScore", "mean", "deviance")

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

  if (any(resample == "cv")) {
    nrepeats <- 1
  }

  if (length(which(standardizes == standardize)) == 0) {
    stop(
      paste(
        "standardize should be one of 'range', 'zScore', 'tScore',",
        "'mean', 'deviance'"
      )
    )
  }

  if (length(standardize) != 1) {
    standardize <- "none"
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

  std.model <- std.train(markers, standardize)
  markersData <- std.model$data

  colnames(markersData) <- c("m1", "m2")
  data <- cbind(status, markersData)

  resample_results <- vector(mode = "list", length = 2)
  names(resample_results) <- c("parameters", "AUC")
  repeated_results <- vector(mode = "list", length = 2)
  names(repeated_results) <- c("parameters", "AUC")

  getFormula <- function(degree1, degree2, interact_ = FALSE) {
    formulaText <- paste0("status ~ m1 + m2")

    if (interact_) {
      formulaText <- paste0(formulaText, " + m1*m2")
    }

    for (j in 2:degree1) {
      formulaText <- paste0(formulaText, " + I(m1", "^", j, ")")
    }
    for (j in 2:degree2) {
      formulaText <- paste0(formulaText, " + I(m2", "^", j, ")")
    }
    return(formulaText)
  }
  getPower <- function(data, degree1, degree2, interact_ = FALSE) {
    space <- data.frame(m1.1 = data$m1)
    for (i in 2:degree1) {
      space[paste0("m1.", i)] <- data$m1^i
    }
    space["m2.1"] <- data$m2
    for (i in 2:degree2) {
      space[paste0("m2.", i)] <- data$m2^i
    }

    if (interact_) {
      space["interact"] <- data$m1 * data$m2
    }
    return(space)
  }
  suppressMessages(suppressWarnings({
    if (method == "polyreg") {
      if (any(resample == "boot")) {
        iters <- caret::createResample(status, niters)

        for (i in (1:niters)) {
          trainMark <- data[iters[[i]], ]
          testMark <- data

          model <-
            glm(
              getFormula(degree1, degree2, include.interact),
              data = trainMark,
              family = binomial(link = "logit")
            )
          comb.score <-
            predict(model, newdata = testMark, type = "response")

          auc_value <- as.numeric(pROC::auc(
            testMark$status,
            as.numeric(comb.score)
          ))

          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- resample_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters, newdata = data, type = "response")
      } else if (any(resample == "cv") ||
        any(resample == "repeatedcv")) {
        for (r in (1:nrepeats)) {
          folds <- caret::createFolds(status, nfolds)

          for (i in (1:nfolds)) {
            trainMark <- data[-folds[[i]], ]
            testMark <- data[folds[[i]], ]


            model <-
              glm(
                getFormula(degree1, degree2, include.interact),
                data = trainMark,
                family = binomial(link = "logit")
              )

            comb.score <-
              predict(model, newdata = testMark, type = "response")

            auc_value <- as.numeric(pROC::auc(
              testMark$status,
              as.numeric(comb.score)
            ))

            resample_results$parameters[[i]] <- model
            resample_results$AUC[[i]] <- auc_value
          }

          max_AUC <- which(resample_results$AUC ==
            max(unlist(resample_results$AUC)))
          if (length(max_AUC) != 1) {
            max_AUC <- max_AUC[[1]]
          }
          repeated_results$parameters[[r]] <-
            resample_results$parameters[[max_AUC]]
          repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        }

        max_AUC <- which(repeated_results$AUC ==
          max(unlist(repeated_results$AUC)))

        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- repeated_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters, newdata = data, type = "response")
      } else {
        model <- glm(
          getFormula(degree1, degree2, include.interact),
          data = data,
          family = binomial(link = "logit")
        )

        comb.score <-
          as.matrix(predict(model, newdata = data, type = "response"))

        parameters <- model
      }
    } else if (method == "ridgereg") {
      if (any(resample == "boot")) {
        iters <- caret::createResample(status, niters)

        for (i in (1:niters)) {
          trainMark <- data[iters[[i]], ]

          testMark <- data


          space <- getPower(trainMark, degree1, degree2, include.interact)
          cv.model <-
            glmnet::cv.glmnet(
              x = as.matrix(space),
              y = trainMark$status,
              alpha = 0,
              family = "binomial"
            )
          model <-
            glmnet::glmnet(
              x = space,
              y = trainMark$status,
              alpha = 0,
              family = "binomial",
              lambda = cv.model$lambda.min
            )

          testspace <- getPower(testMark, degree1, degree2, include.interact)
          dataspace <- getPower(data, degree1, degree2, include.interact)

          comb.score <-
            predict(model, newx = as.matrix(testspace), type = "response")

          auc_value <- as.numeric(pROC::auc(
            testMark$status,
            as.numeric(comb.score)
          ))

          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }

        parameters <- resample_results$parameters[[max_AUC]]
        comb.score <- predict(parameters,
          newx = as.matrix(dataspace),
          type = "response"
        )
      } else if (any(resample == "cv") ||
        any(resample == "repeatedcv")) {
        for (r in (1:nrepeats)) {
          folds <- caret::createFolds(data$status, nfolds)

          for (i in (1:nfolds)) {
            trainMark <- data[-folds[[i]], ]

            testMark <- data


            space <- getPower(trainMark, degree1, degree2, include.interact)

            cv.model <-
              glmnet::cv.glmnet(
                x = as.matrix(space),
                y = trainMark$status,
                alpha = 0,
                family = "binomial"
              )
            model <-
              glmnet::glmnet(
                x = space,
                y = trainMark$status,
                alpha = 0,
                family = "binomial",
                lambda = cv.model$lambda.min
              )


            testspace <- getPower(testMark, degree1, degree2, include.interact)
            dataspace <- getPower(data, degree1, degree2, include.interact)


            comb.score <-
              predict(model,
                newx = as.matrix(testspace),
                type = "response"
              )

            auc_value <- as.numeric(pROC::auc(
              testMark$status,
              as.numeric(comb.score)
            ))

            resample_results$parameters[[i]] <- model
            resample_results$AUC[[i]] <- auc_value
          }

          max_AUC <- which(resample_results$AUC ==
            max(unlist(resample_results$AUC)))
          if (length(max_AUC) != 1) {
            max_AUC <- max_AUC[[1]]
          }
          repeated_results$parameters[[r]] <-
            resample_results$parameters[[max_AUC]]
          repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        }

        max_AUC <- which(repeated_results$AUC ==
          max(unlist(repeated_results$AUC)))

        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- repeated_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters,
            newx = as.matrix(dataspace),
            type = "response"
          )
      } else {
        space <- getPower(data, degree1, degree2, include.interact)

        cv.model <-
          glmnet::cv.glmnet(
            x = as.matrix(space),
            y = status,
            alpha = 0,
            family = "binomial"
          )
        model <- glmnet::glmnet(
          x = space,
          y = status,
          alpha = 0,
          family = "binomial",
          lambda = cv.model$lambda.min
        )


        comb.score <-
          predict(model, newx = as.matrix(space), type = "response")

        parameters <- model
      }
    } else if (method == "lassoreg") {
      if (any(resample == "boot")) {
        iters <- caret::createResample(data$status, niters)

        for (i in (1:niters)) {
          trainMark <- data[iters[[i]], ]

          testMark <- data

          space <- getPower(trainMark, degree1, degree2, include.interact)

          cv.model <-
            glmnet::cv.glmnet(
              x = as.matrix(space),
              y = trainMark$status,
              alpha = 1,
              family = "binomial"
            )
          model <-
            glmnet::glmnet(
              x = space,
              y = trainMark$status,
              alpha = 1,
              family = "binomial",
              lambda = cv.model$lambda.min
            )

          testspace <- getPower(testMark, degree1, degree2, include.interact)
          dataspace <- getPower(data, degree1, degree2, include.interact)

          comb.score <-
            predict(model, newx = as.matrix(testspace), type = "response")

          auc_value <- as.numeric(pROC::auc(
            testMark$status,
            as.numeric(comb.score)
          ))

          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- resample_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters,
            newx = as.matrix(dataspace),
            type = "response"
          )
      } else if (any(resample == "cv") ||
        any(resample == "repeatedcv")) {
        for (r in (1:nrepeats)) {
          folds <- caret::createFolds(data$status, nfolds)

          for (i in (1:nfolds)) {
            trainMark <- data[-folds[[i]], ]

            testMark <- data[folds[[i]], ]


            space <- getPower(trainMark, degree1, degree2, include.interact)

            cv.model <-
              glmnet::cv.glmnet(
                x = as.matrix(space),
                y = trainMark$status,
                alpha = 1,
                family = "binomial"
              )
            model <-
              glmnet::glmnet(
                x = space,
                y = trainMark$status,
                alpha = 1,
                family = "binomial",
                lambda = cv.model$lambda.min
              )


            testspace <- getPower(testMark, degree1, degree2, include.interact)
            dataspace <- getPower(data, degree1, degree2, include.interact)

            comb.score <-
              predict(model,
                newx = as.matrix(testspace),
                type = "response"
              )

            auc_value <- as.numeric(pROC::auc(
              testMark$status,
              as.numeric(comb.score)
            ))

            resample_results$parameters[[i]] <- model
            resample_results$AUC[[i]] <- auc_value
          }

          max_AUC <- which(resample_results$AUC ==
            max(unlist(resample_results$AUC)))
          if (length(max_AUC) != 1) {
            max_AUC <- max_AUC[[1]]
          }
          repeated_results$parameters[[r]] <-
            resample_results$parameters[[max_AUC]]
          repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        }

        max_AUC <- which(repeated_results$AUC ==
          max(unlist(repeated_results$AUC)))

        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- repeated_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters,
            newx = as.matrix(dataspace),
            type = "response"
          )
      } else {
        space <- getPower(data, degree1, degree2, include.interact)

        cv.model <-
          glmnet::cv.glmnet(
            x = as.matrix(space),
            y = status,
            alpha = 1,
            family = "binomial"
          )
        model <- glmnet::glmnet(
          x = space,
          y = status,
          alpha = 1,
          family = "binomial",
          lambda = cv.model$lambda.min
        )


        comb.score <-
          predict(model, newx = as.matrix(space), type = "response")

        parameters <- model
      }
    } else if (method == "elasticreg") {
      if (any(resample == "boot")) {
        iters <- caret::createResample(data$status, niters)

        for (i in (1:niters)) {
          trainMark <- data[iters[[i]], ]
          testMark <- data

          space <- getPower(trainMark, degree1, degree2, include.interact)


          cv.model <-
            glmnet::cv.glmnet(
              x = as.matrix(space),
              y = trainMark$status,
              alpha = alpha,
              family = "binomial"
            )
          model <-
            glmnet::glmnet(
              x = space,
              y = trainMark$status,
              alpha = alpha,
              family = "binomial",
              lambda = cv.model$lambda.min
            )

          testspace <- getPower(testMark, degree1, degree2, include.interact)
          dataspace <- getPower(data, degree1, degree2, include.interact)

          comb.score <-
            predict(model, newx = as.matrix(testspace), type = "response")
          auc_value <-
            as.numeric(pROC::auc(testMark$status, as.numeric(comb.score)))

          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- resample_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters,
            newx = as.matrix(dataspace),
            type = "response"
          )
      } else if (any(resample == "cv") ||
        any(resample == "repeatedcv")) {
        for (r in (1:nrepeats)) {
          folds <- caret::createFolds(data$status, nfolds)

          for (i in (1:nfolds)) {
            trainMark <- data[-folds[[i]], ]
            testMark <- data[folds[[i]], ]

            space <- getPower(trainMark, degree1, degree2, include.interact)

            cv.model <-
              glmnet::cv.glmnet(
                x = as.matrix(space),
                y = trainMark$status,
                alpha = alpha,
                family = "binomial"
              )
            model <-
              glmnet::glmnet(
                x = space,
                y = trainMark$status,
                alpha = alpha,
                family = "binomial",
                lambda = cv.model$lambda.min
              )


            testspace <- getPower(testMark, degree1, degree2, include.interact)
            dataspace <- getPower(data, degree1, degree2, include.interact)

            comb.score <-
              predict(model,
                newx = as.matrix(testspace),
                type = "response"
              )
            auc_value <- as.numeric(pROC::auc(
              testMark$status,
              as.numeric(comb.score)
            ))

            resample_results$parameters[[i]] <- model
            resample_results$AUC[[i]] <- auc_value
          }

          max_AUC <- which(resample_results$AUC ==
            max(unlist(resample_results$AUC)))
          if (length(max_AUC) != 1) {
            max_AUC <- max_AUC[[1]]
          }
          repeated_results$parameters[[r]] <-
            resample_results$parameters[[max_AUC]]
          repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        }

        max_AUC <- which(repeated_results$AUC ==
          max(unlist(repeated_results$AUC)))

        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- repeated_results$parameters[[max_AUC]]
        comb.score <-
          as.matrix(predict(
            parameters,
            newx = as.matrix(dataspace),
            type = "response"
          ))
      } else {
        space <- getPower(data, degree1, degree2, include.interact)

        cv.model <-
          glmnet::cv.glmnet(
            x = as.matrix(space),
            y = status,
            alpha = alpha,
            family = "binomial"
          )
        model <-
          glmnet::glmnet(
            x = space,
            y = status,
            alpha = alpha,
            family = "binomial",
            lambda = cv.model$lambda.min
          )

        comb.score <-
          predict(model, newx = as.matrix(space), type = "response")

        parameters <- model
      }
    } else if (method == "splines") {
      if (any(resample == "boot")) {
        iters <- caret::createResample(status, niters)

        for (i in (1:niters)) {
          trainMark <- data[iters[[i]], ]
          testMark <- data

          model <-
            glm(
              status ~ splines::bs(m1, degree = degree1, df = df1) +
                splines::bs(m2, degree = degree2, df = df2),
              data = trainMark,
              family = binomial
            )

          comb.score <-
            predict(model, newdata = testMark, type = "response")
          auc_value <- as.numeric(pROC::auc(
            testMark$status,
            as.numeric(comb.score)
          ))

          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- resample_results$parameters[[max_AUC]]
        comb.score <- predict(parameters,
          newdata = data,
          type = "response"
        )
      } else if (any(resample == "cv") ||
        any(resample == "repeatedcv")) {
        for (r in (1:nrepeats)) {
          folds <- caret::createFolds(status, nfolds)

          for (i in (1:nfolds)) {
            trainMark <- data[-folds[[i]], ]
            testMark <- data[folds[[i]], ]

            model <-
              glm(
                status ~ splines::bs(m1, degree = degree1, df = df1) +
                  splines::bs(m2, degree = degree2, df = df2),
                data = trainMark,
                family = binomial
              )

            comb.score <-
              predict(model, newdata = testMark, type = "response")

            auc_value <- as.numeric(pROC::auc(
              testMark$status,
              as.numeric(comb.score)
            ))

            resample_results$parameters[[i]] <- model
            resample_results$AUC[[i]] <- auc_value
          }

          max_AUC <- which(resample_results$AUC ==
            max(unlist(resample_results$AUC)))
          if (length(max_AUC) != 1) {
            max_AUC <- max_AUC[[1]]
          }
          repeated_results$parameters[[r]] <-
            resample_results$parameters[[max_AUC]]
          repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        }

        max_AUC <- which(repeated_results$AUC ==
          max(unlist(repeated_results$AUC)))

        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- repeated_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters, newdata = data, type = "response")
      } else {
        model <- glm(
          status ~ splines::bs(m1, degree = degree1, df = df1) +
            splines::bs(m2, degree = degree2, df = df2),
          data = data,
          family = binomial
        )
        comb.score <-
          predict(model, newdata = data, type = "response")

        parameters <- model
      }
    } else if (method == "sgam") {
      if (any(resample == "boot")) {
        iters <- caret::createResample(status, niters)

        for (i in (1:niters)) {
          trainMark <- data[iters[[i]], ]
          testMark <- data

          model <- gam::gam(
            status ~ gam::s(m1, df = df1) +
              gam::s(m2, df = df2),
            data = trainMark,
            family = binomial
          )

          comb.score <-
            predict(model, newdata = testMark, type = "response")

          auc_value <- as.numeric(pROC::auc(
            testMark$status,
            as.numeric(comb.score)
          ))

          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- resample_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters, newdata = data, type = "response")
      } else if (any(resample == "cv") ||
        any(resample == "repeatedcv")) {
        for (r in (1:nrepeats)) {
          folds <- caret::createFolds(status, nfolds)

          for (i in (1:nfolds)) {
            trainMark <- data[-folds[[i]], ]
            testMark <- data[folds[[i]], ]

            model <- gam::gam(
              status ~ gam::s(m1, df = df1) +
                gam::s(m2, df = df2),
              data = trainMark,
              family = binomial
            )

            comb.score <-
              predict(model, newdata = testMark, type = "response")

            auc_value <- as.numeric(pROC::auc(
              testMark$status,
              as.numeric(comb.score)
            ))

            resample_results$parameters[[i]] <- model
            resample_results$AUC[[i]] <- auc_value
          }

          max_AUC <- which(resample_results$AUC ==
            max(unlist(resample_results$AUC)))
          if (length(max_AUC) != 1) {
            max_AUC <- max_AUC[[1]]
          }
          repeated_results$parameters[[r]] <-
            resample_results$parameters[[max_AUC]]
          repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        }

        max_AUC <- which(repeated_results$AUC ==
          max(unlist(repeated_results$AUC)))

        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- repeated_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters, newdata = data, type = "response")
      } else {
        model <- gam::gam(
          status ~ gam::s(m1, df = df1) +
            gam::s(m2, df = df2),
          data = data,
          family = binomial
        )

        comb.score <-
          predict(model, newdata = data, type = "response")

        parameters <- model
      }
    } else if (method == "nsgam") {
      if (any(resample == "boot")) {
        iters <- caret::createResample(status, niters)

        for (i in (1:niters)) {
          trainMark <- data[iters[[i]], ]
          testMark <- data

          model <- gam::gam(
            status ~ splines::ns(m1, df = df1) +
              splines::ns(m2, df = df2),
            data = trainMark,
            family = binomial
          )

          comb.score <-
            predict(model, newdata = testMark, type = "response")

          auc_value <- as.numeric(pROC::auc(
            testMark$status,
            as.numeric(comb.score)
          ))

          resample_results$parameters[[i]] <- model
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- resample_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters, newdata = data, type = "response")
      } else if (any(resample == "cv") ||
        any(resample == "repeatedcv")) {
        for (r in (1:nrepeats)) {
          folds <- caret::createFolds(status, nfolds)

          for (i in (1:nfolds)) {
            trainMark <- data[-folds[[i]], ]
            testMark <- data[folds[[i]], ]

            model <- gam::gam(
              status ~ splines::ns(m1, df = df1) +
                splines::ns(m2, df = df2),
              data = trainMark,
              family = binomial
            )

            comb.score <-
              predict(model, newdata = testMark, type = "response")

            auc_value <- as.numeric(pROC::auc(
              testMark$status,
              as.numeric(comb.score)
            ))

            resample_results$parameters[[i]] <- model
            resample_results$AUC[[i]] <- auc_value
          }

          max_AUC <- which(resample_results$AUC ==
            max(unlist(resample_results$AUC)))
          if (length(max_AUC) != 1) {
            max_AUC <- max_AUC[[1]]
          }
          repeated_results$parameters[[r]] <-
            resample_results$parameters[[max_AUC]]
          repeated_results$AUC[[r]] <- resample_results$AUC[[max_AUC]]
        }

        max_AUC <- which(repeated_results$AUC ==
          max(unlist(repeated_results$AUC)))
        if (length(max_AUC) != 1) {
          max_AUC <- max_AUC[[1]]
        }
        parameters <- repeated_results$parameters[[max_AUC]]
        comb.score <-
          predict(parameters, newdata = data, type = "response")
      } else {
        model <- gam::gam(
          status ~ splines::ns(m1, df = df1) +
            splines::ns(m2, df = df2),
          data = data,
          family = binomial
        )

        comb.score <-
          predict(model, newdata = data, type = "response")

        parameters <- model
      }
    }
  }))


  comb.score <- as.matrix(comb.score)
  status <- data$status

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
    CombType = "nonlinComb",
    Method = method,
    Parameters = parameters,
    Degree1 = degree1,
    Degree2 = degree2,
    Classification = status_levels,
    Interact = include.interact,
    Std.model = std.model$std,
    Standardize = standardize
  )

  allres$fit <- model_fit

  if (show.result) {
    print_model <- list(
      CombType = "nonlinComb",
      Method = method,
      rowcount = nrow(markers),
      colcount = ncol(markers),
      classification = status_levels,
      Pre_processing = standardize,
      Resampling = resample,
      niters = niters,
      nfolds = nfolds,
      nrepeats = nrepeats,
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
