#' @title Combine two diagnostic tests with several linear combination methods.
#'
#' @description The \code{linComb} function calculates the combination
#' scores of two diagnostic tests selected among several linear combination
#' methods and standardization options.
#'
#' @param markers a \code{numeric} a numeric data frame that includes two diagnostic tests
#' results
#'
#' @param status a \code{factor} vector that includes the actual disease
#' status of the patients
#'
#' @param event a \code{character} string that indicates the event in the status
#' to be considered as positive event
#'
#' @param method a \code{character} string specifying the method used for
#' combining the markers. \cr
#' \strong{Notations:}
#' Before getting into these methods,
#'  let us first introduce some notations that will be used
#'  throughout this vignette. Let
#' \eqn{D_i,  i = 1, 2, \ldots, n_1}
#' be the marker values of \eqn{ i\text{th}} individual in diseased group, where
#' \eqn{D_i = (D_{i1}, D_{i2})} and
#' \eqn{H_j, j=1,2, \ldots, n_2}
#' be the marker values of \eqn{ j\text{th}} individual in healthy group, where
#' \eqn{H_j = H_{j1}, H_{j2}}.
#' Let
#' \eqn{x_i1 = c(D_{i1}, H_{j1})} be the values of the first marker, and
#' \eqn{x_i2 = c(D_{i2}, H_{j2})} be values of the second marker for the \eqn{ i\text{th}}
#' individual \eqn{ i= 1,2, \ldots, n}. Let
#' \eqn{D_{i,min} = min(D_{i1}, D_{i2}), D_{i,max} = max(D_{i1}, D_{i2}) ,
#'  H_{j,min} = min(H_{j1}, H_{j2}), H_{j,max} = max(H_{j1}, H_{j2}) } and
#'  \eqn{c_i} be  be the resulting combination score for the \eqn{ i\text{th}} individual.
#'
#'
#' The available methods are:
#'
#'
#' \itemize{
#' \item \bold{Logistic Regression} \code{(logistic)}: Combination score obtained
#'  by fitting a logistic regression modelis as follows:
#'  \deqn{ c_i = \left(\frac{e^ {\beta_0 + \beta_1 x_{i1} + \beta_2 x_{i2}}}{1 + e^{\beta_0 + \beta_1 x_{i1} + \beta_2 x_{i2}}}\right)}
#'  A combination score obtained by fitting a logistic regression model typically refers
#'   to the predicted probability or score assigned to each observation
#'   in a dataset based on the logistic regression model’s
#'   fitted values
#'
#' \item \bold{Scoring based on Logistic Regression} \code{(scoring)}: Combination score is obtained using the
#' slope values of the relevant logistic regression model, slope values are rounded to the number of
#' digits taken from the user.
#' \deqn{c_i = \beta_1 x_{i1} + \beta_2 x_{i2}}
#'
#' \item \bold{Pepe & Thompson’s method} \code{(PT)}: The Pepe and Thompson combination score,
#'  developed using their optimal linear combination technique, aims to maximize
#'  the Mann-Whitney statistic in the same way that the Min-max method does. Unlike
#'  the Min-max method, the Pepe and Thomson method takes into account all marker
#'  values instead of just the lowest and maximum values.
#'  \deqn{  maximize\; U(\alpha) = \left(\frac{1}{n_1,n_2}\right) {\sum_{i=1}^{n_1} {\sum_{j=1}^{n_2}}I(D_{i1} + \alpha D_{i2} >= H_{j1} + \alpha H_{j2})}}
#'  \cr
#'  \deqn{c_i = x_{i1} + \alpha x_{i2}}
#'
#'
#' \item \bold{Pepe, Cai & Langton’s method} \code{(PCL)}: Pepe, Cai and Langton combination score
#'  obtained by using AUC as the parameter of a logistic regression model.
#'  \deqn{maximize\; U(\alpha) = \left(\frac{1}{n_1,n_2}\right) {\sum_{i=1}^{n_1} {\sum_{j=1}^{n_2}}I(D_{i1} + \alpha D_{i2} >}}
#'  \deqn{H_{j1} + \alpha H_{j2}) + \left(\frac{1}{2} \right) I(D_{i1} + \alpha D_{i2} = H_{j1} + \alpha H_{j2})}
#'
#'
#'
#'
#' \item \bold{Min-Max method} \code{(minmax)}: This method linearly combines the minimum
#'  and maximum values of the markers by finding a parameter,\eqn{\alpha}  , that
#'  maximizes the Mann-Whitney statistic, an empirical estimate of the ROC area.
#'  \deqn{ maximize\;U( \alpha ) = \left(\frac{1}{n_1,n_2}\right) {\sum_{i=1}^{n_1} {\sum_{j=1}^{n_2}}I(D_{i,max} + \alpha D_{i,min} > H_{j,max} + \alpha H_{j,min})}}
#'  \cr
#'  \deqn{c_i = x_{i,max} + \alpha x_{i,min}}
#'
#'  where \eqn{x_{i,max} = max(x_{i1},x_{i2})} and \eqn{ x_{i,min} = min(x_{i1}, x_{i2})}
#'
#' \item \bold{Su & Liu’s method} \code{(SL)}: The Su and Liu combination score is computed through
#'  Fisher’s discriminant coefficients, which assumes that the underlying
#'  data follow a multivariate normal distribution, and the covariance matrices across
#'  different classes are assumed to be proportional.Assuming that
#'  \eqn{D\sim N(\mu_D,\textstyle \sum_D)}
#'  and
#'  \eqn{H\sim N(\mu_H,\textstyle \sum_H)} represent
#'  the multivariate normal distributions for the diseased and non-diseased groups,
#'  respectively. The Fisher’s coefficients are as follows:
#'
#' \deqn{(\alpha , \beta) = (\textstyle \sum_{D}+\sum_{H})^{\;-1}\mu}
#' \eqn{ \text{where} \mu_=\mu_D - \mu_H. \text{The combination score in this case is:}}
#' \deqn{ c_i = \alpha x_{i1} + \beta x_{i2}}
#'
#'
#' \item \bold{Minimax approach} \code{(minimax)}: Combination score obtained with the Minimax procedure;
#' \eqn{t} parameter is chosen as the value that gives the maximum AUC from the
#'  combination score. Suppose that D follows a multivariate normal distribution
#'  \eqn{D\sim N(\mu_D,\textstyle \sum_D)}, representing diseased group and H follows
#'  a multivariate normal distribution \eqn{H\sim N(\mu_H,\textstyle \sum_H)} , representing the non-diseased group.
#'  Then Fisher’s coefficients are as follows:
#'  \deqn{ (\alpha , \beta) = {[t { \textstyle \sum_{D}} + (1 - t)  \textstyle \sum_{H}] ^ {-1}}{(\mu_D - \mu_H)}}
#'  \deqn{    c_i = b_1 x_1 + b_2 x_2}
#' \item \bold{Todor & Saplacan’s method} \code{(TS)}:Combination score obtained by using
#'  the trigonometric functions of the \eqn{\Theta} value that optimizes the corresponding AUC.
#'  \deqn{ c_i = sin(\theta) x_{i1} + cos(\theta) x_{i2}}
#' }
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
#' \itemize{
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
#' @param ndigits a \code{integer} value to indicate the number of decimal
#' places to be used for rounding in Scoring method (0, default)
#'
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
#' for the roc curve(0.95, default).
#'
#' @param cutoff.method  a \code{character} string determines the cutoff method
#' for the roc curve.
#'
#' @param show.result a \code{logical} string indicating whether the results
#' should be printed to the console.
#'
#' @param \dots further arguments. Currently has no effect on the results.
#'
#' @return A list of \code{numeric} linear combination scores calculated
#' according to the given method and standardization option.
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' # call data
#' data(exampleData1)
#'
#' # define the function parameters
#' markers <- exampleData1[, -1]
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#'
#' score1 <- linComb(
#'   markers = markers, status = status, event = event,
#'   method = "logistic", resample = "none", show.plot = TRUE,
#'   standardize = "none", direction = "<", cutoff.method = "Youden"
#' )
#'
#' # call data
#' data(exampleData2)
#'
#' # define the function parameters
#' markers <- exampleData2[, -c(1:3, 6:7)]
#' status <- factor(exampleData2$Group, levels = c("normals", "carriers"))
#' event <- "carriers"
#'
#' score2 <- linComb(
#'   markers = markers, status = status, event = event,
#'   method = "PT", resample = "none", standardize = "none", direction = "<",
#'   cutoff.method = "Youden", show.result = "TRUE"
#' )
#'
#' score3 <- linComb(
#'   markers = markers, status = status, event = event,
#'   method = "minmax", resample = "none", direction = "<",
#'   cutoff.method = "Youden"
#' )
#'
#' @export

linComb <- function(markers = NULL,
                    status = NULL,
                    event = NULL,
                    method = c(
                      "scoring",
                      "SL",
                      "logistic",
                      "minmax",
                      "PT",
                      "PCL",
                      "minimax",
                      "TS"
                    ),
                    resample = c("none", "cv", "repeatedcv", "boot"),
                    nfolds = 5,
                    nrepeats = 3,
                    niters = 10,
                    standardize = c(
                      "none", "range",
                      "zScore", "tScore", "mean", "deviance"
                    ),
                    ndigits = 0, show.plot = TRUE,
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
      "scoring",
      "SL",
      "logistic",
      "minmax",
      "PT",
      "PCL",
      "minimax",
      "TS"
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
        "method should be one of 'scoring', 'SL', 'logistic', 'minmax',",
        "'PT', 'PCL', 'minimax', 'TS'"
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

  if (method %in% c("minmax", "PCL") &&
    (!standardize == "range")) {
    warning(
      paste(
        "The used combination method requires range standardization.",
        "All biomarker values are standardized to a range between 0 and 1."
      )
    )
    standardize <- "range"
  }
  if (method %in% "PT" &&
    (!standardize == "zScore")) {
    warning(
      paste("The used combination method requires zScore standardization.")
    )
    standardize <- "zScore"
  }

  std.model <- std.train(markers, standardize)
  markers <- std.model$data

  neg.markers <- markers[status != 1, ]
  pos.markers <- markers[status == 1, ]

  resample_results <- vector(mode = "list", length = 3)
  names(resample_results) <- c("parameters", "AUC")
  repeated_results <- vector(mode = "list", length = 3)
  names(repeated_results) <- c("parameters", "AUC")

  suppressMessages(if (method == "scoring") {
    if (any(resample == "boot")) {
      iters <- caret::createResample(status, niters)

      for (i in (1:niters)) {
        trainMark <- markers[iters[[i]], ]
        testMark <- markers

        trainStat <- status[iters[[i]]]
        testStat <- status

        res <- glm(trainStat ~ trainMark[, 1] + trainMark[, 2],
          family = binomial((link <- "logit"))
        )

        round.coef <- round(res$coefficients, digits = ndigits)
        comb.score <-
          as.matrix(testMark) %*% as.matrix(round.coef[-1])

        auc_value <-
          as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

        resample_results$parameters[[i]] <- round.coef
        resample_results$AUC[[i]] <- auc_value
      }

      max_AUC <-
        which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% as.matrix(parameters[-1])
    } else if (any(resample == "cv") ||
      any(resample == "repeatedcv")) {
      for (r in (1:nrepeats)) {
        folds <- caret::createFolds(status, nfolds)

        for (i in (1:nfolds)) {
          trainMark <- markers[-folds[[i]], ]
          testMark <- markers[folds[[i]], ]

          trainStat <- status[-folds[[i]]]
          testStat <- status[folds[[i]]]

          res <- glm(trainStat ~ trainMark[, 1] + trainMark[, 2],
            family = binomial((link <- "logit"))
          )

          round.coef <- round(res$coefficients, digits = ndigits)
          comb.score <-
            as.matrix(testMark) %*% as.matrix(round.coef[-1])

          auc_value <-
            as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

          resample_results$parameters[[i]] <- round.coef
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <-
          resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <-
          resample_results$AUC[[max_AUC[1]]]
      }

      max_AUC <- which(repeated_results$AUC ==
        max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% as.matrix(parameters[-1])
    } else {
      res <- glm(status ~ markers[, 1] + markers[, 2],
        family = binomial((link <- "logit"))
      )

      round.coef <- round(res$coefficients, digits = ndigits)
      parameters <- round.coef
      comb.score <- as.matrix(markers) %*% as.matrix(round.coef[-1])
    }
  } else if (method == "SL") {
    if (any(resample == "boot")) {
      iters <- caret::createResample(status, niters)

      for (i in (1:niters)) {
        trainMark <- markers[iters[[i]], ]
        testMark <- markers

        trainStat <- status[iters[[i]]]
        testStat <- status

        neg.markers <- trainMark[trainStat != 1, ]
        pos.markers <- trainMark[trainStat == 1, ]

        sum.var <- var(pos.markers) + var(neg.markers)
        subs_mean <- colMeans(pos.markers) - colMeans(neg.markers)
        est.coef <- as.numeric(abs(solve(sum.var) %*% subs_mean))
        comb.score <- as.matrix(testMark) %*% est.coef

        auc_value <-
          as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

        names(est.coef) <- c("alpha", "beta")
        resample_results$parameters[[i]] <- est.coef
        resample_results$AUC[[i]] <- auc_value
      }

      max_AUC <- which(resample_results$AUC ==
        max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% parameters
    } else if (any(resample == "cv") ||
      any(resample == "repeatedcv")) {
      for (r in (1:nrepeats)) {
        folds <- caret::createFolds(status, nfolds)

        for (i in (1:nfolds)) {
          trainMark <- markers[-folds[[i]], ]
          testMark <- markers[folds[[i]], ]

          trainStat <- status[-folds[[i]]]
          testStat <- status[folds[[i]]]

          neg.markers <- trainMark[trainStat != 1, ]
          pos.markers <- trainMark[trainStat == 1, ]

          sum.var <- var(pos.markers) + var(neg.markers)
          subs_mean <- colMeans(pos.markers) - colMeans(neg.markers)
          est.coef <- as.numeric(abs(solve(sum.var) %*% subs_mean))
          comb.score <- as.matrix(testMark) %*% est.coef

          auc_value <-
            as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

          names(est.coef) <- c("alpha", "beta")
          resample_results$parameters[[i]] <- est.coef
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <-
          resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <-
          resample_results$AUC[[max_AUC[1]]]
      }

      max_AUC <- which(repeated_results$AUC ==
        max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% parameters
    } else {
      sum.var <- var(pos.markers) + var(neg.markers)
      subs_mean <- colMeans(pos.markers) - colMeans(neg.markers)
      est.coef <- as.numeric(abs(solve(sum.var) %*% subs_mean))
      names(est.coef) <- c("alpha", "beta")
      parameters <- est.coef
      comb.score <- as.matrix(markers) %*% est.coef
    }
  } else if (method == "logistic") {
    if (any(resample == "boot")) {
      markersData <- markers
      colnames(markersData) <- c("m1", "m2")
      data <- cbind(status, markersData)

      iters <- caret::createResample(status, niters)

      for (i in (1:niters)) {
        trainMark <- data[iters[[i]], ]
        testMark <- data

        res <- glm(status ~ m1 + m2,
          family = binomial((link <- "logit")), data = trainMark
        )

        comb.score <- as.matrix(predict(res,
          newdata =
            testMark, type = "response"
        ))

        auc_value <- as.numeric(pROC::auc(
          testMark$status,
          as.numeric(comb.score)
        ))

        resample_results$parameters[[i]] <- res
        resample_results$AUC[[i]] <- auc_value
      }

      max_AUC <- which(resample_results$AUC ==
        max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(predict(parameters,
        newdata = data,
        type = "response"
      ))
    } else if (any(resample == "cv") ||
      any(resample == "repeatedcv")) {
      markersData <- markers
      colnames(markersData) <- c("m1", "m2")
      data <- cbind(status, markersData)

      for (r in (1:nrepeats)) {
        folds <- caret::createFolds(status, nfolds)

        for (i in (1:nfolds)) {
          trainMark <- data[-folds[[i]], ]
          testMark <- data[folds[[i]], ]

          res <- glm(status ~ m1 + m2,
            family = binomial((link <- "logit")),
            data = trainMark
          )

          comb.score <- as.matrix(predict(res,
            newdata = testMark,
            type = "response"
          ))

          auc_value <- as.numeric(pROC::auc(
            testMark$status,
            as.numeric(comb.score)
          ))

          resample_results$parameters[[i]] <- res
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <-
          resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <-
          resample_results$AUC[[max_AUC[1]]]
      }

      max_AUC <- which(repeated_results$AUC ==
        max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(predict(parameters,
        newdata = data,
        type = "response"
      ))
    } else {
      res <- glm(status ~ .,
        data = markers,
        family = binomial((link <- "logit"))
      )
      parameters <- res
      comb.score <- as.matrix(predict(res,
        newdata = markers,
        type = "response"
      ))
    }
  } else if (method == "minmax") {
    if (any(resample == "boot")) {
      iters <- caret::createResample(status, niters)

      for (i in (1:niters)) {
        trainMark <- markers[iters[[i]], ]
        testMark <- markers

        trainStat <- status[iters[[i]]]
        testStat <- status

        neg.markers <- trainMark[trainStat != 1, ]
        pos.markers <- trainMark[trainStat == 1, ]

        init.param <- runif(1, 0, 1)

        opt.func <- optim(
          par = init.param,
          fn = helper_minmax,
          neg.set = neg.markers,
          pos.set = pos.markers,
          method = "Brent",
          lower = 0,
          upper = 1
        )
        lambda <- as.numeric(opt.func$par)

        comb.score <- as.matrix(apply(testMark, 1, max)
        + lambda * apply(testMark, 1, min))

        auc_value <-
          as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

        names(lambda) <- "lambda"
        resample_results$parameters[[i]] <- lambda
        resample_results$AUC[[i]] <- auc_value
      }

      max_AUC <- which(resample_results$AUC ==
        max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(apply(markers, 1, max)
      + parameters * apply(markers, 1, min))
    } else if (any(resample == "cv") ||
      any(resample == "repeatedcv")) {
      for (r in (1:nrepeats)) {
        folds <- caret::createFolds(status, nfolds)

        for (i in (1:nfolds)) {
          trainMark <- markers[-folds[[i]], ]
          testMark <- markers[folds[[i]], ]


          trainStat <- status[-folds[[i]]]
          testStat <- status[folds[[i]]]

          neg.markers <- trainMark[trainStat != 1, ]
          pos.markers <- trainMark[trainStat == 1, ]

          init.param <- runif(1, 0, 1)

          opt.func <- optim(
            par = init.param,
            fn = helper_minmax,
            neg.set = neg.markers,
            pos.set = pos.markers,
            method = "Brent",
            lower = 0,
            upper = 1
          )
          lambda <- as.numeric(opt.func$par)

          comb.score <- as.matrix(apply(testMark, 1, max)
          + lambda * apply(testMark, 1, min))

          auc_value <-
            as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

          names(lambda) <- "lambda"
          resample_results$parameters[[i]] <- lambda
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <-
          resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <-
          resample_results$AUC[[max_AUC[1]]]
      }

      max_AUC <- which(repeated_results$AUC ==
        max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(apply(markers, 1, max)
      + parameters * apply(markers, 1, min))
    } else {
      init.param <- runif(1, 0, 1)

      opt.func <- optim(
        par = init.param,
        fn = helper_minmax,
        neg.set = neg.markers,
        pos.set = pos.markers,
        method = "Brent",
        lower = 0,
        upper = 1
      )
      lambda <- as.numeric(opt.func$par)

      names(lambda) <- "lambda"
      parameters <- lambda

      comb.score <- as.matrix(apply(markers, 1, max)
      + lambda * apply(markers, 1, min))
    }
  } else if (method == "PT") {
    if (any(resample == "boot")) {
      iters <- caret::createResample(status, niters)

      for (i in (1:niters)) {
        trainMark <- markers[iters[[i]], ]
        testMark <- markers

        trainStat <- status[iters[[i]]]
        testStat <- status

        neg.markers <- trainMark[trainStat != 1, ]
        pos.markers <- trainMark[trainStat == 1, ]

        init.param <- runif(1, -1, 1)

        opt.func <- optim(
          par = init.param,
          fn = helper_PT,
          neg.set = neg.markers,
          pos.set = pos.markers,
          method = "Brent",
          lower = -1,
          upper = 1
        )

        lambda <- as.numeric(opt.func$par)
        markers <- as.matrix(markers)

        comb.score <-
          as.matrix(testMark[, 1] + testMark[, 2] * lambda)

        auc_value <-
          as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

        names(lambda) <- "lambda"
        resample_results$parameters[[i]] <- lambda
        resample_results$AUC[[i]] <- auc_value
      }

      max_AUC <- which(resample_results$AUC ==
        max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <-
        as.matrix(markers[, 1] + markers[, 2] * parameters)
    } else if (any(resample == "cv") ||
      any(resample == "repeatedcv")) {
      for (r in (1:nrepeats)) {
        folds <- caret::createFolds(status, nfolds)

        for (i in (1:nfolds)) {
          trainMark <- markers[-folds[[i]], ]
          testMark <- markers[folds[[i]], ]

          trainStat <- status[-folds[[i]]]
          testStat <- status[folds[[i]]]

          neg.markers <- trainMark[trainStat != 1, ]
          pos.markers <- trainMark[trainStat == 1, ]

          init.param <- runif(1, -1, 1)

          opt.func <- optim(
            par = init.param,
            fn = helper_PT,
            neg.set = neg.markers,
            pos.set = pos.markers,
            method = "Brent",
            lower = -1,
            upper = 1
          )

          lambda <- as.numeric(opt.func$par)
          markers <- as.matrix(markers)

          comb.score <-
            as.matrix(testMark[, 1] + testMark[, 2] * lambda)

          auc_value <-
            as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

          names(lambda) <- "lambda"
          resample_results$parameters[[i]] <- lambda
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <-
          resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <-
          resample_results$AUC[[max_AUC[1]]]
      }

      max_AUC <- which(repeated_results$AUC ==
        max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <-
        as.matrix(markers[, 1] + markers[, 2] * parameters)
    } else {
      init.param <- runif(1, -1, 1)

      opt.func <- optim(
        par = init.param,
        fn = helper_PT,
        neg.set = neg.markers,
        pos.set = pos.markers,
        method = "Brent",
        lower = -1,
        upper = 1
      )

      lambda <- as.numeric(opt.func$par)
      markers <- as.matrix(markers)

      names(lambda) <- "lambda"
      parameters <- lambda

      comb.score <- as.matrix(markers[, 1] + markers[, 2] * lambda)
    }
  } else if (method == "PCL") {
    if (any(resample == "boot")) {
      iters <- caret::createResample(status, niters)

      for (i in (1:niters)) {
        trainMark <- markers[iters[[i]], ]
        testMark <- markers

        trainStat <- status[iters[[i]]]
        testStat <- status

        neg.markers <- trainMark[trainStat != 1, ]
        pos.markers <- trainMark[trainStat == 1, ]

        init.param <- runif(1, 0, 1)

        opt.func <- optim(
          par = init.param,
          fn = helper_PCL,
          neg.set = neg.markers,
          pos.set = pos.markers,
          method = "Brent",
          lower = 0,
          upper = 1
        )

        lambda <- as.numeric(opt.func$par)
        markers <- as.matrix(markers)

        comb.score <-
          as.matrix(testMark[, 1] + testMark[, 2] * lambda)

        auc_value <-
          as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

        names(lambda) <- "lambda"
        resample_results$parameters[[i]] <- lambda
        resample_results$AUC[[i]] <- auc_value
      }

      max_AUC <- which(resample_results$AUC ==
        max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <-
        as.matrix(markers[, 1] + markers[, 2] * parameters)
    } else if (any(resample == "cv") ||
      any(resample == "repeatedcv")) {
      for (r in (1:nrepeats)) {
        folds <- caret::createFolds(status, nfolds)

        for (i in (1:nfolds)) {
          trainMark <- markers[-folds[[i]], ]
          testMark <- markers[folds[[i]], ]

          trainStat <- status[-folds[[i]]]
          testStat <- status[folds[[i]]]

          neg.markers <- trainMark[trainStat != 1, ]
          pos.markers <- trainMark[trainStat == 1, ]

          init.param <- runif(1, 0, 1)

          opt.func <- optim(
            par = init.param,
            fn = helper_PCL,
            neg.set = neg.markers,
            pos.set = pos.markers,
            method = "Brent",
            lower = 0,
            upper = 1
          )

          lambda <- as.numeric(opt.func$par)
          markers <- as.matrix(markers)

          comb.score <-
            as.matrix(testMark[, 1] + testMark[, 2] * lambda)

          auc_value <-
            as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

          names(lambda) <- "lambda"
          resample_results$parameters[[i]] <- lambda
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <-
          resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <-
          resample_results$AUC[[max_AUC[1]]]
      }

      max_AUC <- which(repeated_results$AUC ==
        max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <-
        as.matrix(markers[, 1] + markers[, 2] * parameters)
    } else {
      init.param <- runif(1, 0, 1)

      opt.func <- optim(
        par = init.param,
        fn = helper_PCL,
        neg.set = neg.markers,
        pos.set = pos.markers,
        method = "Brent",
        lower = 0,
        upper = 1
      )

      lambda <- as.numeric(opt.func$par)
      markers <- as.matrix(markers)

      names(lambda) <- "lambda"
      parameters <- lambda

      comb.score <- as.matrix(markers[, 1] + markers[, 2] * lambda)
    }
  } else if (method == "minimax") {
    if (any(resample == "boot")) {
      iters <- caret::createResample(status, niters)

      for (i in (1:niters)) {
        trainMark <- markers[iters[[i]], ]
        testMark <- markers

        trainStat <- status[iters[[i]]]
        testStat <- status

        neg.markers <- trainMark[trainStat != 1, ]
        pos.markers <- trainMark[trainStat == 1, ]

        init.param <- runif(1, 0, 1)

        opt.func <- optim(
          par = init.param,
          fn = helper_minimax,
          neg.set = neg.markers,
          pos.set = pos.markers,
          markers = trainMark,
          status = trainStat,
          method = "Brent",
          lower = 0,
          upper = 1
        )
        t <- as.numeric(opt.func$par)

        b.coef <-
          as.numeric((solve(t * var(
            pos.markers
          )) + (1 - t) *
            var(neg.markers)) %*% (colMeans(pos.markers) -
            colMeans(neg.markers)))
        comb.score <- as.matrix(testMark) %*% b.coef

        auc_value <-
          as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

        names(b.coef) <- c("b1", "b2")
        resample_results$parameters[[i]] <- b.coef
        resample_results$AUC[[i]] <- auc_value
      }

      max_AUC <- which(resample_results$AUC ==
        max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% parameters
    } else if (any(resample == "cv") ||
      any(resample == "repeatedcv")) {
      for (r in (1:nrepeats)) {
        folds <- caret::createFolds(status, nfolds)

        for (i in (1:nfolds)) {
          trainMark <- markers[-folds[[i]], ]
          testMark <- markers[folds[[i]], ]

          trainStat <- status[-folds[[i]]]
          testStat <- status[folds[[i]]]

          neg.markers <- trainMark[trainStat != 1, ]
          pos.markers <- trainMark[trainStat == 1, ]

          init.param <- runif(1, 0, 1)

          opt.func <- optim(
            par = init.param,
            fn = helper_minimax,
            neg.set = neg.markers,
            pos.set = pos.markers,
            markers = trainMark,
            status = trainStat,
            method = "Brent",
            lower = 0,
            upper = 1
          )
          t <- as.numeric(opt.func$par)

          b.coef <-
            as.numeric((solve(t * var(
              pos.markers
            )) + (1 - t) *
              var(neg.markers)) %*% (colMeans(pos.markers) -
              colMeans(neg.markers)))
          comb.score <- as.matrix(testMark) %*% b.coef

          auc_value <-
            as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

          names(b.coef) <- c("b1", "b2")
          resample_results$parameters[[i]] <- b.coef
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <-
          resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <-
          resample_results$AUC[[max_AUC[1]]]
      }

      max_AUC <- which(repeated_results$AUC ==
        max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]
      comb.score <- as.matrix(markers) %*% parameters
    } else {
      init.param <- runif(1, 0, 1)

      opt.func <- optim(
        par = init.param,
        fn = helper_minimax,
        neg.set = neg.markers,
        pos.set = pos.markers,
        markers = markers,
        status = status,
        method = "Brent",
        lower = 0,
        upper = 1
      )
      t <- as.numeric(opt.func$par)

      b.coef <-
        (solve(t * var(pos.markers)) + (1 - t) * var(neg.markers)) %*%
        (colMeans(pos.markers) - colMeans(neg.markers))

      parameters <- b.coef

      comb.score <- as.matrix(markers) %*% b.coef
    }
  } else if (method == "TS") {
    if (any(resample == "boot")) {
      iters <- caret::createResample(status, niters)

      for (i in (1:niters)) {
        trainMark <- markers[iters[[i]], ]
        testMark <- markers

        trainStat <- status[iters[[i]]]
        testStat <- status

        init.param <- runif(1, -1.57079633, 1.57079633)

        opt.func <-
          optim(
            par = init.param,
            fn = helper_TS,
            markers = trainMark,
            status = trainStat,
            method = "Brent",
            lower = -1.57079633,
            upper = 1.57079633
          )
        theta <- as.numeric(opt.func$par)

        a1 <- sin(theta)
        a2 <- cos(theta)

        comb.score <-
          as.matrix(a1 * testMark[, 1] + a2 * testMark[, 2])


        auc_value <-
          as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

        trigFuncs <- as.numeric(c(a1, a2))
        names(trigFuncs) <- c("sin(theta)", "cos(theta)")
        resample_results$parameters[[i]] <- trigFuncs
        resample_results$AUC[[i]] <- auc_value
      }

      max_AUC <-
        which(resample_results$AUC == max(unlist(resample_results$AUC)))
      parameters <- resample_results$parameters[[max_AUC[1]]]
      comb.score <-
        as.matrix(parameters[1] * markers[, 1] + parameters[2] *
          markers[, 2])
    } else if (any(resample == "cv") ||
      any(resample == "repeatedcv")) {
      for (r in (1:nrepeats)) {
        folds <- caret::createFolds(status, nfolds)

        for (i in (1:nfolds)) {
          trainMark <- markers[-folds[[i]], ]
          testMark <- markers[folds[[i]], ]

          trainStat <- status[-folds[[i]]]
          testStat <- status[folds[[i]]]

          init.param <- runif(1, -1.57079633, 1.57079633)

          opt.func <- optim(
            par = init.param,
            fn = helper_TS,
            markers = trainMark,
            status = trainStat,
            method = "Brent",
            lower = -1.57079633,
            upper = 1.57079633
          )
          theta <- as.numeric(opt.func$par)

          a1 <- sin(theta)
          a2 <- cos(theta)

          comb.score <-
            as.matrix(a1 * testMark[, 1] + a2 * testMark[, 2])


          auc_value <-
            as.numeric(pROC::auc(testStat, as.numeric(comb.score)))

          trigFuncs <- as.numeric(c(a1, a2))
          names(trigFuncs) <- c("sin(theta)", "cos(theta)")
          resample_results$parameters[[i]] <- trigFuncs
          resample_results$AUC[[i]] <- auc_value
        }

        max_AUC <- which(resample_results$AUC ==
          max(unlist(resample_results$AUC)))
        repeated_results$parameters[[r]] <-
          resample_results$parameters[[max_AUC[1]]]
        repeated_results$AUC[[r]] <-
          resample_results$AUC[[max_AUC[1]]]
      }

      max_AUC <- which(repeated_results$AUC ==
        max(unlist(repeated_results$AUC)))
      parameters <- repeated_results$parameters[[max_AUC[1]]]

      comb.score <- as.matrix(parameters[1] * markers[, 1] +
        parameters[2] * markers[, 2])
    } else {
      init.param <- runif(1, -1.57079633, 1.57079633)

      opt.func <-
        optim(
          par = init.param,
          fn = helper_TS,
          markers = markers,
          status = status,
          method = "Brent",
          lower = -1.57079633,
          upper = 1.57079633
        )
      theta <- as.numeric(opt.func$par)

      a1 <- sin(theta)
      a2 <- cos(theta)

      trigFuncs <- as.numeric(c(a1, a2))
      names(trigFuncs) <- c("sin(theta)", "cos(theta)")
      parameters <- trigFuncs

      comb.score <- as.matrix(a1 * markers[, 1] + a2 * markers[, 2])
    }
  })

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
    CombType = "linComb",
    Method = method,
    Parameters = parameters,
    Std.model = std.model$std,
    Classification = status_levels,
    Standardize = standardize
  )

  allres$fit <- model_fit

  if (show.result) {
    print_model <- list(
      CombType = "linComb",
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

#' @title Helper function for minmax method.
#'
#' @description The \code{helper_minmax} function estimates optimized value of
#' given biomarkers for the minmax method.
#'
#' @param lambda a \code{numeric} parameter that will be estimated in minmax
#' method for the combination score
#'
#' @param neg.set a \code{numeric} data frame that contains the observations
#' with negative status
#'
#' @param pos.set a \code{numeric} data frame that contains the observations
#' with positive status
#'
#' @return A \code{numeric} value for the estimated optimized value
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' # call data
#' data(exampleData1)
#'
#' # define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' neg.set <- markers[status == levels(status)[1], ]
#' pos.set <- markers[status == levels(status)[2], ]
#'
#' lambda <- 0.5
#'
#' stat <- helper_minmax(lambda, neg.set = neg.set, pos.set = pos.set)
#'
#' @export

helper_minmax <- function(lambda, neg.set, pos.set) {
  Xmax <- as.matrix(apply(neg.set, 1, max))
  Xmin <- as.matrix(apply(neg.set, 1, min))

  Ymax <- as.matrix(apply(pos.set, 1, max))
  Ymin <- as.matrix(apply(pos.set, 1, min))

  W.lambda <- 0
  n <- dim(neg.set)[1]
  m <- dim(pos.set)[1]

  for (i in 1:n) {
    for (j in 1:m) {
      W.lambda <- W.lambda + as.numeric(Ymax[j, ] + lambda * Ymin[j, ] >
        Xmax[i, ] + lambda * Xmin[i, ])
    }
  }

  return(-W.lambda / (n * m))
}



#' @title Helper function for PCL method.
#'
#' @description The \code{helper_PCL} function estimates the optimized value of
#' given biomarkers for the PCL method.
#'
#' @param lambda a \code{numeric} parameter that will be estimated in minmax
#' method for the combination score
#'
#' @param neg.set a \code{numeric} data frame that contains the observation with
#' negative status
#'
#' @param pos.set a \code{numeric} data frame that contains the observation with
#' positive status
#'
#' @return A \code{numeric} value for the estimated optimized value
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' # call data
#' data(exampleData1)
#'
#' # define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' neg.set <- markers[status == levels(status)[1], ]
#' pos.set <- markers[status == levels(status)[2], ]
#'
#' lambda <- 0.5
#'
#' stat <- helper_PCL(lambda, neg.set = neg.set, pos.set = pos.set)
#'
#' @export

helper_PCL <- function(lambda, neg.set, pos.set) {
  YD1 <- as.matrix(pos.set[, 1])
  YD2 <- as.matrix(pos.set[, 2])

  YDN1 <- as.matrix(neg.set[, 1])
  YDN2 <- as.matrix(neg.set[, 2])

  W.lambda <- 0

  n <- dim(pos.set)[1]
  m <- dim(neg.set)[1]

  for (i in 1:n) {
    for (j in 1:m) {
      W.lambda <- W.lambda +
        as.numeric(YD1[i, ] + lambda * YD2[i, ] >
          YDN1[j, ] + lambda * YDN2[j, ])
      +

        as.numeric(YD1[i, ] + lambda * YD2[i, ] == YDN1[j, ] +
          lambda * YDN2[j, ]) / 2
    }
  }

  return(-W.lambda / (n * m))
}

#' @title Helper function for minimax method.
#'
#' @description The \code{helper_minimax} function calculates the combination
#' coefficient and optimized value of given biomarkers for the minimax method.
#'
#' @param t a \code{numeric} parameter that will be estimated in minimax
#' method for the combination score
#'
#' @param neg.set a \code{numeric} data frame that contains the observation with
#' negative status
#'
#' @param pos.set a \code{numeric} data frame that contains the observation with
#' positive status
#'
#' @param markers a \code{numeric} data frame that contains the biomarkers
#'
#' @param status a \code{factor} data frame that includes the actual disease
#' status of the patients
#'
#' @return A \code{numeric} Optimized value calculated with combination scores
#' using t
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' # call data
#' data(exampleData1)
#'
#' # define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' neg.set <- markers[status == levels(status)[1], ]
#' pos.set <- markers[status == levels(status)[2], ]
#'
#' t <- 0.5
#'
#' stat <- helper_minimax(t,
#'   neg.set = neg.set, pos.set = pos.set,
#'   markers = markers, status
#' )
#'
#' @importFrom pROC auc
#'
#' @export

helper_minimax <- function(t, neg.set, pos.set, markers, status) {
  b.coef <- (solve(t * var(pos.set)) + (1 - t) * var(neg.set)) %*%
    (colMeans(pos.set) - colMeans(neg.set))
  comb.score <- as.matrix(markers) %*% b.coef
  comb.score <- as.numeric(comb.score)

  value <-
    suppressMessages(as.numeric(pROC::auc(status, comb.score)))

  return(-(value))
}

#' @title Helper function for TS method.
#'
#' @description The \code{helper_TS} function calculates the combination
#' coefficient and optimized value of given biomarkers for the TS method.
#'
#' @param theta a \code{numeric} parameter that will be estimated in TS
#' method for the combination score
#'
#' @param markers a \code{numeric} data frame that contains the biomarkers
#'
#' @param status a \code{factor} data frame that includes the actual disease
#' status of the patients
#'
#' @return A \code{numeric} Optimized value calculated with combination scores
#' using theta
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' # call data
#' data(exampleData1)
#'
#' # define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' t <- 0.5
#'
#' stat <- helper_TS(theta = t, markers = markers, status = status)
#'
#' @importFrom pROC auc
#'
#' @export

helper_TS <- function(theta, markers, status) {
  a1 <- sin(theta)
  a2 <- cos(theta)
  z <- a1 * markers[, 1] + a2 * markers[, 2]

  roc_obj <- suppressMessages(pROC::roc(status, z))
  value <- as.numeric(pROC::auc(roc_obj))

  return(-(value))
}
#' @title Helper function for PT method.
#'
#' @description The \code{helper_PT} function estimates the optimized value of
#' given biomarkers for the PT method.
#'
#' @param lambda a \code{numeric} parameter that will be estimated in minmax
#' method for the combination score
#'
#' @param neg.set a \code{numeric} data frame that contains the observation with
#' negative status
#'
#' @param pos.set a \code{numeric} data frame that contains the observation with
#' positive status
#'
#' @return A \code{numeric} value for the estimated optimized value
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#' # call data
#' data(exampleData1)
#'
#' # define the function parameters
#' markers <- cbind(exampleData1$ddimer, exampleData1$log_leukocyte)
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#'
#' neg.set <- markers[status == levels(status)[1], ]
#' pos.set <- markers[status == levels(status)[2], ]
#'
#' lambda <- 0.5
#'
#' stat <- helper_PT(lambda, neg.set = neg.set, pos.set = pos.set)
#'
#' @export

helper_PT <- function(lambda, neg.set, pos.set) {
  YD1 <- as.matrix(pos.set[, 1])
  YD2 <- as.matrix(pos.set[, 2])

  YDN1 <- as.matrix(neg.set[, 1])
  YDN2 <- as.matrix(neg.set[, 2])

  W.lambda <- 0

  n <- dim(pos.set)[1]
  m <- dim(neg.set)[1]

  for (i in 1:n) {
    for (j in 1:m) {
      W.lambda <- W.lambda +
        as.numeric(YD1[i, ] + lambda * YD2[i, ] >=
          YDN1[j, ] + lambda * YDN2[j, ])
    }
  }

  return(-W.lambda / (n * m))
}
