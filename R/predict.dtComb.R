#' @title Predict combination scores and labels for new data sets using the
#' training model
#'
#' @description The \code{predict.dtComb} is a function that generates predictions
#' for a new dataset of biomarkers using the parameters from the fitted model.
#' The function takes arguments newdata and model. The function's output is the
#' combination scores and labels of object type.
#'
#' @param newdata a \code{numeric} new data set that includes biomarkers that
#'  have not been introduced to the model before.
#'
#' @param object a \code{list} object where the parameters from the training
#' model are saved.
#'
#' @param \dots further arguments. Currently has no effect on the results.
#'
#' @return A \code{data.frame} predicted combination scores (or probabilities)
#' and labels
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#'
#' # call data
#' data(laparoscopy)
#'
#' # define the function parameters
#' markers <- laparoscopy[, -1]
#' status <- factor(laparoscopy$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#'
#' score1 <- linComb(
#'   markers = markers, status = status, event = event,
#'   method = "logistic", resample = "none",
#'   standardize = "none", direction = "<", cutoff.method = "Youden"
#' )
#'
#' comb.score1 <- predict(score1, markers)
#'
#' score2 <- nonlinComb(
#'   markers = markers, status = status, event = "needed", include.interact = TRUE,
#'   method = "polyreg", resample = "repeatedcv", nfolds = 5,
#'   nrepeats = 10, cutoff.method = "Youden", direction = "auto"
#' )
#'
#' comb.score2 <- predict(score2, markers)
#'
#' score3 <- mathComb(
#'   markers = markers, status = status, event = event,
#'   method = "distance", distance = "euclidean", direction = "auto",
#'   standardize = "tScore", cutoff.method = "Youden"
#' )
#'
#' comb.score3 <- predict(score3, markers)
#'
#' @method  predict dtComb
#' @export
#' @importFrom stats 'predict'


predict.dtComb <- function(object, newdata = NULL, ...) {
  model <- object

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  combtype <- model$fit$CombType



  if (combtype != "mlComb") {
    colnames(newdata) <- c("m1", "m2")
    newdata <- std.test(newdata, model)
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

  if (combtype == "linComb") {
    method <- model$fit$Method
    parameters <- model$fit$Parameters

    if (method == "scoring") {
      comb.score <- as.matrix(newdata) %*% as.matrix(model$fit$Parameters[-1])
    } else if (method == "SL") {
      comb.score <- as.matrix(newdata) %*% model$fit$Parameters
    } else if (method == "logistic") {
      coefficientsColNames <- names(model$fit$Parameters$coefficients[2:3])
      colnames(newdata) <- coefficientsColNames

      comb.score <- predict(model$fit$Parameters,
        newdata = newdata,
        type = "response"
      )
    } else if (method == "minmax") {
      comb.score <- as.matrix(apply(newdata, 1, max)
      + model$fit$Parameters * apply(newdata, 1, min))
    } else if (method == "PT") {
      comb.score <- as.matrix(newdata[, 1] + model$fit$Parameters * newdata[, 2])
    } else if (method == "PCL") {
      comb.score <- as.matrix(newdata[, 1] + newdata[, 2] * model$fit$Parameters)
    } else if (method == "minimax") {
      comb.score <- as.matrix(newdata) %*% model$fit$Parameters
    } else {
      comb.score <- as.matrix(model$fit$Parameters[1] * newdata[, 1] + model$fit$Parameters[2] *
        newdata[, 2])
    }
  } else if (combtype == "nonlinComb") {
    method <- model$fit$Method
    parameters <- model$fit$Parameters
    interact <- model$fit$Interact

    if (method %in% c("ridgereg", "lassoreg", "elasticreg")) {
      dataspace <- getPower(newdata, model$fit$Degree1, model$fit$Degree2, interact)


      comb.score <- predict(model$fit$Parameters,
        newx = as.matrix(dataspace), type = "response"
      )
    } else {
      comb.score <- predict(model$fit$Parameters,
        newdata = newdata,
        type = "response"
      )
    }
  } else if (combtype == "mlComb") {
    model_fit <- model$fit$Model


    comb.score <- predict(model_fit, newdata = newdata, type = "prob")
  } else {
    method <- model$fit$Method
    distance <- model$fit$Distance
    transform <- model$fit$Transform
    power.transform <- model$fit$PowerTransform
    max_power <- model$fit$MaxPower

    if (any(transform == "none")) {
      newdata <- newdata
    } else if (any(transform == "log")) {
      newdata <- log(newdata)
    } else if (any(transform == "exp")) {
      newdata <- exp(newdata)
    } else if (any(transform == "sin")) {
      newdata <- sin(newdata)
    } else {
      newdata <- cos(newdata)
    }

    if (method == "add") {
      markers <- newdata^model$fit$MaxPower
      comb.score <- markers[, 1] + markers[, 2]
    } else if (method == "multiply") {
      comb.score <- newdata[, 1] * newdata[, 2]
    } else if (method == "divide") {
      comb.score <- newdata[, 1] / newdata[, 2]
    } else if (method == "subtract") {
      markers <- newdata^model$fit$MaxPower
      comb.score <- markers[, 1] - markers[, 2]
    } else if (method == "distance") {
      if (distance == "euclidean") {
        comb.score <- sqrt(newdata[, 1]^2 + newdata[, 2]^2)
      } else if (distance == "manhattan") {
        comb.score <- newdata[, 1] + newdata[, 2]
      } else if (distance == "chebyshev") {
        comb.score <- apply(newdata, 1, max)
      } else if (distance == "kulczynski_d") {
        comb.score <- (abs(newdata[, 1]) + abs(newdata[, 2])) / 0.00002
      } else if (distance == "lorentzian") {
        comb.score <- log(abs(newdata[, 1]) + 1) + log(abs(newdata[, 2]) + 1)
      } else if (distance == "taneja") {
        epsilon <- 0.00001
        z1 <- (newdata[, 1] + 0.00001) / 2
        z2 <- (newdata[, 2] + 0.00001) / 2
        comb.score <- (z1 / 2) * log(z1 * sqrt(newdata[, 1] * epsilon)) +
          (z2 / 2) * log(z2 * sqrt(newdata[, 2] * epsilon))
      } else if (distance == "kumar-johnson") {
        epsilon <- 0.00001
        z1 <- (((newdata[, 1]^2) - (epsilon^2))^2) /
          2 * ((newdata[, 1] * epsilon)^1.5)
        z2 <- (((newdata[, 2]^2) - (epsilon^2))^2) /
          2 * ((newdata[, 2] * epsilon)^1.5)
        comb.score <- z1 + z2
      } else {
        comb.score <- (newdata[, 1] + newdata[, 2] + apply(newdata, 1, max)) / 2
      }
    } else if (method == "baseinexp") {
      comb.score <- newdata[, 1]^newdata[, 2]
    } else if (method == "expinbase") {
      comb.score <- newdata[, 2]^newdata[, 1]
    }
  }

  if (combtype != "mlComb") {
    comb.score <- as.matrix(comb.score)
    labels <- (comb.score > model$ThresholdCombined)
    labels[labels == TRUE] <- model$fit$Classification[2]
    labels[labels == "FALSE"] <- model$fit$Classification[1]
    comb.score <- data.frame(comb.score, labels)
    colnames(comb.score) <- c("comb.score", "labels")
  }

  return(comb.score)
}
