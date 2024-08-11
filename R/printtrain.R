#' @title Print the summary of linComb, nonlinComb, mlComb and mathComb
#' functions.
#'
#' @description The \code{print_train} function prints the summary statistics of
#' the fitted model
#'
#' @param print_model a \code{list} of parameters taken from the fitted model
#' that includes the combination method, resampling method, pre-processing
#' method, selected optimum parameters and the results of fit.
#'
#' @return \code{No return value} writes a summary of the results to the console.
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz


print_train <- function(print_model) {
  if (print_model$CombType != "mlComb") {
    cat("Method :", print_model$Method, "\n")

    if (print_model$Method == "distance") {
      cat(paste("Distance", print_model$Distance, sep = " : "), "\n")
    }
    cat(paste("Samples", print_model$rowcount, sep = " : "), "\n")
    cat(paste("Markers", print_model$colcount, sep = " : "), "\n")
    cat("Events :", paste(print_model$classification, collapse = ", "), "\n")
    cat(paste("Standardization", print_model$Pre_processing, sep = " : "), "\n")
    cat(paste("Cut points :", print_model$Cutoff_method, collapse = ", "), "\n")

    if (print_model$CombType == "mathComb") {
      cat(paste("Transform", print_model$Transform, sep = " : "), "\n")
      if (print_model$Method == "add" || print_model$Method == "subtract") {
        cat(paste("MaxPower", print_model$MaxPower, sep = " : "), "\n")
      }
    }

    if (print_model$CombType == "linComb" || print_model$CombType == "nonlinComb") {
      if (print_model$Resampling == "boot") {
        cat("Resampling : boot (niters:", paste(print_model$niters, ")", sep = ""))
      } else if (print_model$Resampling == "cv") {
        cat("Resampling : cv (nfolds:", paste(print_model$nfolds, ")", sep = ""))
      } else if (print_model$Resampling == "repeatedcv") {
        cat(
          "Resampling : repeatedcv (nfolds:", print_model$nfolds, ",", "nrepeats:",
          paste(print_model$nrepeats, ")", sep = "")
        )
      }
    }

    kappa.accuracy <- kappa.accuracy(print_model$DiagStatCombined)
    cat("\n")
    cat(
      " Kappa ", "   ", " Accuracy ", "\n", kappa.accuracy$C.kappa, " ",
      kappa.accuracy$Accuracy
    )
  } else {
    print(print_model$Model)
  }

  cat("\n")
  cat("\n")
  cat("Area Under the Curves of markers and combination score : ", "\n")
  print(print_model$AUC_table)
  cat("------------------------------------------------------------", "\n")

  cat("Area Under the Curve comparison of markers and combination score : ", "\n")
  print(print_model$MultComp_table)
  cat("------------------------------------------------------------", "\n")
  cat("Confusion matrix : ", "\n")
  print(print_model$DiagStatCombined)
  cat("------------------------------------------------------------", "\n")

  cat("Cut-off Results :", "\n")
  cat("Optimal cut-off method :", print_model$Cutoff_method, "\n", sep = " ")
  cat("Optimal cut-off point  :", print_model$ThresholdCombined, "\n", sep = " ")
  cat("Optimal criterion      :", print_model$Criterion, "\n", sep = " ")
  cat("------------------------------------------------------------", "\n")
}


#' @title Calculate Cohen's kappa and accuracy.
#'
#' @description The \code{kappa.accuracy} calculates Cohen's kappa and accuracy.
#'
#' @param DiagStatCombined a \code{numeric} table of confusion matrix of the
#' calculated combination score.
#'
#' @return A \code{list} of Cohen's kappa and accuracy values
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz

kappa.accuracy <- function(DiagStatCombined) {
  Outcome.p <- as.numeric(DiagStatCombined$tab$`   Outcome +`)
  Outcome.n <- as.numeric(DiagStatCombined$tab$`   Outcome -`)

  xtab <- as.table(cbind(Outcome.p, Outcome.n))
  xtab <- xtab[-3, ]
  diagonal.counts <- diag(xtab)
  N <- sum(xtab)
  row.marginal.props <- rowSums(xtab) / N
  col.marginal.props <- colSums(xtab) / N

  Po <- sum(diagonal.counts) / N
  Pe <- sum(row.marginal.props * col.marginal.props)
  k <- (Po - Pe) / (1 - Pe)

  accuracy <- sum(diagonal.counts) / N
  res <- list(
    C.kappa = k,
    Accuracy = accuracy
  )
  return(res)
}
