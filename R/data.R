# TODO: Add comment
#
# Author: serra
###############################################################################

#' Example data for the dtComb package
#'
#' A dataset containing the results of diagnostic laparoscopy procedure for 225
#' patients
#'
#' @docType data
#'
#' @usage data(exampleData1)
#'
#' @name exampleData1
#'
#' @format A data frame with 225 rows and 3 variables:
#' \describe{
#'   \item{group}{Indicator if the procedure was needed, values needed and
#'   not_needed}
#'   \item{ddimer}{Biomarker 1, D-Dimer protein level in blood, ng/mL}
#'   \item{log_leukocyte}{Biomarker 2, Logarithm of Leukocyte count in blood,
#'   per mcL}
#' }
#'
#' @examples
#' data(exampleData1)
#' exampleData1$group<-factor(exampleData1$group)
#' gcol <- c("#E69F00", "#56B4E9")
#' plot(exampleData1$ddimer, exampleData1$log_leukocyte,
#'                   col = gcol[as.numeric(exampleData1$group)])
#'
#'
"exampleData1"

