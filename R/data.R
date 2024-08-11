#' Examples data for the dtComb package
#'
#' A data set containing the results of diagnostic laparoscopy procedures for 225
#' patients.
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
#' exampleData1$group <- factor(exampleData1$group)
#' gcol <- c("#E69F00", "#56B4E9")
#' plot(exampleData1$ddimer, exampleData1$log_leukocyte,
#'   col = gcol[as.numeric(exampleData1$group)]
#' )
#'
"exampleData1"

###############################################################################
#'
#' A data set containing the carriers of a rare genetic disorder for 120 samples.
#'
#' @docType data
#'
#' @usage data(exampleData2)
#'
#' @name exampleData2
#'
#' @format A data frame with 120 rows and 5 variables:
#' \describe{
#'   \item{Group}{Indicator if the person was carriers, values carriers and
#'   normals}
#'   \item{m1}{Biomarker 1, 1. measurement blood sample}
#'   \item{m2}{Biomarker 2, 2. measurement blood sample}
#'   \item{m3}{Biomarker 3, 3. measurement blood sample}
#'   \item{m4}{Biomarker 4, 4. measurement blood sample}
#' }
#'
#' @examples
#' data(exampleData2)
#' exampleData2$Group <- factor(exampleData2$Group)
#' gcol <- c("#E69F00", "#56B4E9")
#' plot(exampleData2$m1, exampleData2$m2,
#'   col = gcol[as.numeric(exampleData2$Group)]
#' )
#'
"exampleData2"

###############################################################################
#' A simulation data containing 250 diseased and 250 healthy individuals.
#' @docType data
#'
#' @usage data(exampleData3)
#'
#' @name exampleData3
#'
#' @format A data frame with 500 rows and 3 variables:
#' \describe{
#'   \item{status}{Indicator of one's condition, values healthy and diseased}
#'   \item{marker1}{1. biomarker}
#'   \item{marker2}{2. biomarker}
#' }
#'
#' @examples
#' data(exampleData3)
#' exampleData3$status <- factor(exampleData3$status)
#' gcol <- c("#E69F00", "#56B4E9")
#' plot(exampleData3$marker1, exampleData3$marker2,
#'   col = gcol[as.numeric(exampleData3$status)]
#' )
#'
"exampleData3"

###############################################################################
#' Includes machine learning models used for the mlComb function
#' @docType data
#'
#' @usage data(allMethods)
#'
#' @name allMethods
#'
#' @format A data frame with 113 rows and 2 variables:
#' \describe{
#'   \item{Method}{Valid name for the function}
#'   \item{Model}{Model name}
#' }
#'
#' @examples
#' data(allMethods)
#' allMethods
#'
"allMethods"
