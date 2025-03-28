% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{laparotomy}
\alias{laparotomy}
\title{Diagnostic laparotomy dataset
A data set containing the results of diagnostic laparotomy procedures for 225
patients.}
\format{
A data frame with 225 rows and 3 variables:
\describe{
\item{group}{Indicator whether the procedure was needed; values:
\code{"needed"} or \code{"not_needed"}.}
\item{ddimer}{Biomarker 1, D-Dimer protein level in blood, ng/mL}
\item{log_leukocyte}{Biomarker 2, Logarithm of Leukocyte count in blood,
per mcL}
}
}
\usage{
data(laparotomy)
}
\description{
Examples data for the dtComb package
}
\examples{
data(laparotomy)
laparotomy$group <- factor(laparotomy$group)
gcol <- c("#E69F00", "#56B4E9")
plot(laparotomy$ddimer, laparotomy$log_leukocyte,
  col = gcol[as.numeric(laparotomy$group)]
)

}
\keyword{datasets}
