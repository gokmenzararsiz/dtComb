% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardization.R
\name{std.train}
\alias{std.train}
\title{Standardization according to the chosen method.}
\usage{
std.train(data, standardize = NULL)
}
\arguments{
\item{data}{a \code{numeric} data frame of biomarkers}

\item{standardize}{a \code{character} string indicating the name of the
standardization method. The default option is no standardization applied.
Available options are:
\itemize{
\item \bold{Z-score} \code{(zScore)}: This method scales the data to have a mean
of 0 and a standard deviation of 1. It subtracts the mean and divides by the standard
deviation for each feature. Mathematically,
\deqn{ Z-score = \frac{x - (\overline x)}{sd(x)}}

where \eqn{x} is the value of a marker, \eqn{\overline{x}} is the mean of the marker and \eqn{sd(x)} is the standard deviation of the marker.
\item \bold{T-score} \code{(tScore)}: T-score is commonly used
in data analysis to transform raw scores into a standardized form.
The standard formula for converting a raw score \eqn{x} into a T-score is:
\deqn{T-score = \Biggl(\frac{x - (\overline x)}{sd(x)}\times 10 \Biggl) +50}
where \eqn{x} is the value of a marker, \eqn{\overline{x}} is the mean of the marker
and \eqn{sd(x)} is the standard deviation of the marker.

\item \bold{min_max_scale} \code{(min_max_scale)}: This method transforms data to
a specific scale, between 0 and 1. The formula for this method is:
\deqn{min_max_scale = \frac{x - min(x)}{max(x) - min(x)}}

\item \bold{scale_mean_to_one} \code{(scale_mean_to_one)}: This method scales
the arithmetic mean to 1. The formula for this method is:
\deqn{scale_mean_to_one =  \frac{x}{\overline{x}}}
where \eqn{x} is the value of a marker and \eqn{\overline{x}} is the mean of the marker.

\item \bold{scale_sd_to_one} \code{(scale_sd_to_one)}: This method, which allows for
comparison of individual data points in relation to the overall spread of
the data, scales the standard deviation to 1. The formula for this method is:
\deqn{scale_sd_to_one = \frac{x}{sd(x)}}
where \eqn{x} is the value of a marker and \eqn{sd(x)} is the standard deviation of the marker.
}}
}
\value{
A \code{numeric} data.frame of standardized biomarkers
}
\description{
The \code{std.train} Standardization (min_max_scale, zScore etc.) can be
estimated from the training data and applied to any dataset with the same
variables.
}
\examples{
# call data
data(laparotomy)

# define the function parameters
markers <- laparotomy[, -1]
markers2 <- std.train(markers, "deviance")

}
\author{
Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
}
