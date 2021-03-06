% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/procD_comparisons.R
\name{procD_comparisons}
\alias{procD_comparisons}
\title{Procrustes p-Value calculations}
\usage{
procD_comparisons(lm_matrix, labels, plot_results = TRUE, robust = TRUE)
}
\arguments{
\item{lm_matrix}{A matrix containing the data}

\item{labels}{A factorial list containing the sample labels}

\item{plot_results}{A boolean value defining whether results are to be plotted.
A maximum of 4 groups can be included for results to be plotted.}

\item{robust}{A boolean value defining whether p-values are to be computed
using robust statistical metrics.}
}
\value{
p-Values computed from procrustes distance matrices
}
\description{
A function used to calculate inter and intra-group Procrustes distance p-Values.
}
\examples{
library(geomorph)
library(shapes)
library(GraphGMM)

# example 1 ----------------------------------

# load data
data(plethodon)

# Generalized Procrustes Fit
GPAshape <- procGPA(plethodon$land)

# procrustes distances

procD_comparisons(GPAshape$rotated, plethodon$species)

# example 2 ----------------------------------

# load data
data(apes)

# Generalized Procrustes Fit
GPAshape <- procGPA(apes$x)

# procrustes distances

procD_comparisons(GPAshape$rotated, apes$group, plot_results = FALSE)

}
\seealso{
\code{\link{mahalanobis_dist_matrix}}
}
\author{
Lloyd A. Courtenay
}
