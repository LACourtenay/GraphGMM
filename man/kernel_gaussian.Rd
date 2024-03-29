% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_gaussian.R
\name{kernel_gaussian}
\alias{kernel_gaussian}
\title{The Gaussian Kernel for Non-Linear Transformations}
\usage{
kernel_gaussian(X, alpha = 2)
}
\arguments{
\item{X}{A matrix of data to be transformed}

\item{alpha}{The alpha hyperparameter defining a scaling degree}
}
\value{
Transformed matrix of values
}
\description{
The Gaussian Kernel for Non-Linear Transformations
}
\section{Bibliography}{

Rahimi, A.; Recht, B. (2007) Random features for large-scale kernel
machines, Proceedings of the International Conference of Neural Information
Processing Systems, 20:1-8. DOI: 10.5555/2981562.2981710
}

\examples{
library(shapes)
library(GraphGMM)

# load data
data(apes)

# Generalized Procrustes Fit
GPAshape <- GPA(apes$x)
data <- vector_from_landmarks(GPAshape$coordinates)

# perform non linear transformation

transformed_data <- kernel_gaussian(data)

# plot non-linear pca

pca_plot(transformed_data, apes$group)

}
\author{
Lloyd A. Courtenay
}
