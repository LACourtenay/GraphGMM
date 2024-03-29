% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separate_samples.R
\name{separate_samples}
\alias{separate_samples}
\title{Miscellaneous: A function to separate landmark samples}
\usage{
separate_samples(GPAmatrix, labels, target_sample)
}
\arguments{
\item{GPAmatrix}{A (landmark x dimension x individual) array containing the landmark
data}

\item{labels}{A factorial list containing the sample labels}

\item{target_sample}{A string defining the name of the sample the user wishes to
separate}
}
\value{
A (landmark x dimension x individual) array of landmarks from the defined sample
}
\description{
This function can be used to separate landmarks from a specific sample.
}
\examples{
library(shapes)
library(GraphGMM)

# load data
data(apes)

# Generalized Procrustes Fit
GPAshape <- GPA(apes$x)

# separate the samples under analysis

female_gorilla <- separate_samples(GPAshape$coordinates,
                                   apes$group,
                                   "gorf")
male_gorilla <- separate_samples(GPAshape$coordinates,
                                 apes$group,
                                 "gorm")

central_female <- calc_central_morph(female_gorilla)
central_male <- calc_central_morph(male_gorilla)

plot(central_female, pch = 19, asp = 1)
plot(central_male, pch = 19, asp = 1)

}
