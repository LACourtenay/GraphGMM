% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_edge_list.R
\name{as_edge_list}
\alias{as_edge_list}
\title{Create an Edge List}
\usage{
as_edge_list(edge_matrix)
}
\arguments{
\item{edge_matrix}{A 2 column matrix specifying which landmarks are connected}
}
\value{
An edge list compatible with the igraph and graphGMM library
}
\description{
The present function can be used to convert a set of indices marking which edges connect
which landmark nodes for further processing in the \code{\link{graph_embeddings}} function
}
\examples{
library(shapes)
library(GraphGMM)

# load data
data(apes)

# Generalized Procrustes Fit
GPAshape <- GPA(apes$x)

# calculate central configuration
central_config <- calc_central_morph(GPAshape$coordinates)

# compute graph edges
edges <- triangulate2d(central_config)

# extract edge list
edge_list <- as_edge_list(edges)

# create graph embeddings
graph_object <- graph_embeddings(GPAshape$coordinates, edge_list)

}
\seealso{
\code{\link{graph_embeddings}}
}
