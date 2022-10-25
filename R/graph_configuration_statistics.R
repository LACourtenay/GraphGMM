
#' Calculate Graph Configuration Statistics
#'
#' This function can be used to calculate the properties of a given landmark graph.
#'
#' @param target_central A (landmark x dimension) matrix containing the central landmark
#' configuration the user wishes to analyse
#' @param target_edges An edge list that defines the graph, using
#' \code{\link{as_edge_list}}
#' @param plot A boolean variable defining whether to plot a histogram of the results
#'
#' @section Details:
#' The objective of this function is to calculate a set of graph-theory based statistics
#' that can be used to describe landmark configurations. This function calculates the
#' overall density of the graph as well as its transitivity - how well the landmarks are
#' connected among themselves - as well as landmark centrality values - how important a
#' landmark is in the general context of the landmark configuration.
#'
#' @return \code{graph_clustering_coe} - The graph clustering coefficient
#' @return \code{graph_density} - The density of edges in the graph
#' @return \code{landmark_degree} - A list of degree values corresponding to each landmark,
#' indicating the number of connections a landmark has with other landmarks
#' @return \code{landmark_eigenvalues} - A list of landmark centrality values according to
#' their eigenvalues, identifying landmarks that are important in the global
#' context of the entire landmark configuration
#' @return \code{landmark_betweenness} - A list of landmark centrality values according to
#' their 'betweenness' (the number of times a landmark acts as a bridge to other landmarks
#' within the configuration).
#'
#' @seealso
#' From the present library; \code{\link{separate_samples}}, \code{\link{plot_lmgraph_stats}}
#'
#' From the \code{igraph} library; \code{\link[igraph]{edge_density}},
#' \code{\link[igraph]{transitivity}}, \code{\link[igraph]{centr_degree}},
#' \code{\link[igraph]{centr_betw}}, \code{\link[igraph]{centr_eigen}}
#'
#' @examples
#' library(shapes)
#' library(GraphGMM)
#'
#' # load data
#' data(apes)
#'
#' # Generalized Procrustes Fit
#' GPAshape <- GPA(apes$x)
#'
#' # separate the samples under analysis
#'
#' female_gorilla <- separate_samples(GPAshape$rotated,
#'                                    apes$group,
#'                                    "gorf")
#' male_gorilla <- separate_samples(GPAshape$rotated,
#'                                  apes$group,
#'                                  "gorm")
#'
#' female_central <- calc_central_morph(female_gorilla)
#' male_central <- calc_central_morph(male_gorilla)
#'
#' # calculate edges
#'
#' female_edges <- knn_graph(female_central, k = 3, radius = 75)
#' male_edges <- knn_graph(male_central, k = 3, radius = 75)
#' female_edges <- as_edge_list(female_edges)
#' male_edges <- as_edge_list(male_edges)
#'
#' # Graph configuration statistics
#'
#' graph_configuration_statistics(female_central, female_edges)
#' graph_configuration_statistics(male_central, male_edges)
#'
#' @export

graph_configuration_statistics <- function(target_central,
                                           target_edges,
                                           plot = TRUE) {

  if("array" %!in% class(target_central) | "array" %!in% class(target_edges)) {
    if ("matrix" %!in% class(target_central) | "matrix" %!in% class(target_edges)) {
      stop("Landmark central configuration must be a (landmark x dimension) array")
    }
  }

  if(dim(target_central)[2] != 2 & dim(target_central)[2] != 3) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if (dim(target_edges)[1] == 0) {
    stop("The number of edges is equal to 0. Please provide a valid edge matrix")
  }

  LM_graph <- create_lm_graph(target_central, target_edges)

  # ratio of edges to the number of possible edges
  edge_density <- igraph::edge_density(LM_graph, loops = FALSE)

  # probability that adjacent vertices are connected - clustering coefficient
  transitivity <- igraph::transitivity(LM_graph, type = "global")

  degree <- igraph::centr_degree(LM_graph, mode = "in", normalized = TRUE)$res

  eigenvalues <- igraph::centr_eigen(LM_graph, directed = FALSE, normalized = TRUE)$vector

  betweeness <- igraph::centr_betw(LM_graph, directed = FALSE, normalized = TRUE)$res

  if (plot == TRUE) {
    par(mfrow = c(1,3), mar = c(5.1, 5, 4.1, 2.))
    hist(degree,
         #breaks = igraph::vcount(LM_graph)/15,
         main = "Landmark Degree Distribution",
         col = "dark gray",
         xlab = "", ylab = "")
    mtext(side = 1, line = 3, "Degree", cex = 1, font = 2)
    mtext(side = 2, line = 3, "Frequency", cex = 1, font = 2)
    hist(eigenvalues,
         #breaks = igraph::vcount(LM_graph)/15,
         main = "Landmark Eigenvector Distribution",
         col = "dark gray",
         xlab = "", ylab = "")
    mtext(side = 1, line = 3, "Eigenvalue", cex = 1, font = 2)
    mtext(side = 2, line = 3, "Frequency", cex = 1, font = 2)
    hist(betweeness,
         #breaks = igraph::vcount(LM_graph)/15,
         main = "Landmark Betweeness Distribution",
         col = "dark gray",
         xlab = "", ylab = "")
    mtext(side = 1, line = 3, "Betweeness", cex = 1, font = 2)
    mtext(side = 2, line = 3, "Frequency", cex = 1, font = 2)
    par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
  }

  return(list(graph_clustering_coe = transitivity,
              graph_density = edge_density,
              landmark_degree = degree,
              landmark_eigenvalues = eigenvalues,
              landmark_betweenness = betweeness))

}
