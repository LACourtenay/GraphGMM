
#' Create Landmark Graph
#'
#' This function can be used to create a landmark graph for further more detailed
#' network analyses of graph theory applications using the \code{igraph} library
#'
#' @param target_central A (landmark x dimension) matrix containing the central landmark
#' configuration the user wishes to analyse
#' @param target_edges An edge list that defines the graph, using
#' \code{\link{as_edge_list}}
#'
#' @section Details:
#' The objective of this function is to create an igraph object for more user specific
#' analyses (that may not be included in the GraphGMM library)
#'
#' @return An object of class igraph
#'
#' @seealso
#' \code{igraph} library
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
#' # calculate central configuration
#' central_config <- calc_central_morph(GPAshape$coordinates)
#'
#' # compute graph edges
#' edges <- triangulate2d(central_config)
#'
#' # extract edge list
#' edge_list <- as_edge_list(edges)
#'
#' # create igraph object
#' lm_graph <- create_lm_graph(central_config, edge_list)
#'
#' @export

create_lm_graph <- function(target_central, target_edges) {

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

  graph_df <- get_graph_df(target_central)

  LM_graph <- igraph::graph_from_data_frame(vertices = graph_df,
                                            d = as.matrix(target_edges),
                                            directed = FALSE)

  return(LM_graph)

}
