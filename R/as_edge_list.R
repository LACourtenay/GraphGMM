
#' Create an Edge List
#'
#' The present function can be used to convert a set of indices marking which edges connect
#' which landmark nodes for further processing in the \code{\link{graph_embeddings}} function
#'
#' @param edge_matrix A 2 column matrix specifying which landmarks are connected
#'
#' @return An edge list compatible with the igraph and graphGMM library
#'
#' @seealso \code{\link{graph_embeddings}}
#'
#' @examples
#' library(geomorph)
#' library(shapes)
#' library(GraphGMM)
#'
#' # load data
#' data(apes)
#'
#' # Generalized Procrustes Fit
#' GPAshape <- procGPA(apes$x)
#'
#' # calculate central configuration
#' central_config <- calc_central_morph(GPAshape$rotated)
#'
#' # compute graph edges
#' edges <- triangulate2d(central_config)
#'
#' # extract edge list
#' edge_list <- as_edge_list(edges)
#'
#' # create graph embeddings
#' graph_object <- graph_embeddings(GPAshape$rotated, edge_list)
#' 
#' @export

as_edge_list <- function(edge_matrix) {

  if(dim(edge_matrix)[2] != 2) {
    stop("Edge matrices must have a maximum of 2 columns")
  }

  if(dim(edge_matrix)[1] == 0) {
    stop("The number of edges is equal to 0. Please provide a valid edge matrix")
  }

  colnames(edge_matrix) <- c("Source", "Target")
  for(row in 1:nrow(edge_matrix)) {
    for(col in 1:ncol(edge_matrix)) {
      edge_matrix[row, col] <- paste("LM", edge_matrix[row, col], sep = "")
    }
  }
  return(edge_matrix)

}
