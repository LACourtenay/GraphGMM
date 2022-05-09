
#' 2D Delauney Triangulation
#'
#' A delauney triangulation function to compute the spatial relationships between landmarks
#'
#' @param central_config A (landmark x dimension) matrix containing the
#' central configuration that will be used to calculate the graph.
#'
#' @section Bibliography:
#' Delaunay, B. (1934) Sur la sphère vide.
#' Bulletin l’Academie de Sciences l’URSS Classe des Sciences Mathematiques et Naturelles
#' 6:793–800.
#'
#' @return A list of edges.
#'
#' @seealso \code{\link[RTriangle]{triangulate}} of the \code{RTriangle} library
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
#' plot_landmark_graph(central_config, edges)
#'
#' @export

triangulate2d <- function(central_config) {

  if("array" %!in% class(central_config)) {
    if ("matrix" %!in% class(central_config)) {
      stop("Landmark central configuration must be a (landmark x dimension) array")
    }
  }

  if(dim(central_config)[2] != 2) {
    stop("This function is only for the triangulation of 2D data")
  }

  tp <- RTriangle::triangulate(RTriangle::pslg(central_config))
  edges <- tp$E

  return(edges)

}
