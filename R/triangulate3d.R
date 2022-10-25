
#' 3D Ball Pivoting Triangulation
#'
#' A ball-pivoting triangulation function to compute the spatial relationships between landmarks
#'
#' @param central_config A (landmark x dimension) matrix containing the
#' central configuration that will be used to calculate the graph.
#' @param pivot_ball_radius The radius of the ball rolling over the 3D landmarks
#' @param pivot_ball_clustering The clustering radius of the ball.
#' @param pivot_ball_angle The angle threshold of the ball.
#'
#' @section Note:
#' The variables clustering and angle for the pivot ball can be used to avoid
#' the creation of too small triangles.
#'
#' @section Bibliography:
#' Bernardini, F.; Mittleman, J.; Rushmeier, H.; Silva, C.; Taubin, G. (1999)
#' The ball-pivoting algorithm for surface reconstruction.
#' IEEE Transaction on Visualization and Computer Graphics, 5(4):349-359
#'
#' @return A list of edges.
#'
#' @seealso \code{\link[Rvcg]{vcgBallPivoting}} of the \code{Rvcg} library
#'
#' @examples
#' library(shapes)
#' library(GraphGMM)
#'
#' # load data
#' data(macf.dat)
#'
#' # Generalized Procrustes Fit
#' GPAshape <- GPA(macf.dat)
#'
#' # calculate central configuration
#' central_config <- calc_central_morph(GPAshape$coordinates)
#'
#' # compute graph edges
#' edges <- triangulate3d(central_config)
#'
#' plot_landmark_graph(central_config, edges)
#'
#' @export

triangulate3d <- function(central_config, pivot_ball_radius = 0,
                          pivot_ball_clustering = 0.2,
                          pivot_ball_angle = pi/2) {

  if("array" %!in% class(central_config)) {
    if ("matrix" %!in% class(central_config)) {
      stop("Landmark central configuration must be a (landmark x dimension) array")
    }
  }

  if(!is.na(dim(central_config)[3])) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if(dim(central_config)[2] != 3) {
    stop("This function is only for the triangulation of 3D data")
  }

  if(pivot_ball_clustering <= 0 | pivot_ball_angle <= 0) {
    stop("Pivot ball clustering and angle parameters must be greater than 0")
  }

  if(pivot_ball_clustering < 0) {
    stop("Pivot ball radius must be equal to or greater than 0")
  }

  # using the ball pivoting algorithm

  mesh <- Rvcg::vcgBallPivoting(central_config,
                                radius = pivot_ball_radius,
                                clustering = pivot_ball_clustering,
                                angle = pivot_ball_angle)
  face_edges <- t(mesh$it)

  edges <- array(numeric(), dim = c(0, 2))

  if (dim(face_edges)[1] < 2) {
    stop("Pivot ball properties are too strict - no edges detected")
  }

  for (edge in 1:dim(face_edges)[1]) {
    edges <- abind::abind(
      edges,
      c(face_edges[edge, 1], face_edges[edge, 2]),
      along = 1
    )
    edges <- abind::abind(
      edges,
      c(face_edges[edge, 1], face_edges[edge, 3]),
      along = 1
    )
    edges <- abind::abind(
      edges,
      c(face_edges[edge, 2], face_edges[edge, 3]),
      along = 1
    )
  }

  edges <- edges[complete.cases(edges),]
  edges <- edges[!duplicated(edges),]

  return(edges)

}
