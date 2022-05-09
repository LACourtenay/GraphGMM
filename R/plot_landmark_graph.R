
#' Plot a landmark graph
#'
#' The present function can be used to visualise the computed graphs on a given landmark
#' configuration.
#'
#' @param central_config A (landmark x dimension) matrix containing the
#' central configuration that will be used to calculate the graph.
#' @param edges A matrix defining the edges for the graph.
#' @param point_size A numeric value defining the size of the points in the visualisation
#' @param point_colour A string defining the colour of points.
#' @param line_size A numeric value defining the thickness of edge lines
#' @param line_colour A string defining the colour of edges.
#'
#' @return A plot of the landmark graph
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
#' edges <- knn_graph(central_config)
#'
#' plot_landmark_graph(central_config, edges)
#'
#' @export

plot_landmark_graph <- function(central_config, edges,
                                point_size = NULL,
                                point_colour = NULL,
                                line_size = NULL,
                                line_colour = NULL) {

  if("array" %!in% class(central_config)) {
    if ("matrix" %!in% class(central_config)) {
      stop("Landmark central configuration must be a (landmark x dimension) array")
    }
  }

  if(!is.na(dim(central_config)[3])) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if("array" %!in% class(edges)) {
    if ("matrix" %!in% class(edges)) {
      stop("Edges must be an array")
    }
  }

  if(dim(central_config)[2] != 2 & dim(central_config)[2] != 3) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if(dim(edges)[2] != 2) {
    stop("Edge matrix must be an array with 2 columns")
  }

  if(dim(central_config)[2] == 2) {
    if(is.null(point_size)) {
      point_size = 1
    }
    if(is.null(line_size)) {
      line_size = 1
    }
  } else {
    if(is.null(point_size)) {
      point_size = 10
    }
    if(is.null(line_size)) {
      line_size = 1
    }
  }
  if(is.null(point_colour)) {
    point_colour = "black"
  }
  if(is.null(line_colour)) {
    line_colour = "red"
  }

  if (dim(central_config)[2] == 3) {
    rgl::open3d()
    rgl::points3d(
      central_config[,1], central_config[,2], central_config[,3],
      size = point_size,
      col = point_colour
    )
    for (i in 1:nrow(edges)) {
      rgl::segments3d(
        x = as.vector(c(
          central_config[edges[i,1], 1], central_config[edges[i,2], 1]
        )),
        y = as.vector(c(
          central_config[edges[i,1],2], central_config[edges[i,2],2]
        )),
        z = as.vector(c(
          central_config[edges[i,1],3], central_config[edges[i,2],3]
        )),
        col = line_colour,
        lwd = line_size
      )
    }
  } else {
    plot(central_config, asp = 1, pch = 19, col = NA, xlab = "", ylab = "")
    for(link in 1:nrow(edges)) {
      segments(
        central_config[edges[link, 1], 1],
        central_config[edges[link, 1], 2],
        central_config[edges[link, 2], 1],
        central_config[edges[link, 2], 2],
        col = line_colour,
        lwd = line_size
      )
    }
    points(central_config[,1], central_config[,2], cex = point_size, pch = 19,
           col = point_colour)
  }

}
