
#' Thin Plate Spline - Shape Predictor
#'
#' The present function can be used to calculate and visualise morphological variation given
#' a set of Procrustes superimposed coordinates and a feature space, based on the concept
#' of Thin Plate Splines (Bookstein, 1989)
#'
#' @param gpa_coordinates A (landmark x dimension x individuals) array (or tensor) containing
#' all landmark coordinates after Generalized Procrustes Analyses.
#' @param pca_data The PC Scores from the Principal Component Analysis
#' @param pcscore The PC Score(s) that will be used for shape prediction
#' @param pccoord The coordinates used for shape prediction
#' @param type The type of visualisation; 'graph' visualisation or 'surface' for 3D
#' surfaces.
#' @param edges A edge matrix used for 'graph' type visualisation
#' @param robust A boolean value indicating whether robust metrics will be incorporated into the
#' prediction of landmarks.
#'
#' @section Details:
#' This function generates a plot of the shape or form differences for a given point
#' in feature space (typically PCA). Using the \code{morphological_predictor} function,
#' the function first calculates the landmark coordinates of the point of interest,
#' followed by the visualisation of these changes using either a graph or a surface plot.
#'
#' @section Bibliography:
#' Bookstein, F.L. (1989) Principal Warps: thin plate spline and the decomposition of
#' deformations, Transactions on Pattern Analysis in Machine Intelligence.
#' 11(6):567-585
#'
#' @return A plotted visualisation of the predicted shape
#' @return The predicted landmark coordinates
#'
#' @author Lloyd A. Courtenay
#'
#' @seealso \code{morphological_predictor}, or \code{\link[geomorph]{shape.predictor}}
#' and \code{\link[geomorph]{plotRefToTarget}}
#' from the \code{geomorph} library
#'
#' @examples
#' library(shapes)
#' library(GraphGMM)
#'
#' # 2D example --------------------------------------
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
#' # create graph embeddings
#' graph_object <- graph_embeddings(GPAshape$coordinates, edge_list,
#'                                  num_convolutions = 2)
#'
#' pca <- pca_plot(graph_object$similarity_vector, apes$group)
#'
#' pca$pca_plot # visualise pca_plot
#'
#' # visualise pc1 extremities
#'
#' tps_visualisation(GPAshape$coordinates, pca$pc_scores, pcscore = 1,
#'                   pccoord = max(pca$pc_scores[,1]), type = "graph",
#'                   edges = edges)
#' tps_visualisation(GPAshape$coordinates, pca$pc_scores, pcscore = 1,
#'                   pccoord = min(pca$pc_scores[,1]), type = "graph",
#'                   edges = edges)
#'
#' # visualise pc2 extremities
#'
#' tps_visualisation(GPAshape$coordinates, pca$pc_scores, pcscore = 2,
#'                   pccoord = max(pca$pc_scores[,2]), type = "graph",
#'                   edges = edges)
#' tps_visualisation(GPAshape$coordinates, pca$pc_scores, pcscore = 2,
#'                   pccoord = min(pca$pc_scores[,2]), type = "graph",
#'                   edges = edges)
#'
#' # 3D example --------------------------------------
#'
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
#' # extract edge list
#' edge_list <- as_edge_list(edges)
#'
#' # create graph embeddings
#' graph_object <- graph_embeddings(GPAshape$coordinates, edge_list,
#'                                  num_convolutions = 2)
#'
#' pca <- pca_plot(graph_object$similarity_vector)
#'
#' pca$pca_plot # visualise pca_plot
#'
#' # visualise pc1 extremities as surface plot
#'
#' tps_visualisation(GPAshape$coordinates, pca$pc_scores, pcscore = 1,
#'                   pccoord = max(pca$pc_scores[,1]), type = "surface")
#' tps_visualisation(GPAshape$coordinates, pca$pc_scores, pcscore = 1,
#'                   pccoord = min(pca$pc_scores[,1]), type = "surface")
#'
#' # visualise pc1 extremities as graph plot
#'
#' tps_visualisation(GPAshape$coordinates, pca$pc_scores, pcscore = 1,
#'                   pccoord = max(pca$pc_scores[,1]), type = "graph",
#'                   edges = edges)
#' tps_visualisation(GPAshape$coordinates, pca$pc_scores, pcscore = 1,
#'                   pccoord = min(pca$pc_scores[,1]), type = "graph",
#'                   edges = edges)
#'
#' @export

tps_visualisation <- function(gpa_coordinates,
                              pca_data,
                              pcscore,
                              pccoord,
                              type,
                              edges = NULL,
                              reference_shape = NULL,
                              robust = FALSE) {

  if (missing(gpa_coordinates)) {
    stop(
      "No gpa_coordinates have been provided"
    )
  }

  if (missing(pca_data)) {
    stop(
      "No pca_data have been provided"
    )
  }

  if ("array" %!in% class(gpa_coordinates) | is.na(dim(gpa_coordinates)[3])) {
    stop("gpa_coordinates must be a (landmark x dimension x individual) array")
  }

  if (dim(gpa_coordinates)[2] != 2 & dim(gpa_coordinates)[2] != 3) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if ("matrix" %!in% class(pca_data)) {
    stop("pca_data must be a matrix")
  }

  if (dim(pca_data)[1] != dim(gpa_coordinates)[3]) {
    stop(
      "PCA data must be a matrix with the same number of rows as number of individuals in gpa_coordinates"
    )
  }

  if (missing(pcscore)) {
    stop(
      "A pcscore or set of two pcscores must be defined."
    )
  }

  if (missing(pccoord)) {
    stop(
      "A pccoord or set of two pccoords must be defined."
    )
  }

  if (missing(type)) {
    stop("The user must supply a visualisation type, between 'graph' or 'surface'")
  }

  possible_types <- c("graph", "surface")

  if (type %!in% possible_types) {
    stop(
      "visualisation type must be either 'graph' or 'surface'"
    )
  }

  if (type == possible_types[1] & is.null(edges)) {
    stop(
      "For 'graph' type visualisations, a valid edge list must be supplied"
    )
  }

  if (type == possible_types[2] & dim(gpa_coordinates)[2] != 3) {
    stop(
      "'surface' type visualisations are only available for 3D data"
    )
  }

  if (!is.null(edges)) {
    if("array" %!in% class(edges)) {
      if ("matrix" %!in% class(edges)) {
        stop("Edges must be an array or a matrix")
      }
    }
    if(dim(edges)[2] != 2) {
      stop("Edge matrix must be an array with 2 columns")
    }
  }

  if (length(pcscore) == 1) {
    if (!is.numeric(pcscore) |
        pcscore <= 0 |
        pcscore %% 1 != 0) {
      stop(
        "pcscore variable must be a single or set of numeric integer values specifying the pc score(s) to be visualised"
      )
    }
  } else if (length(pcscore) == 2) {
    if (length(pccoord) != 2) {
      stop(
        "If two pcscores have been defined then two pccoord must be defined"
      )
    }
    if (!is.numeric(pcscore[1]) | !is.numeric(pcscore[2]) |
        pcscore[1] <= 0 | !is.numeric(pcscore[2]) |
        pcscore[1] %% 1 != 0 | pcscore[2] %% 1 != 0) {
      stop(
        "pcscore variable must be a single numeric integer value specifying the pc score to be visualised"
      )
    }
  } else {
    stop(
      "Invalid pcscore parameter"
    )
  }

  if (!is.logical(robust)) {
    stop("The robust parameter only accepts boolean TRUE or FALSE values")
  }

  tps <- morphological_predictor(
    gpa_coordinates,
    pca_data[, pcscore],
    pccoord,
    robust = robust
  )

  if (type == possible_types[1]) {

    plot_landmark_graph(matrix(tps, nrow = dim(tps)[1], ncol = dim(tps)[2]), edges)

  } else if (type == possible_types[2]) {

    rgl::open3d(); rgl::shade3d(Rvcg::vcgBallPivoting(tps), col = "grey")

  }

  return(tps)

}
