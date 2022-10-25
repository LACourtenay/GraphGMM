
#' Morphological Predictor
#'
#' The present function can be used to predict the landmark coordinates of an individual
#' given a set of morphological variables
#'
#' @param gpa_coordinates A (landmark x dimension x individuals) array (or tensor) containing
#' all landmark coordinates after Generalized Procrustes Analyses.
#' @param morph_descriptors A matrix of morphological descriptors that will be used to construct
#' the multivariate linear function that will be used to predict landmarks. e.g. a matrix
#' of PC scores.
#' @param target_pred A vector of morphological variables that will be used to predict their corresponding
#' landmarks.
#' @param robust A boolean value indicating whether robust metrics will be incorporated into the
#' prediction of landmarks.
#'
#' @section Details:
#' This function uses multivariate linear regression functions to construct a mathematical function
#' linking a set of morphological variables with landmark coordinates. Using this function,
#' the algorithm will then predict the landmark coordinates of a specific point in a
#' given feature space.
#'
#' @return A matrix of predicted landmark coordinates
#'
#' @seealso \code{tps_visualisation}
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
#' # predict the landmarks of a given point
#'
#' tps_example <- morphological_predictor(GPAshape$coordinates,
#'                                        pca$pc_scores[,1:2],
#'                                        target_pred = c(0.01, 0.01))
#'
#' # create simple plot of the predicted landmarks
#'
#' plot(tps_example, asp = 1, pch = 19)
#'
#' @export

morphological_predictor <- function(gpa_coordinates, morph_descriptors,
                                    target_pred, robust = FALSE) {

  if (missing(gpa_coordinates)) {
    stop("No GPA coordinates have been provided")
  }

  if (missing(morph_descriptors)) {
    stop("No morphological descriptors have been provided")
  }

  if (missing(target_pred)) {
    stop("No target morphological descriptors have been provided for a prediction")
  }

  if (!is.array(gpa_coordinates) | is.na(dim(gpa_coordinates)[3])) {
    stop("Landmark configuration must be a (landmark x dimension x individual) array")
  }

  if (!is.array(morph_descriptors) & !is.numeric(morph_descriptors)) {
    stop("Morphological descriptors must be an array or a vector of morphological variables, e.g. PC scores")
  }

  if (!is.logical(robust)) {
    stop("The robust parameter only accepts boolean TRUE or FALSE values")
  }

  p <- dim(gpa_coordinates)[1]
  k <- dim(gpa_coordinates)[2]
  n <- dim(gpa_coordinates)[3]

  if (is.array(morph_descriptors)) {

    if (length(target_pred) != ncol(morph_descriptors)) {
      stop(
        paste0(
          "The target predictor must be a vector of the same length as the number of columns in the",
          " morphological descriptor array"
        )
      )
    }

    if (n != nrow(morph_descriptors)) {
      stop("There must be the same number of morphological descriptors as individuals in the landmark array")
    }

  } else {

    if(is.numeric(morph_descriptors) & length(target_pred) != 1) {
      stop(
        paste0(
          "If only a single vector of morphological descriptors is provided then the",
          " target predictor must be a single numeric value"
        )
      )
    }

    if(n != length(morph_descriptors)) {
      stop("There must be the same number of morphological descriptors as individuals in the landmark array")
    }

  }

  lm_vector <- as.matrix(vector_from_landmarks(gpa_coordinates))
  morph_descriptors <- as.matrix(morph_descriptors)
  target_pred <- as.vector(target_pred)

  formula <- lm_vector ~ morph_descriptors
  formula <- update(formula, ~.+0)

  data_morph <- data.frame(lm_vector = lm_vector, morph_descriptors = morph_descriptors)

  fit_model <- lm(formula, data = data_morph)
  linear_coefficients <- fit_model$coefficients
  predicted_morphology <- crossprod(target_pred, linear_coefficients)

  if (robust == TRUE) {
    predicted_morphology <- matrix(predicted_morphology, p, k, byrow = TRUE) +
      calc_central_morph(gpa_coordinates, method = "median")
  } else {
    predicted_morphology <- matrix(predicted_morphology, p, k, byrow = TRUE) +
      calc_central_morph(gpa_coordinates, method = "mean")
  }

  return(predicted_morphology)

}
