
#' k-Nearest Neighbour Calculation
#'
#' A nearest neighbour function to compute the spatial relationships between landmarks
#'
#' @param central_config A (landmark x dimension) matrix containing the
#' central configuration that will be used to calculate the graph.
#' @param k An integer defining the number of Nearest Neighbours.
#' @param radius A numeric value defining the radius around a point
#' to be considered a neighbour. If a radius is not provided, then all points can
#' be considered neighbours.
#' @param verbose A boolean value indicating whether a progress bar should be printed across
#' screen
#'
#' @section Details:
#' This function calculates the nearest neighbours of each landmark.
#'
#' @return A list of edges.
#'
#' @author Lloyd A. Courtenay
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
#' edges <- knn_graph(central_config)
#'
#' plot_landmark_graph(central_config, edges)
#'
#' @export

knn_graph <- function(central_config, k = 3, radius = NULL,
                      verbose = TRUE) {

  if("array" %!in% class(central_config)) {
    if ("matrix" %!in% class(central_config)) {
      stop("Landmark central configuration must be a (landmark x dimension) array")
    }
  }

  if(dim(central_config)[2] != 2 & dim(central_config)[2] != 3) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if(!is.numeric(k) | k <= 0 | k %% 1 != 0) {
    stop("k must be a non-negative integer, e.g. 2")
  }

  if(k > dim(central_config)[1]) {
    stop("k values must be lower than the total number of landmarks")
  }

  if(!is.null(radius)) {
    if (!is.numeric(radius) & radius <= 0) {
      stop(
        "Radius parameter must be a positive numeric value."
      )
    }
  }

  if (verbose == TRUE) {
    cat("\nCalculating 3D distances\n")
    pb_distances <- txtProgressBar(min = 0, max = dim(central_config)[1],
                                   style = 3, width = 100, char = "=")
  }

  distances <- c()
  for (lm_from in 1:dim(central_config)[1]) {
    for (lm_to in 1:dim(central_config)[1]) {
      distances <- c(
        distances,
        philentropy::euclidean(central_config[lm_from,], central_config[lm_to,],
                               testNA = FALSE)
      )
    }

    if (verbose == TRUE) {
      setTxtProgressBar(pb_distances, lm_from)
    }

  }
  distance_matrix <- matrix(distances, ncol = dim(central_config)[1])
  name_matrix <- matrix(nrow = dim(central_config)[1], ncol = dim(central_config)[1])
  for (lmi in 1:dim(central_config)[1]) {
    for (lmj in 1:dim(central_config)[1]) {
      name_matrix[lmi, lmj] <- paste(lmi, ",", lmj, sep = "")
    }
  }
  single_direction_distances <- distance_matrix[lower.tri(distance_matrix)]
  name_matrix <- name_matrix[lower.tri(name_matrix)]
  name_vectors <- matrix(as.numeric(unlist(strsplit(name_matrix, ","))),
                         ncol = 2, byrow = TRUE)
  distance_table <- data.frame(from = name_vectors[,2], to = name_vectors[,1], dist = single_direction_distances)

  if (verbose == TRUE) {
    cat("\nCalculating nearest neighbours\n")
    pb_graph <- txtProgressBar(min = 0, max = dim(central_config)[1],
                               style = 3, width = 100, char = "=")
  }

  edges <- array(numeric(), dim = c(0, 2))
  for (i in 1:dim(central_config)[1]) {
    target_lm <- distance_table[distance_table$from == i,]
    target_lm <- target_lm[order(target_lm$dist, decreasing = FALSE),]
    if (is.null(radius)) {
      lm_nn <- target_lm[1:k,]
      for(neighbours in 1:nrow(lm_nn)) {
        edges <- abind::abind(edges, lm_nn[,1:2], along = 1)
      }
    } else {
      lm_nn <- target_lm[1:k,]
      lm_nn <- lm_nn[lm_nn$dist < radius, ]
      for(neighbours in 1:nrow(lm_nn)) {
        edges <- abind::abind(edges, lm_nn[,1:2], along = 1)
      }
    }

    if (verbose == TRUE) {
      setTxtProgressBar(pb_graph, lm_from)
    }

  }

  edges <- edges[complete.cases(edges),]
  edges <- edges[!duplicated(edges),]

  return(edges)

}
