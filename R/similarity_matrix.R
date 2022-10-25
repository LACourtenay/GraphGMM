
#' Calculate Similarity Matrix
#'
#' The present function can be used to calculate a similarity matrix from a given landmark
#' configuration.
#'
#' @param landmark_configuration A (landmark x dimension x individuals)
#' array (or tensor) containing landmark coordinates
#' @param method The method to be used to define the similarity matrix
#' ("cosine", "euclid" or "chebyshev").
#' @param unravel A boolean defining whether to turn matrices into similarity vectors
#' @param verbose A boolean value indicating whether a progress bar should be printed across
#' screen
#'
#' @section Description:
#' Landmark configurations can be converted
#' into a similarity matrix, calculating the structural similarity of each landmark in
#' relation to all landmarks within the configuration. For this purpose, three different
#' similarity matrices can be calculated (Cosine, Euclidean and Chebyshev),
#' however optimal results (and the most recommendable
#' approach), is obtained using the cosine similarity function.
#'
#' @section Bibliography:
#' Leskovec, J. (2019) Graph Node Embedding Algorithms, Stanford University CS224W:
#' Machine Learning with Graphs.
#'
#' @return similarity vectors or matrices for each individual
#'
#' @seealso \code{\link{graph_embeddings}}
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
#' lm_similarities <- similarity_matrix(GPAshape$coordinates)
#'
#' pca_plot(lm_similarities, apes$group)
#'
#' @export

similarity_matrix <- function(landmark_configuration, method = "cosine", unravel = TRUE,
                              verbose = TRUE) {

  if (!is.array(landmark_configuration)) {
    stop("Landmark configuration must be a (landmark x dimension x individual) array")
  }

  possible_metrics <- c("cosine", "euclid", "chebyshev")

  if (method %!in% possible_metrics) {
    stop("Incorrect metric")
  }

  lm_vector <- array(
    numeric(),
    dim = c(
      0,
      ((dim(landmark_configuration)[1] * dim(landmark_configuration)[1]) - dim(landmark_configuration)[1])/2
    )
  )

  similarity_whole_matrix <- array(
    numeric(),
    dim = c(dim(landmark_configuration)[1], dim(landmark_configuration)[1],
            0)
  )

  if (verbose == TRUE) {
    pb <- txtProgressBar(min = 0, max = dim(landmark_configuration)[3],
                         style = 3, width = 100, char = "=")
  }

  for (i in 1:dim(landmark_configuration)[3]) {
    dist_matrix <- matrix(numeric(), nrow = dim(landmark_configuration)[1],
                          ncol = dim(landmark_configuration)[1])
    for (lmi in 1:dim(landmark_configuration)[1]) {
      for (lmj in 1:dim(landmark_configuration)[1]) {

        if (method == "cosine") {
          dist_matrix[lmi, lmj] <- cosine_distance(landmark_configuration[lmi,,i],
                                                   landmark_configuration[lmj,,i])
        } else if (method == "euclid") {
          dist_matrix[lmi, lmj] <- philentropy::minkowski(landmark_configuration[lmi,,i],
                                                          landmark_configuration[lmj,,i],
                                                          n = 2, testNA = FALSE)[[1]]
        } else {
          dist_matrix[lmi, lmj] <- philentropy::chebyshev(landmark_configuration[lmi,,i],
                                                          landmark_configuration[lmj,,i],
                                                          testNA = FALSE)[[1]]
        }

      }

    }

    lm_vector <- abind::abind(lm_vector,
                              dist_matrix[lower.tri(dist_matrix)],
                              along = 1)

    similarity_whole_matrix <- abind::abind(
      similarity_whole_matrix,
      dist_matrix,
      along = 3
    )

    if (verbose == TRUE) {
      setTxtProgressBar(pb, i)
    }

  }

  name_matrix <- matrix(nrow = dim(landmark_configuration)[1],
                        ncol = dim(landmark_configuration)[1])

  for (lmi in 1:nrow(landmark_configuration)) {
    for (lmj in 1:nrow(landmark_configuration)) {
      name_matrix[lmi, lmj] <- paste("LM", lmi, "LM", lmj, sep = "")
    }
  }

  lm_vector <- as.data.frame(lm_vector)
  colnames(lm_vector) <- name_matrix[lower.tri(name_matrix)]

  if (unravel == TRUE) {
    return(lm_vector)
  } else {
    return(similarity_whole_matrix)
  }

}
