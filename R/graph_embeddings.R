
#' Calculate Graph Embeddings
#'
#' The present function can be used to calculate embeddings for each landmark configuration,
#' producing a final output of each embedded individual, as well as the similarity matrix.
#'
#' @param landmarks A (num_landmark x num_dimension x num_individuals) array (or tensor) containing
#' superimposed landmark coordinates
#' @param edge_matrix A matrix defining the edges for the graph, either defined by the
#' user or by the \code{calc_edge_matrix} function.
#' @param num_convolutions An integer defining the number of convolutions to be used
#' to calculate the final embedding. num_convolutions can not be more than the diameter
#' of the graph.
#' @param dist_method The method to be used to define the final similarity matrix (
#' "cosine", "euclid" or "chebyshev").
#'
#' @section Details:
#' This function is the central component of the graphGMM library, used to project landmark
#' configurations into a new embedded feature space, and calculate the similarity matrix
#' that represents the structural properties of the entire configuration.
#'
#' The embedding component of this function aims to project landmarks into a new 2 or 3 dimensional
#' feature space, so that nodes in the graph that are structurally similiar are embedded closer
#' together. In this sense, similarity found in the new embedded landmark configuration,
#' should approximate the similarity observed in biological or geometric structure.
#'
#' For embedding, the present algorithm uses graph theory and the concept of a graph convolution,
#' typically used in graph based representation learning. Each landmark configuration
#' is thus represented as an undirected graph \eqn{G}, where every landmark is a node
#' \eqn{LM_{v} \in V}, and is connected to neighbour \eqn{LM_{u}} via an edge. \eqn{G}
#' can thus be represented as an adjacency matrix \eqn{A \in \mathbb{R}^{n \times n}},
#' with a feature matrix \eqn{X \in \mathbb{R}^{n \times d}} of procrustes landmark
#' coordinates, where \eqn{n} is the number of landmarks, and \eqn{d} is the number of dimensions.
#'
#' The embedding process thus consists in a message passing and neighbourhood aggregation
#' scheme, whereby each landmark is transformed to be a function of itself and its neighbours.
#' For this purpose, the convolution formula proposed by Kipf and Welling (2017) was adapted and
#' implemented, using a standardised version of the adjacency matrix transformed with landmark
#' centrality degrees, as well as self loops, to perform a single message-passing convolution on the graph.
#'
#' To ensure that the embedding includes information from landmarks further away across the graph,
#' then the depth of the embedding can be increased using a higher number of convolutions.
#' Nevertheless, the present function will not compute if the number of convolutions
#' is higher than the natural diameter of the graph.
#'
#' Once the landmarks have been embedded, the transformed configuration is converted
#' into a similarity matrix, calculating the structural similarity of each landmark in
#' relation to all landmarks within the configuration. For this purpose, three different
#' similarity matrices can be calculated (Cosine, Euclidean and Chebyshev),
#' however optimal results (and the most recommendable
#' approach), is obtained using the cosine similarity function.
#'
#' @section Bibliography:
#' Kipf, T.N.; Welling, M. (2017) Semi-supervised classification with graph
#' convolutional networks. International Conference on Learning Representations. arXiv:
#' 1609.02907v4
#'
#' Leskovec, J. (2019) Graph Node Embedding Algorithms, Stanford University CS224W:
#' Machine Learning with Graphs.
#'
#' Wang, Y.; Sun, Y.; Lui, Z.; Sarma, S.E.; Bronstein, M.M.; Solomon, J.M. (2019)
#' Dynamic Graph CNN for Learning on Point Clouds, ACM Transaction on Graphics,
#' 1(1):1-13, arXiv: 1801.07829v2
#'
#' Bronstein, M.M.; Bruna, J.; Cohen, T.; Velickovic, P. (2021) Geometric Deep Learning
#' Grids, Groups, Graphs, Geodesics, and Gauges, arXiv: 2104.13478v2
#'
#' Courtenay, L.A.; Aramendi, J.; González-Aguilera, D. (In Prep) A Graph
#' Based Geometric Morphometric approach to the analysis of primate radii:
#' A new mathematical model for the processing of landmark data.
#'
#' @return \code{similarity_vector} - similarity vectors for each individual
#' @return \code{landmark_embeddings} - an array of the embedded landmarks
#'
#' @author Lloyd A. Courtenay
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
#' graph_object <- graph_embeddings(GPAshape$rotated, edge_list,
#'                                  num_convolutions = 2)
#'
#' @export

graph_embeddings <- function(landmarks, edge_matrix, num_convolutions = 2,
                            dist_method = "cosine") {

  if (!is.array(landmarks) | is.na(dim(landmarks)[3])) {
    stop("Landmark configuration must be a (landmark x dimension x individual) array")
  }

  #if(dim(central_config)[2] != 2 & dim(central_config)[2] != 3) {
  #  stop("Landmark central configuration must be a (landmark x dimension) array")
  #}

  if (!is.matrix(edge_matrix)) {
    stop("Edge matrix must be a (landmark x 2) array")
  }

  if(!is.numeric(num_convolutions) | num_convolutions <= 0 | num_convolutions %% 1 != 0) {
    stop("The number of convolutions must be a non-negative integer, e.g. 2")
  }

  possible_metrics <- c("cosine", "euclid", "chebyshev")

  if (dist_method %!in% possible_metrics) {
    stop("Please select a method between cosine, euclid or chebyshev")
  }

  lm_vector <- array(
    numeric(),
    dim = c(
      0,
      ((dim(landmarks)[1] * dim(landmarks)[1]) - dim(landmarks)[1])/2
    )
  )

  embedded_landmarks <- array(
    numeric(),
    dim = dim(landmarks)
  )

  # for loop for the number of convolutions, fill an array with the respected embedded feature matrix
  cat("\nEmbedding data...")
  pb <- txtProgressBar(min = 0, max = dim(landmarks)[3],
                       style = 3, width = 100, char = "=")

  for (i in 1:dim(landmarks)[3]) {

    # create graph vertex attributes

    graph_df <- get_graph_df(landmarks[,,i])

    # define graph

    LM_graph <- igraph::graph_from_data_frame(vertices = graph_df,
                                              d = as.matrix(edge_matrix),
                                              directed = FALSE)

    if (num_convolutions > igraph::diameter(LM_graph)) {
      stop(paste(
        "\nThe number of convolutions cannot be larger than the diameter of the landmark graph.",
        "\nThe current landmark configuration has a diameter of: ",
        igraph::diameter(LM_graph), sep = ""
      ))
    }

    # compute adjacency matrix, feature vector, identity matrix, degree matrix

    adjacency_matrix_A <- as.matrix(igraph::as_adjacency_matrix(LM_graph))
    feature_vector_X <- as.matrix(graph_df[,2:ncol(graph_df)])
    identity_matrix_I <- diag(dim(adjacency_matrix_A)[1])
    new_adj_matrix_A_tilde <- adjacency_matrix_A + identity_matrix_I
    deg <- igraph::degree(igraph::graph.adjacency(new_adj_matrix_A_tilde,
                                          mode = "undirected"))
    inverse_D_tilde <- sqrt(diag(deg^-1))

    feature_embedding <- feature_vector_X
    for (conv in 1:num_convolutions) {
      feature_embedding <- (inverse_D_tilde %*%
                              new_adj_matrix_A_tilde %*%
                              inverse_D_tilde) %*% feature_embedding
    }

    # calculate similarity distances between embedded landmarks

    dist_matrix <- matrix(numeric(), nrow = dim(landmarks)[1],
                          ncol = dim(landmarks)[1])
    for (lmi in 1:nrow(feature_embedding)) {
      for (lmj in 1:nrow(feature_embedding)) {

        if (dist_method == "cosine") {
          dist_matrix[lmi, lmj] <- lsa::cosine(feature_embedding[lmi,],
                                               feature_embedding[lmj,])[[1]]
        }

        if (dist_method == "euclid") {
          dist_matrix[lmi, lmj] <- philentropy::minkowski(feature_embedding[lmi,],
                                                          feature_embedding[lmj,],
                                                          n = 2, testNA = FALSE)[[1]]
        }

        if (dist_method == "chebyshev") {
          dist_matrix[lmi, lmj] <- philentropy::chebyshev(feature_embedding[lmi,],
                                                          feature_embedding[lmj,],
                                                          testNA = FALSE)[[1]]
        }

      }
    }

    # extract embedded vector for each individual

    lm_vector <- abind::abind(lm_vector,
                              dist_matrix[lower.tri(dist_matrix)],
                              along = 1)

    embedded_landmarks[,,i] <- feature_embedding

    setTxtProgressBar(pb, i)

  }

  name_matrix <- matrix(nrow = dim(landmarks)[1], ncol = dim(landmarks)[1])
  for (lmi in 1:nrow(feature_embedding)) {
    for (lmj in 1:nrow(feature_embedding)) {
      name_matrix[lmi, lmj] <- paste("LM", lmi, "LM", lmj, sep = "") # is cosine appr.
    }
  }
  lm_vector <- as.data.frame(lm_vector)
  colnames(lm_vector) <- name_matrix[lower.tri(name_matrix)]

  return(list(
    similarity_vector = lm_vector,
    landmark_embeddings = embedded_landmarks
  ))

}
