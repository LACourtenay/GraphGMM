
#' The Polynomial Kernel for Non-Linear Transformations
#'
#' @param X A matrix of data to be transformed
#' @param degree The degree hyperparameter defining a scaling degree
#'
#' @section Bibliography:
#' Bishop, C. (2006) Pattern Recognition and Machine Learning. Singapore: Springer.
#'
#' Cortes, C.; Vapnik, V. (1995) Support Vector Networks, Machine Learning, 20:273-297
#' DOI: 10.1007/BF00994018
#'
#' @return Transformed matrix of values
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
#' data <- vector_from_landmarks(GPAshape$rotated)
#'
#' # perform non linear transformation
#'
#' transformed_data <- kernel_poly(data)
#'
#' # plot non-linear pca
#'
#' pca_plot(transformed_data, apes$group)
#'
#' @export

kernel_poly <- function(X, degree = 2) {

  if (!is.data.frame(X) & !is.matrix(X)) {
    stop("Input data must be of type data.frame or matrix")
  }

  if (degree <= 0 | !is.numeric(degree)) {
    stop("Degree must be a positive numeric value")
  }

  X_tilde = as.matrix(X)

  for (i in 1:dim(X_tilde)[2]) {
    dimension <- X[,i]
    for (j in 1:dim(X_tilde)[1]) {
      X_tilde[j, i] <- norm(as.matrix(dimension) * dimension[j], "F")^degree
    }
  }

  results <- as.data.frame(X_tilde)
  rownames(results) <- rownames(X)
  colnames(results) <- colnames(X)

  return(results)

}
