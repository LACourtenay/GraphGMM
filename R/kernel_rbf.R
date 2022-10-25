
#' The Radial Basis Function Kernel for Non-Linear Transformations
#'
#' @param X A matrix of data to be transformed
#' @param gamma The gamma hyperparameter defining a scaling degree
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
#' library(shapes)
#' library(GraphGMM)
#'
#' # load data
#' data(apes)
#'
#' # Generalized Procrustes Fit
#' GPAshape <- GPA(apes$x)
#' data <- vector_from_landmarks(GPAshape$rotated)
#'
#' # perform non linear transformation
#'
#' transformed_data <- kernel_rbf(data)
#'
#' # plot non-linear pca
#'
#' pca_plot(transformed_data, apes$group)
#'
#' @export

kernel_rbf <- function(X, gamma = 1.0) {

  if (!is.data.frame(X) & !is.matrix(X)) {
    stop("Input data must be of type data.frame or matrix")
  }

  if (gamma <= 0 | !is.numeric(gamma)) {
    stop("Gamma must be a positive numeric value")
  }

  X_tilde = as.matrix(X)

  for (i in 1:dim(X_tilde)[2]) {
    dimension <- X[,i]
    for (j in 1:dim(X_tilde)[1]) {
      X_tilde[j, i] <- norm(as.matrix(dimension) - dimension[j], "F")^2
    }
  }

  X_tilde <- exp(-gamma * X_tilde)

  results <- as.data.frame(X_tilde)
  rownames(results) <- rownames(X)
  colnames(results) <- colnames(X)

  return(results)

}
