
#' The Cauchy Kernel for Non-Linear Transformations
#'
#' @param X A matrix of data to be transformed
#' @param sigma The sigma hyperparameter defining a scaling degree
#'
#' @section Bibliography:
#' Rahimi, A.; Recht, B. (2007) Random features for large-scale kernel
#' machines, Proceedings of the International Conference of Neural Information
#' Processing Systems, 20:1-8. DOI: 10.5555/2981562.2981710
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
#' transformed_data <- kernel_cauchy(data)
#'
#' # plot non-linear pca
#'
#' pca_plot(transformed_data, apes$group)
#'
#' @export

kernel_cauchy <- function(X, sigma = 1) {

  if (!is.data.frame(X) & !is.matrix(X)) {
    stop("Input data must be of type data.frame or matrix")
  }

  if (sigma <= 0 | !is.numeric(sigma)) {
    stop("Degree must be a positive numeric value")
  }

  X_tilde = as.matrix(X)

  for (i in 1:dim(X_tilde)[2]) {
    dimension <- X[,i]
    for (j in 1:dim(X_tilde)[1]) {
      X_tilde[j, i] <- norm(abs(as.matrix(dimension) - dimension[j]), "F")
    }
  }

  X_tilde <- (2 * (sigma * pi)) / (1 + X_tilde^2)

  results <- as.data.frame(X_tilde)
  rownames(results) <- rownames(X)
  colnames(results) <- colnames(X)

  return(results)

}
