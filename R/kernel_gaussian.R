
#' The Gaussian Kernel for Non-Linear Transformations
#'
#' @param X A matrix of data to be transformed
#' @param alpha The alpha hyperparameter defining a scaling degree
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
#' transformed_data <- kernel_gaussian(data)
#'
#' # plot non-linear pca
#'
#' pca_plot(transformed_data, apes$group)
#'
#' @export

kernel_gaussian <- function(X, alpha = 2) {

  if (!is.data.frame(X) & !is.matrix(X)) {
    stop("Input data must be of type data.frame or matrix")
  }

  if (alpha <= 0 | !is.numeric(alpha)) {
    stop("Degree must be a positive numeric value")
  }

  X_tilde = as.matrix(X)

  for (i in 1:dim(X_tilde)[2]) {
    dimension <- X[,i]
    for (j in 1:dim(X_tilde)[1]) {
      X_tilde[j, i] <- -norm(abs(as.matrix(dimension) - dimension[j]), "F")
    }
  }

  X_tilde <- exp(as.matrix(X_tilde) / 2 * alpha^2)

  results <- as.data.frame(X_tilde)
  rownames(results) <- rownames(X)
  colnames(results) <- colnames(X)

  return(results)

}
