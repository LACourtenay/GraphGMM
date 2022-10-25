
#' The Spline Kernel for Non-Linear Transformations
#'
#' @param X A matrix of data to be transformed
#'
#' @section Bibliography:
#' Gunn, S. (1998) Support Vector Machines for Classification and Regression,
#' ISIS Technical Report. http://www.svms.org/tutorials/Gunn1998.pdf
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
#' data <- vector_from_landmarks(GPAshape$coordinates)
#'
#' # perform non linear transformation
#'
#' transformed_data <- kernel_spline(data)
#'
#' # plot non-linear pca
#'
#' pca_plot(transformed_data, apes$group)
#'
#' @export

kernel_spline <- function(X) {

  if (!is.data.frame(X) & !is.matrix(X)) {
    stop("Input data must be of type data.frame or matrix")
  }

  X_tilde = as.matrix(X)

  for (i in 1:dim(X_tilde)[2]) {
    dimension <- as.matrix(X[,i])
    for (j in 1:dim(X_tilde)[1]) {

      XY <- dimension * dimension[j]
      X_plus_Y <- dimension + dimension[j]
      X_min <- dimension
      for (m in 1:length(dimension)) {
        if (X_min[m] > dimension[j]) {
          X_min[m] = dimension[j]
        }
      }
    }
    X_tilde[,i] <- 1 + XY + (XY * X_min) - ((X_plus_Y/2) * X_min^2) + ((1/3) * X_min^3)
  }

  results <- as.data.frame(X_tilde)
  rownames(results) <- rownames(X)
  colnames(results) <- colnames(X)

  return(results)

}
