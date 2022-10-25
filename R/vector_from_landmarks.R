
#' Miscellaneous: A function to convert a landmark array into a vector
#'
#' This function can be used to convert a landmark array into a vector for dimensionality
#' reduction (PCA, non-linear PCA, or tSNE)
#'
#' @param landmark_array A (landmark x dimension x individual) array containing the landmark
#' data
#'
#' @return A (individual x (landmark * dimension)) matrix of landmarks
#'
#' @seealso
#' \code{\link{kernel_pca}}, \code{\link{kernel_pca_biplot}},
#' \code{\link{pca_plot}}, \code{\link{pca_biplot}}
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
#' data <- vector_from_landmarks(GPAshape$rotated)
#'
#' # plot non-linear pca
#'
#' pca_plot(data, apes$group)
#'
#' @export

vector_from_landmarks <- function(landmark_array) {

  if ("array" %!in% class(landmark_array)) {
    stop("Landmarks must be in a (landmark x dimension x individual) array")
  }

  if (is.na(dim(landmark_array)[3])) {
    stop("Landmarks must be in a (landmark x dimension x individual) array")
  }

  columns <- c()
  if (dim(landmark_array)[2] == 2) {
    for(i in 1:dim(landmark_array)[1]) {
      columns <- c(columns, paste("X", i, sep = ""))
      columns <- c(columns, paste("Y", i, sep = ""))
    }
  } else if (dim(landmark_array)[2] == 3) {
    for(i in 1:dim(landmark_array)[1]) {
      columns <- c(columns, paste("X", i, sep = ""))
      columns <- c(columns, paste("Y", i, sep = ""))
      columns <- c(columns, paste("Z", i, sep = ""))
    }
  } else {
    stop("Function only accepts 2D or 3D landmark configurations")
  }

  flat_proc <- c()
  for (i in 1:dim(landmark_array)[3]) {
    flat <- flatten_matrix(landmark_array[,,i])
    flat_proc <- rbind(flat_proc, flat)
  }
  flat_proc <- as.data.frame(flat_proc)
  colnames(flat_proc) <- columns
  rownames(flat_proc) <- seq(1:dim(landmark_array)[3])

  return(flat_proc)

}
