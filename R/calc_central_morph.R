
#' Mean or Median Central Landmark Configuration
#'
#' The present function can be used to calculate the central configuration of a set of
#' landmark coordinates. To comply with robust statistical parameters, this function provides
#' the option of both the median and mean landmark configuration
#'
#' @param array A (landmark x dimension x individual) array of landmark coordinates
#' @param method Parameter to define the calculation of either the "mean" or the "median" central configuration
#'
#' @return A matrix containing the mean or median central configuraion
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
#' @export

calc_central_morph <- function(array, method = "median") {

  if(!is.array(array)) {
    stop("Landmark configuration must be a (landmark x dimension x individual) array")
  }

  mshape <- array(numeric(), c(dim(array)[1], dim(array)[2]))

  if(method == "mean") {
    for(i in 1:dim(array)[1]) {
      mshape[i,1] = sum(array[i,1,]) / dim(array)[3]
      mshape[i,2] = sum(array[i,2,]) / dim(array)[3]

      if (dim(array)[2] == 3) {
        mshape[i,3] = sum(array[i,3,]) / dim(array)[3]
      }

    }
  } else if (method == "median") {
    for(i in 1:dim(array)[1]) {

      X = array[i,1,]
      Y = array[i,2,]

      if (dim(array)[2] == 3) {
        Z = array[i,3,]; Z= sort(Z)
      }

      X = sort(X); Y = sort(Y)
      mshape[i,1] = X[ceiling(dim(array)[3] * 0.5)]
      mshape[i,2] = Y[ceiling(dim(array)[3] * 0.5)]

      if (dim(array)[2] == 3) {
        mshape[i,3] = Z[ceiling(dim(array)[3] * 0.5)]
      }

    }
  } else {
    stop("Please select either 'mean' or 'median' method")
  }

  return(mshape)

}
