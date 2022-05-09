
#' Orthogonal Projection of Procrustes or Kendall Coordinates
#'
#' The present function performs a orthogonal projection of coordinates
#' from non-linear Kendall space or the Procrustes hemisphere into their corresponding
#' position in tangent space
#'
#' @param kendall_coordinates A (num_landmark x num_dimension x num_individuals)
#' array (or tensor) containing kendall or Procrustes coordinates.
#' @param scale A boolean value indicating whether scaling should be used in the
#' calculation of the pole (central configuration).
#' @param robust A boolean value indicating whether a robust statistics should be used.
#'
#' @section Bibliography:
#'
#' Kendall, D.G. (1984) Shape maniforlds, procrustean metrics, and complex projective
#' spaces, Bulletin of the London Mathematical Society. 16:81-121
#'
#' Rohlf, J.F.; Slice, D.E. (1990) Extension of the Procrustes method for the
#' optimal superimposition of landmarks, Systematic Biology. 39:40-59
#'
#' Bookstein, F.L. (1986) Size and shape spaces for landmark data in two dimensions,
#' Statistical Science. 4(2):181-242
#'
#' Goodall, C. (1991) Procrustes Methods in the Statistical Analysis of Shape,
#' Journal of the Royal Statistical Society B. 53(2):285-339
#'
#' Rohlf, J.F. (1996) Morphometric Spaces, Shape Components and the Effects of
#' Linear Transformations. In: L.F. Marcus, M. Corti, A. Loy, G.J.P. Naylor, D.E.
#' Slice (Eds.) Advances in Morphometrics. The Netherlands: Springer. 117-129
#'
#' Rohlf, J.F. (1999) Shape Statistics: Procrustes Superimpositions and Tangent
#' Spaces, Journal of Classification, 16:197-223
#'
#' Slice, D. (2001) Landmark coordinates aligned by Procrustes Analysis do not
#' lie in Kendall’s shape space, Systematic Biology. 50(1):141-149.
#' DOI: 10.1080/10635150119110
#'
#' @return The non-euclidean coordinates projected orthogonally into tanjent space
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' library(shapes)
#' data(apes)
#'
#' # shape analysis
#'
#' GPAshape <- GPA(apes$x, proc_method = "LS")
#'
#' orthogonal_projection(GPAshape$procrustes_coordinates)
#'
#' @export

orthogonal_projection <- function(kendall_coordinates, scale = TRUE,
                                  robust = FALSE) {

  if (!is.array(kendall_coordinates) | is.na(dim(kendall_coordinates)[3])) {
    stop("Kendall_coordinates must be a (landmark x dimension x individual) array in
         some form of Non-Euclidean shape space")
  }

  if (!is.logical(robust)) {
    stop("The robust parameter only accepts boolean TRUE or FALSE values")
  }

  p <- dim(kendall_coordinates)[1]
  k <- dim(kendall_coordinates)[2]
  n <- dim(kendall_coordinates)[3]

  if (robust == TRUE) {
    central_config <- calc_central_morph(kendall_coordinates, method = "median")
  } else {
    central_config <- calc_central_morph(kendall_coordinates, method = "mean")
  }

  if (scale == TRUE) {
    central_centroid <- find_centroid(central_config, robust = robust)
    central_centroid_size <- centroid_size(central_config, central_centroid)
    central_shape <- central_config / central_centroid_size
    vector_centre <- as.vector(central_shape)
  } else {
    vector_centre <- as.vector(central_config)
  }

  xm <- as.matrix(rep(1, n)) %*% vector_centre
  identity_matrix <- diag(1, k * p)
  vectorized_lm_matrix <- matrix(NA, n, k * p)

  for (individual in 1:n) {
    vectorized_lm_matrix[individual,] <- as.vector(kendall_coordinates[,,individual])
  }
  xp <- vectorized_lm_matrix %*% (
    identity_matrix - (vector_centre %*% t(vector_centre))
  )
  xp <- xp + xm
  return(array(t(xp), dim = c(p, k, n)))

}
