

#' Stereographic Projection of Procrustes or Kendall Coordinates
#'
#' The present function performs a stereographic projection of coordinates
#' from non-linear Kendall space or the Procrustes hemisphere into their corresponding
#' position in tangent space
#'
#' @param kendall_coordinates A (num_landmark x num_dimension x num_individuals)
#' array (or tensor) containing kendall or Procrustes coordinates
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
#' @return The non-euclidean coordinates projected stereographically into tanjent space
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
#' GPAshape <- GPA(apes$x, scale = TRUE)
#'
#' GPAshape <- GPA(apes$x, proc_method = "LS")
#'
#' stereographic_projection(GPAshape$procrustes_coordinates)
#'
#' @export

stereographic_projection <- function(kendall_coordinates,
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
    pole <- calc_central_morph(kendall_coordinates, method = "median")
  } else {
    pole <- calc_central_morph(kendall_coordinates, method = "mean")
  }

  stereographic_coordinates <- array(NA, dim = c(p, k, n))

  for (individual in 1:n) {
    rho <- 2 * asin(
      (sqrt(sum((kendall_coordinates[,,individual] - pole)^2))) / 2
    )
    stereographic_coordinates[,,individual] <- kendall_coordinates[,,individual] /
      (cos(rho))
  }

  return(stereographic_coordinates)

}


