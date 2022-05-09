
#' Generalised Resistant Fit (GRF)
#'
#' The present function performs a full Generalised Resistant Analysis for the
#' registration of landmark configurations into a superimposed configuration.
#' This function by default performs superimposition using translation, rotation,
#' and scaling (shape analysis), however has also been adapted for analyses without
#' scaling (form analysis).
#'
#' @param lm_array A (num_landmark x num_dimension x num_individuals) array (or tensor) containing
#' raw landmark coordinates
#' @param scale A boolean value indicating whether analyses should be performed in
#' shape (scale = TRUE) or form (scale = FALSE) space.
#' @param robust A boolean value indicating whether a robust statistics should be used.
#' @param tolerance_limit A numeric value defining the tolerance acceptable for
#' algorithm convergence
#' @param projection Either select "orthogonal" or "stereographic" for this parameter.
#' This defines whether Procrustes coordinates are projected into tanjent space
#' using an orthogonal or stereographic projection.
#'
#' @section Details:
#' This function calculates a full Generalised Resistant fit on a set of 2D or 3D landmark
#' coordinates. The final product is a set of projected superimposed coordinates. GRF
#' was developed as a means of reducing the effect a single landmark may have on the
#' entire configuration (Siegel and Benson, 1982; Siegel and Pinkerton, 1982;
#' Rohlf and Slice, 1990; Slice, 1996; Dryden and Walker, 1999). Nevertheless,
#' this method is computationally expensive,
#' and in some cases do not guarantee convergence.
#'
#' @section Note:
#'
#' This is only an experimental function, we recommend using traditional GPA.
#'
#' If this function produces an error, it will be due to a lack of convergence
#' in the algorithm.
#'
#' @section Bibliography (and additional references):
#'
#' Siegel, A.F.; Benson, R.H. (1982) A robust comparison of biological shapes
#' Biometrics. 38(2):341-350.
#'
#' Siegel, A.F.; Pinkerton, J.R. (1982) Robust comparison of three dimensional shapes
#' with an application to protein molecule configurations, Technical Report,
#' Department of Statistics, Princeton University. 224
#'
#' Rohlf, J.F.; Slice, D.E. (1990) Extension of the Procrustes method for the
#' optimal superimposition of landmarks, Systematic Biology. 39:40-59
#'
#' Slice, D.E. (1996) Three-dimensional, generalized resistant fitting and the
#' comparison of least-squares and resistant fit residuals. In: L.F. Marcus,
#' M. Corti, A. Loy, G.J.P. Naylor, D.E. Slice (Eds.) Advances in Morphometrics,
#' New York: Plenum Press. 179-199.
#'
#' Dryden, I.L.; Walker, G. (1999) Highly Resistant Regression and Object Matching,
#' Biometrics, 55:820-825.
#'
#' Claude, J. (2008) Morphometrics with R. The Netherlands: Springer
#'
#' Courtenay, L.A.; Aramendi, J.; González-Aguilera, D. (In Prep) A Graph
#' Based Geometric Morphometric approach to the analysis of primate radii:
#' A new mathematical model for the processing of landmark data.
#'
#' @return The superimposed Procrustes coordinates
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' library(shapes)
#' data(apes)
#'
#' # for issues with computational time we set the tolerance to 1
#'
#' GRFshape <- generalised_resistant_fit(apes$x,
#'                                       tolerance_limit = 1)
#' pca_plot(vector_from_landmarks(GRFshape), apes$group)
#'
#' @export

generalised_resistant_fit <- function(lm_array, scale = TRUE, robust = FALSE,
                                  tolerance_limit = 0.005,
                                  projection = c("orthogonal", "stereographic")) {

  if (!is.logical(scale)) {
    stop("The scale parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.logical(robust)) {
    stop("The robust parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.array(lm_array) | is.na(dim(lm_array)[3])) {
    stop("Landmark configuration must be a (landmark x dimension x individual) array")
  }

  if(!is.numeric(tolerance_limit) | tolerance_limit <= 0) {
    stop("tolerance_limit must be a positive numeric value")
  }

  projection <- match.arg(projection)

  p <- dim(lm_array)[1]
  k <- dim(lm_array)[2]
  n <- dim(lm_array)[3]

  preliminary_coordinates <- GPA(
    lm_array, scale = scale, robust = robust,
    projection = projection
  )$coordinates

  if (any(is.na(preliminary_coordinates))) {
    stop("NAs were produced when calculating stereographic projections, consider using orthogonal projections")
  }

  D <- transformed_array <- array(NA, dim = c(p, k, n))

  for (individual in 1:n) {
    transformed_array[,, individual] <- preliminary_coordinates[,, individual] /
      resistant_median_size(preliminary_coordinates[,,individual])
  }

  reference_matrix <- apply(transformed_array, 1:2, median)
  reference_matrix <- reference_matrix / resistant_median_size(reference_matrix)
  tolerance <- 10

  if (k == 2) {

    while (tolerance > tolerance_limit) {

      for (individual in 1:n) {
        pretransformed_array <- as.matrix(
          dist(reference_matrix)) / as.matrix(dist(transformed_array[,,individual])
          )
        beta <- median(apply(pretransformed_array, 2, median, na.rm = TRUE))
        angles_array <- resistant_rotation_matrix(
          transformed_array[,,individual], reference_matrix
        )
        phi <- median(apply(angles_array, 2, median, na.rm = TRUE))
        Gamma <- matrix(
          c(cos(phi), -sin(phi), sin(phi), cos(phi)), 2, 2
        )
        alpha <- reference_matrix - beta * transformed_array[,,individual] %*% Gamma
        D[,,individual] <- beta * transformed_array[,,individual] %*%
          Gamma + matrix(
            apply(alpha, 2, median), p, k, byrow = TRUE
          )
        new_reference <- apply(D, 1:2, median)
      }

      tolerance <- median(
        sqrt(apply((new_reference - reference_matrix)^2, 1, sum))
      )
      reference_matrix <- new_reference <- new_reference / resistant_median_size(new_reference)
      transformed_array <- D

    }

  } else if (k == 3) {

    while (tolerance > tolerance_limit) {

      for (individual in 1:n) {
        pretransformed_array <- as.matrix(
          dist(reference_matrix)) / as.matrix(dist(transformed_array[,,individual])
          )
        beta <- median(apply(pretransformed_array, 2, median, na.rm = TRUE))
        Gamma <- resistant_3D_rotation(
          transformed_array[,,individual], reference_matrix
        )$Gamma
        alpha <- reference_matrix - beta * transformed_array[,,individual] %*% Gamma
        D[,,individual] <- beta * transformed_array[,,individual] %*%
          Gamma + matrix(
            apply(alpha, 2, median), p, k, byrow = TRUE
          )
        new_reference <- apply(D, 1:2, median)
      }

      tolerance <- median(
        sqrt(apply((new_reference - reference_matrix)^2, 1, sum))
      )
      reference_matrix <- new_reference <- new_reference / resistant_median_size(new_reference)
      transformed_array <- D

    }

  }

  return(transformed_array)

}
