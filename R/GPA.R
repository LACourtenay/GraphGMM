
#' Generalised Procrustes Analysis (GPA)
#'
#' The present function performs a full Generalised Procrustes Analysis for the
#' registration of landmark configurations into a superimposed configuration.
#' This function by default performs superimposition using translation, rotation,
#' and scaling (shape analysis), however has also been adapted for analyses without
#' scaling (form analysis). GPA is performed using a Least Square (LS) approach.
#'
#' @param lm_array A (num_landmark x num_dimension x num_individuals) array (or tensor) containing
#' raw landmark coordinates
#' @param scale A boolean value indicating whether analyses should be performed in
#' shape (scale = TRUE) or form (scale = FALSE) space.
#' @param robust A boolean value indicating whether a Generalised Robust Fit (GRF)
#' should be performed.
#' @param proc_method Either select "optimalLS" or "LS" for this parameter. See "Details"
#' for information on the difference between the two approaches
#' @param projection Either select "orthogonal" or "stereographic" for this parameter.
#' This defines whether Procrustes coordinates are projected into tanjent space
#' using an orthogonal or stereographic projection. This is only used if the proc_method
#' is set to "LS".
#'
#' @section Details:
#' This function calculates a full Generalised Procrustes fit on a set of 2D or 3D landmark
#' coordinates. The final product is either a set of superimposed coordinates on the Procrustes
#' (hyper)hemisphere (Slice, 2001), or their projection into tanject space
#' according to either an orthogonal (default) or a stereographic projection. If proc_method = "optimalLS",
#' then the coordinates are already provided as projected coordinates.
#'
#' The Least-Square algorithm employed in the present study has been implemented using
#' two different techniques. The vanilla "LS" approach employs the classic superimposition
#' techniques described by Gower (1975), while "optimalLS" employs the use of a more
#' powerful superimposition algorithm that ensures are more optimal convergence
#' (Rohlf and Slice, 1990). The "optimalLS" approach also includes some computational
#' optimisations used by Ian Dryden in the "shapes" library (Dryden and Mardia, 2016).
#'
#' The robust parameter is used for what Courtenay et al. (In Prep) describe as a
#' Generalised Robust Fit (GRoF), which uses robsut statistical measures (Höhle and Höhle, 2009)
#' for the calculation of both the centroid and the central configuration
#' (used as the reference for LS optimisation).
#'
#' @section Notes:
#'
#' At present GRoF is in early stages of development, and the GPA function at present
#' works best when using robust = FALSE and proc_method = "optimalLS". All other parameters
#' with the exception of scale work similiar to the \code{procGPA} function of the
#' \code{shapes} library.
#'
#' Optimal results are therefore obtained using:
#'
#' GPA(x, scale = TRUE, robust = FALSE, proc_method = "optimalLS")
#'
#' or
#'
#' GPA(x, scale = FALSE, robust = FALSE, proc_method = "optimalLS")
#'
#' The GRoF function (x, robust = TRUE) is therefore only meant to be used
#' for experimental purposes (in this version of GraphGMM)
#'
#' If users wish to use a more stable version of GPA, we strongly recommend the
#' \code{procGPA} function of the \code{shapes} library
#'
#' @section Bibliography (and additional references):
#'
#' Gower, J.C. (1975) Generalized Procrustes Analysis, Psychometrika. 40:33-50
#'
#' Kendall, D.G. (1984) Shape maniforlds, procrustean metrics, and complex projective
#' spaces, Bulletin of the London Mathematical Society. 16:81-121
#'
#' Bookstein, F.L. (1986) Size and shape spaces for landmark data in two dimensions,
#' Statistical Science. 4(2):181-242
#'
#' Rohlf, J.F.; Slice, D.E. (1990) Extension of the Procrustes method for the
#' optimal superimposition of landmarks, Systematic Biology. 39:40-59
#'
#' Goodall, C. (1991) Procrustes Methods in the Statistical Analysis of Shape,
#' Journal of the Royal Statistical Society B. 53(2):285-339
#'
#' Slice, D. (2001) Landmark coordinates aligned by Procrustes Analysis do not
#' lie in Kendall’s shape space, Systematic Biology. 50(1):141-149.
#' DOI: 10.1080/10635150119110
#'
#' Claude, J. (2008) Morphometrics with R. The Netherlands: Springer
#'
#' Höhle, J.; Höhle, M. (2009) Accuracy assessment of digital evelation models
#' by means of robust statistical methods, ISPRS Journal of Photogrammetry and
#' Remote Sensing. 64:398-406
#'
#' Dryden, I.L.; Mardia, K.V. (2016) Statistical Shape Analysis: With Applications
#' in R. West Sussex: Wiley
#'
#' Courtenay, L.A.; Aramendi, J.; González-Aguilera, D. (In Prep) A Graph
#' Based Geometric Morphometric approach to the analysis of primate radii:
#' A new mathematical model for the processing of landmark data.
#'
#' @return \code{coordinates} - The Superimposed Coordinates
#' @return \code{centroid_sizes} - A list of centroid sizes for each individual
#'
#' @author Lloyd A. Courtenay (with elements inspired by Ian Dryden's "shapes" package,
#' and Claude, 2008)
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
#' pca_plot(vector_from_landmarks(GPAshape$coordinates), apes$group)
#'
#' # form analysis
#'
#' GPAform <- GPA(apes$x, scale = FALSE)
#'
#' pca_plot(vector_from_landmarks(GPAform$coordinates), apes$group)
#'
#' @export

GPA <- function(lm_array,
                scale = TRUE,
                robust = FALSE,
                proc_method = c("optimalLS", "LS"),
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

  # if no proc_method or projection is defined, then select optimalLS, otherwise use the
  # method defined by the user

  proc_method <- match.arg(proc_method)
  projection <- match.arg(projection)

  # match.args should catch errors

  p <- dim(lm_array)[1]
  k <- dim(lm_array)[2]
  n <- dim(lm_array)[3]

  if (proc_method == "optimalLS") {

    if (robust == FALSE) {
      lm_array <- dryden_pretransform(lm_array)
    } else {
      for (individual in 1:n){
        lm_array[,,individual] <- kendall_translation(lm_array[,,individual], robust = TRUE)
      }
    }

    transformed_array <- dryden_rotation(lm_array, robust = robust)

    # i think this can be performed outside the proc_method conditional statement

    centroid_sizes <- c()
    for (individual in 1:n) {
      centroid_sizes <- c(centroid_sizes, centroid_size(
        lm_array[,,individual],
        find_centroid(lm_array[,,individual],
                      robust = robust)
      ))
    }

    if (robust == TRUE) {
      rho_coefficients <- apply(lm_array, 3, function(x) {
        riemannianD(x, calc_central_morph(lm_array, method = "median"), robust = TRUE)
      })
    } else {
      rho_coefficients <- apply(lm_array, 3, function(x) {
        riemannianD(x, calc_central_morph(lm_array, method = "mean"), robust = FALSE)
      })
    }

    if (scale == TRUE) {

      x1 <- dryden_dif(lm_array,
                       robust = robust)

      x2 <- dryden_dif(transformed_array$rotated,
                       robust = robust)

      if (robust == FALSE) {
        transformed_array <- dryden_rotation(
          dryden_scale(transformed_array$rotated)
        )
      } else {

        intermediate_array <- transformed_array$rotated

        for (individual in 1:n) {
          cs_ind <- centroid_size(
            intermediate_array[,,individual],
            find_centroid(intermediate_array[,,individual],
                          robust = TRUE)
          )
          intermediate_array[,,individual] <- intermediate_array[,,individual] / cs_ind
        }
        transformed_array <- dryden_rotation(intermediate_array)

      }

      x1 <- x2
      x2 <- dryden_dif(transformed_array$rotated,
                       robust = robust)

      rho <- x1 - x2

      dryden_i = 1

      if (rho > 1e-05) {
        while (rho > 1e-05) {

          x1 <- x2
          dryden_i = dryden_i + 1

          if (robust == FALSE) {
            transformed_array <- dryden_rotation(
              dryden_scale(transformed_array$rotated)
            )
          } else {

            intermediate_array <- transformed_array$rotated

            for (individual in 1:n) {
              cs_ind <- centroid_size(
                intermediate_array[,,individual],
                find_centroid(intermediate_array[,,individual],
                              robust = TRUE)
              )
              intermediate_array[,,individual] <- intermediate_array[,,individual] / cs_ind
            }
            transformed_array <- dryden_rotation(intermediate_array)

          }

          x2 <- dryden_dif(transformed_array$rotated,
                           robust = robust)
          rho <- x1 - x2

        }
      }

      if (robust == TRUE) {
        central_config <- calc_central_morph(transformed_array$rotated, method = "median")
      } else {
        central_config <- calc_central_morph(transformed_array$rotated, method = "mean")
      }

      tan_partial <- matrix(0, p * k - k, n)
      identity_matrix <- diag(rep(1, (k * p - k)))
      Gamma <- as.vector(kendall_preshape(central_config))

      for (individual in 1:n) {
        tan_partial[, individual] <- (identity_matrix - Gamma %*% t(Gamma)) %*%
          as.vector(kendall_preshape(transformed_array$rotated[,,individual]))
      }

      tan <- transformed_array$rotated[, 1, ] - central_config[, 1]
      for (dimension in 2:k) {
        tan <- rbind(tan,
                     transformed_array$rotated[, dimension, ] - central_config[, dimension])
      }

    } else {

      x1 <- transformed_array$rotated
      if (robust == TRUE) {
        central_config <- calc_central_morph(transformed_array$rotated, method = "median")
      } else {
        central_config <- calc_central_morph(transformed_array$rotated, method = "mean")
      }

      tan_partial <- matrix(0, p * k - k, n)
      identity_matrix <- diag(rep(1, (k * p - k)))
      Gamma <- as.vector(kendall_preshape(central_config))

      for (individual in 1:n) {
        tan_partial[, individual] <- (identity_matrix - Gamma %*% t(Gamma)) %*%
          as.vector(kendall_preshape(transformed_array$rotated[,,individual]))
      }

      tan <- transformed_array$rotated[, 1, ] - central_config[, 1]
      for (dimension in 2:k) {
        tan <- rbind(tan,
                     transformed_array$rotated[, dimension, ] - central_config[, dimension])
      }

    }

    return_coordinates <- list(
      coordinates = transformed_array$rotated,
      centroid_sizes = centroid_sizes
    )

  } else if (proc_method == "LS") {

    scaled_configs <- rotated_configs <- array(NA, dim = c(p, k, n))
    centroid_sizes <- numeric(n)

    for (individual in 1:n) {

      centroid <- find_centroid(lm_array[,,individual],
                                robust = robust)
      CS <- centroid_size(lm_array[,,individual],
                          centroid)

      if (scale == TRUE) {
        scaled_configs[,,individual] <- kendall_translation(
          lm_array[,,individual] / CS,
          robust = robust
        )
      } else {
        scaled_configs[,,individual] <- kendall_translation(
          lm_array[,,individual],
          robust = robust
        )
      }

      centroid_sizes[individual] <- CS

    }

    original_loss_dist_matrix <- dist(t(matrix(scaled_configs, k * p, n)))
    loss <- sum(original_loss_dist_matrix)

    while(abs(loss) > 1e-04) {
      for(individual in 1:n) {

        if (robust == TRUE) {
          central_config <- calc_central_morph(scaled_configs[,,-individual], method = "median")
        } else {
          central_config <- calc_central_morph(scaled_configs[,,-individual], method = "mean")
        }

        if (scale == TRUE) {

          rotated_configs[,,individual] <- kendall_full_superimposition(
            scaled_configs[,,individual], central_config, robust = robust
          )$config_1_prima

        } else {

          rotated_configs[,,individual] <- kendall_full_superimposition(
            scaled_configs[,,individual], central_config, scale = FALSE,
            robust = robust
          )$config_1_prima

        }
      }

      epoch_loss_dist_matrix <- dist(t(matrix(scaled_configs, k * p, n)))
      loss <- sum(original_loss_dist_matrix) - sum(epoch_loss_dist_matrix)
      original_loss_dist_matrix <- epoch_loss_dist_matrix
      scaled_configs <- rotated_configs

    }

    kendall_coordinates <- rotated_configs

    if (projection == "orthogonal") {
      projected_coordinates <- orthogonal_projection(kendall_coordinates, robust = robust)
    } else {
      projected_coordinates <- suppressWarnings(
        stereographic_projection(kendall_coordinates, robust = robust)
      )

      if (any(is.na(projected_coordinates))) {
        warning("Stereographic projections produced NA values, consider using Orthogonal projections")
      }

    }

    return_coordinates <- list(
      procrustes_coordinates = kendall_coordinates,
      projected_coordinates = projected_coordinates,
      centroid_sizes = centroid_sizes
    )

  }

  #interlandmark_distance <- function(config_1, config_2) {
  #  return(sqrt(
  #    apply((config_1 - config_2)^2, 1, sum)
  #  ))
  #}
  #
  #pd <- asin(dist(t(matrix(kendall_coordinates, k * p, n))))
  #projected_coordinates <- orthogonal_projection(kendall_coordinates, robust = robust)
  #tau <- matrix(NA, n, n); for (individual_i in 1:n) {
  #  for (individual_j in 1:n) {
  #    tau[individual_i, individual_j] <- sqrt(
  #      sum(
  #        interlandmark_distance(
  #          projected_coordinates[,,individual_i],
  #          projected_coordinates[,,individual_j]
  #        )^2
  #      )
  #    )
  #  }
  #}
  #tau <- as.dist(t(tau))
  #
  #if ((cor(pd, tau))^2 > 0.9) {
  #  projected_coordinates <- tangent_projection(kendall_coordinates, robust = robus)
  #}

  return(return_coordinates)

}





