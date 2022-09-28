
#' Non-Linear Principal Component Analysis
#'
#' This function provides an extension of \code{\link{pca_plot}} to create a
#' non-linear PCA using a selection of kernels.
#'
#' @param data A \code{data.frame} containing the information to plot
#' @param kernel A string defining the type of kernel to be used;
#' "rbf" (Radial Basis Function), "poly" (Polynomial), "laplace"
#' (Laplace), "gauss" (Gaussian), "cauchy" (Cauchy), "spline" (Spline)
#' @param hyperparam A positive numeric value defining the hyperparameter to be
#' applied to each type of kernel (unavailable for the Spline function, and not
#' recommended for the Cauchy function)
#' @param labels A factorial list containing the sample labels
#' @param axes A two-value vector defining the axes to be plotted (e.g. c(1,2))
#' @param label_colours A string vector of colours to be used in the plot.
#' @param centre A boolean variable specifying whether data is to be centred first
#' @param scale A boolean variable specifying whether data is to be scaled first
#' @param point_size A numeric value defining the siz eof plot points
#' @param CI_ellipse A boolean specifying whether Confidence Ellipses are to be
#' included
#' @param CI_level A numeric value between 0 and 1 defining the % of confidence
#' for CI ellipses
#' @param CI_size A numeric value defining the thikness of CI ellipses lines
#' @param Chull A boolean specifying whether convex hulls are to be included
#' @param Ch_size A numeric value defining the thickness of Convex Hull lines
#' @param main A string specifying the plot title (if included)
#'
#' @section Details:
#' This function plots a Non-Linear Princpial Component Analysis using a
#' prior transformation to the data according to a kernel provided by the user.
#'
#' @return \code{eigenvalues} - Eigenvalues for each PC score
#' @return \code{variance} - Explained variance of each PC Score
#' @return \code{pc_scores} - The PC scores
#' @return \code{pca_plot} - a ggplot object containing the pca plot
#'
#' @author Lloyd A. Courtenay
#'
#' @seealso
#' \code{\link{kernel_cauchy}}, \code{\link{kernel_gaussian}},
#' \code{\link{kernel_laplace}}, \code{\link{kernel_poly}},
#' \code{\link{kernel_rbf}}, \code{\link{kernel_spline}},
#' \code{\link{kernel_pca_biplot}},
#' \code{\link{pca_plot}}, \code{\link{pca_biplot}}
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
#' # plot non-linear pca
#'
#' kernel_pca(data, apes$group, kernel = "rbf", hyperparam = 0.001)
#'
#' @export

kernel_pca <- function(data,
                       kernel,
                       hyperparam = 0,
                       labels = NULL,
                       axes = c(1,2),
                       label_colours = NULL,
                       centre = TRUE, scale = FALSE,
                       point_size = 2, CI_ellipse = FALSE,
                       CI_level = 0.95, CI_size = 1,
                       Chull = FALSE, Ch_size = 1,
                       main = NULL) {

  if (!is.data.frame(data)) {
    stop("Input data must be of type data.frame")
  }

  if (length(axes) != 2) {
    stop("Axes requires a list of 2 pc scores to visualize, e.g. c(1,2) for PC1 and PC2")
  }

  if (axes[1] > ncol(data) |
      axes[2] > ncol(data) |
      axes[1] <= 0 |
      axes[2] <= 0 |
      !is.numeric(axes[1]) |
      !is.numeric(axes[2])) {
    stop("Invalid PCA axes selected.")
  }

  if (!is.logical(centre)) {
    stop("The centre parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.logical(scale)) {
    stop("The scale parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.numeric(point_size) | point_size < 0) {
    stop("point_size must be a positive value")
  }

  if (!is.logical(CI_ellipse)) {
    stop("The CI.ellipse parameter only accepts boolean TRUE or FALSE values")
  }

  if (CI_level >= 1 | CI_level <= 0) {
    stop("Confidence Interval levels must be set between 0 and 1")
  }

  if (!is.numeric(CI_size) | CI_size < 0) {
    stop("CI_size must be a positive value")
  }

  if (!is.null(label_colours)) {
    colour_bool <- check_colours(label_colours)
    if (FALSE %in% colour_bool) {
      stop("Invalid colour provided for label colours")
    }
  }

  if (!is.logical(Chull)) {
    stop("The Chull parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.numeric(Ch_size) | Ch_size < 0) {
    stop("CI_size must be a positive value")
  }

  if (missing(kernel)) {
    stop("No kernel has been provided.")
  }

  possible_kernels <- c("rbf", "poly", "laplace", "gauss", "cauchy", "spline")
  if (kernel %!in% possible_kernels) {
    stop(
      "Invalid kernel provided. Select from 'rbf', 'poly', 'laplace', 'gauss', 'cauchy' or 'spline'."
    )
  }

  if (kernel == "rbf") {
    if (hyperparam == 0) {
      X <- kernel_rbf(data)
    } else {
      X <- kernel_rbf(data, gamma = hyperparam)
    }
  } else if (kernel == "poly") {
    if (hyperparam == 0) {
      X <- kernel_poly(data)
    } else {
      X <- kernel_poly(data, degree = hyperparam)
    }
  } else if (kernel == "laplace") {
    if (hyperparam == 0) {
      X <- kernel_laplace(data)
    } else {
      X <- kernel_laplace(data, alpha = hyperparam)
    }
  } else if (kernel == "gauss") {
    if (hyperparam == 0) {
      X <- kernel_gaussian(data)
    } else {
      X <- kernel_gaussian(data, alpha = hyperparam)
    }
  } else if (kernel == "cauchy") {
    if (hyperparam != 0) {
      warning("For the cauchy kernel, the hyperparameter will be ignored.")
    }
    X <- kernel_cauchy(data)
  } else {
    if (hyperparam != 0) {
      warning("For the spline kernel, the hyperparameter will be ignored.")
    }
    X <- kernel_spline(data)
  }

  kernel_pca_data <- pca_plot(X,
                              labels = labels,
                              axes = axes,
                              label_colours = label_colours,
                              centre = centre,
                              scale = scale,
                              point_size = point_size,
                              CI_ellipse = CI_ellipse,
                              CI_level = CI_level,
                              CI_size = CI_size,
                              Chull = Chull,
                              Ch_size = Ch_size,
                              main = main)

  return(kernel_pca_data)

}
