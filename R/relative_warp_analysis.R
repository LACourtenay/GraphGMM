
#' Relative Warp Analysis (RWA)
#'
#' This function is used to calculate the bending energy derived from multiple Thin
#' Plate Spline interpolation functions. This bending energy can then be used to describe
#' patterns in affine and non-affine variation, while using an alpha value to emphasize
#' both local and global morphological changes. The eigendecomposition of Partial
#' Warp Scores is known as Relative Warp Analysis.
#'
#' @param proc_coords A (num_landmark x num_dimension x num_individuals) array (or tensor) containing
#' raw landmark coordinates
#' @param labels A factorial list containing the sample labels for RWA plot
#' @param robust A boolean value indicating whether robust statistics should be used.
#' @param alpha A non-zero numeric value defining the scaling effect the bending matrix
#' has on RWA results. alpha = 1 will highlight large scale differences (global), while alpha = -1 will
#' highlight small scale variations (local).
#' @param projection Whether a projection of Procrustes coordinates to tangent space is needed
#' @param proj_scale Whether an additional scaling procedure is need in the calculation of
#' the pole (central configuration) during projection onto tangent space
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
#'
#' Thin Plate Spline (TPS) interpolation functions are used to describe a smooth function
#' from one specimen to another. This function can be used to plot visual deformations
#' between two different shapes, or used to calculate multiple components of morphological
#' variation. When TPS functions are decomposed into their seperate components,
#' they can be used to describe a series of affine (uniform) and non-affine (not uniform) deformations
#' between the two shapes.
#'
#' Affine transformations can include translation, scaling, and shearing
#' of the overall configuration, without being localised to any particular region of the
#' configuration (subset of landmarks). Non-affine transformations, on the other hand,
#' refer to the expansion, compression or bending of local regions.
#'
#' Affine transformations are described by Uniform shape components.
#'
#' Non-affine transformations are described by principal warps
#'
#' Relative Warp Scores are the eigendecomposition of the weighted principal warps
#'
#' @section Selection of alpha value:
#'
#' Please do not select an alpha value of 0. Scaling partial warp scores to 0
#' has been found to produce almost identical results to a standard PCA on
#' coordinate data. See Rohlf (1993, 1996, 1999) for more details.
#'
#' @section Additional Notes:
#'
#' Many R libraries perform RWA using their own configuration of GPA. The present
#' library allows the user to predefine the superimposition procedure they wish to
#' use. Thus the user must provide the relative_warp_analysis function with
#' coordinates already superimposed in their desired feature space.
#'
#' @section Bibliography (and additional references):
#'
#' Bookstein, F.L. (1989) Principal warps: thin-plate splines and the decomposition
#' of deformations. IEEE Transactions on Pattern Analysis and Machine Intelligence. 11:567-585
#'
#' Bookstein, F.L. (1991)  Morphometric tools for landmark data: geometry and
#' biology. Cambridge: Cambridge University Press.
#'
#' Rohlf, F.J. (1993) Relative warp analysis and an example of its application to
#' mosquito wings. In: L.F. Marcus, E. Bello, A. García-Valdecasas (Eds.)
#' Contributions to Morphometrics. Madrid: Monografías del Museo Nacional de
#' Ciencias Naturales. 8:131-159
#'
#' Rohlf, F.J. (1996) Morphometric spaces, shape components, and the effects of
#' linear transformations. In: L.F. Marcus, M. Corti, A. Loy, G.J.P. Naylor, D.E.
#' Slice (Eds.) Advances in Morphometrics. The Netherlands: Springer. 117-129
#'
#' Rohlf, F.J., Loy, A., Corti, M. (1996) Morphometric analysis of old world
#' talpidae (Mammalia, Insectivora) using partial-warp scores, Systematic Biology.
#' 45(3):344-362
#'
#' Rohlf, F.J. (1999) Shape statistics: Procrustes superimposition and tangent spaces,
#' Journal of Classification. 16:197-223
#'
#' Rohlf, F.J.; Bookstein, F.L. (2003) Computing the uniform component of shape
#' variation. Systematic Biology. 52(1):66-69
#'
#' @return \code{relative_warp_scores} - PC Scores obtained from Relative Warp Analysis
#' @return \code{nonaffine_variation_exp} - The percentage of non-affine
#' variation explained by each relative warp score.
#' @return \code{nonaffine_variation} - descriptors of nonaffine variation (the Partial Warp Scores)
#' @return \code{uniform_scores} - Uniform Scores
#' @return \code{affine_variation} - descriptors of affine variation (the Uniform Component)
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
#' relative_warp_analysis(GPAshape$coordinates, apes$group,
#'                        alpha = -1)
#' relative_warp_analysis(GPAshape$coordinates, apes$group,
#'                        alpha = 1)
#'
#' @export

relative_warp_analysis <- function(proc_coords,
                                   labels = NULL,
                                   robust = FALSE,
                                   alpha = 1,
                                   projection = FALSE,
                                   proj_scale = FALSE,
                                   axes = c(1,2),
                                   label_colours = NULL,
                                   centre = TRUE, scale = FALSE,
                                   point_size = 2, CI_ellipse = FALSE,
                                   CI_level = 0.95, CI_size = 1,
                                   Chull = FALSE, Ch_size = 1,
                                   main = NULL) {

  if (!is.array(proc_coords) | is.na(dim(proc_coords)[3])) {
    stop("proc_coords must be an unprojected (landmark x dimension x individual) array of Procrustes coordinates")
  }

  if (!is.logical(robust)) {
    stop("The robust parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.logical(proj_scale)) {
    stop("The proj_scale parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.logical(projection)) {
    stop("The projection parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.numeric(alpha)) {
    stop("The alpha value should be a numeric value")
  }

  if (alpha == 0) {
    stop("alpha should not be set to 0, see help(relative_warp_analysis)")
  }

  rwa_data <- RWA(proc_coords, robust = robust, alpha = alpha, projection = projection,
                  proj_scale = proj_scale)

  if (length(axes) != 2) {
    stop("Axes requires a list of 2 pc scores to visualize, e.g. c(1,2) for PC1 and PC2")
  }

  if (axes[1] > ncol(rwa_data$relative_warp_scores) |
      axes[2] > ncol(rwa_data$relative_warp_scores) |
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

  if(!is.null(labels)) {

    pca_data <- data.frame(
      x = rwa_data$relative_warp_scores[,axes[1]],
      y = rwa_data$relative_warp_scores[,axes[2]],
      Sample = as.factor(labels)
    )

    base_plot <- ggplot2::ggplot(data = pca_data,
                                 ggplot2::aes(x = x, y = y,
                                              colour = Sample))

    if (!is.null(label_colours)) {

      if (length(label_colours) < length(levels(as.factor(labels)))) {
        stop("Insufficient label_colours arguments provided.")
      }

      plot_colours <- ggplot2::scale_color_manual(values = label_colours)

    } else {
      plot_colours <- ggplot2::scale_color_manual(values = c(
        "black","red","blue","orange","darkgreen","violetred"
      ))
    }

  } else {
    pca_data <- data.frame(x = rwa_data$relative_warp_scores[,axes[1]],
                           y = rwa_data$relative_warp_scores[,axes[2]])
    base_plot <- ggplot2::ggplot(data = pca_data,
                                 ggplot2::aes(x = x, y = y))
    plot_colours <- ggplot2::scale_color_manual(values = c("black"))

  }

  rwa_plot <- base_plot +
    ggplot2::geom_point(stat = "identity", size = point_size) +
    plot_colours

  if (CI_ellipse == TRUE) {
    conf_interval <- ggplot2::stat_ellipse(size = CI_size,
                                           level = CI_level)
    rwa_plot <- rwa_plot + conf_interval
  }

  if (Chull == TRUE) {
    hulls <- pca_data %>%
      dplyr::group_by(Sample) %>%
      dplyr::slice(chull(x, y))
    conf_interval <- ggplot2::geom_polygon(
      data = hulls,
      alpha = 0,
      size = Ch_size
    )
    rwa_plot <- rwa_plot + conf_interval
  }

  rwa_plot <- rwa_plot +
    ggplot2::xlab(paste("PC", axes[1], " (",
                        round(
                          rwa_data$nonaffine_variation_exp[
                            axes[1], 2
                          ] * 100, 2
                        ),"%)", sep = "")) +
    ggplot2::ylab(paste("PC", axes[2], " (",
                        round(
                          rwa_data$nonaffine_variation_exp[
                            axes[2], 2
                          ] * 100, 2
                        ),"%)", sep = "")) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 18, face = "bold"),
      legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
      legend.box.background = ggplot2::element_rect(colour = "black"),
      legend.position = "bottom"
    ) +
    ggplot2::geom_vline(xintercept = 0,
                        colour = "black",
                        size = 0.5,
                        linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0,
                        colour = "black",
                        linetype = "dashed",
                        size = 0.5)

  if (!is.null(main)) {
    rwa_plot <- rwa_plot + ggplot2::ggtitle(main)
  }

  rwa_data$rwa_plot <- rwa_plot

  return(rwa_data)

}

RWA <- function(proc_coords, robust = FALSE, alpha = 1, projection = FALSE,
                proj_scale = FALSE) {

  uniform_scores <- affine_variation <- nonaffine_variation <- NULL
  p <- dim(proc_coords)[2]
  k <- dim(proc_coords)[1]

  if (projection == TRUE) {
    if (proj_scale == TRUE) {
      orthogonal_coordinates <- orthogonal_projection(proc_coords,
                                                      robust = robust,
                                                      scale = TRUE)
    } else {
      orthogonal_coordinates <- orthogonal_projection(proc_coords,
                                                      robust = robust,
                                                      scale = FALSE)
    }
    orthogonal_vector <- vector_from_landmarks(orthogonal_coordinates)
  } else {
    orthogonal_vector <- vector_from_landmarks(proc_coords)
  }



  if (robust == TRUE) {
    central_config <- calc_central_morph(proc_coords, method = "median")
  } else {
    central_config <- calc_central_morph(proc_coords, method = "mean")
  }

  # Be = Bending Energy Matrix

  # I cannot yet work out how to efficiently calculate the submatrix of L, therefore
  # i have used the Morpho libraries CreateL function for this purpose

  Be <- Morpho::CreateL(central_config, output = "Lsubk")$Lsubk
  scale_vecs <- scale(orthogonal_vector, scale = FALSE)

  eigen_Be <- svd(Be)
  eigen_Be$d <- Re(eigen_Be$d)
  eigen_Be$v <- Re(eigen_Be$v)

  tolerance <- which(eigen_Be$d < 1e-10)

  diagonal_inverse <- diagonal_Be <- eigen_Be$d * 0
  diagonal_inverse[-tolerance] <- eigen_Be$d[-tolerance]^(alpha/2)
  diagonal_Be[-tolerance] <- eigen_Be$d[-tolerance]^(-alpha/2)
  IM <- diag(rep(1, p))

  Be2 <- IM %x% (eigen_Be$v %*% diag(diagonal_Be) %*% t(eigen_Be$v))

  covariance_structure <- Be2 %*% t(scale_vecs)
  eigen_covariance_structure <- svd(covariance_structure)
  eigen_covariance_structure$d <- (
    eigen_covariance_structure$d / sqrt(nrow(scale_vecs) - 1)
  )^2

  tolerance <- which(eigen_covariance_structure$d > 1e-10)

  rel_warp_scores <- as.matrix(
    t(
      t(eigen_covariance_structure$u[,tolerance]) %*% Be2 %*% t(scale_vecs)
    )
  )[,tolerance]

  nonaffine_variation <- IM %x% eigen_Be$v
  nonaffine_variation <- as.matrix(
    (nonaffine_variation) %*% diag(rep(diagonal_inverse, p)) %*%
      t(nonaffine_variation) %*% eigen_covariance_structure$u[,tolerance]
  )

  E <- eigen_Be$v[, -which(eigen_Be$d < 1e-10)]
  Lambda <- diag(rep(1, k))
  L <- Lambda - E %*% solve(crossprod(E)) %*% t(E)

  eigen_bend <- svd(scale_vecs %*% (IM %x% L))
  bend_v <- min(ncol(eigen_bend$v), (p + 0.5 * p * (p-1)-1))

  LS <- eigen_bend$u %*% diag(eigen_bend$d)
  useLS <- min(ncol(LS),(p + 0.5 * p * (p-1)-1))

  uniform_scores <- LS[,1:useLS]
  affine_variation <- eigen_bend$v[,1:bend_v]

  eig_cs <- eigen_covariance_structure$d[eigen_covariance_structure$d > 1e-10]
  eig_cs_sum <- sum(eig_cs)
  eig_cs_var <- eig_cs / eig_cs_sum
  eig_cs_cumsum <- cumsum(eig_cs_var)

  nonaffine_variation_warpscore <- data.frame(
    eigenvalues = eig_cs,
    explained_variance = eig_cs_var,
    cumulative_variance = eig_cs_cumsum
  )

  output <- list(
    relative_warp_scores = rel_warp_scores,
    nonaffine_variation_exp = nonaffine_variation_warpscore,
    nonaffine_variation = nonaffine_variation,
    uniform_scores = uniform_scores,
    affine_variation = affine_variation
  )

  return(output)

}
