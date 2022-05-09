
#' Principal Component Analysis
#'
#' This function calculates and plots a PCA.
#'
#' @param data A \code{data.frame} containing the information to plot
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
#' @return \code{variance} - Explained variance of each PC Score
#' @return \code{pc_scores} - The PC scores
#' @return \code{pca_plot} - a ggplot object containing the pca plot
#'
#' @seealso
#' \code{\link{kernel_pca}},
#' \code{\link{pca_biplot}},
#' \code{\link{kernel_pca_biplot}}
#' 
#' @author Lloyd A. Courtenay
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
#' # plot pca
#'
#' pca_plot(data, apes$group, CI_ellipse = TRUE)
#'
#' @export

pca_plot <- function(data,
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

  #if(!is.null(CI_ellipse) & !is.null(Chull)) {
  #  stop("Please select either CI_ellipse or Chull, unable to plot both.")
  #}

  if (!is.numeric(Ch_size) | Ch_size < 0) {
    stop("CI_size must be a positive value")
  }

  #PC <- prcomp(as.matrix(data), centre = centre, scale = scale)
  PC <- single_value_decomposition(as.matrix(data), centre = centre, scale = scale)

  expl_var <- PC$sdev^2 / sum(PC$sdev^2)

  if (!is.null(labels)) {

    pca_data <- data.frame(
      x = PC$pc_scores[,axes[1]], y = PC$pc_scores[,axes[2]],
      Sample = as.factor(labels)
    )
    base_plot <- ggplot2::ggplot(data = pca_data,
                                 ggplot2::aes(x = x, y = y, colour = Sample))

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

  }

  if (is.null(labels)) {
    pca_data <- data.frame(
      x = PC$pc_scores[,axes[1]], y = PC$pc_scores[,axes[2]]
    )
    base_plot <- ggplot2::ggplot(data = pca_data,
                                 ggplot2::aes(x = x, y = y))
    plot_colours <- ggplot2::scale_color_manual(values = c("black"))
  }

  pca_plot <- base_plot +
    ggplot2::geom_point(stat = "identity", size = point_size) +
    plot_colours

  if (CI_ellipse == TRUE) {
    conf_interval <- ggplot2::stat_ellipse(size = CI_size,
                                           level = CI_level)
    pca_plot <- pca_plot + conf_interval
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
    pca_plot <- pca_plot + conf_interval
  }

  pca_plot <- pca_plot +
    ggplot2::xlab(paste("PC", axes[1], " (",
                        round(expl_var[axes[1]] * 100,2),"%)", sep = "")) +
    ggplot2::ylab(paste("PC", axes[2], " (",
                        round(expl_var[axes[2]] * 100,2),"%)", sep = "")) +
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
    pca_plot <- pca_plot + ggplot2::ggtitle(main)
  }

  return(list(
    variance = expl_var,
    pc_scores = PC$pc_scores,
    pca_plot = pca_plot
  ))

}
