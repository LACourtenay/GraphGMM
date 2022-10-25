
#' Principal Component Analysis Biplot
#'
#' This function calculates and plots a PCA biplot.
#'
#' @param data A \code{data.frame} containing the information to plot
#' @param labels A factorial list containing the sample labels
#' @param axes A two-value vector defining the axes to be plotted (e.g. c(1,2))
#' @param label_colours A string vector of colours to be used in the plot.
#' @param centre A boolean variable specifying whether data is to be centred first
#' @param scale A boolean variable specifying whether data is to be scaled first
#' @param n_variables An integer defining the number of variables to be included
#' in the biplot
#' @param arrow_size An integer defining the size of the arrows in the biplot
#' @param arrow_colour A string defining the colour of biplot arrows
#' @param point_size A numeric value defining the siz eof plot points
#' @param CI_ellipse A boolean specifying whether Confidence Ellipses are to be
#' included
#' @param CI_level A numeric value between 0 and 1 defining the % of confidence
#' for CI ellipses
#' @param CI_size A numeric value defining the thikness of CI ellipses lines
#' @param Chull A boolean specifying whether convex hulls are to be included
#' @param Ch_size A numeric value defining the thickness of Convex Hull lines
#' @param legend_position A string defining the position of the legend
#' @param legend_inset A numeric value defining the shift in position of the legend
#' @param main A string specifying the plot title (if included)
#'
#' @return PCA biplot
#' @return \code{variable_contribution} - a list of values defining the weight
#' each variable has on the biplot.
#'
#' @seealso
#' \code{\link{kernel_pca}},
#' \code{\link{pca_plot}},
#' \code{\link{kernel_pca_biplot}}
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
#' data <- vector_from_landmarks(GPAshape$rotated)
#'
#' # plot pca
#'
#' pca_biplot(data, apes$group, arrow_size = 2)
#'
#' @export

pca_biplot <- function(data,
                       labels = NULL,
                       axes = c(1,2),
                       label_colours = NULL,
                       centre = TRUE, scale = FALSE,
                       n_variables = 5, arrow_size = 1.5,
                       arrow_colour = "#CD5C5C",
                       point_size = 1, CI_ellipse = FALSE,
                       CI_level = 0.95, CI_size = 1,
                       Chull = FALSE, Ch_size = 1,
                       legend_position = "topright", legend_inset = 0.02,
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

  if (!is.numeric(Ch_size) | Ch_size <= 0) {
    stop("CI_size must be a positive value")
  }

  possible_legend_positions <- c(
    "topleft", "top", "topright", "left", "centre", "right",
    "bottomleft", "bottom", "bottomright"
  )

  if (legend_position %!in% possible_legend_positions) {
    stop("Incorrect string value for legend position")
  }

  if (!is.numeric(arrow_size) | arrow_size <= 0) {
    stop("Arrow size must be a positive numeric value")
  }

  if (arrow_colour != "#CD5C5C") {
    colour_bool <- check_colours(arrow_colour)
    if (FALSE %in% colour_bool) {
      stop("Invalid colour provided for arrow colours")
    }
  }

  PC <- single_value_decomposition(as.matrix(data), centre = centre, scale = scale)

  expl_var <- PC$sdev^2 / sum(PC$sdev^2)

  contribution_computation <- function(var_cos2, comp_cos2) {var_cos2 * 100 / comp_cos2}

  comp_cos2 <- apply(PC$biplot_var_scores^2, 2, sum)
  variable_contribution <- t(
    apply(
      PC$biplot_var_scores^2, 1, contribution_computation, comp_cos2
    )
  )

  rownames(variable_contribution) <- colnames(data)

  order_variable_coords <- PC$biplot_var_scores
  rownames(order_variable_coords) <- colnames(data)

  order_variable_coords <- order_variable_coords[
    order(variable_contribution[,axes[1]], variable_contribution[,axes[2]],
          decreasing = TRUE),
  ]

  variable_importance_table <- cbind(order_variable_coords[,axes[1]],
                                     order_variable_coords[,axes[2]])
  rownames(variable_importance_table) <- rownames(order_variable_coords)

  arrow_values <- order_variable_coords[1:n_variables, c(axes[1], axes[2])]
  colnames(arrow_values) <- c("x", "y")
  origin <- data.frame(var_name = rownames(order_variable_coords)[1:n_variables]) # error here
  arrow_values <- cbind(arrow_values, origin)

  n <- nrow(PC$biplot_pc_scores)
  p <- nrow(PC$biplot_var_scores)

  var_labs <- rownames(arrow_values)

  unsigned.range <- function(x) {c(-abs(min(x, na.rm = TRUE)),
                                   abs(max(x, na.rm = TRUE)))}
  rangx1 <- unsigned.range(PC$biplot_pc_scores[, axes[1]])
  rangx2 <- unsigned.range(PC$biplot_pc_scores[, axes[2]])
  rangy1 <- unsigned.range(PC$biplot_var_scores[, axes[1]])
  rangy2 <- unsigned.range(PC$biplot_var_scores[, axes[2]])

  xlim = ylim = rangx1 = rangx2 = range(rangx1, rangx2)

  ratio = max(rangy1/rangx1, rangy2/rangx2)

  xlim2 = xlim * ratio
  ylim2 = ylim * ratio

  par(pty = "s")

  if (is.null(labels)) {

    plot(PC$biplot_pc_scores[,axes[1]],
         PC$biplot_pc_scores[,axes[2]],
         pch = 19, cex = point_size,
         xlim = xlim, ylim = ylim,
         xlab = paste("PC", axes[1], " (",
                      round(expl_var[axes[1]] * 100,2),"%)", sep = ""),
         ylab = paste("PC", axes[2], " (",
                      round(expl_var[axes[2]] * 100,2),"%)", sep = ""))
    par(new = TRUE)
    plot(PC$biplot_var_scores[,axes[1]],
         PC$biplot_var_scores[,axes[2]], axes = FALSE, type = "n",
         xlim = xlim2,
         ylim = ylim2,
         xlab = "", ylab = "")
    #axis(3)
    #axis(4)
    text(arrow_values[,1], arrow_values[,2], arrow_values$var_name,
         col = arrow_colour, font = 2)
    arrows(0, 0, arrow_values[,1] * 0.8, arrow_values[,2] * 0.8,
           lwd = arrow_size,
           length = 0.1, col = arrow_colour)

  } else {

    if (!is.null(label_colours)) {

      if (length(label_colours) < length(levels(as.factor(labels)))) {
        stop("Insufficient label_colours arguments provided.")
      }

      plot_colours <- label_colours

    } else {
      plot_colours <- c(
        "black","red","blue","orange","darkgreen","violetred"
      )
    }

    label_levels <- levels(as.factor(labels))

    plot(PC$biplot_pc_scores[labels == label_levels[1],axes[1]],
         PC$biplot_pc_scores[labels == label_levels[1],axes[2]],
         pch = 19, cex = point_size,
         col = plot_colours[1],
         xlim = xlim, ylim = ylim,
         xlab = paste("PC", axes[1], " (",
                      round(expl_var[axes[1]] * 100,2),"%)", sep = ""),
         ylab = paste("PC", axes[2], " (",
                      round(expl_var[axes[2]] * 100,2),"%)", sep = ""))
    if (CI_ellipse == TRUE) {
      ci <- car::dataEllipse(PC$biplot_pc_scores[labels == label_levels[1],
                                                 c(axes[1], axes[2])],
                             levels = CI_level, draw = FALSE)
      lines(ci, lwd = CI_size, col = plot_colours[1])
    }
    if (Chull == TRUE) {
      chull_plot(PC$biplot_pc_scores[labels == label_levels[1],
                                     axes[1]],
                 PC$biplot_pc_scores[labels == label_levels[1],
                                     axes[2]],
                 colour = plot_colours[1],
                 lwd = Ch_size)
    }
    for (i in 2:length(label_levels)) {
      points(PC$biplot_pc_scores[labels == label_levels[i],axes[1]],
             PC$biplot_pc_scores[labels == label_levels[i],axes[2]],
             pch = 19, cex = point_size,
             col = plot_colours[i])
      if (CI_ellipse == TRUE) {
        ci <- car::dataEllipse(PC$biplot_pc_scores[labels == label_levels[i],
                                                   c(axes[1], axes[2])],
                               levels = CI_level, draw = FALSE)
        lines(ci, lwd = CI_size, col = plot_colours[i])
      }
      if (Chull == TRUE) {
        chull_plot(PC$biplot_pc_scores[labels == label_levels[i],
                                       axes[1]],
                   PC$biplot_pc_scores[labels == label_levels[i],
                                       axes[2]],
                   colour = plot_colours[i],
                   lwd = Ch_size)
      }
    }
    par(new = TRUE)
    plot(PC$biplot_var_scores[,axes[1]],
         PC$biplot_var_scores[,axes[2]], axes = FALSE, type = "n",
         xlim = xlim2,
         ylim = ylim2,
         xlab = "", ylab = "")
    #axis(3)
    #axis(4)
    text(arrow_values[,1], arrow_values[,2], arrow_values$var_name,
         col = arrow_colour, font = 2)
    arrows(0, 0, arrow_values[,1] * 0.8, arrow_values[,2] * 0.8,
           lwd = arrow_size,
           length = 0.1, col = arrow_colour)

    legend(legend_position,
           legend = label_levels,
           pch = 19,
           col = plot_colours, inset = legend_inset)

  }

  if (!is.null(main)){
    title(main = main)
  }

  par(pty = "m")

  return(variable_importance_table)

}
