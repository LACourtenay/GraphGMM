
#' t-Distributed Stochastic Neighbour Embedding
#'
#' This function calculates and plots the results of a tSNE, a non-linear
#' dimensionality reduction algorithm.
#'
#' @param data A \code{data.frame} containing the information to plot
#' @param n_iterations An integer defining the number of iterations to be performed
#' @param perplexity A numeric value defining the perplexity hyperparameter. The
#' default is set to the sample size to the power of 1/2.
#' @param labels A factorial list containing the sample labels
#' @param label_colours A string vector of colours to be used in the plot.
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
#' @return \code{tSNE_scores} - The 2 dimensional tSNE scores.
#' @return \code{tSNE_plot} - a ggplot object containing the tSNE plot
#'
#' @section Note:
#' tSNE is a computationally expensive algorithm, and does not work well (if at all) on
#' similarity matrices. t-SNE should therefore be used on landmark embeddings or directly
#' on landmark coordinates.
#'
#' @section Bibliography:
#' Hinton, G.E.; Roweis, S.T. (2003) Stochastic Neighbor Embedding, Advances in Neural
#' Information Processing Systems. 857-864
#'
#' @seealso
#' \code{\link[Rtsne]{Rtsne}} of the \code{Rtsne} library,
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
#' data <- vector_from_landmarks(GPAshape$coordinates)
#'
#' # plot pca
#'
#' tsne_plot(data, apes$group, CI_ellipse = TRUE)
#'
#' @export

tsne_plot <- function(data, labels = NULL,
                      n_iterations = 1000, perplexity = 0,
                      label_colours = NULL,
                      point_size = 2, CI_ellipse = FALSE,
                      CI_level = 0.95, CI_size = 1,
                      Chull = FALSE, Ch_size = 1,
                      main = NULL) {

  if (!is.data.frame(data)) {
    stop("Input data must be of type data.frame")
  }

  if (!is.numeric(point_size) | point_size < 0) {
    stop("point_size must be a positive value")
  }

  if (!is.logical(CI_ellipse)) {
    stop("The CI_ellipse parameter only accepts boolean TRUE or FALSE values")
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

  if (perplexity == 0) {
    perplexity = dim(data)[1] ** (1/2)
  }

  tsne_results <- Rtsne::Rtsne(
    data, dims = 2,
    perplexity = perplexity,
    verbose = TRUE,
    max_iter = n_iterations
  )

  # present posibility of using convex hulls

  if (!is.null(labels)) {

    dim_red <- data.frame(Sample = as.factor(labels),
                          x = tsne_results$Y[,1],
                          y = tsne_results$Y[,2])
    base_plot <- ggplot2::ggplot(data = dim_red,
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
    dim_red <- data.frame(x = tsne_results$Y[,1],
                          y = tsne_results$Y[,2])
    base_plot <- ggplot2::ggplot(data = dim_red,
                                 ggplot2::aes(x = x, y = y))
    plot_colours <- ggplot2::scale_color_manual(values = c("black"))
  }

  tsne_plot <- base_plot +
    ggplot2::geom_point(stat = "identity", size = point_size) +
    plot_colours

  if (CI_ellipse == TRUE) {
    conf_interval <- ggplot2::stat_ellipse(size = CI_size,
                                           level = CI_level)
    tsne_plot <- tsne_plot + conf_interval
  }

  if (Chull == TRUE) {
    hulls <- dim_red %>%
      dplyr::group_by(Sample) %>%
      dplyr::slice(chull(x, y))
    conf_interval <- ggplot2::geom_polygon(
      data = hulls,
      alpha = 0,
      size = Ch_size
    )
    tsne_plot <- tsne_plot + conf_interval
  }

  tsne_plot <- tsne_plot +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
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
                   legend.position = "bottom") +
    ggplot2::geom_vline(xintercept = 0, colour = "black", size = 0.5, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linetype = "dashed", size = 0.5)

  if (!is.null(main)) {
    tsne_plot <- tsne_plot + ggplot2::ggtitle(main)
  }

  return(list(tSNE_scores = dim_red,
              tSNE_plot = tsne_plot))

}
