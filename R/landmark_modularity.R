
#' Calculate Landmark Modules
#'
#' A function used to calculate and plot the landmark modules for a given configuration
#'
#' @param target_central A (landmark x dimension) matrix containing the central landmark
#' configuration the user wishes to analyse
#' @param target_edges An edge list that defines the graph, using
#' \code{\link{as_edge_list}}
#' @param dimensions A two-value vector defining the dimensions to plot in the case
#' of using a 2D visualisation (e.g. c("x","y"))
#' @param point_size A numeric value defining the size of plot points
#' @param Chull A boolean specifying whether convex hulls are to be included
#' @param Ch_size A numeric value defining the thickness of Convex Hull lines
#' @param Ch_alpha A numeric value defining how transparent the convex hulls are
#' @param plot_type "2D" or "3D" defining the type of visualisation
#' @param main A string specifying the plot title (if included)
#'
#' @section Details:
#' This function uses the Louvian algorithm to calculate neighbourhood structure
#' (modularity) within the landmark graph.
#'
#' @return \code{modularity_results} - a table of landmarks and their corresponding
#' module
#' @return \code{visualisation} - a ggplot object visualising the modules (for 2D
#' visualisation only)
#'
#' @section Bibliography:
#' Blondel, V.D.; Guillaume, J.L.; Lambiotte, R.; Lefebvre, E. (2008)
#' Fast unfolding of communities in large networks. Journal of Statistical
#' Mechanics. arXiv: 0803.0476v2
#'
#' @seealso
#' \code{\link[igraph]{cluster_louvain}} from the Igraph library
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
#' # separate the samples under analysis
#'
#' female_gorilla <- separate_samples(GPAshape$coordinates,
#'                                    apes$group,
#'                                    "gorf")
#' male_gorilla <- separate_samples(GPAshape$coordinates,
#'                                  apes$group,
#'                                  "gorm")
#'
#' female_central <- calc_central_morph(female_gorilla)
#' male_central <- calc_central_morph(male_gorilla)
#'
#' # calculate edges
#'
#' female_edges <- knn_graph(female_central, k = 3, radius = 75)
#' male_edges <- knn_graph(male_central, k = 3, radius = 75)
#' female_edges <- as_edge_list(female_edges)
#' male_edges <- as_edge_list(male_edges)
#'
#' # Graph configuration statistics
#'
#' landmark_modularity(female_central, female_edges, Chull = FALSE,
#'                     point_size = 5)
#' landmark_modularity(male_central, male_edges, Chull = FALSE,
#'                     point_size = 5)
#' 
#' @export

landmark_modularity <- function(target_central, target_edges,
                                dimensions = c("x","y"),
                                point_size = 2,
                                Chull = TRUE,
                                Ch_size = 1,
                                Ch_alpha = 0,
                                plot_type = "2D",
                                main = NULL) {

  if("array" %!in% class(target_central) | "array" %!in% class(target_edges)) {
    if ("matrix" %!in% class(target_central) | "matrix" %!in% class(target_edges)) {
      stop("Landmark central configuration must be a (landmark x dimension) array")
    }
  }

  if(dim(target_central)[2] != 2 & dim(target_central)[2] != 3) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if (length(dimensions) != 2) {
    stop("Dimensions requires a list of 2 pc scores to visualize, e.g. c(1,2) for PC1 and PC2")
  }

  if (!is.character(dimensions)) {
    stop("Dimensions must be specified as 'x', 'y' or 'z'")
  }

  possible_plots <- c("2D", "3D")
  if (plot_type %!in% possible_plots) {
    stop("User must define plot_type as either 2D or 3D")
  }

  if (plot_type == "3D" & dim(target_central)[2] != 3) {
    stop("plot_type can only be 3D if the landmark configuration being used is in three dimensions")
  }

  graph_df <- get_graph_df(target_central)

  LM_graph <- igraph::graph_from_data_frame(vertices = graph_df,
                                            d = as.matrix(target_edges),
                                            directed = FALSE)

  iso <- which(igraph::degree(LM_graph) == 0)
  g2 <- igraph::delete.vertices(LM_graph, iso)

  # ceb <- igraph::cluster_edge_betweenness(igraph::as.undirected(g2))
  ceb <- igraph::cluster_louvain(igraph::as.undirected(g2)) # same as gephi

  groups <- data.frame(id = character(), group = factor())
  for (i in 1:length(ceb)) {
    new <- data.frame(id = as.character(ceb[i][[1]]), group = as.factor(i))
    groups <- rbind(groups, new)
  }; groups$group = as.factor(groups$group)

  modularity_data <- graph_df %>%
    dplyr::mutate(module = as.character(groups$group[match(
      x = graph_df$id, table = groups$id
    )])) %>%
    dplyr::mutate(id = as.character(id),
                  module = as.factor(dplyr::coalesce(module, "0")))

  if (plot_type == "2D") {

    if (dim(target_central)[2] != 2) {
      if (TRUE %!in% (dimensions == "x")) {
        modularity_data <- modularity_data %>% dplyr::select(-x)
      }
      if (TRUE %!in% (dimensions == "y")) {
        modularity_data <- modularity_data %>% dplyr::select(-y)
      }
      if (TRUE %!in% (dimensions == "z")) {
        modularity_data <- modularity_data %>% dplyr::select(-z)
      }
    }

    plot_mod_data <- modularity_data %>% dplyr::select(-id)
    if (dimensions[1] == "y") {
      plot_mod_data <- plot_mod_data %>% dplyr::mutate(
        x2 = y, y2 = x
      ) %>% dplyr::select(x2, y2, module)
    }
    colnames(plot_mod_data) <- c("x", "y", "module")

    gg_visualisation <- ggplot2::ggplot() +
      ggplot2::geom_point(data = plot_mod_data[plot_mod_data$module != 0,],
                          ggplot2::aes(x, y, col = module),
                          size = point_size) +
      ggplot2::geom_point(data = plot_mod_data[plot_mod_data$module == 0,],
                          ggplot2::aes(x, y),
                          col = "black",
                          size = point_size) +
      ggplot2::theme_bw() +
      ggplot2::coord_equal() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", size = 20)
      )

    if (Chull == TRUE) {

      hulls <- plot_mod_data %>%
        dplyr::filter(module != 0) %>%
        dplyr::group_by(module) %>%
        dplyr::slice(chull(x, y))

      gg_visualisation <- gg_visualisation +
        ggplot2::geom_polygon(
          data = hulls,
          ggplot2::aes(x, y, col = module, fill = module),
          alpha = Ch_alpha,
          size = Ch_size
        )

    }

    modularity_results <- modularity_data %>% dplyr::select(id, module)

    if (!is.null(main)) {
      gg_visualisation <- gg_visualisation + ggplot2::ggtitle(main)
    }

    return(list(modularity_results = modularity_results,
                visualisation = gg_visualisation))

  } else {

    if (Chull == TRUE) {

    }

    cols <- RColorBrewer::brewer.pal(n = length(levels(modularity_data$module)) - 1,
                                     name = "Spectral")
    rgl::open3d(); for (i in 1:length(levels(modularity_data$module))) {
      lm_group <- modularity_data[
        modularity_data$module == levels(modularity_data$module)[i],
      ]
      if(levels(modularity_data$module)[i] == 0) {
        rgl::spheres3d(lm_group[,2],
                       lm_group[,3],
                       lm_group[,4],
                       r = point_size,
                       col = "black")
      } else {
        rgl::spheres3d(lm_group[,2],
                       lm_group[,3],
                       lm_group[,4],
                       r = point_size,
                       col = cols[i])
      }

    }

    modularity_results <- modularity_data %>% dplyr::select(id, module)

    return(modularity_results)

  }

}
