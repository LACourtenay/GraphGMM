
#' Calculate Graph Configuration Statistics
#'
#' This function can be used to visualise the calculated properties of a given landmark graph.
#'
#' @param target_central A (landmark x dimension) matrix containing the central landmark
#' configuration the user wishes to analyse
#' @param variable A vector containing the variables to be visualized, as extracted
#' from \code{\link{graph_configuration_statistics}}
#' @param scale_factor A numeric value defining how to scale the statistical
#' metrics
#' @param log_values A boolean variable defining whether the log function of
#' statistical variables should be calculated.
#' @param colour_grad A boolean variable defining whether visualisations should
#' include colour gradients as well.
#' @param colour_values A vector defining two colours that will be used to calculate
#' the colour gradients.
#' @param n_col_ranks A numeric value defining the number of levels in the colour
#' gradients.
#' @param surface A 3D surface model (ASCII .ply file) containing the surface the user
#' wishes to plot the landmarks on (note colour_grad is ignored if a surface is provided)
#'
#' @return A visualisation of the landmark graph
#'
#' @seealso
#' \code{\link{separate_samples}}, \code{\link{graph_configuration_statistics}}
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
#' female_results <- graph_configuration_statistics(female_central, female_edges)
#' male_results <- graph_configuration_statistics(male_central, male_edges)
#'
#' # Visualise results
#'
#' plot_lmgraph_stats(female_central, female_results$landmark_degree,
#'                    colour_values = c("red","blue"),
#'                    n_col_ranks = 5)
#' plot_lmgraph_stats(female_central, female_results$landmark_eigenvalues,
#'                    scale = 5,
#'                    colour_values = c("red","blue"),
#'                    n_col_ranks = 5)
#' plot_lmgraph_stats(female_central, female_results$landmark_betweenness,
#'                    scale = 0.5,
#'                    colour_values = c("red","blue"),
#'                    n_col_ranks = 5)
#'
#' plot_lmgraph_stats(male_central, female_results$landmark_degree,
#'                    colour_values = c("red","blue"),
#'                    n_col_ranks = 5)
#' plot_lmgraph_stats(male_central, female_results$landmark_eigenvalues,
#'                    scale = 5,
#'                    colour_values = c("red","blue"),
#'                    n_col_ranks = 5)
#' plot_lmgraph_stats(male_central, female_results$landmark_betweenness,
#'                    scale = 0.5,
#'                    colour_values = c("red","blue"),
#'                    n_col_ranks = 5)
#'
#' @export

plot_lmgraph_stats <- function(target_central, variable, scale_factor = 1,
                               log_values = FALSE, colour_grad = FALSE,
                               colour_values = c("black", "red"),
                               n_col_ranks = 10, surface = NULL) {

  if("array" %!in% class(target_central)) {
    if ("matrix" %!in% class(target_central)) {
      stop("Landmark central configuration must be a (landmark x dimension) array")
    }
  }

  if(dim(target_central)[2] != 2 & dim(target_central)[2] != 3) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if(!is.vector(variable)) {
    stop("The variable to be visualised must be a vector")
  }

  if (length(variable) != dim(target_central)[1]) {
    stop("The variable vector to be visualised must be the same length as the number of landmarks")
  }

  if (log_values == TRUE & 0 %in% variable) {
    stop("log_values cannot be equal to TRUE if the vector to be visualised contains 0 values")
  }

  if (!is.logical(log_values)) {
    stop("Thelog_values parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.logical(colour_grad)) {
    stop("The colour_grad parameter only accepts boolean TRUE or FALSE values")
  }

  if (length(colour_values) != 2) {
    stop("Only two colour values can be supplied for the colour gradients")
  }

  colour_bool <- check_colours(colour_values)
  if (FALSE %in% colour_bool) {
    stop("Invalid colour provided for colour_values")
  }

  if (!is.null(surface)) {
    if ("mesh3d" %!in% class(surface)) {
      stop("Surface must be a 3D mesh, preferably loaded using the read.ply function in geomorph")
    }
  }

  variable = variable * scale_factor

  if (log_values == TRUE) {
    variable = log(variable)
  }

  if (colour_grad == TRUE) {
    ranks <- as.factor(as.numeric(cut(
      variable, n_col_ranks
    )))
    palette_function <- colorRampPalette(c(colour_values[1], colour_values[2]))
    colours <- palette_function(n_col_ranks)
    #colours <- paletteer::paletteer_c(palette = "viridis::inferno", n = 13)
    #colours <- paletteer::paletteer_c(palette = "grDevices::heat.colors", n = 10, direction = -1)
  }

  if (dim(target_central)[2] == 2) {

    plot(target_central[,1], target_central[,2], col = NA, asp = 1)
    if (colour_grad == TRUE) {
      for (i in 1:dim(target_central)[1]) {
        if (variable[i] != 0) {
          points(target_central[i,1], target_central[i,2],
                 cex = variable[i], pch = 19,
                 col = colours[ranks[i]])
        }
      }
    } else {
      for (i in 1:dim(target_central)[1]) {
        if (variable[i] != 0) {
          points(target_central[i,1], target_central[i,2],
                 cex = variable[i], pch = 19)
        }
      }
    }


  } else {

    rgl::open3d()

    if (!is.null(surface)) {
      if (!is.null(surface) & colour_grad == TRUE) {
        warning("colour_grad is not applicable to plots where surfaces are included.")
      }
      rgl::shade3d(surface, col = "grey")
      for (i in 1:dim(target_central)[1]) {
          rgl::spheres3d(
            target_central[i,1], target_central[i,2], target_central[i,3],
            radius = variable[i]
          )
        }
    } else {
      if (colour_grad == TRUE) {
        for (i in 1:dim(target_central)[1]) {
          if (variable[i] != 0) {
            rgl::points3d(target_central[i,1], target_central[i,2], target_central[i,3],
                          size = variable[i], col = colours[ranks[i]])
          }
        }
      } else {
        for (i in 1:dim(target_central)[1]) {
          if (variable[i] != 0) {
            rgl::points3d(target_central[i,1], target_central[i,2], target_central[i,3],
                          size = variable[i])
          }
        }
      }
    }

    

  }

}
