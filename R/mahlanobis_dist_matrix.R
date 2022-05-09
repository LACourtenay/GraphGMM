
#' Mahalanobis Distances and p-Values
#'
#' A function used to calculate inter and intra-group mahalanobis distances. p-Values
#' are computed to calculate statistical differences between samples
#'
#' @param data A matrix containing the data
#' @param labels A factorial list containing the sample labels
#' @param plot_results A boolean value defining whether results are to be plotted.
#' A maximum of 4 groups can be included for results to be plotted.
#' @param robust A boolean value defining whether p-values are to be computed
#' using robust statistical metrics.
#'
#' @return \code{distances} - A distance matrix
#' @return \code{p_values} - p-Values computed from the distance matrix
#'
#' @seealso
#' \code{\link{procD_comparisons}}
#' 
#' @author Lloyd A. Courtenay
#'
#' @examples
#' library(geomorph)
#' library(shapes)
#' library(GraphGMM)
#'
#' # example 1 ----------------------------------
#'
#' # load data
#' data(plethodon)
#'
#' # Generalized Procrustes Fit
#' GPAshape <- procGPA(plethodon$land)
#' data <- vector_from_landmarks(GPAshape$rotated)
#'
#' # pca
#'
#' pca <- pca_plot(data, plethodon$species)
#'
#' pca$pca_plot # visualise pca
#' pc_scores <- pca$pc_scores # extract pc_scores
#'
#' # mahalanobis distances
#'
#' mahalanobis_dist_matrix(pc_scores, plethodon$species)
#'
#' # example 2 ----------------------------------
#'
#' # load data
#' data(apes)
#'
#' # Generalized Procrustes Fit
#' GPAshape <- procGPA(apes$x)
#' data <- vector_from_landmarks(GPAshape$rotated)
#'
#' # pca
#'
#' pca <- pca_plot(data, apes$group)
#'
#' pca$pca_plot # visualise pca
#' pc_scores <- pca$pc_scores # extract pc_scores
#'
#' # mahalanobis distances
#'
#' mahalanobis_dist_matrix(pc_scores, apes$group, plot_results = FALSE)
#'
#' @export

mahalanobis_dist_matrix <- function(data, labels,
                                   plot_results = TRUE,
                                   robust = TRUE) {

  if(!is.matrix(data)) {
    stop("Input must be a matrix")
  }

  if(missing(data) | missing(labels)) {
    stop("Both a matrix containing multivariate distributions and a set of labels need to be specified for this function")
  }

  level_list = levels(as.factor(labels))

  if (plot_results == TRUE & length(level_list) > 4) {
    stop(
      "plot_results can only be TRUE when comparing 4 samples due to issues with R plot margins"
    )
  }

  data = as.matrix(data)
  level_list = levels(as.factor(labels))

  dist_matrix = array(numeric(), dim = c(length(level_list), length(level_list)),
                      dimnames = list(level_list, level_list))
  p_value_matrix = array(numeric(), dim = c(length(level_list), length(level_list)),
                         dimnames = list(level_list, level_list))

  if (plot_results != FALSE) {
    dev.new(); par(mfrow = c(length(level_list), length(level_list)))
  }

  pb <- txtProgressBar(min = 0, max = length(level_list),
                       style = 3, width = 100, char = "=")

  for (i in 1:length(level_list)) {

    reference <- data[labels == level_list[i],]

    reference_distribution <- data.frame(
      distances = mahalanobis_distances_to_target(reference,
                                                 reference),
      Sample = as.factor(level_list[i])
    )

    for (j in 1:length(level_list)) {

      target <- data[labels == level_list[j],]

      target_distribution <- data.frame(
        distances = mahalanobis_distances_to_target(
          target, reference
        ),
        Sample = as.factor(level_list[j])
      )

      compare <- rbind(reference_distribution, target_distribution)

      if (plot_results != FALSE) {
        sm::sm.density.compare(compare$distances, compare$Sample,
                               lwd = 2, col = c("black","red"))
        title(main = paste(
          levels(compare$Sample)[1], " vs ",
          if (is.na(levels(compare$Sample)[2])) {
            levels(compare$Sample)[1]
          } else {
            levels(compare$Sample)[2]
          }
        ))
      }

      if (levels(target_distribution$Sample) == levels(reference_distribution$Sample)) {
        p_val = 0
        dist = 0
      } else {

        index_1 <- table(compare$Sample)[1]
        index_2 <- table(compare$Sample)[2]

        if (robust == TRUE) {
          p_val <- kruskal.test(compare$distances, compare$Sample)$"p.value"
          dist <- median(compare[index_1:index_2, 1])
        } else {
          p_val <- summary(aov(compare$distances ~ compare$Sample))[[1]]$`Pr(>F)`[1]
          dist <- mean(compare[index_1:index_2, 1])
        }

      }

      dist_matrix[i, j] = dist
      p_value_matrix[i, j] = p_val

    }

    setTxtProgressBar(pb, i)

  }

  if (plot_results != FALSE) {
    par(mfrow = c(1,1))
  }

  cat("\n")

  return(list(distances = as.dist(t(dist_matrix)),
              p_values = as.dist(t(p_value_matrix))))
}
