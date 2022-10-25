
#' Procrustes p-Value calculations
#'
#' A function used to calculate inter and intra-group Procrustes distance p-Values.
#'
#' @param lm_matrix A matrix containing the data
#' @param labels A factorial list containing the sample labels
#' @param plot_results A boolean value defining whether results are to be plotted.
#' A maximum of 4 groups can be included for results to be plotted.
#' @param robust A boolean value defining whether p-values are to be computed
#' using robust statistical metrics.
#'
#' @return p-Values computed from procrustes distance matrices
#'
#' @seealso
#' \code{\link{mahalanobis_dist_matrix}}
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
#' GPAshape <- GPA(plethodon$land)
#'
#' # procrustes distances
#'
#' procD_comparisons(GPAshape$rotated, plethodon$species)
#'
#' # example 2 ----------------------------------
#'
#' # load data
#' data(apes)
#'
#' # Generalized Procrustes Fit
#' GPAshape <- GPA(apes$x)
#'
#' # procrustes distances
#'
#' procD_comparisons(GPAshape$rotated, apes$group, plot_results = FALSE)
#'
#' @export

procD_comparisons <- function(lm_matrix, labels, plot_results = TRUE,
                                 robust = TRUE) {

  if(!is.array(lm_matrix) | is.na(dim(lm_matrix)[3])) {
    stop("Input must be a (n x n x individual) array")
  }

  if(missing(lm_matrix) | missing(labels)) {
    stop("Both a landmark array and set of labels need to be specified for this function")
  }

  level_list = levels(as.factor(labels))

  if (plot_results == TRUE & length(level_list) > 4) {
    stop(
      "plot_results can only be TRUE when comparing 4 samples due to issues with R plot margins"
    )
  }

  matrices = lm_matrix
  p_value_matrix = array(numeric(), dim = c(length(level_list), length(level_list)),
                         dimnames = list(level_list, level_list))
  reference_matrices = array(numeric(), dim = c(dim(matrices)[1], dim(matrices)[2], 0))

  for (sample in level_list) {
    target_sample <- separate_samples(matrices, as.factor(labels), sample)
    if (robust == TRUE) {
      mshape <- calc_central_morph(target_sample, method = "median")
    } else {
      mshape <- calc_central_morph(target_sample, method = "mean")
    }
    reference_matrices <- abind::abind(reference_matrices, mshape, along = 3)
  }

  if (plot_results != FALSE) {
    dev.new(); par(mfrow = c(length(level_list), length(level_list)))
  }

  pb <- txtProgressBar(min = 0, max = dim(reference_matrices)[3],
                       style = 3, width = 100, char = "=")

  for (i in 1:dim(reference_matrices)[3]) {
    from_target <- c()
    for (ind in 1:dim(matrices)[3]) {
      from_target <- c(from_target, shapes::procdist(matrices[,,ind], reference_matrices[,,i]))
    }
    target_distribution <- data.frame(
      proc_dist = from_target[as.factor(labels) == level_list[i]],
      Sample = as.factor(level_list[i])
    )
    for(j in 1:dim(reference_matrices)[3]) {
      compare_distribution <- data.frame(
        proc_dist = from_target[as.factor(labels) == level_list[j]],
        Sample = as.factor(level_list[j])
      )
      compare <- rbind(target_distribution, compare_distribution)

      if (plot_results != FALSE) {
        sm::sm.density.compare(compare$proc_dist, compare$Sample,
                               lwd = 2, col = c("black", "red"))
        title(main = paste(
          levels(compare$Sample)[1], " vs ",
          if (is.na(levels(compare$Sample)[2])) {
            levels(compare$Sample)[1]
          } else {
            levels(compare$Sample)[2]
          }
        ))
      }

      if (levels(target_distribution$Sample) == levels(compare_distribution$Sample)) {
        p_val = 0
        dist = 0
      } else {

        index_1 <- table(compare$Sample)[1]
        index_2 <- table(compare$Sample)[2]

        if (robust == TRUE) {
          p_val <- kruskal.test(compare$proc_dist, compare$Sample)$"p.value"
          dist = median(compare[index_1:index_2,1])
        } else {
          p_val = summary(aov(compare$proc_dist ~ compare$Sample))[[1]]$`Pr(>F)`[1]
          dist = mean(compare[index_1:index_2,1])
        }

      }

      #dist_matrix[i, j] = dist
      p_value_matrix[i, j] = p_val

    }
    setTxtProgressBar(pb, i)
  }

  if (plot_results != FALSE) {
    par(mfrow = c(1,1))
  }

  cat("\n")

  return(as.dist(t(p_value_matrix)))
}
