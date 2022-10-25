
#' Calculate the Optimal Number of PC Scores
#'
#' The present function can be used to calculate the optimal number of PC Scores
#' required for a study
#'
#' @param data A \code{data.frame} or \code{matrix} containing the data that
#' will be used for Principal
#' Component Analysis
#' @param N_permutations An integer defining the number of permutations to be
#' used in the calculation.
#' @param centre A boolean variable specifying whether data is to be centred first
#' @param scale A boolean variable specifying whether data is to be scaled first
#'
#' @section Details:
#' The calculation of the optimal number of PC scores is performed by calculating the
#' observed variance explained by each PC with the permuted variance.
#'
#' @return A plot presenting the results from the permuted calculation of the
#' optimal number PC scores
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
#' # calculate optimal number of PC scores
#' data <- vector_from_landmarks(GPAshape$rotated)
#'
#' calculate_optimal_pc_scores(data, N_permutations = 10000)
#'  
#' @export

calculate_optimal_pc_scores <- function(data, N_permutations = 10000,
                                        centre = TRUE, scale = FALSE) {

  if(!is.matrix(data) & !is.data.frame(data)) {
    stop("Data parameter must be a data.frame or a matrix of values")
  }

  if(N_permutations <= 0 | N_permutations %% 1 != 0) {
    stop("N_permutations parameter must be a positive, non-zero integer")
  }

  if (!is.logical(centre)) {
    stop("The centre parameter only accepts boolean TRUE or FALSE values")
  }

  if (!is.logical(scale)) {
    stop("The scale parameter only accepts boolean TRUE or FALSE values")
  }

  pb <- txtProgressBar(min = 0, max = N_permutations,
                       style = 3, width = 100, char = "=")

  PC <- prcomp(as.matrix(data), center = centre, scale = scale)

  expl_var <- PC$sdev^2 / sum(PC$sdev^2)

  expl_var_perm <- matrix(NA, ncol = length(PC$sdev), nrow = N_permutations)

  for (i in 1:N_permutations) {

    expr_perm <- apply(as.matrix(data), 2, sample)
    PC_perm <- prcomp(as.matrix(expr_perm), center = centre, scale = scale)
    expl_var_perm[i, ] <- PC_perm$sdev^2 / sum(PC_perm$sdev^2)

    setTxtProgressBar(pb, i)

  }

  par(mfrow = c(1,2))

  plot(expl_var[1:length(expl_var)] ~ seq(1:length(expl_var)), ylab = "Explained Variance",
       col = "darkgrey", type = "o", xlab = "PCs")

  lines(colMeans(expl_var_perm)[1:length(expl_var)] ~ seq(1:length(expl_var)), col = "red")

  legend("topright", c("Explained by PCs", "Explained by Chance"),
         fill = c("darkgrey","red"), inset = 0.02)

  pval <- apply(t(expl_var_perm) >= expl_var, 1, sum) / 10000

  plot(pval[1:length(expl_var)] ~ seq(1:length(expl_var)), col = "darkred", type = "o",
       xlab = "PCs", ylab = "p-Value")

  optPC <- head(which(pval > 0.005), 1)-1

  if (length(optPC) == 0) {
    optPC <- length(expl_var)
  }

  mtext(paste0("Optimal Number of PCs = ", optPC, " (", round(sum(expl_var[1:optPC]) * 100, 2), "% Var.)"))

  par(mfrow = c(1,1))

}
