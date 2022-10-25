
#' Miscellaneous: A function to separate landmark samples
#'
#' This function can be used to separate landmarks from a specific sample.
#'
#' @param GPAmatrix A (landmark x dimension x individual) array containing the landmark
#' data
#' @param labels A factorial list containing the sample labels
#' @param target_sample A string defining the name of the sample the user wishes to
#' separate
#'
#' @return A (landmark x dimension x individual) array of landmarks from the defined sample
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
#' central_female <- calc_central_morph(female_gorilla)
#' central_male <- calc_central_morph(male_gorilla)
#'
#' plot(central_female, pch = 19, asp = 1)
#' plot(central_male, pch = 19, asp = 1)
#'
#' @export

separate_samples <- function(GPAmatrix, labels, target_sample) {

  if (!is.array(GPAmatrix) | is.na(dim(GPAmatrix)[3])) {
    stop("Landmark configuration must be a (landmark x dimension x individual) array")
  }

  if (missing(labels) | missing(target_sample)) {
    stop("labels and a target sample must be specified")
  }

  labels <- as.factor(labels)
  index <- labels == target_sample

  flat_matrix<-c(); for (i in 1:length(labels)) {
    flat <- flatten_matrix(GPAmatrix[,,i])
    flat_matrix<-rbind(flat_matrix, flat)
  }; flat_matrix<-as.data.frame(flat_matrix)

  sample <- flat_matrix[index,]

  dimensions <- dim(GPAmatrix)
  coord_tensor <- array(numeric(), c(dimensions[1], dimensions[2], 0))

  for (i in 1:nrow(sample)) {
    mat <- matrix(sample[i, ], nrow = dimensions[1], byrow = TRUE)
    coord_tensor <- abind::abind(coord_tensor, mat, along = 3)
  }

  return(coord_tensor)

}
