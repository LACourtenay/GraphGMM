
#' Miscellaneous: Convert an array of landmarks into a morphologika file
#'
#' A miscellaneous function that can be used to format landmarks into a morphologika
#' file
#'
#' @param file_name A string specifying the name of the file that will be produced
#' @param landmarks A (landmark x dimension x individual)  or (landmark x dimension)
#' array containing landmark data
#' @param labels A factorial list specifying the label that corresponds to each individual
#'
#' @return Produces a morphologika file.
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#' library(shapes)
#' library(GraphGMM)
#'
#' # load data
#' data(macf.dat)
#' data(macm.dat)
#'
#' library(abind)
#'
#' mac.dat <- abind(macf.dat, macm.dat, along = 3)
#'
#' f_list <- rep("F", dim(macf.dat)[3])
#' m_list <- rep("M", dim(macm.dat)[3])
#' label_list <- as.factor(c(f_list, m_list))
#'
#' write_morphologika_file("mac_data", mac.dat,
#'                         label_list)
#'
#' @export

write_morphologika_file <- function(file_name,
                                     landmarks,
                                     labels) {

  if (!is.character(file_name)) {
    stop("file_name must be a valid string")
  }

  if (!is.array(landmarks)) {
    stop("Landmark configuration must be a (landmark x dimension x individual) array")
  }

  if (missing(labels) | missing(landmarks) | missing(file_name)) {
    stop("The function requires a file_name, landmark array and factorial vector of labels.")
  }

  n_individuals <- dim(landmarks)[3]
  n_dimensions <- dim(landmarks)[2]
  n_landmarks <- dim(landmarks)[1]

  samples <- data.frame(sample_names = levels(labels),
                        number_samples = NA)
  for (i in 1:nrow(samples)) {
    samples[i,2] <- summary(labels)[[i]]
  }

  output <- paste0("[individuals]\n", n_individuals)
  output <- paste0(output, "\n[landmarks]\n", n_landmarks)
  output <- paste0(output, "\n[dimensions]\n", n_dimensions)
  output <- paste0(output, "\n[names]\n")

  name_list <- c()
  sample_list <- c()
  for (sample_name in 1:nrow(samples)) {
    for (sample_number in 1:samples[sample_name,2]) {
      name_list <- c(name_list, paste0(samples[sample_name,1], sample_number))
      sample_list <- c(sample_list, samples[sample_name,1])
    }
  }

  output <- paste0(output, name_list[1])

  for (i in 2:length(name_list)) {
    output <- paste0(output, "\n", name_list[i])
  }

  output <- paste0(output, "\n[labels]\nSample\n[labelvalues]\n", sample_list[1])

  for (i in 2:length(sample_list)) {
    output <- paste0(output, "\n", sample_list[i])
  }

  output <- paste0(output, "\n[rawpoints]\n")

  for (individual in 1:n_individuals) {
    configuration <- ""
    for (landmark in 1:n_landmarks) {
      for (dimension in 1:n_dimensions) {
        if (dimension == 1) {
          configuration <- paste0(
            configuration, landmarks[landmark,dimension,individual]
          )
        } else {
          configuration <- paste0(
            configuration, " ", landmarks[landmark,dimension,individual]
          )
        }
      }
      configuration <- paste0(
        configuration, "\n"
      )
    }
    output <- paste0(output, "'", name_list[individual], "\n", configuration, "\n")
  }

  suff <- substr(file_name, nchar(file_name) - 3, nchar(file_name))

  if (suff != ".txt") {
    file_name <- paste0(file_name, ".txt")
  }

  cat(output, file = file_name)

}
