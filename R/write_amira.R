#
# #' Miscellaneous: Write Amira files from landmarks
# #'
# #' A miscellaneous function that can be used to convert landmark arrays into amira
# #' files
# #'
# #' @param file_name A string specifying the name of the file that will be produced
# #' @param landmarks A (landmark x dimension x individual)  or (landmark x dimension)
# #' array containing landmark data
# #'
# #' @param Note:
# #' If an array is provided with more than one individual (landmark x dimension x individual),
# #' a multitude of amira files will be created, duplicating the file_name multiple times;
# #' e.g. \code{write_amira("landmarks", landmarks)} will produce landmarks_1.txt,
# #' landmarks_2.txt, landmarks_3.txt etc.
# #'
# #' @return Produces an amira ASCII file.
# #'
# #' @author Lloyd A. Courtenay
# #'
# #' @examples
# #' library(geomorph)
# #' library(shapes)
# #' library(GraphGMM)
# #'
# #' # load data
# #' data(apes)
# #'
# #' # Generalized Procrustes Fit
# #' GPAshape <- procGPA(apes$x)
# #'
# #' # calculate central configuration
# #' central_config <- calc_central_morph(GPAshape$rotated)
# #'
# #' write_amira("central_landmark_configuration", central_config)
# #'

write_amira <- function(file_name,
                             landmarks) {

  n_dimensions <- dim(landmarks)[2]
  n_landmarks <- dim(landmarks)[1]

  if (n_dimensions != 3) {
    stop("Amira files can only be created for 3D landmark data")
  }

  if(is.na(dim(landmarks)[3])) {

    single_amira_file(file_name, landmarks, n_landmarks)

  } else {
    for (individual in 1:dim(landmarks)[3]) {
      ind_file_name <- paste0(file_name, "_", individual)
      single_amira_file(ind_file_name, landmarks[,,individual], n_landmarks)
    }
  }

}
