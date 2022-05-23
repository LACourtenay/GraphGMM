
`%!in%` = Negate(`%in%`)

check_colours <- function(col_list) {
  sapply(col_list, function (values) {
    tryCatch(is.matrix(col2rgb(values)),
             error = function(e) FALSE)
  })
}

add_alpha <- function(col, alpha=1){ # Function for setting the alpha of colours
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

print_cum_var <- function(percents) {
  cummulative = 0
  for (i in 1:length(percents)) {
    cummulative <- cummulative + percents[i]
    cat(paste("PC", i, ": ", round(cummulative, 2), "%\n"))
  }
}

flatten_matrix <- function(lm_matrix) {

  lm_vector <- t(lm_matrix)
  dim(lm_vector) <- NULL
  return(lm_vector)

}

cosine_distance <- function(x, y) {

  return(crossprod(x, y) / sqrt(crossprod(x) * crossprod(y)))

}

single_value_decomposition <- function(X, centre = TRUE, scale = FALSE) {

  X <- as.matrix(X)

  tX <- scale(X, center = centre, scale = scale)

  eigenvectors_svd <- svd(tX)$v # same as prcomp rotation
  scores_svd <- svd(tX)$u %*% diag(svd(tX)$d) # pc scores
  eigenvalues_svd <- svd(tX)$d^2 / (nrow(X) - 1) # same as prcomp sdev^2
  loadings <- eigenvectors_svd %*% diag(prcomp(tX)$sdev)

  sdev <- prcomp(tX)$sdev

  lam <- prcomp(tX)$sdev
  n <- nrow(scores_svd)
  lam <- lam * sqrt(n) # lambda

  biplot_scores <- t(t(scores_svd) / lam)
  biplot_variables <- t(t(eigenvectors_svd) * lam)

  return(list(
    eigenvectors = eigenvectors_svd,
    eigenvalues = eigenvalues_svd,
    pc_scores = scores_svd,
    loadings = loadings,
    sdev = sdev,
    biplot_pc_scores = biplot_scores,
    biplot_var_scores = biplot_variables
  ))

}

chull_plot <- function(x, y, colour, lwd) {
  pts <- chull(x = x, y = y)
  pts <- c(pts, pts[1])
  lines(x[pts], y[pts], col = colour, lwd = lwd)
}

single_amira_file <- function(file_name, landmarks, n_landmarks) {

  suff <- substr(file_name, nchar(file_name) - 3, nchar(file_name))

  if (suff != ".txt") {
    file_name <- paste0(file_name, ".txt")
  }

  output <- paste0("# AmiraMesh 3D ASCII 2.0\n\n\ndefine Markers ", n_landmarks)
  output <- paste0(
    output,
    "\n\nParameters {\n\tNumSets 1,\n\tContentType \"LandmarkSet\"\n}",
    "\n\nMarkers { float[3] Coordinates } @1",
    "\n\n# Data section follows\n@1\n"
  )

  for (lm in 1:n_landmarks) {
    line <- paste(
      landmarks[lm,1], landmarks[lm, 2], landmarks[lm, 3]
    )
    output <- paste0(
      output,
      line, " \n"
    )
  }
  output <- paste0(output, "\n")

  cat(output, file = paste(file_name))
}

mahalanobis_distances <- function(X, robust = TRUE) {
  distances <- c()
  if (robust != TRUE) {
    M <- colMeans(X)
  } else {
    M <- robustbase::colMedians(X)
  }
  c_m1 <- MASS::ginv(cov(X))
  for (i in 1:dim(X)[1]) {
    x <- X[i,]
    d2 <- t(x - M) %*% c_m1 %*% (x - M)
    distances <- c(distances, d2)
  }
  return(distances)
}

mahalanobis_d2 <- function(x, m, robust = TRUE) {
  if (robust != TRUE) {
    mu <- colMeans(m)
  } else {
    mu <- robustbase::colMedians(m)
  }
  c_m1 <- MASS::ginv(cov(m))
  d2 <- t(x - mu) %*% c_m1 %*% (x - mu)
  return(d2)
}

mahalanobis_distances_to_target <- function(from, to) {
  distances <- c()
  M <- robustbase::colMedians(to)
  c_m1 <- MASS::ginv(cov(to))
  for (i in 1:dim(from)[1]) {
    x <- from[i,]
    d2 <- t(x - M) %*% c_m1 %*% (x - M)
    distances <- c(distances, d2)
  }
  return(distances)
}

get_graph_df <- function(target_central) {
  if (dim(target_central)[2] == 2) {
    graph_df <- data.frame(id = rep("LM", dim(target_central)[1]))
    graph_df$x <- target_central[,1]
    graph_df$y <- target_central[,2]
    graph_df$id <- as.character(graph_df$id)
    for(id in 1:nrow(graph_df)) {
      graph_df[id,][["id"]] <- paste("LM", id, sep = "")
    }
  }  else {
    graph_df <- data.frame(id = rep("LM", dim(target_central)[1]))
    graph_df$x <- target_central[,1]
    graph_df$y <- target_central[,2]
    graph_df$z <- target_central[,3]
    graph_df$id <- as.character(graph_df$id)
    for(id in 1:nrow(graph_df)) {
      graph_df[id,][["id"]] <- paste("LM", id, sep = "")
    }
  }
  return(graph_df)
}

# experimental functions ---------------------

lmconfig_isomorphism <- function(central1, edges1, central2, edges2) {

  `%!in%` = Negate(`%in%`)

  if("array" %!in% class(central1) |
     "array" %!in% class(edges1) |
     "array" %!in% class(central2) |
     "array" %!in% class(edges2)) {
    stop("Landmark central configurations must be a (landmark x dimension) array")
  }

  if(dim(central1)[2] != 2 & dim(central1)[2] != 3) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if(dim(central2)[2] != 2 & dim(central2)[2] != 3) {
    stop("Landmark central configuration must be a (landmark x dimension) array")
  }

  if (dim(central1)[2] == 2) {
    graph_df_1 <- data.frame(id = rep("LM", dim(central1)[1]))
    graph_df_1$x <- central1[,1]
    graph_df_1$y <- central1[,2]
    graph_df_1$id <- as.character(graph_df_1$id)
    for(id in 1:nrow(graph_df_1)) {
      graph_df_1[id,][["id"]] <- paste("LM", id, sep = "")
    }
  }  else {
    graph_df_1 <- data.frame(id = rep("LM", dim(central1)[1]))
    graph_df_1$x <- central1[,1]
    graph_df_1$y <- central1[,2]
    graph_df_1$z <- central1[,3]
    graph_df_1$id <- as.character(graph_df_1$id)
    for(id in 1:nrow(graph_df_1)) {
      graph_df_1[id,][["id"]] <- paste("LM", id, sep = "")
    }
  }

  if (dim(central2)[2] == 2) {
    graph_df_2 <- data.frame(id = rep("LM", dim(central2)[1]))
    graph_df_2$x <- central2[,1]
    graph_df_2$y <- central2[,2]
    graph_df_2$id <- as.character(graph_df_2$id)
    for(id in 1:nrow(graph_df_2)) {
      graph_df_2[id,][["id"]] <- paste("LM", id, sep = "")
    }
  }  else {
    graph_df_2 <- data.frame(id = rep("LM", dim(central2)[1]))
    graph_df_2$x <- central2[,1]
    graph_df_2$y <- central2[,2]
    graph_df_2$z <- central2[,3]
    graph_df_2$id <- as.character(graph_df_2$id)
    for(id in 1:nrow(graph_df_2)) {
      graph_df_2[id,][["id"]] <- paste("LM", id, sep = "")
    }
  }

  LM_graph_1 <- igraph::graph_from_data_frame(vertices = graph_df_1,
                                              d = as.matrix(edges1),
                                              directed = FALSE)

  LM_graph_2 <- igraph::graph_from_data_frame(vertices = graph_df_2,
                                              d = as.matrix(edges2),
                                              directed = FALSE)

  results <- igraph::isomorphic(LM_graph_1, LM_graph_2,
                                method = c("bliss"))
  #edges_in_common <- length(intersect(igraph::E(LM_graph_1),
  #                                    igraph::E(LM_graph_2)))

  edges_in_common <- nrow(intersect(igraph::as_data_frame(LM_graph_1),
                                    igraph::as_data_frame(LM_graph_2)))

  edges_1 <- length(igraph::E(LM_graph_1))
  edges_2 <- length(igraph::E(LM_graph_2))
  pos_edges <- (edges_1 + edges_2) / 2

  if (results == TRUE) {
    cat(
      paste0(
        "Isomorphism results:\n",
        "The two landmark configurations were found to present isomorphism; \ni.e. configurations are equal.",
        "\n", edges_in_common, " edges were found to be present in both configurations, accounting for ",
        "100% of structural similarities."
      )
    )
  } else {
    cat(
      paste0(
        "Isomorphism results:\n",
        "The two landmark configurations were not found to present isomorphism; \ni.e. configurations are structurally different.",
        "\n", edges_in_common, " edges were found to be present in both configurations, accounting for \n",
        round((edges_in_common / pos_edges) * 100, 2), "% of structural similarities."
      )
    )

  }
  if((edges_in_common / pos_edges) == 0) {
    message(cat("\n\nWARNING: If landmark configurations present dense semilandmark configurations,",
                "0% similarities are common. It may be worth testing for isomorphism on fixed landmarks only."))
  }

}

procD_from_graphSim <- function(lm_matrix, labels, plot_results = TRUE,
                                robust = TRUE, permutations = 1000,
                                method = "cosine") {

  if(!is.array(matrices)) {
    stop("Input must be a (n x n x individual) array")
  }

  if(missing(lm_matrix) | missing(labels)) {
    stop("Both a landmark array and set of labels need to be specified for this function")
  }

  if (permutations <= 0 | permutations %% 1 != 0) {
    stop(
      "The number of permutations must be a non-negative integer"
    )
  }

  level_list = levels(as.factor(labels))
  cat("\n\nCalculating large similarity matrices from graphs")
  matrices = similarity_matrix(lm_matrix, method = "cosine", unravel = FALSE)
  cat("\n\n")

  dist_matrix = array(numeric(), dim = c(length(level_list), length(level_list)),
                      dimnames = list(level_list, level_list))
  p_value_matrix = array(numeric(), dim = c(length(level_list), length(level_list)),
                         dimnames = list(level_list, level_list))
  reference_matrices = array(numeric(), dim = c(dim(matrices)[1], dim(matrices)[2], 0))

  for (sample in level_list) {
    mshape <- array(numeric(), c(dim(matrices)[1], dim(matrices)[2]))
    target_sample <- separate_samples(matrices, as.factor(labels), sample)
    if (robust == TRUE) {
      for (rows in 1:dim(target_sample)[1]) {
        for (cols in 1:dim(target_sample)[2]) {
          mshape[rows,cols] = sum(target_sample[rows,cols,] / dim(target_sample)[3])
        }
      }
    } else {
      for (rows in 1:dim(target_sample)[1]) {
        for (cols in 1:dim(target_sample)[2]) {
          mshape[rows,cols] = median(target_sample[rows,cols,])
        }
      }
    }
    reference_matrices <- abind::abind(reference_matrices, mshape, along = 3)
  }

  if (plot_results != FALSE) {
    dev.new(); par(mfrow = c(length(level_list), length(level_list)))
  }

  cat("Calculating procrustes distances")
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
        #c = 0

        index_1 <- table(compare$Sample)[1]
        index_2 <- table(compare$Sample)[2]

        #if (robust == TRUE) {
        #  o1 <- median(compare[1:index_1, 1])
        #  o2 <- median(compare[index_1:index_2, 1])
        #} else {
        #  o1 <- mean(compare[1:index_1, 1])
        #  o2 <- mean(compare[index_1:index_2, 1])
        #}

        #t_base<-abs(o1-o2)
        #for (perm in 1:permutations) {
        #  perm <- compare[sample(nrow(compare)), 1]
        #  X1 <- perm[1:index_1]
        #  X2 <- perm[index_1:index_2]
        #  if (robust == TRUE) {
        #    m1 <- median(X1)
        #    m2 <- median(X2)
        #  } else {
        #    m1 <- mean(X1)
        #    m2 <- mean(X2)
        #  }
        #  t_perm <- abs(m1 - m2)
        #  if (t_perm > t_base) {
        #    c = c + 1
        #  }

        #}

        #p_val <- c / permutations

        if (robust == TRUE) {
          p_val <- kruskal.test(compare$proc_dist, compare$Sample)$"p.value"
          dist = median(compare[index_1:index_2,1])
        } else {
          p_val = summary(aov(compare$proc_dist ~ compare$Sample))[[1]]$`Pr(>F)`
          dist = mean(compare[index_1:index_2,1])
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

  return(list(distances = as.dist(t(dist_matrix)),
              p_values = as.dist(t(p_value_matrix))))
}
