
find_centroid <- function(configuration, robust = FALSE) {
  if (robust == FALSE) {
    centroid <- apply(configuration, 2, mean)
  } else {
    centroid <- apply(configuration, 2, median)
  }
  return(centroid)
}

euclid <- function(point_1, point_2) {
  return(sqrt(sum((point_1 - point_2)^2)))
}

centroid_size <- function(configuration, centroid) {
  distances <- c()
  for (lm in 1:nrow(configuration)) {
    distances <- c(distances, euclid(configuration[lm,], centroid)^2)
  }
  return(sqrt(sum(distances)))
}

kendall_translation <- function(configuration, robust = FALSE) {
  translated_matrix <- configuration - matrix(
    find_centroid(configuration, robust = robust),
    nrow(configuration), ncol(configuration),
    byrow = TRUE
  )
  return(translated_matrix)
}

helmert <- function(p, type) {

  Helm_M <- matrix(0, p, p + 1)

  if (type == 1) {
    val <- 1
    while(val <= p) {
      val2 <- 1
      while(val2 <= val) {
        Helm_M[val, val2] <- -1 / sqrt(val * (val + 1))
        val2 <- val2 + 1
      }
      Helm_M[val, val + 1] <- val / sqrt(val * (val + 1))
      val <- val + 1
    }
  } else {
    if (p > 0) {
      for (val in seq(1, p)) {
        value = -1 / sqrt(val * (val + 1))
        Helm_M[val, seq(1, val)] = value
        Helm_M[val, val + 1] = - val * value
      }
    }
  }

  return(Helm_M)
}

helmert_M <- function(configuration) {
  Helm_M <- helmert(nrow(configuration))[-1,]
  return(Helm_M %*% configuration)
}

kendall_2D_coord <- function(configuration) {
  helmert_configuration <- helmert_M(configuration)
  kendall_coordinates <- complex(nrow(helmert_configuration),
                                 helmert_configuration[,1],
                                 helmert_configuration[,2])
  kendall_coordinates <- (kendall_coordinates / kendall_coordinates[1])[-1]
  return(kendall_coordinates)
}

kendall_full_superimposition <- function(config_1, config_2,
                                              scale = TRUE, robust = FALSE) {

  k <- ncol(config_1) # k = dimensions

  # center and scale

  centroid_1 <- find_centroid(config_1, robust = robust)
  centroid_2 <- find_centroid(config_2, robust = robust)
  cs1 <- centroid_size(config_1, centroid_1)
  cs2 <- centroid_size(config_2, centroid_2)

  if (scale == TRUE) {
    scaled_config_1 <- config_1 / cs1
    scaled_config_2 <- config_2 / cs2
    translated_matrix_1 <- kendall_translation(scaled_config_1,
                                               robust = robust)
    translated_matrix_2 <- kendall_translation(scaled_config_2,
                                               robust = robust)
  } else {
    translated_matrix_1 <- kendall_translation(config_1,
                                               robust = robust)
    translated_matrix_2 <- kendall_translation(config_2,
                                               robust = robust)
  }

  # remove remaining scale via single value decomposition

  eigendecomp <- svd(t(translated_matrix_2) %*% translated_matrix_1)

  U_prima <- eigendecomp$v
  V <- eigendecomp$u

  if (scale == TRUE) {
    Delta <- eigendecomp$d
  }

  # remove reflection

  sign_det <- sign(det(t(translated_matrix_2) %*% translated_matrix_1))

  if (scale == TRUE) {
    Delta[k] <- sign_det * abs(Delta[k])
  }

  V[,k] <- sign_det * V[,k]

  Gamma <- U_prima %*% t(V)

  # best fit parameter

  if (scale == TRUE) {
    beta <- sum(Delta)
  }

  # superimposed coordinates

  if (scale == TRUE) {
    M1 <- beta * translated_matrix_1 %*% Gamma
  } else {
    M1 <- translated_matrix_1 %*% Gamma
  }

  M2 <- translated_matrix_2

  # parameters

  if (scale == TRUE) {
    proc_d <- sqrt(1 - beta^2) # full procrustes distance
  } else {
    proc_d <- sqrt(apply((M1 - M2)^2, 1, sum))
  }

  rotation <- Gamma # rotation matrix
  if (scale == TRUE) {
    scalar <- beta # scalar value
    rho <- acos(beta)
  } else {
    scalar <- NA
    rho <- NA
  }

  kendall_coordinates <- list(
    config_1_prima = M1,
    config_2_prima = M2,
    ProcD = proc_d,
    Gamma = Gamma,
    scalar = scale,
    rho = rho,
    centroid_sizes = c(cs1, cs2)
  )

  return(kendall_coordinates)

}

riemannianD <- function(matrix1, matrix2, robust = FALSE) {

  dimensions <- dim(matrix1)[2]

  if (sum((matrix1 - matrix2) ** 2) == 0) {
    rho = 0
  } else {
    if (dimensions < 3) {

      rho <- c(
        acos(
          min(
            1, (
              Mod(t(Conj(kendall_preshape(matrix1))) %*% kendall_preshape(matrix2))
            )
          )
        )
      )

    } else {

      k <- dim(matrix1)[2]
      matrix_x <- kendall_preshape(matrix1)
      matrix_y <- kendall_preshape(matrix2)
      Q <- t(matrix_x) %*% matrix_y %*% t(matrix_y) %*% matrix_x
      eigenvalues <- eigen(t(matrix_x) %*% matrix_y)$values

      check <- 1
      for (iteration in 1:k) {
        check <- check * eigenvalues[iteration]
      }

      eigenvalues <- sqrt(abs(eigen(Q, symmetric = TRUE)$values))
      if (Re(check) < 0) {
        eigenvalues[k] <- -eigenvalues[k]
      }
      rho <- acos(min(sum(eigenvalues), 1))

    }
  }

  return(rho)

}

kendall_preshape <- function(target_matrix, robust = FALSE) {

  p <- dim(target_matrix)[1]
  if (p < 100) {
    h_type = 1
  }  else {
    h_type = 2
  }
  h <- helmert(p - 1, h_type)

  preshape_x <- h %*% target_matrix
  centroid <- find_centroid(target_matrix,
                            robust = robust)
  CS <- centroid_size(target_matrix,
                      centroid)
  preshape_x <- preshape_x / CS

  return(preshape_x)

}

dryden_pretransform <- function(lm_array) {

  # function based on code from the "shapes" package

  pretransform_fun_1 <- function(n_value) {
    value_1 <- matrix(1:1, n_value, n_value)
    return(diag(n_value) - (1 / n_value) * value_1)
  }

  pretransform_fun_2 <- function(lm_array_2) {
    n <- dim(lm_array_2)[1]
    return(pretransform_fun_1(n) %*% lm_array_2)
  }

  transformed_array <- apply(lm_array, 3, pretransform_fun_2)
  transformed_array <- array(transformed_array, dim(lm_array))
  return(transformed_array)

}

dryden_rotation <- function(lm_array, robust = FALSE) {

  # function based on code from the "shapes" package

  transformed_array <- list(
    rotated = 0,
    dif = 0,
    r_no = 0,
    inc = 0
  )

  array_to_transform <- lm_array

  dryden_l <- dim(lm_array)[3]
  dryden_d <- 0; dryden_d[1] <- 1e+12; dryden_d[2] <- dryden_dif(array_to_transform,
                                                                 robust = robust)
  dryden_a <- 0; dryden_a[1] <- dryden_d[2]
  dryden_n <- 0
  dryden_s <- dryden_add(array_to_transform)

  tolerance <- dryden_dif(array_to_transform, robust = robust)
  if (tolerance > 1e-05) {
    while(dryden_d[1] - dryden_d[2] > 1e-05) {

      dryden_n <- dryden_n + 1
      dryden_d[1] <- dryden_d[2]

      for (individual in 1:dryden_l) {
        previous_step <- array_to_transform[,,individual]
        array_to_transform[,,individual] <- previous_step %*% kendall_rotation(
          ((1 / (dryden_l - 1)) * (dryden_s - previous_step)), previous_step
        )
        dryden_s <- dryden_s - previous_step + array_to_transform[,,individual]
      }

      dryden_d[2] <- dryden_dif(array_to_transform)
      dryden_a[dryden_n + 1] <- dryden_d[2]

    }
  }

  transformed_array$rotated <- array_to_transform
  transformed_array$dif <- dryden_a
  transformed_array$r_no <- dryden_n
  transformed_array$inc <- dryden_d[1] - dryden_d[2]

  return(transformed_array)

}

kendall_rotation <- function(matrix1, matrix2) {

  k <- dim(matrix1)[2]

  eigendecomp <- svd(t(matrix1) %*% matrix2)

  U_prima <- eigendecomp$v
  V <- eigendecomp$u
  Delta <- eigendecomp$d

  # remove reflection

  sign_det <- sign(det(t(matrix1) %*% matrix2))

  Delta[k] <- sign_det * abs(Delta[k])
  V[,k] <- sign_det * V[,k]

  Gamma <- U_prima %*% t(V)

  return(Gamma)

}

# function based on code from the "shapes" package

dryden_add <- function(lm_array) {

  sum_value <- 0
  for (individual in 1:dim(lm_array)[3]) {
    sum_value <- sum_value + lm_array[,,individual]
  }
  return(sum_value)

}

# function based on code from the "shapes" package

dryden_dif <- function(lm_array, robust = FALSE) {

  # centroid_size
  centroid <- find_centroid(lm_array, robust = robust)
  CS <- centroid_size((dryden_add(lm_array) / dim(lm_array)[3]), centroid)
  if (robust == TRUE) {
    X <- sweep(lm_array, c(1,2), apply(lm_array, c(1,2), median))
  } else {
    X <- sweep(lm_array, c(1,2), apply(lm_array, c(1,2), mean))
  }

  enorm <- sqrt(sum((as.vector(X) / CS)^2))

  return(enorm^2 / dim(lm_array)[3])


}

# function based on code from the "shapes" package

dryden_scale <- function(lm_array) {

  coeff <- dryden_b_coefficients(lm_array)
  iteration <- rep(dim(lm_array)[1] * dim(lm_array)[2], dim(lm_array)[3])
  iteration_sequence <- rep(coeff, iteration)

  scaled_data <- array(as.vector(lm_array) * iteration_sequence, dim(lm_array))

  return(scaled_data)

}

# function based on code from the "shapes" package

dryden_b_coefficients <- function(lm_array) {

  dryden_vec_function <- function(X) {
    return(
      matrix(X, dim(X)[1] * dim(X)[2], dim(X)[3])
    )
  }

  dryden_h <- 0
  transformed_array <- lm_array

  n <- dim(lm_array)[3]

  transformations <- apply(
    transformed_array, 3,
    function(x) {
      sqrt(sum(x^2))
    }
  )

  sum_transformations <- 0
  sum_transformations <- sum(transformations)

  o_matrix <- t(dryden_vec_function(transformed_array))
  o_k <- dim(o_matrix)[2]
  o_n <- dim(o_matrix)[1]
  if (o_n > o_k) {

    dryden_q <- rep(0, o_n)
    for (iteration in 1:n) {
      dryden_q[iteration] <- var(o_matrix[iteration, ]) * (n - 1) / n
      o_matrix[iteration,] <- o_matrix[iteration, ] - mean(o_matrix[iteration, ])
    }

    o_matrix <- diag(sqrt(1 / dryden_q)) %*% o_matrix
    n <- o_k
    L_matrix <- t(o_matrix) %*% o_matrix / n

    eigenvalues <- eigen(L_matrix, symmetric = TRUE)

    U <- eigenvalues$vectors
    lambda <- eigenvalues$values

    V <- o_matrix %*% U

    dryden_v <- rep(0, n)
    for (iteration in 1:n) {
      dryden_v[iteration] <- sqrt(t(V[, iteration]) %*% V[, iteration])
      V[, iteration] <- V[, iteration] / dryden_v[iteration]
    }

    delta <- sqrt(abs(lambda / n)) ** dryden_v
    delta <- delta[order(delta, decreasing = TRUE)]
    V <- V[, order(delta, decreasing = TRUE)]

    coefs <- sqrt(sum_transformations / transformations) * V[, 1]

  } else if (o_k >= o_n) {

    dryden_z <- cor(dryden_vec_function(transformed_array))
    coefs <- sqrt(sum_transformations / transformations) * eigen(dryden_z)$vectors[, 1]

  }

  return(abs(coefs))

}

resistant_median_size <- function(lm_array) {

  dist_matrix <- as.matrix(dist(lm_array))
  rms <- median(apply(dist_matrix, 2, median, na.rm = TRUE))
  return(rms)

}

resistant_rotation_matrix <- function(config_1, config_2) {

  p <- dim(config_1)[1]
  rotation_matrix <- matrix(NA, p, p)

  for (i in 1:p) {
    for (j in 1:p) {
      rotation_matrix[i,j] <- Arg(
        complex(1, config_2[i, 1], config_2[i, 2]) -
          complex(1, config_2[j, 1], config_2[j, 2])
      ) - Arg(
        complex(1, config_1[i, 1], config_1[i, 2]) -
          complex(1, config_1[j, 1], config_1[j, 2])
      )
    }
  }

  return((rotation_matrix + pi %% (2 * pi)) - pi)

}

resistant_2D_ordinary_superimposition <- function(config_1, config_2, scale = TRUE) {

  p <- dim(config_1)[1]
  k <- dim(config_1)[2]

  original_fsup_coords <- kendall_full_superimposition(config_1,
                                                       config_2,
                                                       scale = scale,
                                                       robust = FALSE)
  x1 <- original_fsup_coords$config_1_prima
  x2 <- original_fsup_coords$config_2_prima

  x <- as.matrix(dist(x2)) / as.matrix(dist(x1))
  beta <- median(apply(x, 2, median, na.rm = TRUE))

  rotation_matrix <- resistant_rotation_matrix(x1, x2)
  phi <- median(apply(rotation_matrix, 2, median, na.rm = TRUE))

  Gamma <- matrix(c(cos(phi),
                    -sin(phi),
                    sin(phi),
                    cos(phi)), 2, 2)

  alpha <- x2 - beta * x1 %*% Gamma
  alpha <- matrix(apply(alpha, 2, median), p, k, byrow = TRUE)

  x1 <- beta * x1 %*% Gamma + alpha

  return(list(
    config_1_prima = x1,
    config_2_prima = x2
  ))

}

resistant_3D_rotation <- function(config_1, config_2) {

  p <- dim(config_1)[1]
  k <- dim(config_1)[2]

  identity_matrix <- diag(1, 3)
  angles_array <- array(NA, dim = c(p, p, p, 3, 3))

  for (index_i in 1:p) {
    for (index_j in (1:p)[-index_i]) {
      for (index_k in (1:p)[-c(index_i, index_j)]) {

        M1i <- matrix(config_1[index_i,])
        M1j <- matrix(config_1[index_j,])
        M1k <- matrix(config_1[index_k,])

        M2i <- matrix(config_2[index_i,])
        M2j <- matrix(config_2[index_j,])
        M2k <- matrix(config_2[index_k,])

        M1ij <- (M1j - M1i) / sqrt(sum((M1j - M1i)^2))
        M2ij <- (M2j - M2i) / sqrt(sum((M2j - M2i)^2))

        M1ijk <- (
          (M1k - M1i) - as.numeric(
            t(M1k - M1i) %*% M1ij
          ) * M1ij
        ) / sqrt(
          sum(
            ((M1k - M1i) - as.numeric(
              t(M1k - M1i) %*% M1ij
            ) * M1ij)^2
          )
        )

        M2ijk <- (
          (M2k - M2i) - as.numeric(
            t(M2k - M2i) %*% M2ij
          ) * M2ij
        ) / sqrt(
          sum(
            ((M2k - M2i) - as.numeric(
              t(M2k - M2i) %*% M2ij
            ) * M2ij)^2
          )
        )

        M1s <- matrix(
          c(
            M1ij[2] * M1ijk[3] - M1ij[3] * M1ijk[2],
            M1ij[3] * M1ijk[1] - M1ij[1] * M1ijk[3],
            M1ij[1] * M1ijk[2] - M1ij[2] * M1ijk[1]
          )
        )

        M2s <- matrix(
          c(
            M2ij[2] * M2ijk[3] - M2ij[3] * M2ijk[2],
            M2ij[3] * M2ijk[1] - M2ij[1] * M2ijk[3],
            M2ij[1] * M2ijk[2] - M2ij[2] * M2ijk[1]
          )
        )

        Y1 <- cbind(M1ij, M1ijk, M1s)
        Y2 <- cbind(M2ij, M2ijk, M2s)

        Gamma <- t(Y1) %*% Y2

        angles_array[index_i, index_j, index_k, , ] <- S <- solve(
          Gamma + identity_matrix
        ) %*% (
          Gamma - identity_matrix
        )

      }
    }
  }

  S <- apply(
    apply(
      apply(
        angles_array, c(1, 2, 4, 5), median, na.rm = TRUE
      ), c(1, 3, 4), median, na.rm = TRUE
    ), 2:3, median, na.rm = TRUE
  )
  Gamma <- (identity_matrix + S) %*% solve(identity_matrix - S)
  return(list(
    reference = config_2,
    target = config_1 %*% t(Gamma),
    Gamma = Gamma
  ))

}

resistant_3D_superimposition <- function(config_1, config_2, scale = TRUE) {

  p <- dim(config_1)[1]
  k <- dim(config_1)[2]

  original_fsup_coords <- kendall_full_superimposition(config_1,
                                                       config_2,
                                                       scale = scale,
                                                       robust = FALSE)

  x1 <- original_fsup_coords$config_1_prima
  x2 <- original_fsup_coords$config_2_prima

  x <- as.matrix(dist(x2)) / as.matrix(dist(x1))

  reference_matrix <- t(x)
  reference_matrix[row(x) < col(x)] <- x[row(x) < col(x)]
  reference_matrix <- as.matrix(dist(x2)) / as.matrix(dist(x1))

  beta <- median(apply(reference_matrix, 2, median, na.rm = TRUE))

  Gamma <- resistant_3D_rotation(x1, x2)$Gamma
  alpha <- x2 - beta * x1 %*% Gamma
  alpha <- matrix(apply(alpha, 2, median), p, k, byrow = TRUE)

  return(list(
    preshape_1 = beta * x1 %*% Gamma + alpha,
    preshape_2 = x2
  ))

}
