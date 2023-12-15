############################################################################
###                                                                      ###
###              Multivariate t Mixture for Incomplete Data              ###
###                                                                      ###
############################################################################

MtM_incomplete_data <- function(
    X,
    G,
    max_iter       = 20,
    epsilon        = 0.01,
    init_method    = c("kmedoids", "kmeans", "hierarchical", "manual"),
    clusters       = NULL,
    outlier_cutoff = 0.95,
    progress       = TRUE
) {

  #----------------------#
  #    Input Checking    #
  #----------------------#

  if (outlier_cutoff <= 0 | outlier_cutoff >= 1) {
    stop('outlier_cutoff must be in (0, 1)')
  }

  #------------------------------------#
  #    Objects for the EM Algorithm    #
  #------------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  do <- rowSums(!is.na(X))
  R  <- is.na(X)
  M  <- unique(R)
  np <- nrow(M)

  Im <- vector('list', np)    # which observations with missing pattern j

  for (j in 1:np) {
    Im[[j]] <- which( apply(R, 1, function(r) all(r == M[j, ]) ) )
  }

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  w_tilde     <- matrix(NA, nrow = n, ncol = G)
  zw_tilde    <- matrix(NA, nrow = n, ncol = G)
  X_tilde     <- array(rep(X, G), dim = c(n, d, G))
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  py    <- rep(NA, G)
  mu    <- matrix(NA, nrow = G, ncol = d)
  Sigma <- array(NA, dim = c(d, d, G))
  df    <- rep(10, G)

  dens   <- matrix(NA, nrow = n, ncol = G)
  iter   <- 0
  loglik <- NULL

  #--------------------------------#
  #    Parameter Initialization    #
  #--------------------------------#

  init_method <- match.arg(init_method)

  if (progress) {
    cat('\nInitialization:', init_method, '\n')
  }

  X_imp <- mean_impute(X)

  if (G == 1) {

    max_iter <- 1

    pars <- cluster_pars(X = X_imp, clusters = rep(1, n))

    py    <- 1
    mu    <- pars$mu
    Sigma <- pars$Sigma

  } else {

    init <- initialize_clusters(
      X           = X_imp,
      G           = G,
      init_method = init_method,
      clusters    = clusters
    )

    py    <- init$pi
    mu    <- init$mu
    Sigma <- init$Sigma

  }

  #---------------------#
  #  The EM algorithm   #
  #---------------------#

  if (progress) {
    cat('\nModel Fitting:\n')
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3, width = 75, char = "=")
  }

  while (iter < max_iter & getall(loglik) > epsilon) {

    #++++ E-step ++++#

    for (g in 1:G) {
      for (j in 1:np) {

        m    <- M[j, ]                         # missing pattern j
        o    <- !m                             # observed pattern j
        Xo_j <- X[Im[[j]], o, drop = FALSE]    # observations with missing pattern j

        mu_o     <- mu[g, o]
        Sigma_oo <- Sigma[o, o, g]

        z[Im[[j]], g]       <- mvtnorm::dmvt(Xo_j, delta = mu_o, sigma = as.matrix(Sigma_oo), df = df[g], log = FALSE)
        w_tilde[Im[[j]], g] <- (df[g] + sum(o)) / (df[g] + mahalanobis(Xo_j, mu_o, Sigma_oo, tol = 1e-20))

      }
    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    for (g in 1:G) {

      zw_tilde[, g] <- z_tilde[, g] * w_tilde[, g]

      for (j in 1:np) {

        m <- M[j, ]    # missing pattern j

        if (any(m)) {

          o <- !m        # observed pattern j

          mu_m <- mu[g, m]
          mu_o <- mu[g, o]

          Sigma_oo     <- Sigma[o, o, g]
          Sigma_mo     <- Sigma[m, o, g]
          Sigma_oo_inv <- mnormt::pd.solve(Sigma_oo)

          for (i in Im[[j]]) {
            xi <- X[i, ]

            x_ig_tilde       <- mu_m + Sigma_mo %*% Sigma_oo_inv %*% (xi[o] - mu_o)
            X_tilde[i, m, g] <- x_ig_tilde
          }

        }

      }

    }

    #++++ M-step: pi ++++#

    N  <- colSums(z_tilde)
    py <- N / n

    #++++ M-step: mu ++++#

    for (g in 1:G) {
      mu_num  <- colSums(z_tilde[, g] * w_tilde[, g] * X_tilde[, , g])
      mu_den  <- sum(z_tilde[, g] * w_tilde[, g])
      mu[g, ] <- mu_num / mu_den
    }

    #++++ M-step: Prepare Sigma tilde ++++#

    for (g in 1:G) {
      for (j in 1:np) {

        m <- M[j, ]    # missing pattern j

        if (any(m)) {

          o <- !m        # observed pattern j

          mu_m <- mu[g, m]
          mu_o <- mu[g, o]

          Sigma_oo     <- Sigma[o, o, g]
          Sigma_om     <- Sigma[o, m, g]
          Sigma_mo     <- Sigma[m, o, g]
          Sigma_mm     <- Sigma[m, m, g]
          Sigma_oo_inv <- mnormt::pd.solve(Sigma_oo)

          S_mm <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% Sigma_om

          for (i in Im[[j]]) {
            xi <- X[i, ]

            Sigma_tilde[o, o, i, g] <- tcrossprod(xi[o] - mu_o)
            Sigma_tilde[o, m, i, g] <- tcrossprod(xi[o] - mu_o, X_tilde[i, m, g] - mu_m)
            Sigma_tilde[m, o, i, g] <- t(Sigma_tilde[o, m, i, g])
            Sigma_tilde[m, m, i, g] <- tcrossprod(X_tilde[i, m, g] - mu_m) + S_mm / w_tilde[i, g]
          }

        } else {

          X_centrd <- sweep(X[Im[[j]], ], 2, mu[g, ], '-')
          cr_prods <- apply(X_centrd, 1, tcrossprod)

          Sigma_tilde[, , Im[[j]], g] <- array(
            data = unlist(cr_prods),
            dim  = c(d, d, length(Im[[j]]))
          )

        }

      }
    }


    for (g in 1:G) {

      #++++ M-step: Sigma ++++#

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(zw_tilde[slc_ind, g] * Sigma_tilde[, , ,g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ M-step: Degrees of Freedom ++++#

      root <- tryCatch({

        uniroot(function(a) {

          A <- -digamma(a / 2) + log(a / 2) + 1
          B <-  sum( z_tilde[, g] * (log(w_tilde[, g]) - w_tilde[, g]) ) / N[g]
          C <-  sum( z_tilde[, g] * (digamma((df[g] + do) / 2) - log((df[g] + do) / 2)) ) / N[g]

          return(A + B + C)

        }, lower = 0.001, upper = 200)$root

      }, error = function(e) return(NULL))

      df[g] <- ifelse(!is.null(root), root, df[g])

    }

    #++++ Observed Log-Likelihood ++++#

    for (g in 1:G) {
      for (j in 1:np) {

        m    <- M[j, ]                         # missing pattern j
        o    <- !m                             # observed pattern j
        Xo_j <- X[Im[[j]], o, drop = FALSE]    # observations with missing pattern j

        mu_o     <- mu[g, o]
        Sigma_oo <- Sigma[o, o, g]

        dens[Im[[j]], g] <- mvtnorm::dmvt(Xo_j, delta = mu_o, sigma = as.matrix(Sigma_oo), df = df[g], log = FALSE)

      }
    }

    lik                   <- dens %*% py
    lik[lik <= 10^(-323)] <- 10^(-323)
    final_loglik          <- sum(log(lik))
    loglik                <- c(loglik, final_loglik)

    #++++ Update Progress ++++#

    iter <- iter + 1

    if (progress) {
      setTxtProgressBar(pb, iter)
    }

  }

  if (progress) {
    close(pb)
    if (iter < max_iter) {
      cat('\nConvergence was reached before', max_iter, 'iterations\n')
    }
  }

  #---------------------------#
  #    Cluster Memberships    #
  #---------------------------#

  clusters <- apply(z_tilde, 1, which.max)

  #------------------#
  #    Imputation    #
  #------------------#

  X_imputed <- X
  complete  <- complete.cases(X)

  for (i in which(!complete)) {
    X_imputed[i, ] <- X_tilde[i, , clusters[i]]
  }

  #-------------------------#
  #    Outlier Detection    #
  #-------------------------#

  delta <- matrix(NA, nrow = n, ncol = G)

  for (g in 1:G) {
    delta[, g] <- mahalanobis(X_imputed, center = mu[g, ], cov = Sigma[, , g], tol = 1e-20)
  }

  cluster_matr <- clusters_to_matrix(clusters, G)
  delta        <- rowSums(delta * cluster_matr)
  outliers     <- 1 - pchisq(delta, df = d) < 1 - outlier_cutoff

  #------------------------#
  #  Number of parameters  #
  #------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G * d,
    Sigma = G * d * (d + 1) / 2,
    df    = G
  )
  npar$total <- Reduce('+', npar)

  #----------------------------#
  #    Information Criteria    #
  #----------------------------#

  AIC <- -2 * final_loglik + 2 * npar$total
  BIC <- -2 * final_loglik + npar$total * log(n)

  KIC  <- -2 * final_loglik + 3 * (npar$total + 1)
  KICc <- -2 * final_loglik + 2 * (npar$total + 1) * n/(n-npar$total -2) - n * digamma((n-npar$total)/2) + n * log(n/2)

  AIC3 <- -2 * final_loglik + 3 * npar$total
  CAIC <- -2 * final_loglik + npar$total * (1 + log(n))
  AICc <- -2 * final_loglik + 2 * npar$total * n/(n - npar$total - 1)

  ent <- apply(z_tilde, 1, max)
  ICL <- BIC - sum(ent * log(ent))

  AWE <- -2 * (final_loglik + sum(ent * log(ent))) + 2 * npar$total * (3/2 + log(n))
  CLC <- -2 * final_loglik + 2 * sum(ent * log(ent))

  #----------------------#
  #    Prepare Output    #
  #----------------------#

  c_names <- paste('comp', 1:G, sep = '')
  v_names <- colnames(X)

  if (is.null(v_names)) {
    v_names <- 1:d
  }

  names(py)       <- c_names
  rownames(mu)    <- c_names
  colnames(mu)    <- v_names
  dimnames(Sigma) <- list(v_names, v_names, c_names)
  names(df)       <- c_names

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
  }

  outputs <- list(
    model         = 't_incomplete_data',
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    df            = df,
    z_tilde       = z_tilde,
    clusters      = clusters,
    outliers      = outliers,
    data          = X_imputed,
    complete      = complete,
    npar          = npar,
    max_iter      = max_iter,
    iter_stop     = iter,
    final_loglik  = final_loglik,
    loglik        = loglik,
    AIC           = AIC,
    BIC           = BIC,
    KIC           = KIC,
    KICc          = KICc,
    AIC3          = AIC3,
    CAIC          = CAIC,
    AICc          = AICc,
    ent           = ent,
    ICL           = ICL,
    AWE           = AWE,
    CLC           = CLC,
    init_method   = init_method
  )
  class(outputs) <- 'MixtureMissing'

  return(outputs)

}

####################################################################
###                                                              ###
###           Multivariate t Mixture for Complete Data           ###
###                                                              ###
####################################################################

MtM_complete_data <- function(
    X,
    G,
    max_iter       = 20,
    epsilon        = 0.01,
    init_method    = c("kmedoids", "kmeans", "hierarchical", "manual"),
    clusters       = NULL,
    equal_prop     = FALSE,
    identity_cov   = FALSE,
    outlier_cutoff = 0.95,
    progress       = TRUE
) {

  #----------------------#
  #    Input Checking    #
  #----------------------#

  if (outlier_cutoff <= 0 | outlier_cutoff >= 1) {
    stop('outlier_cutoff must be in (0, 1)')
  }

  #------------------------------------#
  #    Objects for the EM Algorithm    #
  #------------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  w_tilde     <- matrix(NA, nrow = n, ncol = G)
  zw_tilde    <- matrix(NA, nrow = n, ncol = G)
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  py    <- rep(NA, G)
  mu    <- matrix(NA, nrow = G, ncol = d)
  Sigma <- array(NA, dim = c(d, d, G))
  df    <- rep(10, G)

  dens   <- matrix(NA, nrow = n, ncol = G)
  iter   <- 0
  loglik <- NULL

  #--------------------------------#
  #    Parameter Initialization    #
  #--------------------------------#

  init_method <- match.arg(init_method)

  if (progress) {
    cat('\nInitialization:', init_method, '\n')
  }

  if (G == 1) {

    max_iter <- 1

    pars <- cluster_pars(X = X, clusters = rep(1, n))

    py    <- 1
    mu    <- pars$mu
    Sigma <- pars$Sigma

  } else {

    init <- initialize_clusters(
      X           = X,
      G           = G,
      init_method = init_method,
      clusters    = clusters
    )

    py    <- init$pi
    mu    <- init$mu
    Sigma <- init$Sigma

  }

  #------------------------#
  #    The EM Algorithm    #
  #------------------------#

  if (progress) {
    cat('\nModel Fitting:\n')
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3, width = 75, char = "=")
  }

  while (iter < max_iter & getall(loglik) > epsilon) {

    #++++ E-step ++++#

    for (g in 1:G) {
      z[, g]       <- mvtnorm::dmvt(X, delta = mu[g, ], sigma = as.matrix(Sigma[, , g]), df = df[g], log = FALSE)
      w_tilde[, g] <- (df[g] + d) / (df[g] + mahalanobis(X, mu[g, ], Sigma[, , g], tol = 1e-20))
    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    for (g in 1:G) {
      zw_tilde[, g] <- z_tilde[, g] * w_tilde[, g]
    }

    #++++ M-Step: pi ++++#

    N  <- colSums(z_tilde)
    py <- N / n

    #++++ M-Step: mu ++++#

    for (g in 1:G) {
      mu_num  <- colSums(z_tilde[, g] * w_tilde[, g] * X)
      mu_den  <- sum(z_tilde[, g] * w_tilde[, g])
      mu[g, ] <- mu_num / mu_den
    }

    #++++ M-Step: Prepare Sigma tilde ++++#

    for (g in 1:G) {
      X_centered           <- sweep(X, 2, mu[g, ])
      X_centered_crossprod <- apply(X_centered, 1, tcrossprod)
      Sigma_tilde[, , , g] <- array(X_centered_crossprod, dim = c(d, d, n))
    }

    for (g in 1:G) {

      #++++ M-Step: Sigma ++++#

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(zw_tilde[slc_ind, g] * Sigma_tilde[, , , g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ M-step: Degree of Freedom ++++#

      root <- tryCatch({

        uniroot(function(a) {

          A <- -digamma(a / 2) + log(a / 2) + 1
          B <-  sum(z_tilde[, g] * (log(w_tilde[, g]) - w_tilde[, g])) / N[g]
          C <-  digamma((df[g] + d) / 2) - log((df[g] + d) / 2)

          return(A + B + C)

        }, lower = 2, upper = 200)$root

      }, error = function(e) return(NULL))

      df[g] <- ifelse(!is.null(root), root, df[g])

    }

    #++++ Observed Log-likelihood ++++#

    for (g in 1:G) {
      dens[, g] <- mvtnorm::dmvt(X, delta = mu[g, ], sigma = as.matrix(Sigma[, , g]), df = df[g], log = FALSE)
    }

    lik                   <- dens %*% py
    lik[lik <= 10^(-323)] <- 10^(-323)
    final_loglik          <- sum(log(lik))
    loglik                <- c(loglik, final_loglik)

    #++++ Update progress ++++#

    iter <- iter + 1

    if (progress) {
      setTxtProgressBar(pb, iter)
    }

  }

  if (progress) {
    close(pb)
    if (iter < max_iter) {
      cat('\nConvergence was reached before', max_iter, 'iterations\n')
    }
  }

  #---------------------------#
  #    Cluster Memberships    #
  #---------------------------#

  clusters <- apply(z_tilde, 1, which.max)

  #-------------------------#
  #    Outlier Detection    #
  #-------------------------#

  delta <- matrix(NA, nrow = n, ncol = G)

  for (g in 1:G) {
    delta[, g] <- mahalanobis(X, center = mu[g, ], cov = Sigma[, , g], tol = 1e-20)
  }

  cluster_matr <- clusters_to_matrix(clusters, G)
  delta        <- rowSums(delta * cluster_matr)
  outliers     <- 1 - pchisq(delta, df = d) < 1 - outlier_cutoff

  #------------------------#
  #  Number of parameters  #
  #------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G * d,
    Sigma = G * d * (d + 1) / 2,
    df    = G
  )
  npar$total <- Reduce('+', npar)

  #----------------------------#
  #    Information Criteria    #
  #----------------------------#

  AIC <- -2 * final_loglik + 2 * npar$total
  BIC <- -2 * final_loglik + npar$total * log(n)

  KIC  <- -2 * final_loglik + 3 * (npar$total + 1)
  KICc <- -2 * final_loglik + 2 * (npar$total + 1) * n/(n-npar$total -2) - n * digamma((n-npar$total)/2) + n * log(n/2)

  AIC3 <- -2 * final_loglik + 3 * npar$total
  CAIC <- -2 * final_loglik + npar$total * (1 + log(n))
  AICc <- -2 * final_loglik + 2 * npar$total * n/(n - npar$total - 1)

  ent <- apply(z_tilde, 1, max)
  ICL <- BIC - sum(ent * log(ent))

  AWE <- -2 * (final_loglik + sum(ent * log(ent))) + 2 * npar$total * (3/2 + log(n))
  CLC <- -2 * final_loglik + 2 * sum(ent * log(ent))

  #----------------------#
  #    Prepare Output    #
  #----------------------#

  c_names <- paste('comp', 1:G, sep = '')
  v_names <- colnames(X)

  if (is.null(v_names)) {
    v_names <- 1:d
  }

  names(py)       <- c_names
  rownames(mu)    <- c_names
  colnames(mu)    <- v_names
  dimnames(Sigma) <- list(v_names, v_names, c_names)
  names(df)       <- c_names

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
  }

  output <- list(
    model         = 't_complete_data',
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    df            = df,
    z_tilde       = z_tilde,
    clusters      = clusters,
    outliers      = outliers,
    data          = X,
    complete      = complete.cases(X),
    npar          = npar,
    max_iter      = max_iter,
    iter_stop     = iter,
    final_loglik  = final_loglik,
    loglik        = loglik,
    AIC           = AIC,
    BIC           = BIC,
    KIC           = KIC,
    KICc          = KICc,
    AIC3          = AIC3,
    CAIC          = CAIC,
    AICc          = AICc,
    ent           = ent,
    ICL           = ICL,
    AWE           = AWE,
    CLC           = CLC,
    init_method   = init_method
  )
  class(output) <- 'MixtureMissing'

  return(output)

}
