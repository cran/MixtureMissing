###########################################################################
###                                                                     ###
###           Multivariate Skew-t Mixture for Incomplete Data           ###
###                                                                     ###
###########################################################################

MStM_incomplete_data <- function(
    X,
    G,
    model        = c('St', 'SC'),
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "mclust", "manual"),
    clusters     = NULL,
    deriv_ctrl   = list(eps = 1e-8, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-7),
                        r = 6, v = 2, show.details = FALSE),
    progress     = TRUE
) {

  #------------------------------------#
  #    Objects for the EM Algorithm    #
  #------------------------------------#

  model <- match.arg(model)

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

  z           <- matrix(NA_real_, nrow = n, ncol = G)
  z_tilde     <- matrix(NA_real_, nrow = n, ncol = G)
  X_tilde     <- array(rep(X, G), dim = c(n, d, G))
  X_hat       <- array(rep(X, G), dim = c(n, d, G))
  Sigma_tilde <- array(NA_real_, dim = c(d, d, n, G))

  a     <- matrix(NA_real_, nrow = n, ncol = G)
  b     <- matrix(NA_real_, nrow = n, ncol = G)
  c     <- matrix(NA_real_, nrow = n, ncol = G)
  N     <- rep(NA_real_, G)
  a_bar <- rep(NA_real_, G)
  b_bar <- rep(NA_real_, G)
  c_bar <- rep(NA_real_, G)
  X_bar <- matrix(NA_real_, nrow = G, ncol = d)

  py    <- rep(NA_real_, G)
  mu    <- matrix(NA_real_, nrow = G, ncol = d)
  Sigma <- array(NA_real_, dim = c(d, d, G))
  beta <- matrix(0.01, nrow = G, ncol = d)

  if (model == 'St') {
    df <- rep(10, G)
  }

  if (model == 'SC') {
    df <- rep(1, G)
  }

  log_dens <- matrix(NA_real_, nrow = n, ncol = G)
  iter     <- 0
  loglik   <- NULL

  #--------------------------------#
  #    Parameter Initialization    #
  #--------------------------------#

  init_method <- match.arg(init_method)

  X_imp <- mean_impute(X)

  if (G == 1) {

    max_iter <- 1

    pars <- cluster_pars(
      X        = X_imp,
      clusters = rep(1, n)
    )

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
      for (j in 1:np) {

        m    <- M[j, ]                         # missing pattern j
        o    <- !m                             # observed pattern j
        Xo_j <- X[Im[[j]], o, drop = FALSE]    # observations with missing pattern j

        mu_o         <- mu[g, o]
        Sigma_oo     <- Sigma[o, o, g]
        Sigma_oo_inv <- solve(Sigma_oo)
        beta_o       <- beta[g, o]

        z[Im[[j]], g] <- dSt(Xo_j, mu = mu_o, Sigma = Sigma_oo, beta = beta_o, df = df[g])

        psi <- c( t(beta_o) %*% Sigma_oo_inv %*% beta_o )
        chi <- df[g] + mahalanobis(Xo_j, center = mu_o, cov = Sigma_oo)

        s1 <- sqrt(psi * chi)
        s2 <- sqrt(chi / psi)

        bessel_num <- besselK(s1, nu = -df[g]/2 - sum(o)/2 + 1, expon.scaled = TRUE)
        bessel_den <- besselK(s1, nu = -df[g]/2 - sum(o)/2, expon.scaled = TRUE)
        bessel_den[bessel_den < 10^-323] <- 10^-323

        a[Im[[j]], g] <- s2 * (bessel_num / bessel_den)
        b[Im[[j]], g] <- (df[g] + sum(o)) / chi + (bessel_num / bessel_den) / s2
        c[Im[[j]], g] <- log(s2) + numDeriv::grad(log_besselK, x = rep(-df[g]/2 - sum(o)/2, nrow(Xo_j)), y = s1,
                                                  method = 'Richardson', method.args = deriv_ctrl)

      }
    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    N <- colSums(z_tilde)

    a_bar <- colSums(z_tilde * a) / N
    b_bar <- colSums(z_tilde * b) / N
    c_bar <- colSums(z_tilde * c) / N

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
          Sigma_oo_inv <- solve(Sigma_oo)

          beta_o <- beta[g, o]
          beta_m <- beta[g, m]

          for (i in Im[[j]]) {

            xi <- X[i, ]

            mu_m_o    <- mu_m + Sigma_mo %*% Sigma_oo_inv %*% (xi[o] - mu_o)
            Sigma_m_o <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% Sigma_om
            beta_m_o  <- beta_m - Sigma_mo %*% Sigma_oo_inv %*% beta_o

            X_hat[i, m, g]   <- mu_m_o + a[i, g] * beta_m_o
            X_tilde[i, m, g] <- b[i, g] * mu_m_o + beta_m_o

            Sigma_tilde[m, m, i, g] <- Sigma_m_o + b[i, g] * tcrossprod(mu_m_o) + mu_m_o %*% t(beta_m_o) + beta_m_o %*% t(mu_m_o)
            Sigma_tilde[m, m, i, g] <- Sigma_tilde[m, m, i, g] + a[i, g] * tcrossprod(beta_m_o)

          }

        }

      }
    }

    for (g in 1:G) {
      X_bar[g, ] <- colSums(z_tilde[, g] * X_hat[, , g]) / N[g]
    }

    #++++ M-step: pi ++++#

    py <- N / n

    #++++ M-step: mu (location) and beta (skewness) ++++#

    for (g in 1:G) {

      num_mu <- colSums( z_tilde[, g] * ( (!R) * (a_bar[g] * b[, g] - 1) * X_tilde[, , g] + R * (a_bar[g] * X_tilde[, , g] - X_hat[, , g]) ) )
      den_mu <- sum( z_tilde[, g] * (a_bar[g] * b[, g] - 1) )

      mu[g, ] <- num_mu / den_mu

      num_beta <- colSums( z_tilde[, g] * ( (!R) * (b_bar[g] - b[, g]) * X_tilde[, , g] + R * (b_bar[g] * X_hat[, , g] - X_tilde[, , g]) ) )
      den_beta <- den_mu

      beta[g, ] <- num_beta / den_beta

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
          Sigma_oo_inv <- solve(Sigma_oo)

          beta_o <- beta[g, o]
          beta_m <- beta[g, m]

          for (i in Im[[j]]) {

            xi <- X[i, ]

            Sigma_tilde[o, o, i, g] <- b[i, g] * tcrossprod(xi[o] - mu_o)
            Sigma_tilde[o, m, i, g] <- tcrossprod(xi[o] - mu_o, X_tilde[i, m, g] - b[i, g] * mu_m)
            Sigma_tilde[m, o, i, g] <- t(Sigma_tilde[o, m, i, g])

            Sigma_tilde[m, m, i, g] <- Sigma_tilde[m, m, i, g] - X_tilde[i, m, g] %*% t(mu_m) - mu_m %*% t(X_tilde[i, m, g])
            Sigma_tilde[m, m, i, g] <- Sigma_tilde[m, m, i, g] + b[i, g] * mu_m %*% t(mu_m)

          }

        } else {

          X_centrd <- sweep(X[Im[[j]], ], 2, mu[g, ], '-')
          cr_prods <- apply(X_centrd, 1, tcrossprod)

          S_tilde <- array(
            data = unlist(cr_prods),
            dim  = c(d, d, length(Im[[j]]))
          )

          slc_ind                     <- slice.index(S_tilde, 3)
          Sigma_tilde[, , Im[[j]], g] <- b[Im[[j]], g][slc_ind] * S_tilde

        }

      }
    }

    for (g in 1:G) {

      #++++ M-step: Sigma (dispersion) ++++#

      slc_ind      <- slice.index(Sigma_tilde[, , , g], 3)
      Sigma[, , g] <- rowSums(z_tilde[slc_ind, g] * Sigma_tilde[, , , g], dims = 2) / N[g]
      Sigma[, , g] <- Sigma[, , g] - beta[g, ] %*% t(X_bar[g, ] - mu[g, ]) - (X_bar[g, ] - mu[g, ]) %*% t(beta[g, ])
      Sigma[, , g] <- Sigma[, , g] + a_bar[g] * tcrossprod(beta[g, ])

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ M-step: df (degree of freedom) ++++#

      if (model == 'St') {
        root <- tryCatch(
          uniroot(df_MSt_func, ws_term = sum(z_tilde[, g] * (c[, g] + b[, g])) / N[g],
                  lower = 2, upper = 200)$root,
          error = function(e) { return(NULL) }
        )

        df[g] <- ifelse(!is.null(root), root, df[g])
      }

      if (model == 'SC') {
        df[g] <- 1
      }

    }

    #++++ Observed Log-Likelihood ++++#

    for (g in 1:G) {
      for (j in 1:np) {

        m    <- M[j, ]                         # missing pattern j
        o    <- !m                             # observed pattern j
        Xo_j <- X[Im[[j]], o, drop = FALSE]    # observations with missing pattern j

        mu_o     <- mu[g, o]
        Sigma_oo <- Sigma[o, o, g]
        beta_o   <- beta[g, o]

        log_dens[Im[[j]], g] <- dSt(Xo_j, mu = mu_o, Sigma = Sigma_oo, beta = beta_o, df = df[g], log = TRUE)

      }
    }

    log_py_dens  <- sweep(log_dens, 2, log(py), FUN = '+')
    final_loglik <- sum( apply(log_py_dens, 1, log_sum_exp) )
    loglik       <- c(loglik, final_loglik)

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

  #------------------#
  #    Imputation    #
  #------------------#

  X_imputed <- X
  complete  <- complete.cases(X)

  for (i in which(!complete)) {
    X_imputed[i, ] <- X_tilde[i, , clusters[i]]
  }

  #------------------------#
  #  Number of parameters  #
  #------------------------#

  npar <- list(
    pi     = G - 1,
    mu     = G * d,
    Sigma  = G * d * (d + 1) / 2,
    beta   = G * d,
    df     = G
  )

  if (model == 'SC') {
    npar$df <- NULL
  }

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
  rownames(beta)  <- c_names
  colnames(beta)  <- v_names
  names(df)       <- c_names

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
    beta  <- beta[1, ]
  }

  output <- list(
    model         = paste(model, '_incomplete_data', sep = ''),
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    beta          = beta,
    df            = df,
    z_tilde       = z_tilde,
    clusters      = clusters,
    data          = X_imputed,
    complete      = !is.na(X),
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

  if (model == 'SC') {
    output$df <- NULL
  }

  class(output) <- 'MixtureMissing'

  return(output)
}

MStM_complete_data <- function(
    X,
    G,
    model        = c('St', 'SC', 'C'),
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "mclust", "manual"),
    clusters     = NULL,
    deriv_ctrl   = list(eps = 1e-8, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-7),
                        r = 6, v = 2, show.details = FALSE),
    progress     = TRUE
) {

  #------------------------------------#
  #    Objects for the EM Algorithm    #
  #------------------------------------#

  model <- match.arg(model)

  n <- nrow(X)
  d <- ncol(X)

  z           <- matrix(NA_real_, nrow = n, ncol = G)
  z_tilde     <- matrix(NA_real_, nrow = n, ncol = G)
  Sigma_tilde <- array(NA_real_, dim = c(d, d, n, G))

  a     <- matrix(NA_real_, nrow = n, ncol = G)
  b     <- matrix(NA_real_, nrow = n, ncol = G)
  c     <- matrix(NA_real_, nrow = n, ncol = G)
  N     <- rep(NA_real_, G)
  a_bar <- rep(NA_real_, G)
  b_bar <- rep(NA_real_, G)
  c_bar <- rep(NA_real_, G)
  X_bar <- matrix(NA_real_, nrow = G, ncol = d)

  py    <- rep(NA_real_, G)
  mu    <- matrix(NA_real_, nrow = G, ncol = d)
  Sigma <- array(NA_real_, dim = c(d, d, G))
  beta <- matrix(0.01, nrow = G, ncol = d)

  if (model == 'St') {
    df <- rep(10, G)
  }

  if (model == 'SC') {
    df <- rep(1, G)
  }

  log_dens <- matrix(NA_real_, nrow = n, ncol = G)
  iter     <- 0
  loglik   <- NULL

  #--------------------------------#
  #    Parameter Initialization    #
  #--------------------------------#

  init_method <- match.arg(init_method)

  if (G == 1) {

    max_iter <- 1

    pars <- cluster_pars(
      X        = X,
      clusters = rep(1, n)
    )

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
      z[, g] <- dSt(X, mu = mu[g, ], Sigma = Sigma[, , g], beta = beta[g, ], df = df[g])

      Sigma_inv <- solve(Sigma[, , g])

      psi <- c( t(beta[g, ]) %*% Sigma_inv %*% beta[g, ] )
      chi <- df[g] + mahalanobis(X, center = mu[g, ], cov = Sigma[, , g])

      s1 <- sqrt(psi * chi)
      s2 <- sqrt(chi / psi)

      bessel_num <- besselK(s1, nu = -df[g]/2 - d/2 + 1, expon.scaled = TRUE)
      bessel_den <- besselK(s1, nu = -df[g]/2 - d/2, expon.scaled = TRUE)
      bessel_den[bessel_den < 10^-323] <- 10^-323

      a[, g] <- s2 * (bessel_num / bessel_den)
      b[, g] <- (df[g] + d) / chi + (bessel_num / bessel_den) / s2
      c[, g] <- log(s2) + numDeriv::grad(log_besselK, x = rep(-df[g]/2 - d/2, n), y = s1,
                                         method = 'Richardson', method.args = deriv_ctrl)

    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    N <- colSums(z_tilde)

    a_bar <- colSums(z_tilde * a) / N
    b_bar <- colSums(z_tilde * b) / N
    c_bar <- colSums(z_tilde * c) / N

    for (g in 1:G) {
      X_bar[g, ] <- colSums(z_tilde[, g] * X) / N[g]
    }

    #++++ M-step: pi ++++#

    py <- N / n

    for (g in 1:G) {

      #++++ M-step: mu (location) and beta (skewness) ++++#

      num_mu <- colSums( z_tilde[, g] * X * (a_bar[g] * b[, g] - 1) )
      den_mu <- sum( z_tilde[, g] * (a_bar[g] * b[, g] - 1) )

      mu[g, ] <- num_mu / den_mu

      num_beta <- colSums( z_tilde[, g] * X * (b_bar[g] - b[, g]) )
      den_beta <- den_mu

      beta[g, ] <- num_beta / den_beta

      #++++ M-step: Sigma (dispersion) ++++#

      X_centrd         <- sweep(X, 2, mu[g, ])
      X_centrd_crsprod <- apply(X_centrd, 1, tcrossprod)
      Sigma_tilde      <- array(X_centrd_crsprod, dim = c(d, d, n))

      slc_ind      <- slice.index(Sigma_tilde, 3)
      Sigma[, , g] <- rowSums(z_tilde[slc_ind, g] * b[slc_ind, g] * Sigma_tilde, dims = 2) / N[g]
      Sigma[, , g] <- Sigma[, , g] - beta[g, ] %*% t(X_bar[g, ] - mu[g, ]) - (X_bar[g, ] - mu[g, ]) %*% t(beta[g, ])
      Sigma[, , g] <- Sigma[, , g] + a_bar[g] * tcrossprod(beta[g, ])

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ M-step: df (degree of freedom) ++++#

      if (model == 'St') {
        root <- tryCatch(
          uniroot(df_MSt_func, ws_term = sum(z_tilde[, g] * (c[, g] + b[, g])) / N[g],
                  lower = 2, upper = 200)$root,
          error = function(e) { return(NULL) }
        )

        df[g] <- ifelse(!is.null(root), root, df[g])
      }

      if (model == 'SC') {
        df[g] <- 1
      }

    }

    #++++ Observed Log-likelihood ++++#

    for (g in 1:G) {
      log_dens[, g] <- dSt(X, mu = mu[g, ], Sigma = Sigma[, , g], beta = beta[g, ], df = df[g], log = TRUE)
    }

    log_py_dens  <- sweep(log_dens, 2, log(py), FUN = '+')
    final_loglik <- sum( apply(log_py_dens, 1, log_sum_exp) )
    loglik       <- c(loglik, final_loglik)

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

  #----------------------------#
  #    Number of Parameters    #
  #----------------------------#

  npar <- list(
    pi     = G - 1,
    mu     = G * d,
    Sigma  = G * d * (d + 1) / 2,
    beta   = G * d,
    df     = G
  )

  if (model == 'SC') {
    npar$df <- NULL
  }

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
  rownames(beta)  <- c_names
  colnames(beta)  <- v_names
  names(df)       <- c_names

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
    beta  <- beta[1, ]
  }

  output <- list(
    model         = paste(model, '_complete_data', sep = ''),
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    beta          = beta,
    df            = df,
    z_tilde       = z_tilde,
    clusters      = clusters,
    data          = X,
    complete      = !is.na(X),
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

  if (model == 'SC') {
    output$df <- NULL
  }

  class(output) <- 'MixtureMissing'

  return(output)

}

###########################################################################
###                                                                     ###
###                  Objective Function to Optimize nu                  ###
###                                                                     ###
###########################################################################

df_MSt_func <- function(x, ws_term) {

  log(x / 2) + 1 - digamma(x / 2) - ws_term

}


###########################################################################
###                                                                     ###
###        Density Function for Multivariate Skew-t Distribution        ###
###                                                                     ###
###########################################################################

dSt <- function(
    X,
    mu    = rep(0, d),    # location
    Sigma = diag(d),      # dispersion
    beta  = rep(0, d),    # skewness
    df    = 10,           # degree of freedom
    log   = FALSE
) {

  #----------------------#
  #    Input Checking    #
  #----------------------#

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (is.vector(X)) {
    X <- matrix(X, nrow = 1, ncol = length(X))
  }

  if (!is.matrix(X)) {
    stop('X must be a vector or matrix')
  }

  if (is.vector(Sigma)) {
    if (length(Sigma) == 1) {
      Sigma <- matrix(Sigma, nrow = 1, ncol = 1)
    }
  }

  n <- nrow(X)
  d <- ncol(X)

  if (length(mu) != d) {
    stop('mu must be a vector of length d')
  }

  if (nrow(Sigma) != d | ncol(Sigma) != d) {
    stop('Sigma must be a d x d matrix')
  }

  if (length(beta) != d) {
    stop('beta must be a vector of length d')
  }

  X_centrd  <- sweep(X, 2, mu, '-')
  Sigma_inv <- solve(Sigma)

  psi <- c( t(beta) %*% Sigma_inv %*% beta )
  chi <- df + mahalanobis(X, center = mu, cov = Sigma_inv, inverted = TRUE)

  s1 <- sqrt(psi * chi)

  res <- (-df/2 - d/2) * log(s1) + (df/2 + d/2) * log(psi)
  res <- res + (df / 2) * log(df) + log(besselK(s1, nu = -df/2 - d/2, expon.scaled = TRUE)) - s1

  if (is.nan(log(det(Sigma)))) {
    Sigma <- diag(d)
  }

  res <- res - d/2 * (log(2) + log(pi)) - 1/2 * log(det(Sigma))
  res <- res - lgamma(df / 2) - (df/2 - 1) * log(2)
  res <- res + X_centrd %*% Sigma_inv %*% beta

  if (!log) {
    res <- c(exp(res))
    res[res <= 10^(-323)] <- 10^(-323)
  }

  return(res)

  # lvx <- (df / 2) * log(df) + (-df/2 - d/2) * log(s1)
  # lvx <- lvx + log(besselK(s1, nu = -df/2 - d/2, expon.scaled = TRUE)) - s1
  # lvx <- lvx + X_centrd %*% Sigma_inv %*% beta
  #
  # lv <- -1/2 * log(det(Sigma)) - d/2 * (log(2) + log(pi))
  # lv <- lv - log(gamma(df / 2)) - (df/2 - 1) * log(2)
  # lv <- lv + (df/2 + d/2) * log(psi)
  #
  # return(
  #   c( exp(lvx + lv) )
  # )

}
