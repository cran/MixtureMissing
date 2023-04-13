#######################################################
###                                                 ###
###           Multivariate Normal Mixture           ###
###                                                 ###
#######################################################

#' Multivariate Normal Mixture (MNM)
#'
#' Carries out model-based clustering using a multivariate normal
#' mixture (MNM). The function will determine itself if the data set is
#' complete or incomplete and fit the appropriate model accordingly. In the incomplete
#' case, the data set must be at least bivariate, and missing values are assumed to
#' be missing at random (MAR).
#'
#' @param X An \eqn{n} by \eqn{d} matrix or data frame where \eqn{n} is the number of
#'   observations and \eqn{d} is the number of variables.
#' @param G The number of clusters, which must be at least 1. If \code{G = 1}, then
#'   both \code{init_method} and \code{clusters} are ignored.
#' @param max_iter (optional) A numeric value giving the maximum number of
#'   iterations each EM algorithm is allowed to use; 20 by default.
#' @param epsilon (optional) A number specifying the epsilon value for the
#'   Aitken-based stopping criterion used in the EM algorithm: 0.01 by default.
#' @param init_method (optional) A string specifying the method to initialize
#'   the EM algorithm. "kmedoids" clustering is used by default. Alternative
#'   methods include "kmeans", "hierarchical", "manual", "emEM",
#'   and "RndEM". When "manual" is chosen, a vector \code{clusters} of
#'   length \eqn{n} must be specified.
#' @param clusters (optional) A vector of length \eqn{n} that specifies the initial
#'   cluster memberships of the user when \code{init_method} is set to "manual".
#'   This argument is NULL by default, so that it is ignored whenever other given
#'   initialization methods are chosen.
#' @param impute (optional) A logical value indicating whether missing values should
#'   be imputed for initialization. It is FALSE by default, in which only complete
#'   observations are used for obtaining initial parameters. When it is TRUE, imputation
#'   varies depending on the initialization method selected. For "emEM" and "RndEM",
#'   after observations are randomly assigned cluster memberships, missing values
#'   are replaced by the corresponding cluster means. For other heuristic methods,
#'   mean imputation is applied on the whole data set as a pre-processing step.
#' @param equal_prop (optional) A logical value indicating whether mixing
#'   proportions should be equal when initialized with emEM or RndEM; FALSE by
#'   default.
#' @param identity_cov (optional) A logical value indicating whether covariance
#'   matrices should be set to identity matrices when initialized with emEM or RndEM;
#'   FALSE by default.
#' @param progress (optional) A logical value indicating whether the
#'   fitting progress should be displayed; TRUE by default.
#' @param n_run (optional) Number of random sets to consider for initialization
#'   if \code{init_method = "emEM"} or \code{init_method = "RndEM"}; 100 by default.
#' @param n_short (optional) Number of iterations in each run of the short EM
#'   phase if \code{init_method = "emEM"}. It is ignored when another initialization
#'   method is used. When \code{init_method = "emEM"}, emEM reduces to RndEM. It is
#'   NULL by default.
#' @param short_eps (optional) The epsilon value for the Aitken-based stopping criterion
#'   used the short EM phase. The value is ignored if \code{n_short} is specified (not NULL).
#'   By default, it is 0.1.
#'
#' @return An object of class \code{MixtureMissing} with:
#'   \item{model}{The model used to fit the data set}
#'   \item{pi}{Mixing proportions.}
#'   \item{mu}{Component mean vectors.}
#'   \item{Sigma}{Component covariance matrices.}
#'   \item{z_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
#'     probabilities that the corresponding observation belongs to each cluster.}
#'   \item{clusters}{A numeric vector of length \eqn{n} indicating cluster
#'     memberships determined by the model.}
#'   \item{data}{The original data set if it is complete; otherwise, this is
#'     the data set with missing values imputed by appropriate expectations.}
#'   \item{complete}{A logical vector of length \eqn{n} indicating which observation(s)
#'     have no missing values.}
#'   \item{npar}{The breakdown of the number of parameters to estimate.}
#'   \item{max_iter}{Maximum number of iterations allowed in the EM algorithm.}
#'   \item{iter_stop}{The actual number of iterations needed when fitting the
#'     data set.}
#'   \item{final_loglik}{The final value of log-likelihood.}
#'   \item{loglik}{All the values of log-likelihood.}
#'   \item{AIC}{Akaike information criterion.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{KIC}{Kullback information criterion.}
#'   \item{KICc}{Corrected Kullback information criterion.}
#'   \item{AIC3}{Modified AIC.}
#'   \item{CAIC}{Bozdogan's consistent AIC.}
#'   \item{AICc}{Small-sample version of AIC.}
#'   \item{ent}{Entropy}
#'   \item{ICL}{Integrated Completed Likelihood criterion.}
#'   \item{AWE}{Approximate weight of evidence.}
#'   \item{CLC}{Classification likelihood criterion.}
#'   \item{init_method}{The initialization method used in model fitting.}
#'   \item{n_run}{Number of random sets considered for initialization if emEM or RndEM is used.}
#'   \item{n_short}{Number of iterations used in each run of the short EM phase.}
#'   \item{short_eps}{The epsilon value for the Aitken-based stopping criterion used the short EM phase.}
#'
#' @references
#' Wolfe, J. H. (1965). A computer program for the maximum likelihood analysis
#'   of types. Technical report, NAVAL PERSONNEL RESEARCH ACTIVITY SAN DIEGO United States. \cr \cr
#' Ghahramani, Z. and Jordan, M. I. (1995).  Learning from incomplete data.
#'
#' @examples
#' data('nm_5_noise_close_100')
#'
#' #++++ With no missing values ++++#
#'
#' X <- nm_5_noise_close_100[, 1:2]
#' mod <- MNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#'
#' summary(mod)
#' plot(mod)
#'
#' #++++ With missing values ++++#
#'
#' set.seed(1234)
#'
#' X <- hide_values(nm_5_noise_close_100[, 1:2], prop_cases = 0.1)
#' mod <- MNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#'
#' summary(mod)
#' plot(mod)

#' @import mvtnorm
#' @export
MNM <- function(
    X,
    G,
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "manual", "emEM", "RndEM"),
    clusters     = NULL,
    impute       = FALSE,
    equal_prop   = FALSE,
    identity_cov = FALSE,
    progress     = TRUE,
    n_run        = 100,
    n_short      = NULL,
    short_eps    = 0.1
) {

  #----------------------#
  #    Input Checking    #
  #----------------------#

  if (G < 1) {
    stop('Number of clusters G must be at least 1')
  }

  if (G %% 1 != 0) {
    stop('Number of clusters G must be an integer')
  }

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X)) {
    X <- matrix(X, nrow = length(X), ncol = 1)
  }

  if (!is.numeric(X)) {
    stop('X must be a numeric matrix, data frame or vector')
  }

  #---------------------#
  #    Model Fitting    #
  #---------------------#

  if (any(is.na(X))) {

    if (ncol(X) < 2) {
      stop('If X contains NAs, X must be at least bivariate')
    }

    mod <- MNM_incomplete_data(
      X              = X,
      G              = G,
      max_iter       = max_iter,
      epsilon        = epsilon,
      init_method    = init_method,
      clusters       = clusters,
      impute         = impute,
      equal_prop     = equal_prop,
      identity_cov   = identity_cov,
      progress       = progress,
      n_run          = n_run,
      n_short        = n_short,
      short_eps      = short_eps
    )

  } else {

    mod <- MNM_complete_data(
      X              = X,
      G              = G,
      max_iter       = max_iter,
      epsilon        = epsilon,
      init_method    = init_method,
      clusters       = clusters,
      equal_prop     = equal_prop,
      identity_cov   = identity_cov,
      progress       = progress,
      n_run          = n_run,
      n_short        = n_short,
      short_eps      = short_eps
    )

  }

  return(mod)
}

###########################################################################
###                                                                     ###
###           Multivariate Normal Mixture for Incomplete Data           ###
###                                                                     ###
###########################################################################

MNM_incomplete_data <- function(
    X,
    G,
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "manual", "emEM", "RndEM"),
    clusters     = NULL,
    impute       = FALSE,
    equal_prop   = FALSE,
    identity_cov = FALSE,
    progress     = TRUE,
    n_run        = 100,
    n_short      = NULL,
    short_eps    = 0.1
) {

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
  X_tilde     <- array(rep(X, G), dim = c(n, d, G))
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  py    <- rep(NA, G)
  mu    <- matrix(NA, nrow = G, ncol = d)
  Sigma <- array(NA, dim = c(d, d, G))

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

  X_imp <- X

  if (G == 1 | !(init_method %in% c('emEM', 'RndEM')) ) {
    n_run     <- 0
    n_short   <- 0
    short_eps <- -Inf
  }

  if (G == 1) {

    max_iter <- 1

    if (impute) {
      X_imp <- mean_impute(X)
    }

    pars <- cluster_pars(X = X_imp, clusters = rep(1, n))

    py    <- 1
    mu    <- pars$mu
    Sigma <- pars$Sigma

  } else {

    if (init_method %in% c('emEM', 'RndEM')) {

      if (init_method == 'RndEM') {
        n_short   <- 1
        short_eps <- Inf
      }

      best_loglik <- -Inf
      best_run    <- NULL

      if (progress) {
        pb <- txtProgressBar(min = 0, max = n_run, style = 3, width = 75, char = "=")
      }

      if (equal_prop) {
        py <- rep(1/G, G)
      }

      if (identity_cov) {
        Sigma <- array(diag(d), dim = c(d, d, G))
      }

      for (r in 1:n_run) {

        run <- tryCatch({

          if (!equal_prop) {
            repeat {
              py <- runif(G)
              py <- py / sum(py)
              if (min(py) > 0.05) break
            }
          }

          clusters <- sample(1:G, n, prob = py, replace = TRUE)

          if (impute) {
            X_imp <- cluster_impute(X, clusters)
          }

          pars <- cluster_pars(X = X_imp, clusters = clusters)

          mu <- pars$mu

          if (!identity_cov) {
            Sigma <- pars$Sigma
          }

          EM_MNM_incomplete(
            X       = X,
            G       = G,
            py      = py,
            mu      = mu,
            Sigma   = Sigma,
            eps     = short_eps,
            n_iter  = n_short
          )

        }, error = function(e) {
          # warning('Run ', r, ' was ignored\n')
          return(NULL)
        })

        if ( !is.null(run) ) {
          if ( length(run$final_loglik) == 1 ) {
            if (run$final_loglik > best_loglik) {
              best_run    <- run
              best_loglik <- run$final_loglik
            }
          }
        }

        if (progress) {
          setTxtProgressBar(pb, r)
        }

      }

      if (progress) {
        close(pb)
      }

      py    <- best_run$py
      mu    <- best_run$mu
      Sigma <- best_run$Sigma

    } else {

      if (impute) {
        X_imp <- mean_impute(X)
      }

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

        z[Im[[j]], g] <- mvtnorm::dmvnorm(Xo_j, mean = mu_o, sigma = as.matrix(Sigma_oo))

      }
    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    for (g in 1:G) {
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
      mu_num  <- colSums(z_tilde[, g] * X_tilde[, , g])
      mu[g, ] <- mu_num / N[g]
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
            Sigma_tilde[m, m, i, g] <- tcrossprod(X_tilde[i, m, g] - mu_m) + S_mm
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

    #++++ M-step: Sigma ++++#

    for (g in 1:G) {

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(z_tilde[slc_ind, g] * Sigma_tilde[, , , g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
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

        dens[Im[[j]], g] <- mvtnorm::dmvnorm(Xo_j, mean = mu_o, sigma = as.matrix(Sigma_oo))

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

  #------------------------#
  #  Number of parameters  #
  #------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G * d,
    Sigma = G * d * (d + 1) / 2
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

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
  }

  output <- list(
    model         = 'MNM_incomplete_data',
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    z_tilde       = z_tilde,
    clusters      = clusters,
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
    init_method   = init_method,
    n_run         = n_run,
    n_short       = n_short,
    short_eps     = short_eps
  )
  class(output) <- 'MixtureMissing'

  return(output)

}

############################################################################
###                                                                      ###
###           Short EM Iterations for MNM with Incomplete Data           ###
###                                                                      ###
############################################################################

EM_MNM_incomplete <- function(
    X      = X,
    G      = G,
    py     = py,
    mu     = mu,
    Sigma  = Sigma,
    eps    = 0.1,
    n_iter = NULL
) {

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
  X_tilde     <- array(rep(X, G), dim = c(n, d, G))
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  dens   <- matrix(NA, nrow = n, ncol = G)
  iter   <- 0
  loglik <- NULL

  #------------------------#
  #    The EM Algorithm    #
  #------------------------#

  if (is.null(n_iter)) {
    n_iter <- Inf
  } else {
    eps <- -Inf
  }

  while (iter < n_iter & getall(loglik) > eps) {

    #++++ E-step ++++#

    for (g in 1:G) {
      for (j in 1:np) {

        m    <- M[j, ]                         # missing pattern j
        o    <- !m                             # observed pattern j
        Xo_j <- X[Im[[j]], o, drop = FALSE]    # observations with missing pattern j

        mu_o     <- mu[g, o]
        Sigma_oo <- Sigma[o, o, g]

        z[Im[[j]], g] <- mvtnorm::dmvnorm(Xo_j, mean = mu_o, sigma = as.matrix(Sigma_oo))

      }
    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    for (g in 1:G) {
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
      mu_num  <- colSums(z_tilde[, g] * X_tilde[, , g])
      mu[g, ] <- mu_num / N[g]
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
            Sigma_tilde[m, m, i, g] <- tcrossprod(X_tilde[i, m, g] - mu_m) + S_mm
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

    #++++ M-step: Sigma ++++#

    for (g in 1:G) {

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(z_tilde[slc_ind, g] * Sigma_tilde[, , , g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
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

        dens[Im[[j]], g] <- mvtnorm::dmvnorm(Xo_j, mean = mu_o, sigma = as.matrix(Sigma_oo))

      }
    }

    lik                   <- dens %*% py
    lik[lik <= 10^(-323)] <- 10^(-323)
    final_loglik          <- sum(log(lik))
    loglik                <- c(loglik, final_loglik)

    #++++ Update Progress ++++#

    iter <- iter + 1

  }

  #----------------------#
  #    Prepare Output    #
  #----------------------#

  output <- list(
    py           = py,
    mu           = mu,
    Sigma        = Sigma,
    z_tilde      = z_tilde,
    loglik       = loglik,
    final_loglik = final_loglik
  )

  return(output)

}

###########################################################################
###                                                                     ###
###            Multivariate Normal Mixture for Complete Data            ###
###                                                                     ###
###########################################################################

MNM_complete_data <- function(
    X,
    G,
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "manual", "emEM", "RndEM"),
    clusters     = NULL,
    equal_prop   = FALSE,
    identity_cov = FALSE,
    progress     = TRUE,
    n_run        = 100,
    n_short      = NULL,
    short_eps    = 0.1
) {

  #------------------------------------#
  #    Objects for the EM Algorithm    #
  #------------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  py    <- rep(NA, G)
  mu    <- matrix(NA, nrow = G, ncol = d)
  Sigma <- array(NA, dim = c(d, d, G))

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

  if (G == 1 | !(init_method %in% c('emEM', 'RndEM')) ) {
    n_run     <- 0
    n_short   <- 0
    short_eps <- -Inf
  }

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

    if (init_method %in% c('emEM', 'RndEM')) {

      if (init_method == 'RndEM') {
        n_short   <- 1
        short_eps <- Inf
      }

      best_loglik <- -Inf
      best_run    <- NULL

      if (progress) {
        pb <- txtProgressBar(min = 0, max = n_run, style = 3, width = 75, char = "=")
      }

      if (equal_prop) {
        py <- rep(1/G, G)
      }

      if (identity_cov) {
        Sigma <- array(diag(d), dim = c(d, d, G))
      }

      for (r in 1:n_run) {

        run <- tryCatch({

          if (!equal_prop) {
            repeat {
              py <- runif(G)
              py <- py / sum(py)
              if (min(py) > 0.05) break
            }
          }

          pars <- cluster_pars(X = X, clusters = sample(1:G, n, prob = py, replace = TRUE))

          mu <- pars$mu

          if (!identity_cov) {
            Sigma <- pars$Sigma
          }

          EM_MNM_complete(
            X       = X,
            G       = G,
            py      = py,
            mu      = mu,
            Sigma   = Sigma,
            eps     = short_eps,
            n_iter  = n_short
          )

        }, error = function(e) {
          # warning('Run ', r, ' was ignored\n')
          return(NULL)
        })

        if ( !is.null(run) ) {
          if ( length(run$final_loglik) == 1 ) {
            if (run$final_loglik > best_loglik) {
              best_run    <- run
              best_loglik <- run$final_loglik
            }
          }
        }

        if (progress) {
          setTxtProgressBar(pb, r)
        }

      }

      if (progress) {
        close(pb)
      }

      py    <- best_run$py
      mu    <- best_run$mu
      Sigma <- best_run$Sigma

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
      z[, g] <- mvtnorm::dmvnorm(X, mean = mu[g, ], sigma = as.matrix(Sigma[, , g]))
    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    #++++ M-Step: pi ++++#

    N  <- colSums(z_tilde)
    py <- N / n

    #++++ M-Step: mu ++++#

    for (g in 1:G) {
      mu_num  <- colSums(z_tilde[, g] * X)
      mu[g, ] <- mu_num / N[g]
    }

    #++++ M-Step: Prepare Sigma tilde ++++#

    for (g in 1:G) {
      X_centered           <- sweep(X, 2, mu[g, ])
      X_centered_crossprod <- apply(X_centered, 1, tcrossprod)
      Sigma_tilde[, , , g] <- array(X_centered_crossprod, dim = c(d, d, n))
    }

    #++++ M-Step: Sigma ++++#

    for (g in 1:G) {

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(z_tilde[slc_ind, g] * Sigma_tilde[, , , g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

    }

    #++++ Observed Log-likelihood ++++#

    for (g in 1:G) {
      dens[, g] <- mvtnorm::dmvnorm(X, mean = mu[g, ], sigma = as.matrix(Sigma[, , g]))
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

  #------------------------#
  #  Number of parameters  #
  #------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G * d,
    Sigma = G * d * (d + 1) / 2
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

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
  }

  output <- list(
    model         = 'MNM_complete_data',
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    z_tilde       = z_tilde,
    clusters      = clusters,
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
    init_method   = init_method,
    n_run         = n_run,
    n_short       = n_short,
    short_eps     = short_eps
  )
  class(output) <- 'MixtureMissing'

  return(output)

}


############################################################################
###                                                                      ###
###            Short EM Iterations for MNM with Complete Data            ###
###                                                                      ###
############################################################################

EM_MNM_complete <- function(
    X      = X,
    G      = G,
    py     = py,
    mu     = mu,
    Sigma  = Sigma,
    eps    = 0.1,
    n_iter = NULL
) {

  #------------------------------------#
  #    Objects for the EM Algorithm    #
  #------------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  dens   <- matrix(NA, nrow = n, ncol = G)
  iter   <- 0
  loglik <- NULL

  #------------------------#
  #    The EM Algorithm    #
  #------------------------#

  if (is.null(n_iter)) {
    n_iter <- Inf
  } else {
    eps <- -Inf
  }

  while (iter < n_iter & getall(loglik) > eps) {

    #++++ E-step ++++#

    for (g in 1:G) {
      z[, g] <- mvtnorm::dmvnorm(X, mean = mu[g, ], sigma = as.matrix(Sigma[, , g]))
    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    #++++ M-Step: pi ++++#

    N  <- colSums(z_tilde)
    py <- N / n

    #++++ M-Step: mu ++++#

    for (g in 1:G) {
      mu_num  <- colSums(z_tilde[, g] * X)
      mu[g, ] <- mu_num / N[g]
    }

    #++++ M-Step: Prepare Sigma tilde ++++#

    for (g in 1:G) {
      X_centered           <- sweep(X, 2, mu[g, ])
      X_centered_crossprod <- apply(X_centered, 1, tcrossprod)
      Sigma_tilde[, , , g] <- array(X_centered_crossprod, dim = c(d, d, n))
    }

    #++++ M-Step: Sigma ++++#

    for (g in 1:G) {

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(z_tilde[slc_ind, g] * Sigma_tilde[, , , g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

    }

    #++++ Observed Log-likelihood ++++#

    for (g in 1:G) {
      dens[, g] <- mvtnorm::dmvnorm(X, mean = mu[g, ], sigma = as.matrix(Sigma[, , g]))
    }

    lik                   <- dens %*% py
    lik[lik <= 10^(-323)] <- 10^(-323)
    final_loglik          <- sum(log(lik))
    loglik                <- c(loglik, final_loglik)

    #++++ Update progress ++++#

    iter <- iter + 1

  }

  #----------------------#
  #    Prepare Output    #
  #----------------------#

  output <- list(
    py           = py,
    mu           = mu,
    Sigma        = Sigma,
    z_tilde      = z_tilde,
    loglik       = loglik,
    final_loglik = final_loglik
  )

  return(output)

}
