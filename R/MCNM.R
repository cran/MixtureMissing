####################################################################
###                                                              ###
###           Multivariate Contaminated Normal Mixture           ###
###                                                              ###
####################################################################

#' Multivariate Contaminated Normal Mixture (MCNM)
#'
#' Carries out model-based clustering using a multivariate contaminated normal
#' mixture (MCNM). The function will determine itself if the data set is
#' complete or incomplete and fit the appropriate model accordingly. In the incomplete
#' case, the data set must be at least bivariate, and missing values are assumed to
#' be missing at random (MAR).
#'
#' @param X An \eqn{n} x \eqn{d} matrix or data frame where \eqn{n} is the number of
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
#' @param clusters (optional) A numeric vector of length \eqn{n} that specifies the initial
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
#' @param eta_min (optional) A numeric value close to 1 to the right specifying
#'   the minimum value of eta; 1.001 by default.
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
#'   \item{alpha}{Component proportions of good observations.}
#'   \item{eta}{Component degrees of contamination.}
#'   \item{z_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
#'     probabilities that the corresponding observation belongs to each cluster.}
#'   \item{v_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
#'     probabilities that the corresponding observation is good with respect
#'     to each cluster.}
#'   \item{clusters}{A numeric vector of length \eqn{n} indicating cluster
#'     memberships determined by the model.}
#'   \item{outliers}{A logical vector of length \eqn{n} indicating observations that are ouliers.}
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
#' Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537. \cr \cr
#' Tong, H. and, Tortora, C., 2022. Model-based clustering and outlier detection
#'   with missing data. \emph{Advances in Data Analysis and Classification}.
#'
#' @examples
#'
#' data('nm_5_noise_close_100')
#'
#' #++++ With no missing values ++++#
#'
#' X <- nm_5_noise_close_100[, 1:2]
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#'
#' summary(mod)
#' plot(mod)
#'
#' #++++ With missing values ++++#
#'
#' set.seed(1234)
#'
#' X <- hide_values(nm_5_noise_close_100[, 1:2], prop_cases = 0.1)
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#'
#' summary(mod)
#' plot(mod)
#'
#' @importFrom stats complete.cases cov cutree dist dnorm hclust kmeans
#'   mahalanobis pchisq rmultinom runif var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
MCNM <- function(
    X,
    G,
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "manual", "emEM", "RndEM"),
    clusters     = NULL,
    impute       = FALSE,
    equal_prop   = FALSE,
    identity_cov = FALSE,
    eta_min      = 1.001,
    progress     = TRUE,
    n_run        = 100,
    n_short      = NULL,
    short_eps    = 0.1
) {

  #----------------------#
  #    Input checking    #
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

    mod <- MCNM_incomplete_data(
      X            = X,
      G            = G,
      max_iter     = max_iter,
      epsilon      = epsilon,
      init_method  = init_method,
      clusters     = clusters,
      impute       = impute,
      equal_prop   = equal_prop,
      identity_cov = identity_cov,
      eta_min      = eta_min,
      progress     = progress,
      n_run        = n_run,
      n_short      = n_short,
      short_eps    = short_eps
    )

  } else {

    mod <- MCNM_complete_data(
      X            = X,
      G            = G,
      max_iter     = max_iter,
      epsilon      = epsilon,
      init_method  = init_method,
      clusters     = clusters,
      equal_prop   = equal_prop,
      identity_cov = identity_cov,
      eta_min      = eta_min,
      progress     = progress,
      n_run        = n_run,
      n_short      = n_short,
      short_eps    = short_eps
    )

  }

  return(mod)
}

############################################################################
###                                                                      ###
###     Multivariate Contaminated Normal Mixture for Incomplete Data     ###
###                                                                      ###
############################################################################

MCNM_incomplete_data <- function(
    X,
    G,
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical",
                     "manual", "emEM", "RndEM"),
    clusters     = NULL,
    impute       = FALSE,
    equal_prop   = FALSE,
    identity_cov = FALSE,
    eta_min      = 1.001,
    progress     = TRUE,
    n_run        = 100,
    n_short      = NULL,
    short_eps    = 0.1
) {

  #-------------------------------------#
  #    Objects for the ECM Algorithm    #
  #-------------------------------------#

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
  v_tilde     <- matrix(NA, nrow = n, ncol = G)
  w_tilde     <- matrix(NA, nrow = n, ncol = G)
  zw_tilde    <- matrix(NA, nrow = n, ncol = G)
  X_tilde     <- array(rep(X, G), dim = c(n, d, G))
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  py    <- rep(NA, G)
  mu    <- matrix(NA, nrow = G, ncol = d)
  Sigma <- array(NA, dim = c(d, d, G))
  alpha <- rep(0.6, G)
  eta   <- rep(1.4, G)

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

          ECM_MCNM_incomplete(
            X       = X,
            G       = G,
            py      = py,
            mu      = mu,
            Sigma   = Sigma,
            alpha   = alpha,
            eta     = eta,
            eta_min = eta_min,
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
      alpha <- best_run$alpha
      eta   <- best_run$eta

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

  #-------------------------#
  #    The ECM Algorithm    #
  #-------------------------#

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

        z[Im[[j]], g]       <- dCN(Xo_j, mu = mu_o, Sigma = Sigma_oo, alpha = alpha[g], eta = eta[g])
        v_tilde[Im[[j]], g] <- alpha[g] * mvtnorm::dmvnorm(Xo_j, mean = mu_o, sigma = as.matrix(Sigma_oo)) / z[Im[[j]], g]

      }
    }

    v_tilde[is.nan(v_tilde)] <- 0

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    for (g in 1:G) {

      w_tilde[, g]  <- v_tilde[, g] + (1 - v_tilde[, g]) / eta[g]
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

    #++++ CM-Step 1: pi and alpha ++++#

    N  <- colSums(z_tilde)
    py <- N / n

    alpha              <- colSums(z_tilde * v_tilde) / N
    alpha[alpha < 0.5] <- 0.5
    alpha[alpha > 1]   <- 1

    #++++ CM-Step 1: mu ++++#

    for (g in 1:G) {

      mu_num  <- colSums(z_tilde[, g] * w_tilde[, g] * X_tilde[, , g])
      mu_den  <- sum(z_tilde[, g] * w_tilde[, g])
      mu[g, ] <- mu_num / mu_den

    }

    #++++ CM-Step 1: Prepare Sigma tilde ++++#

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

      #++++ CM-Step 1: Sigma ++++#

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(zw_tilde[slc_ind, g] * Sigma_tilde[, , ,g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ CM-Step 2: eta ++++#

      delta_o <- rep(NA, n)

      for (j in 1:np) {

        m <- M[j, ]    # missing pattern j
        o <- !m        # observed pattern j

        delta_o[Im[[j]]] <- mahalanobis(X[Im[[j]], o, drop = FALSE], mu[g, o], Sigma[o, o, g], tol = 1e-20)

      }

      eta_num <- sum(z_tilde[, g] * (1 - v_tilde[, g]) * delta_o)
      eta_den <- sum(do * z_tilde[, g] * (1 - v_tilde[, g]))
      eta[g]  <- eta_num / eta_den
      eta[g]  <- max(eta[g], eta_min)

    }

    #++++ Observed Log-Likelihood ++++#

    for (g in 1:G) {
      for (j in 1:np) {

        m    <- M[j, ]                         # missing pattern j
        o    <- !m                             # observed pattern j
        Xo_j <- X[Im[[j]], o, drop = FALSE]    # observations with missing pattern j

        mu_o     <- mu[g, o]
        Sigma_oo <- Sigma[o, o, g]

        dens[Im[[j]], g] <- dCN(Xo_j, mu = mu_o, Sigma = as.matrix(Sigma_oo), alpha = alpha[g], eta = eta[g])

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

  #-------------------------#
  #    Outlier Detection    #
  #-------------------------#

  cluster_matr <- clusters_to_matrix(clusters, G)
  outliers     <- rowSums(v_tilde * cluster_matr) < 0.5

  #------------------#
  #    Imputation    #
  #------------------#

  X_imputed <- X
  complete  <- complete.cases(X)

  for (i in which(!complete)) {
    X_imputed[i, ] <- X_tilde[i, , clusters[i]]
  }

  #----------------------------#
  #    Number of Parameters    #
  #----------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G * d,
    Sigma = G * d * (d + 1) / 2,
    alpha = G,
    eta   = G
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
    model         = 'MCNM_incomplete_data',
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    alpha         = alpha,
    eta           = eta,
    z_tilde       = z_tilde,
    v_tilde       = v_tilde,
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
###          Short ECM Iterations for MCNM with Incomplete Data          ###
###                                                                      ###
############################################################################

ECM_MCNM_incomplete <- function(
    X,
    G,
    py,
    mu,
    Sigma,
    alpha,
    eta,
    eta_min = 1.001,
    eps     = 0.1,
    n_iter  = NULL
) {

  #-------------------------------------#
  #    Objects for the ECM Algorithm    #
  #-------------------------------------#

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
  v_tilde     <- matrix(NA, nrow = n, ncol = G)
  w_tilde     <- matrix(NA, nrow = n, ncol = G)
  zw_tilde    <- matrix(NA, nrow = n, ncol = G)
  X_tilde     <- array(rep(X, G), dim = c(n, d, G))
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  dens   <- matrix(NA, nrow = n, ncol = G)
  iter   <- 0
  loglik <- NULL

  #-------------------------#
  #    The ECM Algorithm    #
  #-------------------------#

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

        z[Im[[j]], g]       <- dCN(Xo_j, mu = mu_o, Sigma = Sigma_oo, alpha = alpha[g], eta = eta[g])
        v_tilde[Im[[j]], g] <- alpha[g] * mvtnorm::dmvnorm(Xo_j, mean = mu_o, sigma = as.matrix(Sigma_oo)) / z[Im[[j]], g]

      }
    }

    v_tilde[is.nan(v_tilde)] <- 0

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    for (g in 1:G) {

      w_tilde[, g]  <- v_tilde[, g] + (1 - v_tilde[, g]) / eta[g]
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

    #++++ CM-Step 1: pi and alpha ++++#

    N  <- colSums(z_tilde)
    py <- N / n

    alpha              <- colSums(z_tilde * v_tilde) / N
    alpha[alpha < 0.5] <- 0.5
    alpha[alpha > 1]   <- 1

    #++++ CM-Step 1: mu ++++#

    for (g in 1:G) {
      mu_num  <- colSums(z_tilde[, g] * w_tilde[, g] * X_tilde[, , g])
      mu_den  <- sum(z_tilde[, g] * w_tilde[, g])
      mu[g, ] <- mu_num / mu_den
    }

    #++++ CM-Step 1: Prepare Sigma tilde ++++#

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

      #++++ CM-Step 1: Sigma ++++#

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(zw_tilde[slc_ind, g] * Sigma_tilde[, , ,g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ CM-Step 2: eta ++++#

      delta_o <- rep(NA, n)

      for (j in 1:np) {

        m <- M[j, ]    # missing pattern j
        o <- !m        # observed pattern j

        delta_o[Im[[j]]] <- mahalanobis(X[Im[[j]], o, drop = FALSE], mu[g, o], Sigma[o, o, g], tol = 1e-20)

      }

      eta_num <- sum(z_tilde[, g] * (1 - v_tilde[, g]) * delta_o)
      eta_den <- sum(do * z_tilde[, g] * (1 - v_tilde[, g]))
      eta[g]  <- eta_num / eta_den
      eta[g]  <- max(eta[g], eta_min)

    }

    #++++ Observed Log-Likelihood ++++#

    for (g in 1:G) {
      for (j in 1:np) {

        m    <- M[j, ]                         # missing pattern j
        o    <- !m                             # observed pattern j
        Xo_j <- X[Im[[j]], o, drop = FALSE]    # observations with missing pattern j

        mu_o     <- mu[g, o]
        Sigma_oo <- Sigma[o, o, g]

        dens[Im[[j]], g] <- dCN(Xo_j, mu = mu_o, Sigma = as.matrix(Sigma_oo), alpha = alpha[g], eta = eta[g])

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
    alpha        = alpha,
    eta          = eta,
    z_tilde      = z_tilde,
    v_tilde      = v_tilde,
    loglik       = loglik,
    final_loglik = final_loglik
  )

  return(output)

}

######################################################################################
###                                                                                ###
###           Multivariate Contaminated Normal Mixture for Complete Data           ###
###                                                                                ###
######################################################################################

MCNM_complete_data <- function(
    X,
    G,
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "manual", "emEM", "RndEM"),
    clusters     = NULL,
    equal_prop   = FALSE,
    identity_cov = FALSE,
    eta_min      = 1.001,
    progress     = TRUE,
    n_run        = 100,
    n_short      = NULL,
    short_eps    = 0.1
) {

  #-------------------------------------#
  #    Objects for the ECM Algorithm    #
  #-------------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  v_tilde     <- matrix(NA, nrow = n, ncol = G)
  w_tilde     <- matrix(NA, nrow = n, ncol = G)
  zw_tilde    <- matrix(NA, nrow = n, ncol = G)
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  py    <- rep(NA, G)
  mu    <- matrix(NA, nrow = G, ncol = d)
  Sigma <- array(NA, dim = c(d, d, G))
  alpha <- rep(0.6, G)
  eta   <- rep(1.4, G)

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

    pars <- cluster_pars(X = X, clusters = rep(1, n))

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

          ECM_MCNM_complete(
            X       = X,
            G       = G,
            py      = py,
            mu      = mu,
            Sigma   = Sigma,
            alpha   = alpha,
            eta     = eta,
            eta_min = eta_min,
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
      alpha <- best_run$alpha
      eta   <- best_run$eta

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

  #-------------------------#
  #    The ECM Algorithm    #
  #-------------------------#

  if (progress) {
    cat('\nModel Fitting:\n')
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3, width = 75, char = "=")
  }

  while (iter < max_iter & getall(loglik) > epsilon) {

    #++++ E-step ++++#

    for (g in 1:G) {
      z[, g]       <- dCN(X, mu = mu[g, ], Sigma = Sigma[, , g], alpha = alpha[g], eta = eta[g])
      v_tilde[, g] <- alpha[g] * mvtnorm::dmvnorm(X, mean = mu[g, ], sigma = as.matrix(Sigma[, , g])) / z[, g]
      w_tilde[, g] <- v_tilde[, g] + (1 - v_tilde[, g]) / eta[g]
    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    for (g in 1:G) {
      zw_tilde[, g] <- z_tilde[, g] * w_tilde[, g]
    }

    #++++ CM-step 1: pi and alpha ++++#

    N  <- colSums(z_tilde)
    py <- N / n

    alpha              <- colSums(z_tilde * v_tilde) / N
    alpha[alpha < 0.5] <- 0.5
    alpha[alpha > 1]   <- 1

    #++++ CM-step 1: mu ++++#

    for (g in 1:G) {
      mu_num  <- colSums(z_tilde[, g] * w_tilde[, g] * X)
      mu_den  <- sum(z_tilde[, g] * w_tilde[, g])
      mu[g, ] <- mu_num / mu_den
    }

    #++++ CM-step 1: Prepare Sigma tilde ++++#

    for (g in 1:G) {
      X_centered           <- sweep(X, 2, mu[g, ])
      X_centered_crossprod <- apply(X_centered, 1, tcrossprod)
      Sigma_tilde[, , , g] <- array(X_centered_crossprod, dim = c(d, d, n))
    }

    for (g in 1:G) {

      #++++ CM-step 1: Sigma ++++#

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(zw_tilde[slc_ind, g] * Sigma_tilde[, , , g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ CM-step 2: eta ++++#

      delta <- mahalanobis(X, mu[g, ], Sigma[, , g], tol = 1e-20)

      eta_num <- sum(z_tilde[, g] * (1 - v_tilde[, g]) * delta)
      eta_den <- d * sum(z_tilde[, g] * (1 - v_tilde[, g]))
      eta[g]  <- eta_num / eta_den
      eta[g]  <- max(eta[g], eta_min)

    }

    #++++ Observed Log-Likelihood ++++#

    for (g in 1:G) {
      dens[, g] <- dCN(X, mu = mu[g, ], Sigma = Sigma[, , g], alpha = alpha[g], eta = eta[g])
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

  #-------------------------#
  #    Outlier Detection    #
  #-------------------------#

  cluster_matr <- clusters_to_matrix(clusters, G)
  outliers     <- rowSums(v_tilde * cluster_matr) < 0.5

  #----------------------------#
  #    Number of Parameters    #
  #----------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G * d,
    Sigma = G * d * (d + 1) / 2,
    alpha = G,
    eta   = G
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
    model         = 'MCNM_complete_data',
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    alpha         = alpha,
    eta           = eta,
    z_tilde       = z_tilde,
    v_tilde       = v_tilde,
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
###           Short ECM Iterations for MCNM with Complete Data           ###
###                                                                      ###
############################################################################

ECM_MCNM_complete <- function(
    X,
    G,
    py,
    mu,
    Sigma,
    alpha,
    eta,
    eta_min = 1.001,
    eps     = 0.1,
    n_iter  = NULL
) {

  #-------------------------------------#
  #    Objects for the ECM Algorithm    #
  #-------------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  v_tilde     <- matrix(NA, nrow = n, ncol = G)
  w_tilde     <- matrix(NA, nrow = n, ncol = G)
  zw_tilde    <- matrix(NA, nrow = n, ncol = G)
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  dens   <- matrix(NA, nrow = n, ncol = G)
  iter   <- 0
  loglik <- NULL

  #-------------------------#
  #    The ECM Algorithm    #
  #-------------------------#

  if (is.null(n_iter)) {
    n_iter <- Inf
  } else {
    eps <- -Inf
  }

  while (iter < n_iter & getall(loglik) > eps) {

    #++++ E-step ++++#

    for (g in 1:G) {
      z[, g]       <- dCN(X, mu = mu[g, ], Sigma = Sigma[, , g], alpha = alpha[g], eta = eta[g])
      v_tilde[, g] <- alpha[g] * mvtnorm::dmvnorm(X, mean = mu[g, ], sigma = as.matrix(Sigma[, , g])) / z[, g]
      w_tilde[, g] <- v_tilde[, g] + (1 - v_tilde[, g]) / eta[g]
    }

    z_tilde <- sweep(z, 2, py, '*')
    z_tilde <- z_tilde / rowSums(z_tilde)

    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G

    for (g in 1:G) {
      zw_tilde[, g] <- z_tilde[, g] * w_tilde[, g]
    }

    #++++ CM-step 1: pi and alpha ++++#

    N  <- colSums(z_tilde)
    py <- N / n

    alpha              <- colSums(z_tilde * v_tilde) / N
    alpha[alpha < 0.5] <- 0.5
    alpha[alpha > 1]   <- 1

    #++++ CM-step 1: mu ++++#

    for (g in 1:G) {
      mu_num  <- colSums(z_tilde[, g] * w_tilde[, g] * X)
      mu_den  <- sum(z_tilde[, g] * w_tilde[, g])
      mu[g, ] <- mu_num / mu_den
    }

    #++++ CM-step 1: Prepare Sigma tilde ++++#

    for (g in 1:G) {
      X_centered           <- sweep(X, 2, mu[g, ])
      X_centered_crossprod <- apply(X_centered, 1, tcrossprod)
      Sigma_tilde[, , , g] <- array(X_centered_crossprod, dim = c(d, d, n))
    }

    for (g in 1:G) {

      #++++ CM-step 1: Sigma ++++#

      slc_ind      <- slice.index(Sigma_tilde[, , , g, drop = FALSE], 3)
      Sigma_num    <- rowSums(zw_tilde[slc_ind, g] * Sigma_tilde[, , , g, drop = FALSE], dims = 2)
      Sigma[, , g] <- Sigma_num / N[g]

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ CM-step 2: eta ++++#

      delta <- mahalanobis(X, mu[g, ], Sigma[, , g], tol = 1e-20)

      eta_num <- sum(z_tilde[, g] * (1 - v_tilde[, g]) * delta)
      eta_den <- d * sum(z_tilde[, g] * (1 - v_tilde[, g]))
      eta[g]  <- eta_num / eta_den
      eta[g]  <- max(eta[g], eta_min)

    }

    #++++ Observed Log-likelihood ++++#

    for (g in 1:G) {
      dens[, g] <- dCN(X, mu = mu[g, ], Sigma = Sigma[, , g], alpha = alpha[g], eta = eta[g])
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
    alpha        = alpha,
    eta          = eta,
    z_tilde      = z_tilde,
    v_tilde      = v_tilde,
    loglik       = loglik,
    final_loglik = final_loglik
  )

  return(output)

}

###########################################################################
###                                                                     ###
###        Density Function for Contaminated Distribution               ###
###                                                                     ###
###########################################################################

dCN <- function(
    X,
    mu    = rep(0, d),    # location
    Sigma = diag(d),      # dispersion
    alpha = 0.99,         # proportion of good observations
    eta   = 1.01          # degree of contamination
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

  if (alpha < 0 | alpha > 1) {
    stop('alpha must be between 0 and 1')
  }

  if (eta <= 0) {
    stop('eta must be greater than 0')
  }

  good_norm <- exp( mvtnorm::dmvnorm(X, mu, Sigma, log = TRUE) )
  bad_norm  <- exp( mvtnorm::dmvnorm(X, mu, eta * Sigma, log = TRUE) )

  return(
    alpha * good_norm + (1 - alpha) * bad_norm
  )

}

