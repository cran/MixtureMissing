###########################################################################
###                                                                     ###
###                     Multivariate Skew-t Mixture                     ###
###                                                                     ###
###########################################################################

#' Multivariate Skew-\emph{t} Mixture (MStM)
#'
#' Carries out model-based clustering using a multivariate Skew-\emph{t} mixture (MStM).
#' The function will determine itself if the data set is complete or incomplete and
#' fit the appropriate model accordingly. In the incomplete case, the data set
#' must be at least bivariate, and missing values are assumed to be missing at random (MAR).
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
#' @param df0 (optional) Starting value of component degrees of freedom; 10 by
#'   default.
#' @param deriv_ctrl (optional) A list containing arguments to control the numerical
#'   procedures for calculating the first and second derivatives. Some values are
#'   suggested by default. Refer to functions \code{grad} and \code{hessian} under
#'   the package \code{numDeriv} for more information.
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
#'   \item{mu}{Component mean vectors (location).}
#'   \item{Sigma}{Component covariance matrices (dispersion).}
#'   \item{alpha}{Component skewness vectors.}
#'   \item{df}{Component degrees of freedom.}
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
#' Lin, T.-I. (2010). Robust mixture modeling using multivariate skew \emph{t} distributions.
#'   \emph{Statistics and Computing}, 20(3):343–356.
#' Wang, W.-L. and Lin, T.-I. (2015). Robust model-based clustering via mixtures of
#'   skew-\emph{t} distributions with missing information. \emph{Advances in Data Analysis and Classification},
#'   9(4):423–445.
#' Wei, Y., Tang, Y., and McNicholas, P. D. (2019). Mixtures of generalized hyperbolic
#'   distributions and mixtures of skew-\emph{t} distributions for model-based clustering
#'    with incomplete data. \emph{Computational Statistics & Data Analysis}, 130:18–41.
#'
#' @examples
#'
#' data('bankruptcy')
#'
#' #++++ With no missing values ++++#
#'
#' X <- bankruptcy[, 2:3]
#' mod <- MStM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#'
#' summary(mod)
#' plot(mod)
#'
#' #++++ With missing values ++++#
#'
#' set.seed(1234)
#'
#' X <- hide_values(bankruptcy[, 2:3], prop_cases = 0.1)
#' mod <- MStM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#'
#' summary(mod)
#' plot(mod)
#'
#' @import numDeriv Bessel
#' @importFrom stats complete.cases cov cutree dist dnorm hclust kmeans
#'   mahalanobis pchisq rmultinom runif var uniroot
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
MStM <- function(
    X,
    G,
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "manual", "emEM", "RndEM"),
    clusters     = NULL,
    impute       = FALSE,
    equal_prop   = FALSE,
    identity_cov = FALSE,
    df0          = rep(10, G),
    deriv_ctrl   = list(eps = 1e-8, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-7),
                        r = 6, v = 2, show.details = FALSE),
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

  #---------------------#
  #    Model Fitting    #
  #---------------------#

  if (any(is.na(X))) {

    if (ncol(X) < 2) {
      stop('If X contains NAs, X must be at least bivariate')
    }

    mod <- MStM_incomplete_data(
      X            = X,
      G            = G,
      max_iter     = max_iter,
      epsilon      = epsilon,
      init_method  = init_method,
      clusters     = clusters,
      impute       = impute,
      equal_prop   = equal_prop,
      identity_cov = identity_cov,
      df0          = df0,
      deriv_ctrl   = deriv_ctrl,
      progress     = progress,
      n_run        = n_run,
      n_short      = n_short,
      short_eps    = short_eps
    )

  } else {

    mod <- MStM_complete_data(
      X            = X,
      G            = G,
      max_iter     = max_iter,
      epsilon      = epsilon,
      init_method  = init_method,
      clusters     = clusters,
      equal_prop   = equal_prop,
      identity_cov = identity_cov,
      df0          = df0,
      deriv_ctrl   = deriv_ctrl,
      progress     = progress,
      n_run        = n_run,
      n_short      = n_short,
      short_eps    = short_eps
    )

  }

  return(mod)

}

###########################################################################
###                                                                     ###
###           Multivariate Skew-t Mixture for Incomplete Data           ###
###                                                                     ###
###########################################################################

MStM_incomplete_data <- function(
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
    df0          = rep(10, G),
    deriv_ctrl   = list(eps = 1e-8, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-7),
                        r = 6, v = 2, show.details = FALSE),
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
  X_hat       <- array(rep(X, G), dim = c(n, d, G))
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  a     <- matrix(NA, nrow = n, ncol = G)
  b     <- matrix(NA, nrow = n, ncol = G)
  c     <- matrix(NA, nrow = n, ncol = G)
  N     <- rep(NA, G)
  a_bar <- rep(NA, G)
  b_bar <- rep(NA, G)
  c_bar <- rep(NA, G)
  X_bar <- matrix(NA, nrow = G, ncol = d)

  py     <- rep(NA, G)
  mu     <- matrix(NA, nrow = G, ncol = d)
  Sigma  <- array(NA, dim = c(d, d, G))
  alpha  <- matrix(0.01, nrow = G, ncol = d)
  df     <- df0

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

    pars <- cluster_pars(
      X        = X_imp,
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

          clusters <- sample(1:G, n, prob = py, replace = TRUE)

          if (impute) {
            X_imp <- cluster_impute(X, clusters)
          }

          pars <- cluster_pars(X = X_imp, clusters = clusters)

          mu <- pars$mu

          if (!identity_cov) {
            Sigma <- pars$Sigma
          }

          EM_MStM_incomplete(
            X          = X,
            G          = G,
            py         = py,
            mu         = mu,
            Sigma      = Sigma,
            alpha      = alpha,
            df         = df,
            deriv_ctrl = deriv_ctrl,
            eps        = short_eps,
            n_iter     = n_short
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
      df    <- best_run$df

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
        alpha_o      <- alpha[g, o]

        z[Im[[j]], g] <- dSt(Xo_j, mu = mu_o, Sigma = Sigma_oo, alpha = alpha_o, df = df[g])

        psi <- c( t(alpha_o) %*% Sigma_oo_inv %*% alpha_o )
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

          alpha_o <- alpha[g, o]
          alpha_m <- alpha[g, m]

          for (i in Im[[j]]) {

            xi <- X[i, ]

            mu_m_o    <- mu_m + Sigma_mo %*% Sigma_oo_inv %*% (xi[o] - mu_o)
            Sigma_m_o <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% Sigma_om
            alpha_m_o <- alpha_m - Sigma_mo %*% Sigma_oo_inv %*% alpha_o

            X_hat[i, m, g]   <- mu_m_o + a[i, g] * alpha_m_o
            X_tilde[i, m, g] <- b[i, g] * mu_m_o + alpha_m_o

            Sigma_tilde[m, m, i, g] <- Sigma_m_o + b[i, g] * tcrossprod(mu_m_o) + mu_m_o %*% t(alpha_m_o) + alpha_m_o %*% t(mu_m_o)
            Sigma_tilde[m, m, i, g] <- Sigma_tilde[m, m, i, g] + a[i, g] * tcrossprod(alpha_m_o)

          }

        }

      }
    }

    for (g in 1:G) {
      X_bar[g, ] <- colSums(z_tilde[, g] * X_hat[, , g]) / N[g]
    }

    #++++ M-step: pi ++++#

    py <- N / n

    #++++ M-step: mu (location) and alpha (skewness) ++++#

    for (g in 1:G) {

      num_mu <- colSums( z_tilde[, g] * ( (!R) * (a_bar[g] * b[, g] - 1) * X_tilde[, , g] + R * (a_bar[g] * X_tilde[, , g] - X_hat[, , g]) ) )
      den_mu <- sum( z_tilde[, g] * (a_bar[g] * b[, g] - 1) )

      mu[g, ] <- num_mu / den_mu

      num_alpha <- colSums( z_tilde[, g] * ( (!R) * (b_bar[g] - b[, g]) * X_tilde[, , g] + R * (b_bar[g] * X_hat[, , g] - X_tilde[, , g]) ) )
      den_alpha <- den_mu

      alpha[g, ] <- num_alpha / den_alpha

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

          alpha_o <- alpha[g, o]
          alpha_m <- alpha[g, m]

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
      Sigma[, , g] <- Sigma[, , g] - alpha[g, ] %*% t(X_bar[g, ] - mu[g, ]) - (X_bar[g, ] - mu[g, ]) %*% t(alpha[g, ])
      Sigma[, , g] <- Sigma[, , g] + a_bar[g] * tcrossprod(alpha[g, ])

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ M-step: df (degree of freedom) ++++#

      root <- tryCatch(
        uniroot(df_MSt_func, ws_term = sum(z_tilde[, g] * (c[, g] + b[, g])) / N[g],
                lower = 2, upper = 200)$root,
        error = function(e) { return(NULL) }
      )

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
        alpha_o  <- alpha[g, o]

        dens[Im[[j]], g] <- dSt(Xo_j, mu = mu_o, Sigma = Sigma_oo, alpha = alpha_o, df = df[g])

      }
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

  #------------------#
  #    Imputation    #
  #------------------#

  X_imputed <- X
  complete  <- complete.cases(X)

  for (i in which(!complete)) {
    X_imputed[i, ] <- X_hat[i, , clusters[i]]
  }

  #----------------------------#
  #    Number of Parameters    #
  #----------------------------#

  npar <- list(
    pi     = G - 1,
    mu     = G * d,
    Sigma  = G * d * (d + 1) / 2,
    alpha  = G * d,
    df     = G
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
    alpha <- alpha[1, ]
  }

  output <- list(
    model         = 'MStM_incomplete_data',
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    alpha         = alpha,
    df            = df,
    z_tilde       = z_tilde,
    clusters      = clusters,
    data          = X_imputed,
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

###########################################################################
###                                                                     ###
###          Short EM Iterations for MStM with Incomplete Data          ###
###                                                                     ###
###########################################################################

EM_MStM_incomplete <- function(
    X,
    G,
    py,
    mu,
    Sigma,
    alpha,
    df,
    deriv_ctrl = list(eps = 1e-8, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-7),
                      r = 6, v = 2, show.details = FALSE),
    eps        = 0.1,
    n_iter     = NULL
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
  X_hat       <- array(rep(X, G), dim = c(n, d, G))
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  a     <- matrix(NA, nrow = n, ncol = G)
  b     <- matrix(NA, nrow = n, ncol = G)
  c     <- matrix(NA, nrow = n, ncol = G)
  N     <- rep(NA, G)
  a_bar <- rep(NA, G)
  b_bar <- rep(NA, G)
  c_bar <- rep(NA, G)
  X_bar <- matrix(NA, nrow = G, ncol = d)

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

        mu_o         <- mu[g, o]
        Sigma_oo     <- Sigma[o, o, g]
        Sigma_oo_inv <- solve(Sigma_oo)
        alpha_o      <- alpha[g, o]

        z[Im[[j]], g] <- dSt(Xo_j, mu = mu_o, Sigma = Sigma_oo, alpha = alpha_o, df = df[g])

        psi <- c( t(alpha_o) %*% Sigma_oo_inv %*% alpha_o )
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

          alpha_o <- alpha[g, o]
          alpha_m <- alpha[g, m]

          for (i in Im[[j]]) {

            xi <- X[i, ]

            mu_m_o    <- mu_m + Sigma_mo %*% Sigma_oo_inv %*% (xi[o] - mu_o)
            Sigma_m_o <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% Sigma_om
            alpha_m_o <- alpha_m - Sigma_mo %*% Sigma_oo_inv %*% alpha_o

            X_hat[i, m, g]   <- mu_m_o + a[i, g] * alpha_m_o
            X_tilde[i, m, g] <- b[i, g] * mu_m_o + alpha_m_o

            Sigma_tilde[m, m, i, g] <- Sigma_m_o + b[i, g] * tcrossprod(mu_m_o) + mu_m_o %*% t(alpha_m_o) + alpha_m_o %*% t(mu_m_o)
            Sigma_tilde[m, m, i, g] <- Sigma_tilde[m, m, i, g] + a[i, g] * tcrossprod(alpha_m_o)

          }

        }

      }
    }

    for (g in 1:G) {
      X_bar[g, ] <- colSums(z_tilde[, g] * X_hat[, , g]) / N[g]
    }

    #++++ M-step: pi ++++#

    py <- N / n

    #++++ M-step: mu (location) and alpha (skewness) ++++#

    for (g in 1:G) {

      num_mu <- colSums( z_tilde[, g] * ( (!R) * (a_bar[g] * b[, g] - 1) * X_tilde[, , g] + R * (a_bar[g] * X_tilde[, , g] - X_hat[, , g]) ) )
      den_mu <- sum( z_tilde[, g] * (a_bar[g] * b[, g] - 1) )

      mu[g, ] <- num_mu / den_mu

      num_alpha <- colSums( z_tilde[, g] * ( (!R) * (b_bar[g] - b[, g]) * X_tilde[, , g] + R * (b_bar[g] * X_hat[, , g] - X_tilde[, , g]) ) )
      den_alpha <- den_mu

      alpha[g, ] <- num_alpha / den_alpha

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

          alpha_o <- alpha[g, o]
          alpha_m <- alpha[g, m]

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
      Sigma[, , g] <- Sigma[, , g] - alpha[g, ] %*% t(X_bar[g, ] - mu[g, ]) - (X_bar[g, ] - mu[g, ]) %*% t(alpha[g, ])
      Sigma[, , g] <- Sigma[, , g] + a_bar[g] * tcrossprod(alpha[g, ])

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ M-step: df (degree of freedom) ++++#

      root <- tryCatch(
        uniroot(df_MSt_func, ws_term = sum(z_tilde[, g] * (c[, g] + b[, g])) / N[g],
                lower = 2, upper = 200)$root,
        error = function(e) { return(NULL) }
      )

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
        alpha_o  <- alpha[g, o]

        dens[Im[[j]], g] <- dSt(Xo_j, mu = mu_o, Sigma = Sigma_oo, alpha = alpha_o, df = df[g])

      }
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
    df           = df,
    z_tilde      = z_tilde,
    loglik       = loglik,
    final_loglik = final_loglik
  )

  return(output)

}

###########################################################################
###                                                                     ###
###            Multivariate Skew-t Mixture for Complete Data            ###
###                                                                     ###
###########################################################################

MStM_complete_data <- function(
    X,
    G,
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "manual", "emEM", "RndEM"),
    clusters     = NULL,
    equal_prop   = FALSE,
    identity_cov = FALSE,
    df0          = rep(10, G),
    deriv_ctrl   = list(eps = 1e-8, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-7),
                         r = 6, v = 2, show.details = FALSE),
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

  a     <- matrix(NA, nrow = n, ncol = G)
  b     <- matrix(NA, nrow = n, ncol = G)
  c     <- matrix(NA, nrow = n, ncol = G)
  N     <- rep(NA, G)
  a_bar <- rep(NA, G)
  b_bar <- rep(NA, G)
  c_bar <- rep(NA, G)
  X_bar <- matrix(NA, nrow = G, ncol = d)

  py    <- rep(NA, G)
  mu    <- matrix(NA, nrow = G, ncol = d)
  Sigma <- array(NA, dim = c(d, d, G))
  alpha <- matrix(0.01, nrow = G, ncol = d)
  df    <- df0

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

          EM_MStM_complete(
            X          = X,
            G          = G,
            py         = py,
            mu         = mu,
            Sigma      = Sigma,
            alpha      = alpha,
            df         = df,
            deriv_ctrl = deriv_ctrl,
            eps        = short_eps,
            n_iter     = n_short
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
      df    <- best_run$df

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
      z[, g] <- dSt(X, mu = mu[g, ], Sigma = Sigma[, , g], alpha = alpha[g, ], df = df[g])

      Sigma_inv <- solve(Sigma[, , g])

      psi <- c( t(alpha[g, ]) %*% Sigma_inv %*% alpha[g, ] )
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

      #++++ M-step: mu (location) and alpha (skewness) ++++#

      num_mu <- colSums( z_tilde[, g] * X * (a_bar[g] * b[, g] - 1) )
      den_mu <- sum( z_tilde[, g] * (a_bar[g] * b[, g] - 1) )

      mu[g, ] <- num_mu / den_mu

      num_alpha <- colSums( z_tilde[, g] * X * (b_bar[g] - b[, g]) )
      den_alpha <- den_mu

      alpha[g, ] <- num_alpha / den_alpha

      #++++ M-step: Sigma (dispersion) ++++#

      X_centrd         <- sweep(X, 2, mu[g, ])
      X_centrd_crsprod <- apply(X_centrd, 1, tcrossprod)
      Sigma_tilde      <- array(X_centrd_crsprod, dim = c(d, d, n))

      slc_ind      <- slice.index(Sigma_tilde, 3)
      Sigma[, , g] <- rowSums(z_tilde[slc_ind, g] * b[slc_ind, g] * Sigma_tilde, dims = 2) / N[g]
      Sigma[, , g] <- Sigma[, , g] - alpha[g, ] %*% t(X_bar[g, ] - mu[g, ]) - (X_bar[g, ] - mu[g, ]) %*% t(alpha[g, ])
      Sigma[, , g] <- Sigma[, , g] + a_bar[g] * tcrossprod(alpha[g, ])

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ M-step: df (degree of freedom) ++++#

      root <- tryCatch(
        uniroot(df_MSt_func, ws_term = sum(z_tilde[, g] * (c[, g] + b[, g])) / N[g],
                lower = 2, upper = 200)$root,
        error = function(e) { return(NULL) }
      )

      df[g] <- ifelse(!is.null(root), root, df[g])

    }

    #++++ Observed Log-likelihood ++++#

    for (g in 1:G) {
      dens[, g] <- dSt(X, mu = mu[g, ], Sigma = Sigma[, , g], alpha = alpha[g, ], df = df[g])
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

  #----------------------------#
  #    Number of Parameters    #
  #----------------------------#

  npar <- list(
    pi     = G - 1,
    mu     = G * d,
    Sigma  = G * d * (d + 1) / 2,
    alpha  = G * d,
    df     = G
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
    alpha <- alpha[1, ]
  }

  output <- list(
    model         = 'MStM_complete_data',
    pi            = py,
    mu            = mu,
    Sigma         = Sigma,
    alpha         = alpha,
    df            = df,
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

###########################################################################
###                                                                     ###
###           Short EM Iterations for MStM with Complete Data           ###
###                                                                     ###
###########################################################################

EM_MStM_complete <- function(
    X,
    G,
    py,
    mu,
    Sigma,
    alpha,
    df,
    deriv_ctrl = list(eps = 1e-8, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-7),
                      r = 6, v = 2, show.details = FALSE),
    eps        = 0.1,
    n_iter     = NULL
) {

  #------------------------------------#
  #    Objects for the EM Algorithm    #
  #------------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  Sigma_tilde <- array(NA, dim = c(d, d, n, G))

  a     <- matrix(NA, nrow = n, ncol = G)
  b     <- matrix(NA, nrow = n, ncol = G)
  c     <- matrix(NA, nrow = n, ncol = G)
  N     <- rep(NA, G)
  a_bar <- rep(NA, G)
  b_bar <- rep(NA, G)
  c_bar <- rep(NA, G)
  X_bar <- matrix(NA, nrow = G, ncol = d)

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
      z[, g] <- dSt(X, mu = mu[g, ], Sigma = Sigma[, , g], alpha = alpha[g, ], df = df[g])

      Sigma_inv <- solve(Sigma[, , g])

      psi <- c( t(alpha[g, ]) %*% Sigma_inv %*% alpha[g, ] )
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

      #++++ M-step: mu (location) and alpha (skewness) ++++#

      num_mu <- colSums( z_tilde[, g] * X * (a_bar[g] * b[, g] - 1) )
      den_mu <- sum( z_tilde[, g] * (a_bar[g] * b[, g] - 1) )

      mu[g, ] <- num_mu / den_mu

      num_alpha <- colSums( z_tilde[, g] * X * (b_bar[g] - b[, g]) )
      den_alpha <- den_mu

      alpha[g, ] <- num_alpha / den_alpha

      #++++ M-step: Sigma (dispersion) ++++#

      X_centrd         <- sweep(X, 2, mu[g, ])
      X_centrd_crsprod <- apply(X_centrd, 1, tcrossprod)
      Sigma_tilde      <- array(X_centrd_crsprod, dim = c(d, d, n))

      slc_ind      <- slice.index(Sigma_tilde, 3)
      Sigma[, , g] <- rowSums(z_tilde[slc_ind, g] * b[slc_ind, g] * Sigma_tilde, dims = 2) / N[g]
      Sigma[, , g] <- Sigma[, , g] - alpha[g, ] %*% t(X_bar[g, ] - mu[g, ]) - (X_bar[g, ] - mu[g, ]) %*% t(alpha[g, ])
      Sigma[, , g] <- Sigma[, , g] + a_bar[g] * tcrossprod(alpha[g, ])

      if (max(abs(Sigma[, , g] - t(Sigma[, , g]))) > .Machine$double.eps) {
        matr <- Sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        Sigma[, , g] <- matr
      }

      #++++ M-step: df (degree of freedom) ++++#

      root <- tryCatch(
        uniroot(df_MSt_func, ws_term = sum(z_tilde[, g] * (c[, g] + b[, g])) / N[g],
                lower = 2, upper = 200)$root,
        error = function(e) { return(NULL) }
      )

      df[g] <- ifelse(!is.null(root), root, df[g])

    }

    #++++ Observed Log-likelihood ++++#

    for (g in 1:G) {
      dens[, g] <- dSt(X, mu = mu[g, ], Sigma = Sigma[, , g], alpha = alpha[g, ], df = df[g])
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
    df           = df,
    z_tilde      = z_tilde,
    loglik       = loglik,
    final_loglik = final_loglik
  )

  return(output)

}

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
    alpha = rep(0, d),    # skewness
    df    = 10            # degree of freedom
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

  if (length(alpha) != d) {
    stop('alpha must be a vector of length d')
  }

  X_centrd  <- sweep(X, 2, mu, '-')
  Sigma_inv <- solve(Sigma)

  psi <- c( t(alpha) %*% Sigma_inv %*% alpha )
  chi <- df + mahalanobis(X, center = mu, cov = Sigma_inv, inverted = TRUE)

  s1 <- sqrt(psi * chi)

  res <- (-df/2 - d/2) * log(s1) + (df/2 + d/2) * log(psi)
  res <- res + (df / 2) * log(df) + log(besselK(s1, nu = -df/2 - d/2, expon.scaled = TRUE)) - s1

  if (is.nan(log(det(Sigma)))) {
    Sigma <- diag(d)
  }

  res <- res - d/2 * (log(2) + log(pi)) - 1/2 * log(det(Sigma))
  res <- res - lgamma(df / 2) - (df/2 - 1) * log(2)
  res <- res + X_centrd %*% Sigma_inv %*% alpha

  return( c(exp(res)) )

  # lvx <- (df / 2) * log(df) + (-df/2 - d/2) * log(s1)
  # lvx <- lvx + log(besselK(s1, nu = -df/2 - d/2, expon.scaled = TRUE)) - s1
  # lvx <- lvx + X_centrd %*% Sigma_inv %*% alpha
  #
  # lv <- -1/2 * log(det(Sigma)) - d/2 * (log(2) + log(pi))
  # lv <- lv - log(gamma(df / 2)) - (df/2 - 1) * log(2)
  # lv <- lv + (df/2 + d/2) * log(psi)
  #
  # return(
  #   c( exp(lvx + lv) )
  # )

}
