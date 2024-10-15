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
#' @param G An integer vector specifying the numbers of clusters, which must be at least 1.
#' @param criterion A character string indicating the information criterion for model
#'   selection. "BIC" is used by default. See the details section for a list of available
#'   information criteria.
#' @param max_iter (optional) A numeric value giving the maximum number of
#'   iterations each EM algorithm is allowed to use; 20 by default.
#' @param epsilon (optional) A number specifying the epsilon value for the
#'   Aitken-based stopping criterion used in the EM algorithm: 0.01 by default.
#' @param init_method (optional) A string specifying the method to initialize
#'   the EM algorithm. "kmedoids" clustering is used by default. Alternative
#'   methods include "kmeans", "hierarchical", "mclust", and "manual". When "manual" is chosen,
#'   a vector \code{clusters} of length \eqn{n} must be specified. If the data set is
#'   incomplete, missing values will be first filled based on the mean imputation method.
#' @param clusters (optional) A numeric vector of length \eqn{n} that specifies the initial
#'   cluster memberships of the user when \code{init_method} is set to "manual".
#'   This argument is NULL by default, so that it is ignored whenever other given
#'   initialization methods are chosen.
#' @param eta_min (optional) A numeric value close to 1 to the right specifying
#'   the minimum value of eta; 1.001 by default.
#' @param progress (optional) A logical value indicating whether the
#'   fitting progress should be displayed; TRUE by default.
#'
#' @details Available information criteria include
#' \itemize{
#'  \item AIC - Akaike information criterion
#'   \item BIC - Bayesian information criterion
#'   \item KIC - Kullback information criterion
#'   \item KICc - Corrected Kullback information criterion
#'   \item AIC3 - Modified AIC
#'   \item CAIC - Bozdogan's consistent AIC
#'   \item AICc - Small-sample version of AIC
#'   \item ICL - Integrated Completed Likelihood criterion
#'   \item AWE - Approximate weight of evidence
#'   \item CLC - Classification likelihood criterion
#' }
#' @return An object of class \code{MixtureMissing} with:
#'   \item{model}{The model used to fit the data set.}
#'   \item{pi}{Mixing proportions.}
#'   \item{mu}{Component location vectors.}
#'   \item{Sigma}{Component dispersion matrices.}
#'   \item{alpha}{Component proportions of good observations.}
#'   \item{eta}{Component degrees of contamination.}
#'   \item{z_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
#'     probabilities that the corresponding observation belongs to each cluster.}
#'   \item{v_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
#'     probabilities that the corresponding observation is good with respect
#'     to each cluster.}
#'   \item{clusters}{A numeric vector of length \eqn{n} indicating cluster
#'     memberships determined by the model.}
#'   \item{outliers}{A logical vector of length \eqn{n} indicating observations that are outliers.}
#'   \item{data}{The original data set if it is complete; otherwise, this is
#'     the data set with missing values imputed by appropriate expectations.}
#'   \item{complete}{An \eqn{n} by \eqn{d} logical matrix indicating which cells have no missing values.}
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
#'   \item{ent}{Entropy.}
#'   \item{ICL}{Integrated Completed Likelihood criterion.}
#'   \item{AWE}{Approximate weight of evidence.}
#'   \item{CLC}{Classification likelihood criterion.}
#'   \item{init_method}{The initialization method used in model fitting.}
#'
#' @references
#' Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537. \cr \cr
#' Tong, H. and, Tortora, C., 2022. Model-based clustering and outlier detection
#'   with missing data. \emph{Advances in Data Analysis and Classification}.
#'
#' @examples
#'
#' data('auto')
#'
#' #++++ With no missing values ++++#
#'
#' X <- auto[, c('engine_size', 'city_mpg', 'highway_mpg')]
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#'
#' summary(mod)
#' plot(mod)
#'
#' #++++ With missing values ++++#
#'
#' X <- auto[, c('normalized_losses', 'horsepower', 'highway_mpg', 'price')]
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
    criterion    = c('BIC', 'AIC', 'KIC', 'KICc', 'AIC3', 'CAIC', 'AICc', 'ICL', 'AWE', 'CLC'),
    max_iter     = 20,
    epsilon      = 0.01,
    init_method  = c("kmedoids", "kmeans", "hierarchical", "mclust", "manual"),
    clusters     = NULL,
    eta_min      = 1.001,
    progress     = TRUE
) {

  #----------------------#
  #    Input checking    #
  #----------------------#

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X)) {
    X <- matrix(X, nrow = length(X), ncol = 1)
  }

  if (!is.numeric(X)) {
    stop('X must be a numeric matrix, data frame or vector')
  }

  G <- unique(G)

  if (any(G < 1)) {
    stop('Number of clusters G must be at least 1')
  }

  if (any(G %% 1 != 0)) {
    stop('Number of clusters G must be an integer')
  }

  #---------------------#
  #    Model Fitting    #
  #---------------------#

  criterion <- match.arg(criterion)
  best_info <- Inf
  best_mod  <- NULL

  infos        <- rep(NA_real_, length(G))
  names(infos) <- G

  if ( any(is.na(X)) ) {
    if (ncol(X) < 2) {
      stop('If X contains NAs, X must be at least bivariate')
    }

    if (progress) {
      cat('\nMixture: Contaminated Normal (CN)\n')
      cat('Data Set: Incomplete\n')
    }
  } else {
    if (progress) {
      cat('\nMixture: Contaminated Normal (CN)\n')
      cat('Data Set: Complete\n')
    }
  }

  init_method <- match.arg(init_method)

  if (progress) {
    cat('Initialization:', init_method, '\n\n')
  }

  if (length(G) > 1 & init_method == 'manual') {
    stop('Mannual initialization can only be done if length(G) is 1')
  }

  G_vec <- G
  iter  <- 1

  for (G in G_vec) {

    #++++ Fit each model in G ++++#

    if (progress) {
      cat('Fitting G = ', G, sep = '')
    }

    if (any(is.na(X))) {

      mod <- tryCatch({
        MCNM_incomplete_data(
          X            = X,
          G            = G,
          max_iter     = max_iter,
          epsilon      = epsilon,
          init_method  = init_method,
          clusters     = clusters,
          eta_min      = eta_min,
          progress     = progress
        )
      }, error = function(err) { return(NULL) })

    } else {

      mod <- tryCatch({
        MCNM_complete_data(
          X            = X,
          G            = G,
          max_iter     = max_iter,
          epsilon      = epsilon,
          init_method  = init_method,
          clusters     = clusters,
          eta_min      = eta_min,
          progress     = progress
        )
      }, error = function(err) { return(NULL) })

    }    # end if ( any(is.na(X)) )

    #++++ Compare to the current best model ++++#

    if (!is.null(mod)) {
      infos[iter] <- mod[[criterion]]

      if (best_info > infos[iter]) {
        best_info <- infos[iter]
        best_mod  <- mod
      }
    }

    #++++ Update progress ++++#

    if (progress) {
      if (is.null(mod)) {
        cat(' was failed\n', sep = '')
      } else {
        cat(' was successful with ', mod$iter_stop, '/', max_iter, ' iterations\n', sep = '')
      }
    }

    iter <- iter + 1

  }    # end for (G in G_vec)

  if (progress) {
    cat('\n')
  }

  #--------------------------------------------#
  #    Summarize Results and Prepare Output    #
  #--------------------------------------------#

  if (progress) {

    if (length(G_vec) > 1) {
      if (sum(is.na(infos)) == length(G_vec)) {
        cat('The best mixture model cannot be determined\n')
      } else {
        cat('According to ', criterion, ', the best mixture model is based on G = ', G_vec[which.min(infos)], sep = '')
      }

      cat('\nModel rank according to ', criterion, ':', sep = '')
      infos <- sort(infos, na.last = TRUE)
      for(j in 1:length(infos)){
        cat('\n')
        if (!is.na(infos[j])) {
          cat('  ', j, '. G = ', names(infos)[j], ': ', round(infos[j], digits = 4), sep = '')
        } else {
          cat('  ', j, '. G = ', names(infos)[j], ': Failed', sep = '')
        }
      }
    } else {
      if (!is.na(infos)) {
        cat('The fitted mixture model with G = ', G, ' has ', criterion, ' = ', infos, sep = '')
      }
    }

    cat('\n\n')
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
    init_method  = c("kmedoids", "kmeans", "hierarchical", "mclust", "manual"),
    clusters     = NULL,
    eta_min      = 1.001,
    progress     = TRUE
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

  z           <- matrix(NA_real_, nrow = n, ncol = G)
  z_tilde     <- matrix(NA_real_, nrow = n, ncol = G)
  v_tilde     <- matrix(NA_real_, nrow = n, ncol = G)
  w_tilde     <- matrix(NA_real_, nrow = n, ncol = G)
  zw_tilde    <- matrix(NA_real_, nrow = n, ncol = G)
  X_tilde     <- array(rep(X, G), dim = c(n, d, G))
  Sigma_tilde <- array(NA_real_, dim = c(d, d, n, G))

  py    <- rep(NA_real_, G)
  mu    <- matrix(NA_real_, nrow = G, ncol = d)
  Sigma <- array(NA_real_, dim = c(d, d, G))
  alpha <- rep(0.6, G)
  eta   <- rep(1.4, G)

  log_dens <- matrix(NA_real_, nrow = n, ncol = G)
  iter     <- 0
  loglik   <- NULL

  #--------------------------------#
  #    Parameter Initialization    #
  #--------------------------------#

  init_method <- match.arg(init_method)

  X_imp <- X
  X_imp <- mean_impute(X)

  if (G == 1) {

    max_iter <- 1

    pars  <- cluster_pars(X = X_imp, clusters = rep(1, n))

    py    <- 1
    mu    <- pars$mu
    Sigma <- pars$Sigma

  } else {

    init  <- initialize_clusters(
      X           = X_imp,
      G           = G,
      init_method = init_method,
      clusters    = clusters
    )

    py    <- init$pi
    mu    <- init$mu
    Sigma <- init$Sigma

  }

  #-------------------------#
  #    The ECM Algorithm    #
  #-------------------------#

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

      delta_o <- rep(NA_real_, n)

      for (j in 1:np) {

        m <- M[j, ]    # missing pattern j
        o <- !m        # observed pattern j

        delta_o[Im[[j]]] <- mahalanobis(X[Im[[j]], o, drop = FALSE], mu[g, o], Sigma[o, o, g], tol = 1e-20)

      }

      eta_num <- sum(z_tilde[, g] * (1 - v_tilde[, g]) * delta_o)
      eta_den <- sum(do * z_tilde[, g] * (1 - v_tilde[, g]))
      eta[g]  <- eta_num / eta_den
      eta[g]  <- max(eta[g], eta_min)

      if (eta[g] == eta_min) {
        alpha[g] <- 0.999
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

        log_dens[Im[[j]], g] <- log( dCN(Xo_j, mu = mu_o, Sigma = as.matrix(Sigma_oo), alpha = alpha[g], eta = eta[g]) )

      }
    }

    log_py_dens  <- sweep(log_dens, 2, log(py), FUN = '+')
    final_loglik <- sum( apply(log_py_dens, 1, log_sum_exp) )
    loglik       <- c(loglik, final_loglik)

    #++++ Update Progress ++++#

    iter <- iter + 1

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

  c_names <- paste('comp', 1:G, sep = '')
  v_names <- colnames(X)

  if (is.null(v_names)) {
    v_names <- 1:d
  }

  names(py)       <- c_names
  rownames(mu)    <- c_names
  colnames(mu)    <- v_names
  dimnames(Sigma) <- list(v_names, v_names, c_names)
  names(alpha)    <- c_names
  names(eta)      <- c_names

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
  }

  output <- list(
    model         = 'CN_incomplete_data',
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
  class(output) <- 'MixtureMissing'

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
    init_method  = c("kmedoids", "kmeans", "hierarchical", "mclust", "manual"),
    clusters     = NULL,
    eta_min      = 1.001,
    progress     = TRUE
) {

  #-------------------------------------#
  #    Objects for the ECM Algorithm    #
  #-------------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  z           <- matrix(NA_real_, nrow = n, ncol = G)
  z_tilde     <- matrix(NA_real_, nrow = n, ncol = G)
  v_tilde     <- matrix(NA_real_, nrow = n, ncol = G)
  w_tilde     <- matrix(NA_real_, nrow = n, ncol = G)
  zw_tilde    <- matrix(NA_real_, nrow = n, ncol = G)
  Sigma_tilde <- array(NA_real_, dim = c(d, d, n, G))

  py    <- rep(NA_real_, G)
  mu    <- matrix(NA_real_, nrow = G, ncol = d)
  Sigma <- array(NA_real_, dim = c(d, d, G))
  alpha <- rep(0.6, G)
  eta   <- rep(1.4, G)

  log_dens <- matrix(NA_real_, nrow = n, ncol = G)
  iter     <- 0
  loglik   <- NULL

  #--------------------------------#
  #    Parameter Initialization    #
  #--------------------------------#

  init_method <- match.arg(init_method)

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

  #-------------------------#
  #    The ECM Algorithm    #
  #-------------------------#

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

      if (eta[g] == eta_min) {
        alpha[g] <- 0.999
      }

    }

    #++++ Observed Log-Likelihood ++++#

    for (g in 1:G) {
      log_dens[, g] <- log( dCN(X, mu = mu[g, ], Sigma = Sigma[, , g], alpha = alpha[g], eta = eta[g]) )
    }

    log_py_dens  <- sweep(log_dens, 2, log(py), FUN = '+')
    final_loglik <- sum( apply(log_py_dens, 1, log_sum_exp) )
    loglik       <- c(loglik, final_loglik)

    #++++ Update Progress ++++#

    iter <- iter + 1

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

  c_names <- paste('comp', 1:G, sep = '')
  v_names <- colnames(X)

  if (is.null(v_names)) {
    v_names <- 1:d
  }

  names(py)       <- c_names
  rownames(mu)    <- c_names
  colnames(mu)    <- v_names
  dimnames(Sigma) <- list(v_names, v_names, c_names)
  names(alpha)    <- c_names
  names(eta)      <- c_names

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
  }

  output <- list(
    model         = 'CN_complete_data',
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
  class(output) <- 'MixtureMissing'

  return(output)

}

####################################################################
###                                                              ###
###        Density Function for Contaminated Distribution        ###
###                                                              ###
####################################################################

dCN <- function(
    X,
    mu    = rep(0, d),    # location
    Sigma = diag(d),      # dispersion
    alpha = 0.99,         # proportion of good observations
    eta   = 1.01,         # degree of contamination
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

  if (alpha < 0 | alpha > 1) {
    stop('alpha must be between 0 and 1')
  }

  if (eta <= 0) {
    stop('eta must be greater than 0')
  }

  good_norm <- exp( mvtnorm::dmvnorm(X, mu, Sigma, log = TRUE) )
  bad_norm  <- exp( mvtnorm::dmvnorm(X, mu, eta * Sigma, log = TRUE) )

  dens <- alpha * good_norm + (1 - alpha) * bad_norm
  dens[dens <= 10^(-323)] <- 10^(-323)

  return(dens)

}
