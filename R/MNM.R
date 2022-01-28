#######################################################
###                                                 ###
###           Multivariate Normal Mixture           ###
###                                                 ###
#######################################################

#' Multivariate Normal Mixture (MNM)
#'
#' Carries out model-based clustering using a multivariate normal
#' mixture (MNM). The function will determine itself if the data set is
#' complete or incomplete and fit the appropriate model accordingly. When using
#' this function, the data set must be at least bivariate, and missing values
#' are assumed to be missing at random (MAR).
#'
#' @param X An \eqn{n} by \eqn{d} matrix or data frame where \eqn{n} is the number of
#'   observations and \eqn{d} is the number of columns or variables.
#' @param G The number of clusters.
#' @param max_iter (optional) A numeric value giving the maximum number of
#'   iterations each EM algorithm is allowed to use; 20 by default.
#' @param epsilon (optional) A number specifying the epsilon value for the
#'   Aitken-based stopping criterion used in the EM algorithm: 0.01 by default.
#' @param init_method (optional) A string specifying the method to initialize
#'   the EM algorithm. "kmedoids" clustering is used by default. Alternative
#'   methods include "kmeans", "hierarchical", "manual", "soft", "hard". When
#'   "manual" is chosen, a vector \code{manual_clusters} of length \eqn{n} must
#'   be specified.
#' @param equal_prop (optional) A logical value indicating whether mixing
#'   proportions should be equal at initialization of the EM algorithm; FALSE by
#'   default.
#' @param identity_cov (optional) A logical value indicating whether covariance
#'   matrices should be initialized as identity matrices; FALSE by default.
#' @param show_progress (optional) A logical value indicating whether the
#'   fitting progress should be displayed; TRUE by default.
#' @param manual_clusters A vector of length \eqn{n} that specifies the initial
#'   cluster memberships of the user when \code{init_method} is set to "manual".
#'   Both numeric and character vectors are acceptable. This argument is NULL by
#'   default, so that it is ignored whenever other given initialization methods
#'   are chosen.
#'
#' @return An object of class \code{MixtureMissing} with:
#'   \item{model}{The model used to fit the data set}
#'   \item{pi}{Mixing proportions.}
#'   \item{mu}{Component mean vectors.}
#'   \item{sigma}{Component covariance matrices.}
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
#'   \item{final_lik}{The final value of likelihood.}
#'   \item{final_loglik}{The final value of log-likelihood.}
#'   \item{lik}{All the values of likelihood.}
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
  X,                      # numeric data matrix or data frame
  G,                      # number of clusters
  max_iter = 20,          # maximum number of iterations
  epsilon = 0.01,         # epsilon value for the Aitken convergence criterion
  init_method = c("kmedoids", "kmeans", "hierarchical",
                  "manual", "soft", "hard"),
  equal_prop = FALSE,     # if TRUE, initialize mixing proportions with 1/G
  identity_cov = FALSE,   # if TRUE, initialize sigma with identity matrices
  show_progress = TRUE,   # if TRUE, show the fitting progress
  manual_clusters = NULL  # cluster memberships specified by the user if init_method = "manual"
) {

  #-------------------#
  #  Input checking   #
  #-------------------#

  if (is.null(G)) {
    stop('Number of clusters G must be specified')
  }

  if (G < 2) {
    stop('Number of clusters G must be at least 2')
  }

  if (G %% 1 != 0) {
    stop('Number of clusters G must be an integer')
  }

  if(is.null(ncol(X)) | ncol(X) < 2) {
    stop("X must be at least bivariate. For a univariate data set, please use function NM.")
  }

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.numeric(X)) {
    stop('X must be a numeric matrix')
  }

  #--------------------------------------------------------------------------#
  #  Fit the appropriate model according to the presence of missing values   #
  #--------------------------------------------------------------------------#

  if (any(is.na(X))) {
    mod <- MNM_incomplete_data(X, G, max_iter, epsilon, init_method, equal_prop,
                               identity_cov, show_progress, manual_clusters)
  } else {
    mod <- MNM_complete_data(X, G, max_iter, epsilon, init_method, equal_prop,
                             identity_cov, show_progress, manual_clusters)
  }

  return(mod)
}

######################################################################
###                                                                ###
###           Multivariate t Mixture for Incomplete Data           ###
###                                                                ###
######################################################################

MNM_incomplete_data <- function(
  X,                      # numeric data matrix or data frame
  G,                      # number of clusters
  max_iter = 20,          # maximum number of iterations
  epsilon = 0.01,         # epsilon value for the Aitken convergence criterion
  init_method = c("kmedoids", "kmeans", "hierarchical",
                  "manual", "soft", "hard"),
  equal_prop = TRUE,      # if TRUE, initialize mixing proportions with 1/G
  identity_cov = TRUE,        # if TRUE, initialize sigma with identity matrices
  show_progress = TRUE,   # if TRUE, show the fitting progress
  manual_clusters = NULL  # cluster memberships specified by the user if init_method = "manual"
) {

  #-------------------------------------#
  #  Initialization of the parameters   #
  #-------------------------------------#

  init_method <- match.arg(init_method)
  init_clustering <- initialize_clusters(X, G, init_method, manual_clusters)

  n <- nrow(X) # number of observations
  d <- ncol(X) # number of variables/dimensions

  do <- apply(X, 1, function(xi) { # number of observed values per observation
    m <- is.na(xi)
    o <- !m
    sum(o)
  })

  mu <- init_clustering$mu

  if (equal_prop) {
    py <- rep(1/G, G)
  } else {
    py <- init_clustering$pi
  }

  if (identity_cov) {
    sigma <- array(rep(diag(d), G), dim = c(d, d, G))
  } else {
    sigma <- init_clustering$sigma
  }

  #---------------------------------------#
  #  Initialization of E-step quantities  #
  #---------------------------------------#

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  X_tilde     <- array(rep(X, G), dim = c(n, d, G))
  sigma_tilde <- array(NA, dim = c(d, d, n, G))

  #------------------------#
  #  Other initialization  #
  #------------------------#

  iter <- 0    # current iteration
  l    <- NULL # likelihood over iterations
  ll   <- NULL # log-likelihood over iterations

  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3, width = 75, char = "=")
  }

  #---------------------#
  #  The EM algorithm  #
  #---------------------#

  while(iter < max_iter & getall(ll) > epsilon) {

    #++++ E-step ++++#

    for (g in 1:G) {
      for (i in 1:n) {
        xi <- X[i, ]
        m <- is.na(xi)
        o <- !m

        mu_o <- mu[g, o]
        sigma_oo <- as.matrix(sigma[o, o, g])

        z[i, g] <- mvtnorm::dmvnorm(xi[o], mean = mu[g, o], sigma = sigma_oo)
      }
    }

    z_tilde <- z * matrix(rep(py, each = n), nrow = n, ncol = G) / matrix(rep(z %*% py, G), nrow = n, ncol = G)
    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G
    N <- colSums(z_tilde)

    for (g in 1:G) {
      for (i in 1:n) {
        xi <- X[i, ]
        m <- is.na(xi)

        if (any(m)) {
          o <- !m

          mu_m <- mu[g, m]
          mu_o <- mu[g, o]

          sigma_mo <- sigma[m, o, g]
          sigma_om <- sigma[o, m, g]
          sigma_mm <- sigma[m, m, g]
          sigma_oo_inv <- mnormt::pd.solve(sigma[o, o, g])
          # sigma_oo_inv <- solve(sigma[o, o, g])

          x_ig_tilde <- mu_m + sigma_mo %*% sigma_oo_inv %*% (xi[o] - mu_o)
          X_tilde[i, m, g] <- x_ig_tilde

          sigma_tilde[o, o, i, g] <- tcrossprod(xi[o] - mu_o)
          sigma_tilde[o, m, i, g] <- tcrossprod(xi[o] - mu_o, x_ig_tilde - mu_m)
          sigma_tilde[m, o, i, g] <- t(sigma_tilde[o, m, i, g])
          M <- sigma_mm - sigma_mo %*% sigma_oo_inv %*% sigma_om
          sigma_tilde[m, m, i, g] <- tcrossprod(x_ig_tilde - mu_m) + M

        } else {
          sigma_tilde[, , i, g] <- tcrossprod(xi - mu[g, ])
        }
      }
    }

    #++++ M-step: pi ++++#
    py <- N / n

    for (g in 1:G) {

      #++++ M-step: mu and sigma ++++#

      mu[g, ] <- colSums(z_tilde[, g] * X_tilde[, , g]) / N[g]

      zz_tilde <- array(rep(z_tilde[, g], each = d * d), dim = c(d, d, n))
      sigma_den <- apply(zz_tilde * sigma_tilde[, , , g], 1:2, sum)
      sigma[, , g] <- sigma_den / N[g]

      if (isSymmetric.matrix(sigma[, , g]) & max(abs(sigma[, , g] - t(sigma[, , g]))) > .Machine$double.eps) {
        matr <- sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        sigma[, , g] <- matr
      }
    }

    #++++ Obtain observed likelihood and log-likelihood ++++#

    dens <- matrix(NA, nrow = n, ncol = G)
    for (g in 1:G) {
      for (i in 1:n) {
        xi <- X[i, ]
        m <- is.na(xi)
        o <- !m

        mu_o <- mu[g, o]
        sigma_oo <- as.matrix(sigma[o, o, g])

        dens[i, g] <- mvtnorm::dmvnorm(xi[o], mean = mu[g, o], sigma = sigma_oo)
      }
    }

    lik <- dens %*% py
    lik[lik <= 10^(-323)] <- 10^(-323)
    final_lik    <- prod(lik)
    final_loglik <- sum(log(lik))
    l  <- c(l, final_lik)
    ll <- c(ll, final_loglik)

    #++++ Update progress ++++#

    iter <- iter + 1

    if (show_progress) {
      setTxtProgressBar(pb, iter)
      cat(' Iteration', iter, '/', max_iter)
    }
  }

  if (show_progress) {
    close(pb)
  }

  if (iter < max_iter) {
    cat('\nConvergence was reached before', max_iter, 'iterations\n')
  }

  #------------------------------#
  #  Obtain cluster memberships  #
  #------------------------------#

  clusters <- apply(z_tilde, 1, which.max)

  #-------------------------------------------------------#
  #  Impute missing values with appropriate expectations  #
  #-------------------------------------------------------#

  # X_imputed <- X
  # for (g in 1:G) {
  #   X_imputed[clusters == g, ] <- X_tilde[, , g]
  # }

  X_imputed <- X
  for (i in which(!complete.cases(X))) {
    xi <- X[i, ]
    m <- is.na(xi)
    o <- !m

    mu_m <- mu[clusters[i], m]
    mu_o <- mu[clusters[i], o]
    sigma_mo <- sigma[m, o, clusters[i]]
    sigma_om <- sigma[o, m, clusters[i]]
    sigma_mm <- sigma[m, m, clusters[i]]
    sigma_oo_inv <- solve(sigma[o, o, clusters[i]])

    x_ig_tilde <- mu_m + sigma_mo %*% sigma_oo_inv %*% (xi[o] - mu_o)
    X_imputed[i, m] <- x_ig_tilde
  }

  #------------------------#
  #  Number of parameters  #
  #------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G * d,
    sigma = G * d * (d - 1) / 2
  )
  npar$total <- Reduce('+', npar)

  #------------------------#
  #  Information criteria  #
  #------------------------#

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

  #-------------------#
  #  Prepare outputs  #
  #-------------------#

  outputs <- list(
    model        = 'MNM_incomplete_data',
    pi           = py,
    mu           = mu,
    sigma        = sigma,
    z_tilde      = z_tilde,
    clusters     = clusters,
    data         = X_imputed,
    complete     = complete.cases(X),
    npar         = npar,
    max_iter     = max_iter,
    iter_stop    = iter,
    final_lik    = final_lik,
    final_loglik = final_loglik,
    lik          = l,
    loglik       = ll,
    AIC          = AIC,
    BIC          = BIC,
    KIC          = KIC,
    KICc         = KICc,
    AIC3         = AIC3,
    CAIC         = CAIC,
    AICc         = AICc,
    ent          = ent,
    ICL          = ICL,
    AWE          = AWE,
    CLC          = CLC,
    init_method  = init_method
  )
  class(outputs) <- 'MixtureMissing'

  return(outputs)
}


#########################################################################
###                                                                   ###
###           Multivariate Normla Mixture for Complete Data           ###
###                                                                   ###
#########################################################################

MNM_complete_data <- function(
  X,                      # numeric data matrix
  G,                      # number of clusters
  max_iter = 20,          # maximum number of iterations
  epsilon = 0.01,         # epsilon value for the Aitken convergence criterion
  init_method = c("kmedoids", "kmeans", "hierarchical",
                  "manual", "soft", "hard"),
  equal_prop = TRUE,      # if TRUE, initialize mixing proportions with 1/G
  identity_cov = TRUE,    # if TRUE, initialize sigma with identity matrices
  show_progress = TRUE,   # if TRUE, show the fitting progress
  manual_clusters = NULL  # cluster memberships specified by the user if init_method = "manual"
) {

  #-------------------------------------#
  #  Initialization of the parameters   #
  #-------------------------------------#

  init_method <- match.arg(init_method)
  init_clustering <- initialize_clusters(X, G, init_method, manual_clusters)

  n <- nrow(X) # number of observations
  d <- ncol(X) # number of variables/dimensions

  mu <- init_clustering$mu

  if (equal_prop) {
    py <- rep(1/G, G)
  } else {
    py <- init_clustering$pi
  }

  if (identity_cov) {
    sigma <- array(rep(diag(d), G), dim = c(d, d, G))
  } else {
    sigma <- init_clustering$sigma
  }

  #---------------------------------------#
  #  Initialization of E-step quantities  #
  #---------------------------------------#

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  X_tilde     <- array(rep(X, G), dim = c(n, d, G))
  sigma_tilde <- array(NA, dim = c(d, d, n, G))

  #------------------------#
  #  Other initialization  #
  #------------------------#

  iter <- 0    # current iteration
  l    <- NULL # likelihood over iterations
  ll   <- NULL # log-likelihood over iterations

  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3, width = 75, char = "=")
  }

  #---------------------#
  #  The EM algorithm  #
  #---------------------#

  while(iter < max_iter & getall(ll) > epsilon) {

    #++++ E-step ++++#

    for (g in 1:G) {
      z[, g] <- mvtnorm::dmvnorm(X, mean = mu[g, ], sigma = sigma[, , g])

      X_centered           <- sweep(X, 2, mu[g, ])
      X_centered_crossprod <- apply(X_centered, 1, tcrossprod)
      sigma_tilde[, , , g] <- array(X_centered_crossprod, dim = c(d, d, n))
    }

    z_tilde <- z * matrix(rep(py, each = n), nrow = n, ncol = G) / matrix(rep(z %*% py, G), nrow = n, ncol = G)
    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G
    N <- colSums(z_tilde)

    #++++ M-step: pi ++++#
    py <- N / n

    for (g in 1:G) {

      #++++ M-step: mu and sigma ++++#

      mu[g, ] <- colSums(z_tilde[, g] * X_tilde[, , g]) / N[g]

      zz_tilde <- array(rep(z_tilde[, g], each = d * d), dim = c(d, d, n))
      sigma_den <- apply(zz_tilde * sigma_tilde[, , , g], 1:2, sum)
      sigma[, , g] <- sigma_den / N[g]
    }

    #++++ Obtain observed likelihood and log-likelihood ++++#

    dens <- matrix(NA, nrow = n, ncol = G)
    for (g in 1:G) {
      dens[, g] <- mvtnorm::dmvnorm(X, mean = mu[g, ], sigma = sigma[, , g])
    }

    lik <- dens %*% py
    lik[lik <= 10^(-323)] <- 10^(-323)
    final_lik    <- prod(lik)
    final_loglik <- sum(log(lik))
    l  <- c(l, final_lik)
    ll <- c(ll, final_loglik)

    #++++ Update progress ++++#

    iter <- iter + 1

    if (show_progress) {
      setTxtProgressBar(pb, iter)
      cat(' Iteration', iter, '/', max_iter)
    }
  }

  if (show_progress) {
    close(pb)
  }

  if (iter < max_iter) {
    cat('\nConvergence was reached before', max_iter, 'iterations\n')
  }

  #------------------------------#
  #  Obtain cluster memberships  #
  #------------------------------#

  clusters <- apply(z_tilde, 1, which.max)

  #------------------------#
  #  Number of parameters  #
  #------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G * d,
    sigma = G * d * (d - 1) / 2
  )
  npar$total <- Reduce('+', npar)

  #------------------------#
  #  Information criteria  #
  #------------------------#

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

  #-------------------#
  #  Prepare outputs  #
  #-------------------#

  outputs <- list(
    model        = 'MNM_complete_data',
    pi           = py,
    mu           = mu,
    sigma        = sigma,
    z_tilde      = z_tilde,
    clusters     = clusters,
    data         = X,
    complete     = complete.cases(X),
    npar         = npar,
    max_iter     = max_iter,
    iter_stop    = iter,
    final_lik    = final_lik,
    final_loglik = final_loglik,
    lik          = l,
    loglik       = ll,
    AIC          = AIC,
    BIC          = BIC,
    KIC          = KIC,
    KICc         = KICc,
    AIC3         = AIC3,
    CAIC         = CAIC,
    AICc         = AICc,
    ent          = ent,
    ICL          = ICL,
    AWE          = AWE,
    CLC          = CLC,
    init_method  = init_method
  )
  class(outputs) <- 'MixtureMissing'

  return(outputs)

}
