####################################################################
###                                                              ###
###           Multivariate Contaminated Normal Mixture           ###
###                                                              ###
####################################################################

#' Multivariate Contaminated Normal Mixture (MCNM)
#'
#' Carries out model-based clustering using a multivariate contaminated normal
#' mixture (MCNM). The function will determine itself if the data set is
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
#' @param eta_min (optional) A numeric value close to 1 to the right specifying
#'   the minimum value of eta; 1.001 by default.
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
#'   \item{alpha}{Component proportions of good observations.}
#'   \item{eta}{Component degrees of contamination.}
#'   \item{z_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
#'     probabilities that the corresponding observation belongs to each cluster.}
#'   \item{v_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
#'     probabilities that the corresponding observation is outlying with respect
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
#' @import ContaminatedMixt
#' @importFrom stats complete.cases cov cutree dist dnorm hclust kmeans
#'   mahalanobis pchisq rmultinom runif var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
MCNM <- function(
  X,                      # numeric data matrix or data frame
  G,                      # number of clusters
  max_iter = 20,          # maximum number of iterations
  epsilon = 0.01,         # epsilon value for the Aitken convergence criterion
  init_method = c("kmedoids", "kmeans", "hierarchical",
                  "manual", "soft", "hard"),
  equal_prop = FALSE,     # if TRUE, initialize mixing proportions with 1/G
  identity_cov = FALSE,   # if TRUE, initialize sigma with identity matrices
  eta_min = 1.001,        # minimum degree of contamination possible
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
    stop("X must be at least bivariate. For a univariate data set, use function CNM")
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
    mod <- MCNM_incomplete_data(X, G, max_iter, epsilon, init_method, equal_prop,
                                identity_cov, eta_min, show_progress, manual_clusters)
  } else {
    mod <- MCNM_complete_data(X, G, max_iter, epsilon, init_method, equal_prop,
                              identity_cov, eta_min, show_progress, manual_clusters)
  }

  return(mod)
}

########################################################################################
###                                                                                  ###
###           Multivariate Contaminated Normal Mixture for Incomplete Data           ###
###                                                                                  ###
########################################################################################

MCNM_incomplete_data <- function(
  X,                      # numeric data matrix
  G,                      # number of clusters
  max_iter = 20,          # maximum number of iterations
  epsilon = 0.01,         # epsilon value for the Aitken convergence criterion
  init_method = c("kmedoids", "kmeans", "hierarchical",
                  "manual", "soft", "hard"),
  equal_prop = TRUE,      # if TRUE, initialize mixing proportions with 1/G
  identity_cov = TRUE,    # if TRUE, initialize sigma with identity matrices
  eta_min = 1.001,        # minimum degree of contamination possible
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

  alpha <- rep(0.6, G)
  eta   <- rep(1.4, G)

  #---------------------------------------#
  #  Initialization of E-step quantities  #
  #---------------------------------------#

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  v_tilde     <- matrix(NA, nrow = n, ncol = G)
  w_tilde     <- matrix(NA, nrow = n, ncol = G)
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
  #  The ECM algorithm  #
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
        # sigma_oo <- round(sigma_oo, digits = 16)

        z[i, g] <- dCN(xi[o], mu = mu_o, Sigma = sigma_oo, alpha = alpha[g], eta = eta[g])

        v_tilde[i, g] <- alpha[g] * mvtnorm::dmvnorm(xi[o], mean = mu_o, sigma = sigma_oo) / z[i, g]
        if (is.nan(v_tilde[i, g])) {
          v_tilde[i, g] <- 0
        }
      }
    }

    z_tilde <- z * matrix(rep(py, each = n), nrow = n, ncol = G) / matrix(rep(z %*% py, G), nrow = n, ncol = G)
    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G
    N <- colSums(z_tilde)

    for (g in 1:G) {
      w_tilde[, g] <- v_tilde[, g] + (1 - v_tilde[, g]) / eta[g]

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
          sigma_tilde[m, m, i, g] <- tcrossprod(x_ig_tilde - mu_m) + M / w_tilde[i, g]

        } else {
          sigma_tilde[, , i, g] <- tcrossprod(xi - mu[g, ])
        }
      }
    }

    #++++ CM-step 1: pi and alpha ++++#
    py <- N / n
    alpha <- colSums(z_tilde * v_tilde) / N
    alpha[alpha < 0.5] <- 0.5
    alpha[alpha > 1] <- 1

    for (g in 1:G) {

      #++++ CM-step 1: mu and sigma ++++#

      mu_num <- colSums(z_tilde[, g] * w_tilde[, g] * X_tilde[, , g])
      mu_den <- sum(z_tilde[, g] * w_tilde[, g])
      mu[g, ] <- mu_num / mu_den

      zw_tilde     <- array(rep(z_tilde[, g] * w_tilde[, g], each = d * d), dim = c(d, d, n))
      sigma_den    <- apply(zw_tilde * sigma_tilde[, , , g], 1:2, sum)
      sigma[, , g] <- sigma_den / N[g]

      if (isSymmetric.matrix(sigma[, , g]) & max(abs(sigma[, , g] - t(sigma[, , g]))) > .Machine$double.eps) {
        matr <- sigma[, , g]
        matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
        sigma[, , g] <- matr
      }

      # sigma[, , g] <- round(sigma[, , g], digits = 5)
      # max(abs(sigma[, , g] - t(sigma[, , g]))) > .Machine$double.eps

      #++++ CM-step 2: eta ++++#

      delta_o <- apply(X, 1, function(xi) {
        m <- is.na(xi)
        o <- !m
        mahalanobis(xi[o], mu[g, o], cov = sigma[o, o, g])
      })

      eta_num <- sum(z_tilde[, g] * (1 - v_tilde[, g]) * delta_o)
      eta_den <- sum(do * z_tilde[, g] * (1 - v_tilde[, g]))
      eta[g] <- eta_num / eta_den
      eta[g] <- max(eta[g], eta_min)
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

        dens[i, g] <- dCN(xi[o], mu = mu_o, Sigma = sigma_oo, alpha = alpha[g], eta = eta[g])
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

  #-----------------------------#
  #  Perform outlier detection  #
  #-----------------------------#

  outliers <- sapply(1:n, function(i) {
    v_tilde[i, clusters[i]] < 0.5
  })

  #-------------------------------------------------------#
  #  Impute missing values with appropriate expectations  #
  #-------------------------------------------------------#

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
    sigma = G * d * (d - 1) / 2,
    alpha = G,
    eta   = G
  )
  npar$total <- Reduce('+', npar)

  #------------------------#
  #  Information criteria  #
  #------------------------#

  AIC <- -2 * final_loglik + 2 * npar$total
  BIC <- -2 * final_loglik + npar$total * log(n)

  KIC  <- -2 * final_loglik + 3 * (npar$total + 1)
  KICc <- -2 * final_loglik + 2 * (npar$total + 1) * n/(n-npar$total - 2) - n * digamma((n-npar$total)/2) + n * log(n/2)

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
    model        = 'MCNM_incomplete_data',
    pi           = py,
    mu           = mu,
    sigma        = sigma,
    alpha        = alpha,
    eta          = eta,
    z_tilde      = z_tilde,
    v_tilde      = v_tilde,
    clusters     = clusters,
    outliers     = outliers,
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

######################################################################################
###                                                                                ###
###           Multivariate Contaminated Normal Mixture for Complete Data           ###
###                                                                                ###
######################################################################################

MCNM_complete_data <- function(
  X,                      # numeric data matrix
  G,                      # number of clusters
  max_iter = 20,          # maximum number of iterations
  epsilon = 0.01,         # epsilon value for the Aitken convergence criterion
  init_method = c("kmedoids", "kmeans", "hierarchical",
                  "manual", "soft", "hard"),
  equal_prop = FALSE,     # if TRUE, initialize mixing proportions with 1/G
  identity_cov = FALSE,       # if TRUE, initialize sigma with identity matrices
  eta_min = 1.001,        # minimum degree of contamination possible
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

  alpha <- rep(0.6, G)
  eta   <- rep(1.4, G)

  #---------------------------------------#
  #  Initialization of E-step quantities  #
  #---------------------------------------#

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
  v_tilde     <- matrix(NA, nrow = n, ncol = G)
  w_tilde     <- matrix(NA, nrow = n, ncol = G)
  sigma_tilde <-  array(NA, dim = c(d, d, n, G))

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
  #  The ECM algorithm  #
  #---------------------#

  while(iter < max_iter & getall(ll) > epsilon) {

    #++++ E-step ++++#

    for (g in 1:G) {
      z[, g] <- dCN(X, mu = mu[g, ], Sigma = sigma[, , g], alpha = alpha[g], eta = eta[g])
      v_tilde[, g] <- alpha[g] * mvtnorm::dmvnorm(X, mean = mu[g, ], sigma = sigma[, , g]) / z[, g]
      w_tilde[, g] <- v_tilde[, g] + (1 - v_tilde[, g]) / eta[g]

      X_centered           <- sweep(X, 2, mu[g, ])
      X_centered_crossprod <- apply(X_centered, 1, tcrossprod)
      sigma_tilde[, , , g] <- array(X_centered_crossprod, dim = c(d, d, n))
    }

    z_tilde <- z * matrix(rep(py, each = n), nrow = n, ncol = G) / matrix(rep(z %*% py, G), nrow = n, ncol = G)
    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G
    N <- colSums(z_tilde)


    #++++ CM-step 1: pi and alpha ++++#
    py <- N / n
    alpha <- colSums(z_tilde * v_tilde) / N
    alpha[alpha < 0.5] <- 0.5
    alpha[alpha > 1] <- 1

    for (g in 1:G) {

      #++++ CM-step 1: mu and sigma ++++#

      mu_num <- colSums(z_tilde[, g] * w_tilde[, g] * X)
      mu_den <- sum(z_tilde[, g] * w_tilde[, g])
      mu[g, ] <- mu_num / mu_den

      zw_tilde     <- array(rep(z_tilde[, g] * w_tilde[, g], each = d * d), dim = c(d, d, n))
      sigma_den    <- apply(zw_tilde * sigma_tilde[, , , g], 1:2, sum)
      sigma[, , g] <- sigma_den / N[g]

      #++++ CM-step 2: eta ++++#

      delta <- mahalanobis(X, mu[g, ], sigma[, , g])

      eta_num <- sum(z_tilde[, g] * (1 - v_tilde[, g]) * delta)
      eta_den <- d * sum(z_tilde[, g] * (1 - v_tilde[, g]))
      eta[g] <- eta_num / eta_den
      eta[g] <- max(eta[g], eta_min)
    }

    #++++ Obtain observed likelihood and log-likelihood ++++#

    dens <- matrix(NA, nrow = n, ncol = G)
    for (g in 1:G) {
      dens[, g] <- dCN(X, mu = mu[g, ], Sigma = sigma[, , g], alpha = alpha[g], eta = eta[g])
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

  #-----------------------------#
  #  Perform outlier detection  #
  #-----------------------------#

  outliers <- sapply(1:n, function(i) {
    v_tilde[i, clusters[i]] < 0.5
  })

  #------------------------#
  #  Number of parameters  #
  #------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G * d,
    sigma = G * d * (d - 1) / 2,
    alpha = G,
    eta   = G
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
    model        = 'MCNM_complete_data',
    pi           = py,
    mu           = mu,
    sigma        = sigma,
    alpha        = alpha,
    eta          = eta,
    z_tilde      = z_tilde,
    v_tilde      = v_tilde,
    clusters     = clusters,
    outliers     = outliers,
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
