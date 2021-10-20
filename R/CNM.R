###########################################################################
###                                                                     ###
###           Contaminated normal mixture for univariate data           ###
###                                                                     ###
###########################################################################

#' Contaminated Normal Mixture (CNM)
#'
#' Carries out model-based clustering using a contaminated normal
#' mixture (CNM) for complete univariate data set.
#'
#' @param X A vector of \emph{n} observations.
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
#' @param unit_var (optional) A logical value indicating whether variance should
#'  be initialized as 1; FALSE by default.
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
#'   \item{data}{The original data set.}
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
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
#'
#' @examples
#' set.seed(1234)
#'
#' mod <- CNM(iris$Sepal.Length, G = 3, init_method = 'kmedoids', max_iter = 30)
#'
#' plot(mod)
#' summary(mod)
#'
#' @import ContaminatedMixt
#' @export
CNM <- function(
  X,                      # numeric data vector
  G,                      # number of clusters
  max_iter = 20,          # maximum number of iterations
  epsilon = 0.01,         # epsilon value for the Aitken convergence criterion
  init_method = c("kmedoids", "kmeans", "hierarchical",
                  "manual", "soft", "hard"),
  equal_prop = FALSE,     # if TRUE, initialize mixing proportions with 1/G
  unit_var = FALSE,       # if TRUE, initialize sigma with 1
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

  if (!is.vector(X) | !is.numeric(X)) {
    stop('X must be a numeric vector. For a multivariate data set, please use function MCNM')
  }

  if (any(is.na(X))) {
    stop('X must not have missing values')
  }

  #-------------------------------------#
  #  Initialization of the parameters   #
  #-------------------------------------#

  init_method <- match.arg(init_method)
  init_clustering <- initialize_clusters(X, G, init_method, manual_clusters)

  n <- length(X) # number of observations

  mu <- init_clustering$mu

  if (equal_prop) {
    py <- rep(1/G, G)
  } else {
    py <- init_clustering$pi
  }

  if (unit_var) {
    sigma <- rep(1, G)
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
      z[, g] <- dCN(X, mu = mu[g], Sigma = sigma[g], alpha = alpha[g], eta = eta[g])
      v_tilde[, g] <- alpha[g] * dnorm(X, mean = mu[g], sd = sqrt(sigma[g])) / z[, g]
      w_tilde[, g] <- v_tilde[, g] + (1 - v_tilde[, g]) / eta[g]
    }

    z_tilde <- z * matrix(rep(py, each = n), nrow = n, ncol = G) / matrix(rep(z %*% py, G), nrow = n, ncol = G)
    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G
    N <- colSums(z_tilde)


    #++++ CM-step 1: pi, alpha, mu ++++#
    py <- N / n
    alpha <- colSums(z_tilde * v_tilde) / N
    alpha[alpha < 0.5] <- 0.5
    alpha[alpha > 1] <- 1
    mu <- colSums(z_tilde * w_tilde * X) / colSums(z_tilde * w_tilde)

    for (g in 1:G) {

      #++++ CM-step 1: sigma ++++#
      sigma_den <- sum(z_tilde[, g] * w_tilde[, g] * (X - mu[g])^2)
      sigma[g] <- sigma_den / N[g]

      #++++ CM-step 2: eta ++++#

      delta <- (X - mu[g])^2 / sigma[g]

      eta_num <- sum(z_tilde[, g] * (1 - v_tilde[, g]) * delta)
      eta_den <- sum(z_tilde[, g] * (1 - v_tilde[, g]))
      eta[g] <- eta_num / eta_den
      eta[g] <- max(eta[g], eta_min)
    }

    #++++ Obtain observed likelihood and log-likelihood ++++#

    dens <- matrix(NA, nrow = n, ncol = G)
    for (g in 1:G) {
      dens[, g] <- dCN(X, mu = mu[g], Sigma = sigma[g], alpha = alpha[g], eta = eta[g])
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
    cat('Convergence was reached before', max_iter, 'iterations\n')
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
    mu    = G,
    sigma = G,
    alpha = G,
    eta   = G
  )
  npar$total <- Reduce('+', npar)

  #------------------------#
  #  Information criteria  #
  #------------------------#

  AIC <- 2 * final_loglik - 2 * npar$total
  BIC <- 2 * final_loglik - npar$total * log(n)

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
    model        = 'CNM',
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
