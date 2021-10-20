#########################################################
###                                                   ###
###           t mixture for univariate data           ###
###                                                   ###
#########################################################

#' \emph{t} Mixture (\emph{t}M)
#'
#' Carries out model-based clustering using a \emph{t} mixture (\emph{t}M) for
#' complete univariate data set.
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
#' @param df0 (optional) Starting values of the degrees of freedom; 10 for all
#'  clusters by default.
#' @param outlier_cutoff (optional) A percentile for outlier detection; 0.95 by default.
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
#'   \item{df}{Component degrees of freedom.}
#'   \item{z_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
#'     probabilities that the corresponding observation belongs to each cluster.}
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
#' Peel,  D.  and  McLachlan,  G.  J.  (2000).   Robust  mixture  modelling
#'   using the \emph{t}  distribution. \emph{Statistics and computing, 10}(4):339-348.
#'
#' @examples
#' set.seed(1234)
#'
#' mod <- tM(iris$Sepal.Length, G = 3, init_method = 'kmedoids', max_iter = 30)
#'
#' plot(mod)
#' summary(mod)
#'
#' @import mvtnorm
#' @export
tM <- function(
  X,                      # numeric data matrix
  G,                      # number of clusters
  max_iter = 20,          # maximum number of iterations
  epsilon = 0.01,         # epsilon value for the Aitken convergence criterion
  init_method = c("kmedoids", "kmeans", "hierarchical",
                  "manual", "soft", "hard"),
  equal_prop = TRUE,      # if TRUE, initialize mixing proportions with 1/G
  unit_var = FALSE,       # if TRUE, initialize sigma with 1
  df0 = rep(10, G),       # starting values of degrees of freedom
  outlier_cutoff = 0.95,  # a percentile for outlier detection
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
    stop('X must be a numeric vector. For a multivariate data set, please use function MtM.')
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

  if (!is.vector(df0) | !is.numeric(df0)) {
    stop('df0 must be a numeric vector')
  }

  if (length(df0) != G) {
    stop('df0 must have length G')
  }

  if (any(df0 < 0)) {
    stop('df0 must be positive')
  }

  df <- df0

  if (length(outlier_cutoff) > 1 | !is.numeric(outlier_cutoff)) {
    stop('outlier_cutoff must be a number')
  }

  if (outlier_cutoff <= 0 | outlier_cutoff >= 1) {
    stop('outlier_cutoff must be in (0, 1)')
  }

  #---------------------------------------#
  #  Initialization of E-step quantities  #
  #---------------------------------------#

  z           <- matrix(NA, nrow = n, ncol = G)
  z_tilde     <- matrix(NA, nrow = n, ncol = G)
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
  #  The EM algorithm  #
  #---------------------#

  while(iter < max_iter & getall(ll) > epsilon) {

    #++++ E-step ++++#

    for (g in 1:G) {
      z[, g] <- sapply(X, function(xi) {
        mvtnorm::dmvt(xi, delta = mu[g], sigma = as.matrix(sigma[g]), df = df[g], log = FALSE)
      })
      w_tilde[, g] <- (df[g] + 1) / (df[g] + (X - mu[g])^2 / sigma[g])
    }

    z_tilde <- z * matrix(rep(py, each = n), nrow = n, ncol = G) / matrix(rep(z %*% py, G), nrow = n, ncol = G)
    z_tilde[is.infinite(z_tilde) | is.nan(z_tilde)] <- 1/G
    N <- colSums(z_tilde)

    #++++ M-step: pi and mu ++++#
    py <- N / n
    mu <- colSums(z_tilde * w_tilde * X) / colSums(z_tilde * w_tilde)

    for (g in 1:G) {

      #++++ M-step: sigma ++++#

      sigma_den <- sum(z_tilde[, g] * w_tilde[, g] * (X - mu[g])^2)
      sigma[g] <- sigma_den / N[g]

      #++++ M-step 2: degree of freedom ++++#

      root <- uniroot.all(function(a) {
        A <- -digamma(a / 2) + log(a / 2) + 1
        B <-  sum(z_tilde[, g] * (log(w_tilde[, g]) - w_tilde[, g])) / N[g]
        C <-  digamma((df[g] + 1) / 2) - log((df[g] + 1) / 2)

        return(A + B + C)
      }, lower = 2, upper = 200)

      df[g] <- ifelse(length(root) != 0, root, df[g])
    }

    #++++ Obtain observed likelihood and log-likelihood ++++#

    dens <- matrix(NA, nrow = n, ncol = G)
    for (g in 1:G) {
      dens[, g] <- sapply(X, function(xi) {
        mvtnorm::dmvt(xi, delta = mu[g], sigma = as.matrix(sigma[g]), df = df[g], log = FALSE)
      })
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

  delta <- (X - mu[g])^2
  outliers <- 1 - pchisq(delta, df = 1) < 1 - outlier_cutoff

  #------------------------#
  #  Number of parameters  #
  #------------------------#

  npar <- list(
    pi    = G - 1,
    mu    = G,
    sigma = G,
    df    = G
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
    model        = 'tM',
    pi           = py,
    mu           = mu,
    sigma        = sigma,
    df           = df,
    z_tilde      = z_tilde,
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
