###########################################################################
###                                                                     ###
###             Multivariate Generalized Hyperbolic Mixture             ###
###                                                                     ###
###########################################################################

#' Multivariate Generalized Hyperbolic Mixture (MGHM)
#'
#' Carries out model-based clustering using a multivariate generalized hyperbolic
#' mixture (MGHM). The function will determine itself if the data set is
#' complete or incomplete and fit the appropriate model accordingly. In the incomplete
#' case, the data set must be at least bivariate, and missing values are assumed to
#' be missing at random (MAR).
#'
#' @param X An \eqn{n} x \eqn{d} matrix or data frame where \eqn{n} is the number of
#'   observations and \eqn{d} is the number of variables.
#' @param G An integer vector specifying the numbers of clusters, which must be at least 1.
#' @param model A string indicating the mixture model to be fitted; "GH" for generalized
#' hyperbolic by default. See the details section for a list of available distributions.
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
#' @param clusters (optional) A vector of length \eqn{n} that specifies the initial
#'   cluster memberships of the user when \code{init_method} is set to "manual".
#'   Both numeric and character vectors are acceptable. This argument is NULL by
#'   default, so that it is ignored whenever other given initialization methods
#'   are chosen.
#' @param outlier_cutoff (optional) A number between 0 and 1 indicating the
#'   percentile cutoff used for outlier detection. This is only relevant for t mixture.
#' @param deriv_ctrl (optional) A list containing arguments to control the numerical
#'   procedures for calculating the first and second derivatives. Some values are
#'   suggested by default. Refer to functions \code{grad} and \code{hessian} under
#'   the package \code{numDeriv} for more information.
#' @param progress (optional) A logical value indicating whether the
#'   fitting progress should be displayed; TRUE by default.
#'
#' @details Beside the generalized hyperbolic distribution, the function can fit mixture
#'   via its special and limiting cases. Available distributions include
#' \itemize{
#'   \item GH - Generalized Hyperbolic
#'   \item NIG - Normal-Inverse Gaussian
#'   \item SNIG - Symmetric Normal-Inverse Gaussian
#'   \item SC - Skew-Cauchy
#'   \item C - Cauchy
#'   \item St - Skew-\emph{t}
#'   \item t - Student's \emph{t}
#'   \item N - Normal or Gaussian
#'   \item SGH - Symmetric Generalized Hyperbolic
#'   \item HUM- Hyperbolic Univariate Marginals
#'   \item H - Hyperbolic
#'   \item SH - Symmetric Hyperbolic
#' }
#' Available information criteria include
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
#'
#' @return An object of class \code{MixtureMissing} with:
#'   \item{model}{The model used to fit the data set.}
#'   \item{pi}{Mixing proportions.}
#'   \item{mu}{Component location vectors.}
#'   \item{Sigma}{Component dispersion matrices.}
#'   \item{beta}{Component skewness vectors. Only available if \code{model} is GH, NIG, SNIG, SC, SGH, HUM, H, or SH; NULL otherwise.}
#'   \item{lambda}{Component index parameters. Only available if \code{model} is GH, NIG, SNIG, SGH, HUM, H, or SH; NULL otherwise.}
#'   \item{omega}{Component concentration parameters. Only available if \code{model} is GH, NIG, SNIG, SGH, HUM, H, or SH; NULL otherwise.}
#'   \item{df}{Component degrees of freedom. Only available if \code{model} is St or t; NULL otherwise.}
#'   \item{z_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
#'     probabilities that the corresponding observation belongs to each cluster.}
#'   \item{clusters}{A numeric vector of length \eqn{n} indicating cluster
#'     memberships determined by the model.}
#'   \item{outliers}{A logical vector of length \eqn{n} indicating observations that are outliers. Only available if \code{model} is t; NULL otherwise.}
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
#' Browne, R. P. and McNicholas, P. D. (2015). A mixture of generalized hyperbolic distributions.
#'   \emph{Canadian Journal of Statistics}, 43(2):176–198. \cr \cr
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
#' mod <- MGHM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#'
#' summary(mod)
#' plot(mod)
#'
#' #++++ With missing values ++++#
#'
#' set.seed(1234)
#'
#' X <- hide_values(bankruptcy[, 2:3], prop_cases = 0.1)
#' mod <- MGHM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#'
#' summary(mod)
#' plot(mod)
#'
#' @import numDeriv Bessel mclust
#' @importFrom stats complete.cases cov cutree dist dnorm hclust kmeans
#'   mahalanobis pchisq rmultinom runif var uniroot
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
MGHM <- function(
    X,
    G,
    model          = c('GH', 'NIG', 'SNIG', 'SC', 'C', 'St', 't', 'N', 'SGH', 'HUM', 'H', 'SH'),
    criterion      = c('BIC', 'AIC', 'KIC', 'KICc', 'AIC3', 'CAIC', 'AICc', 'ICL', 'AWE', 'CLC'),
    max_iter       = 20,
    epsilon        = 0.01,
    init_method    = c("kmedoids", "kmeans", "hierarchical", "mclust", "manual"),
    clusters       = NULL,
    outlier_cutoff = 0.95,
    deriv_ctrl     = list(eps = 1e-8, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-7),
                          r = 6, v = 2, show.details = FALSE),
    progress       = TRUE
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

  model <- match.arg(model)

  m_codes <- c('GH',
               'NIG', 'SNIG',
               'SC', 'C',
               'St', 't', 'N',
               'SGH', 'HUM',
               'H', 'SH')

  m_names <- c('Generalized Hyperbolic',
               'Normal-Inverse Gaussian', 'Symmetric Normal-Inverse Gaussian',
               'Skew-Cauchy', 'Cauchy',
               'Skew-t', 't', 'Gaussian',
               'Symmetric Generalized Hyperbolic', 'Hyperbolic Univariate Marginals',
               'Hyperbolic', 'Symmetric Hyperbolic')

  if ( any(is.na(X)) ) {
    if (ncol(X) < 2) {
      stop('If X contains NAs, X must be at least bivariate')
    }

    if (progress) {
      cat('\nMixture:', m_names[match(model, m_codes)], paste('(', model, ')', sep = ''), '\n')
      cat('Data Set: Incomplete\n')
    }
  } else {
    if (progress) {
      cat('\nMixture:', m_names[match(model, m_codes)], paste('(', model, ')', sep = ''), '\n')
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

    if ( any(is.na(X)) ) {

      ### Generalized Hyperbolic
      ### Normal-Inverse Gaussian
      ### Skew Normal-Inverse Gaussian
      ### Symmetric Generalized Hyperbolic
      ### Hyperbolic Univarate Marginals
      ### Hyperbolic
      ### Symmetric Hyperbolic

      if (model %in% c('GH', 'NIG', 'SNIG', 'SGH', 'HUM', 'H', 'SH')) {
        mod <- tryCatch({
          MGHM_incomplete_data(
            X           = X,
            G           = G,
            model       = model,
            max_iter    = max_iter,
            epsilon     = epsilon,
            init_method = init_method,
            clusters    = clusters,
            deriv_ctrl  = deriv_ctrl,
            progress    = progress
          )
        }, error = function(err) { return(NULL) })
      }

      ### Skew-t
      ### Skew-Cauchy

      if (model %in% c('St', 'SC')) {
        mod <- tryCatch({
          MStM_incomplete_data(
            X            = X,
            G            = G,
            model        = model,
            max_iter     = max_iter,
            epsilon      = epsilon,
            init_method  = init_method,
            clusters     = clusters,
            deriv_ctrl   = deriv_ctrl,
            progress     = progress
          )
        }, error = function(err) { return(NULL) })
      }

      ### t
      ### Cauchy

      if (model %in% c('t', 'C')) {
        mod <- tryCatch({
          MtM_incomplete_data(
            X              = X,
            G              = G,
            model          = model,
            max_iter       = max_iter,
            epsilon        = epsilon,
            init_method    = init_method,
            clusters       = clusters,
            outlier_cutoff = outlier_cutoff,
            progress       = progress
          )
        }, error = function(err) { return(NULL) })
      }

      ### Normal

      if (model == 'N') {
        mod <- tryCatch({
          MNM_incomplete_data(
            X              = X,
            G              = G,
            max_iter       = max_iter,
            epsilon        = epsilon,
            init_method    = init_method,
            clusters       = clusters,
            progress       = progress
          )
        }, error = function(err) { return(NULL) })
      }

    } else {

      ### Generalized Hyperbolic
      ### Normal-Inverse Gaussian
      ### Skew Normal-Inverse Gaussian
      ### Symmetric Generalized Hyperbolic
      ### Hyperbolic Univarate Marginals
      ### Hyperbolic
      ### Symmetric Hyperbolic

      if (model %in% c('GH', 'NIG', 'SNIG', 'SGH', 'HUM', 'H', 'SH')) {
        mod <- tryCatch({
          MGHM_complete_data(
            X           = X,
            G           = G,
            model       = model,
            max_iter    = max_iter,
            epsilon     = epsilon,
            init_method = init_method,
            clusters    = clusters,
            deriv_ctrl  = deriv_ctrl,
            progress    = progress
          )
        }, error = function(err) { return(NULL) })
      }

      ### Skew-t
      ### Skew-Cauchy

      if (model %in% c('St', 'SC')) {
        mod <- tryCatch({
          MStM_complete_data(
            X            = X,
            G            = G,
            model        = model,
            max_iter     = max_iter,
            epsilon      = epsilon,
            init_method  = init_method,
            clusters     = clusters,
            deriv_ctrl   = deriv_ctrl,
            progress     = progress
          )
        }, error = function(err) { return(NULL) })
      }

      ### t
      ### Cauchy

      if (model %in% c('t', 'C')) {
        mod <- tryCatch({
          MtM_complete_data(
            X              = X,
            G              = G,
            model          = model,
            max_iter       = max_iter,
            epsilon        = epsilon,
            init_method    = init_method,
            clusters       = clusters,
            outlier_cutoff = outlier_cutoff,
            progress       = progress
          )
        }, error = function(err) { return(NULL) })
      }

      ### Normal

      if (model == 'N') {
        mod <- tryCatch({
          MNM_complete_data(
            X              = X,
            G              = G,
            max_iter       = max_iter,
            epsilon        = epsilon,
            init_method    = init_method,
            clusters       = clusters,
            progress       = progress
          )
        }, error = function(err) { return(NULL) })
      }

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

  return(best_mod)

}

###########################################################################
###                                                                     ###
###   Multivariate Generalized Hyperbolic Mixture for Incomplete Data   ###
###                                                                     ###
###########################################################################

MGHM_incomplete_data <- function(
    X,
    G,
    model        = c('GH', 'NIG', 'SNIG', 'SGH', 'HUM', 'H', 'SH'),
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

  py     <- rep(NA_real_, G)
  mu     <- matrix(NA_real_, nrow = G, ncol = d)
  Sigma  <- array(NA_real_, dim = c(d, d, G))
  beta   <- matrix(0, nrow = G, ncol = d)
  omega  <- rep(1, G)

  if (model %in% c('GH', 'NIG', 'SNIG', 'SGH')) {
    lambda <- rep(-1/2, G)
  }

  if (model == 'HUM') {
    lambda <- rep(1, G)
  }

  if (model %in% c('H', 'SH')) {
    lambda <- rep(-(d + 1)/2, G)
  }

  log_dens   <- matrix(NA_real_, nrow = n, ncol = G)
  iter   <- 0
  loglik <- NULL

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

        z[Im[[j]], g] <- dGH(Xo_j, mu = mu_o, Sigma = Sigma_oo, beta = beta_o,
                             lambda = lambda[g], omega = omega[g])

        psi <- c(omega[g] + t(beta_o) %*% Sigma_oo_inv %*% beta_o)
        chi <- omega[g] + mahalanobis(Xo_j, center = mu_o, cov = Sigma_oo)

        s1 <- sqrt(psi * chi)
        s2 <- sqrt(chi / psi)

        bessel_num <- besselK(s1, nu = lambda[g] - sum(o)/2 + 1, expon.scaled = TRUE)
        bessel_den <- besselK(s1, nu = lambda[g] - sum(o)/2, expon.scaled = TRUE)
        bessel_den[bessel_den < 10^-323] <- 10^-323

        a[Im[[j]], g] <- s2 * (bessel_num / bessel_den)
        b[Im[[j]], g] <- -(2 * lambda[g] - sum(o)) / chi + (bessel_num / bessel_den) / s2
        c[Im[[j]], g] <- log(s2) + numDeriv::grad(log_besselK, x = rep(lambda[g] - sum(o)/2, nrow(Xo_j)), y = s1,
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

      if (model %in% c('GH', 'NIG', 'H')) {
        num_mu <- colSums( z_tilde[, g] * ( (!R) * (a_bar[g] * b[, g] - 1) * X_tilde[, , g] + R * (a_bar[g] * X_tilde[, , g] - X_hat[, , g]) ) )
        den_mu <- sum( z_tilde[, g] * (a_bar[g] * b[, g] - 1) )

        mu[g, ] <- num_mu / den_mu

        num_beta <- colSums( z_tilde[, g] * ( (!R) * (b_bar[g] - b[, g]) * X_tilde[, , g] + R * (b_bar[g] * X_hat[, , g] - X_tilde[, , g]) ) )
        den_beta <- den_mu

        beta[g, ] <- num_beta / den_beta
      }

      if (model %in% c('SNIG', 'SGH', 'HUM', 'SH')) {
        num_mu <- colSums(z_tilde[, g] * b[, g] * ((!R) * X_tilde[, , g] + R * (X_tilde[, , g] - X_hat[, , g])/(b[, g] - 1)))
        den_mu <- sum(z_tilde[, g] * b[, g])

        mu[g, ] <- num_mu / den_mu

        beta[g, ] <- 0
      }

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

      #++++ M-step: lambda (index) and omega (concentration) ++++#

      if (model %in% c('GH', 'SGH')) {
        if (c_bar[g] == 0) {
          lambda[g] <- 0
        } else {
          grad_lambda <- numDeriv::grad(log_besselK, x = lambda[g], y = omega[g],
                                        method = 'Richardson', method.args = deriv_ctrl)

          lambda[g] <- c_bar[g] * lambda[g] / grad_lambda
        }
      }

      if (model %in% c('NIG', 'SNIG')) {
        lambda[g] <- -1/2
      }

      if (model == 'HUM') {
        lambda[g] <- 1
      }

      if (model %in% c('H', 'SH')) {
        lambda[g] <- -(d + 1)/2
      }

      grad_omega <- numDeriv::grad(q_func, x = omega[g], lambda = lambda[g],
                                   a_bar = a_bar[g], b_bar = b_bar[g], c_bar = c_bar[g],
                                   method = 'Richardson', method.args = deriv_ctrl)

      hess_omega <- numDeriv::hessian(q_func, x = omega[g], lambda = lambda[g],
                                      a_bar = a_bar[g], b_bar = b_bar[g], c_bar = c_bar[g],
                                      method = 'Richardson', method.args = deriv_ctrl)

      if (omega[g] - grad_omega / hess_omega > 0) {
        omega[g] <- omega[g] - grad_omega / hess_omega
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

        log_dens[Im[[j]], g] <- dGH(Xo_j, mu = mu_o, Sigma = Sigma_oo, beta = beta_o,
                                lambda = lambda[g], omega = omega[g], log = TRUE)

      }
    }

    log_py_dens  <- sweep(log_dens, 2, log(py), FUN = '+')
    final_loglik <- sum( apply(log_py_dens, 1, log_sum_exp) )
    loglik       <- c(loglik, final_loglik)

    #++++ Update progress ++++#

    iter <- iter + 1

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
    beta   = G * d,
    lambda = G,
    omega  = G
  )

  if (model %in% c('NIG', 'H')) {
    npar$lambda <- 0
  }

  if (model %in% c('SNIG', 'HUM', 'SH')) {
    npar$beta   <- 0
    npar$lambda <- 0
  }

  if (model == 'SGH') {
    npar$beta <- 0
  }

  npar$total <- Reduce('+', npar)

  #----------------------------#
  #    Information Criteria    #
  #----------------------------#

  AIC <- -2 * final_loglik + 2 * npar$total
  BIC <- -2 * final_loglik + npar$total * log(n)

  KIC  <- -2 * final_loglik + 3 * (npar$total + 1)
  KICc <- -2 * final_loglik + 2 * (npar$total + 1) * n/(n-npar$total - 2) - n * digamma((n - npar$total)/2) + n * log(n/2)

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
  names(lambda)   <- c_names
  names(omega)    <- c_names

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
    beta  <- beta[1, ]
  }

  output <- list(
    model        = paste(model, '_incomplete_data', sep = ''),
    pi           = py,
    mu           = mu,
    Sigma        = Sigma,
    beta         = beta,
    lambda       = lambda,
    omega        = omega,
    z_tilde      = z_tilde,
    clusters     = clusters,
    data         = X_imputed,
    complete     = !is.na(X),
    npar         = npar,
    max_iter     = max_iter,
    iter_stop    = iter,
    final_loglik = final_loglik,
    loglik       = loglik,
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
  class(output) <- 'MixtureMissing'

  return(output)

}

###########################################################################
###                                                                     ###
###    Multivariate Generalized Hyperbolic Mixture for Complete Data    ###
###                                                                     ###
###########################################################################

MGHM_complete_data <- function(
    X,
    G,
    model        = c('GH', 'NIG', 'SNIG', 'SGH', 'HUM', 'H', 'SH'),
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

  py     <- rep(NA_real_, G)
  mu     <- matrix(NA_real_, nrow = G, ncol = d)
  Sigma  <- array(NA_real_, dim = c(d, d, G))
  beta   <- matrix(0, nrow = G, ncol = d)
  omega  <- rep(1, G)

  if (model %in% c('GH', 'NIG', 'SNIG', 'SGH')) {
    lambda <- rep(-1/2, G)
  }

  if (model == 'HUM') {
    lambda <- rep(1, G)
  }

  if (model %in% c('H', 'SH')) {
    lambda <- rep(-(d + 1)/2, G)
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

  while (iter < max_iter & getall(loglik) > epsilon) {

    #++++ E-step ++++#

    for (g in 1:G) {
      z[, g] <- dGH(X, mu = mu[g, ], Sigma = Sigma[, , g], beta = beta[g, ], lambda = lambda[g], omega = omega[g])

      Sigma_inv <- solve(Sigma[, , g])

      psi <- c(omega[g] + t(beta[g, ]) %*% Sigma_inv %*% beta[g, ])
      chi <- omega[g] + mahalanobis(X, center = mu[g, ], cov = Sigma[, , g])

      s1 <- sqrt(psi * chi)
      s2 <- sqrt(chi / psi)

      bessel_num <- besselK(s1, nu = lambda[g] - d/2 + 1, expon.scaled = TRUE)
      bessel_den <- besselK(s1, nu = lambda[g] - d/2, expon.scaled = TRUE)
      bessel_den[bessel_den < 10^-323] <- 10^-323

      a[, g] <- s2 * (bessel_num / bessel_den)
      b[, g] <- -(2 * lambda[g] - d) / chi + (bessel_num / bessel_den) / s2
      c[, g] <- log(s2) + numDeriv::grad(log_besselK, x = rep(lambda[g] - d/2, n), y = s1,
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

      if (model %in% c('GH', 'NIG', 'H')) {
        num_mu <- colSums( z_tilde[, g] * X * (a_bar[g] * b[, g] - 1) )
        den_mu <- sum( z_tilde[, g] * (a_bar[g] * b[, g] - 1) )

        mu[g, ] <- num_mu / den_mu

        num_beta <- colSums( z_tilde[, g] * X * (b_bar[g] - b[, g]) )
        den_beta <- den_mu

        beta[g, ] <- num_beta / den_beta
      }

      if (model %in% c('SNIG', 'SGH', 'HUM', 'SH')) {
        num_mu <- colSums( z_tilde[, g] * b[, g] * X )
        den_mu <- sum( z_tilde[, g] * b[, g] )

        mu[g, ] <- num_mu / den_mu

        beta[g, ] <- 0
      }

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

      #++++ M-step: lambda (index) and omega (concentration) ++++#

      if (model %in% c('GH', 'SGH')) {
        if (c_bar[g] == 0) {
          lambda[g] <- 0
        } else {
          grad_lambda <- numDeriv::grad(log_besselK, x = lambda[g], y = omega[g],
                                        method = 'Richardson', method.args = deriv_ctrl)

          lambda[g] <- c_bar[g] * lambda[g] / grad_lambda
        }
      }

      if (model %in% c('NIG', 'SNIG')) {
        lambda[g] <- -1/2
      }

      if (model == 'HUM') {
        lambda[g] <- 1
      }

      if (model %in% c('H', 'SH')) {
        lambda[g] <- -(d + 1)/2
      }

      grad_omega <- numDeriv::grad(q_func, x = omega[g], lambda = lambda[g],
                                   a_bar = a_bar[g], b_bar = b_bar[g], c_bar = c_bar[g],
                                   method = 'Richardson', method.args = deriv_ctrl)

      hess_omega <- numDeriv::hessian(q_func, x = omega[g], lambda = lambda[g],
                                      a_bar = a_bar[g], b_bar = b_bar[g], c_bar = c_bar[g],
                                      method = 'Richardson', method.args = deriv_ctrl)

      if (omega[g] - grad_omega / hess_omega > 0) {
        omega[g] <- omega[g] - grad_omega / hess_omega
      }

    }

    #++++ Observed Log-likelihood ++++#

    for (g in 1:G) {
      log_dens[, g] <- dGH(X, mu = mu[g, ], Sigma = Sigma[, , g], beta = beta[g, ], lambda = lambda[g], omega = omega[g], log = TRUE)
    }

    log_py_dens  <- sweep(log_dens, 2, log(py), FUN = '+')
    final_loglik <- sum( apply(log_py_dens, 1, log_sum_exp) )
    loglik       <- c(loglik, final_loglik)

    #++++ Update progress ++++#

    iter <- iter + 1

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
    lambda = G,
    omega  = G
  )

  if (model %in% c('NIG', 'H')) {
    npar$lambda <- 0
  }

  if (model %in% c('SNIG', 'HUM', 'SH')) {
    npar$beta   <- 0
    npar$lambda <- 0
  }

  if (model == 'SGH') {
    npar$beta <- 0
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
  names(lambda)   <- c_names
  names(omega)    <- c_names

  if (G == 1) {
    mu    <- mu[1, ]
    Sigma <- Sigma[, , 1]
    beta  <- beta[1, ]
  }

  output <- list(
    model        = paste(model, '_complete_data', sep = ''),
    pi           = py,
    mu           = mu,
    Sigma        = Sigma,
    beta         = beta,
    lambda       = lambda,
    omega        = omega,
    z_tilde      = z_tilde,
    clusters     = clusters,
    data         = X,
    complete     = !is.na(X),
    npar         = npar,
    max_iter     = max_iter,
    iter_stop    = iter,
    final_loglik = final_loglik,
    loglik       = loglik,
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
  class(output) <- 'MixtureMissing'

  return(output)

}

############################################################################
###                                                                      ###
###       Density Function for Generalized Hyperbolic Distribution       ###
###                                                                      ###
############################################################################

dGH <- function(
    X,
    mu     = rep(0, d),    # location
    Sigma  = diag(d),      # dispersion
    beta   = rep(0, d),    # skewness
    lambda = 0.5,          # index
    omega  = 1,            # concentration
    log    = FALSE
) {

  #----------------------#
  #    Input Checking    #
  #----------------------#

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

  if (omega <= 0) {
    stop('omega must be positive')
  }

  if (is.nan(log(det(Sigma)))) {
    Sigma <- diag(d)
  }

  X_centrd  <- sweep(X, 2, mu, '-')
  Sigma_inv <- solve(Sigma)

  psi <- c(omega + t(beta) %*% Sigma_inv %*% beta)
  chi <- omega + mahalanobis(X, center = mu, cov = Sigma_inv, inverted = TRUE)

  s1 <- sqrt(psi * chi)

  res <- (lambda - d/2) * log(s1) + (d/2 - lambda) * log(psi)
  res <- res + log(besselK(s1, nu = lambda - d/2, expon.scaled = TRUE)) - s1

  res <- res - d/2 * (log(2) + log(pi)) - 1/2 * log(det(Sigma))
  res <- res + omega - log(besselK(omega, nu = lambda, expon.scaled = TRUE))
  res <- res + X_centrd %*% Sigma_inv %*% beta

  if (!log) {
    res <- c(exp(res))
    res[res <= 10^(-323)] <- 10^(-323)
  }

  return(res)

  # lvx <- (lambda - d/2) * log(s1) + X_centrd %*% Sigma_inv %*% beta
  # lvx <- lvx + log(besselK(s1, nu = lambda - d/2, expon.scaled = TRUE)) - s1
  #
  # if (is.nan(log(det(Sigma)))) {
  #   Sigma <- diag(d)
  # }
  #
  # lv <- -1/2 * log(det(Sigma)) - d/2 * (log(2) + log(pi))
  # lv <- lv + omega - log(besselK(omega, nu = lambda, expon.scaled = TRUE))
  # lv <- lv + (d/2 - lambda) * log(psi)
  #
  # return(
  #   c( exp(lvx + lv) )
  # )

}

############################################################################
###                                                                      ###
###                       Log of Bessel K Function                       ###
###                                                                      ###
############################################################################

log_besselK <- function(x, y) {

  val <- log(besselK(y, nu = x, expon.scaled = FALSE))
  sun <- is.infinite(val)

  val[sun] <- Bessel::besselK.nuAsym(x = y[sun], nu = abs(x[sun]), k.max = 4, log = TRUE, expon.scaled = FALSE)

  val

}

###########################################################################
###                                                                     ###
###           Objective Function to Optimize lambda and omega           ###
###                                                                     ###
###########################################################################

q_func <- function(x, lambda, a_bar, b_bar, c_bar) {
  -log_besselK(x = lambda, y = x) + (lambda - 1) * c_bar - x * (a_bar + b_bar) / 2
}

# q_func <- function(x, lambda, a_bar, b_bar, c_bar) {
#   -log(besselK(x, nu = lambda)) + (lambda - 1) * c_bar - x * (a_bar + b_bar) / 2
# }

Rlambda <- function (x, lambda = NULL)  {

  v1  <- besselK(x, nu = lambda + 1, expon.scaled = FALSE)
  v0  <- besselK(x, nu = lambda, expon.scaled = FALSE)
  val <- v1/v0

  set <- is.infinite(v1) | is.infinite(v0) | v0 == 0 | v1 == 0

  if (any(set)) {

    lv1 <- Bessel::besselK.nuAsym(x = x[set], nu = abs(lambda[set] + 1), k.max = 4, expon.scaled = FALSE, log = TRUE)
    lv0 <- Bessel::besselK.nuAsym(x = x[set], nu = abs(lambda[set]), k.max = 4, expon.scaled = FALSE, log = TRUE)

    val[set] <- exp(lv1 - lv0)

  }

  return(val)
}

