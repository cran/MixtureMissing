###########################################################################
###                                                                     ###
###                       Mixture Model Selection                       ###
###                                                                     ###
###########################################################################

#' Mixture Model Selection
#'
#' Fit mixtures via various distributions and decide the best model based on a given
#' information criterion. The distributions include multivariate contaminated normal,
#' multivariate generalized hyperbolic, special and limiting cases of multivariate
#' generalized hyperbolic.
#'
#' @param X An \eqn{n} x \eqn{d} matrix or data frame where \eqn{n} is the number of
#'   observations and \eqn{d} is the number of variables.
#' @param G The number of clusters, which must be at least 1. If \code{G = 1}, then
#'   both \code{init_method} and \code{clusters} are ignored.
#' @param model A vector of character strings indicating the mixture model(s) to be fitted.
#'   See the details section for a list of available distributions. However, all distributions
#'   will be considered by default.
#' @param criterion A character string indicating the information criterion for model
#'   selection. "BIC" is used by default. See the details section for a list of available
#'   information criteria.
#' @param max_iter (optional) A numeric value giving the maximum number of
#'   iterations each EM algorithm is allowed to use; 20 by default.
#' @param epsilon (optional) A number specifying the epsilon value for the
#'   Aitken-based stopping criterion used in the EM algorithm: 0.01 by default.
#' @param init_method (optional) A string specifying the method to initialize
#'   the EM algorithm. "kmedoids" clustering is used by default. Alternative
#'   methods include "kmeans", "hierarchical", and "manual". When "manual" is chosen,
#'   a vector \code{clusters} of length \eqn{n} must be specified. If the data set is
#'   incomplete, missing values will be first filled based on the mean imputation method.
#' @param clusters (optional) A vector of length \eqn{n} that specifies the initial
#'   cluster memberships of the user when \code{init_method} is set to "manual".
#'   Both numeric and character vectors are acceptable. This argument is NULL by
#'   default, so that it is ignored whenever other given initialization methods
#'   are chosen.
#' @param eta_min (optional) A numeric value close to 1 to the right specifying
#'   the minimum value of eta; 1.001 by default. This is only relevant for CN mixture
#' @param outlier_cutoff (optional) A number between 0 and 1 indicating the
#'   percentile cutoff used for outlier detection. This is only relevant for t mixture.
#' @param deriv_ctrl (optional) A list containing arguments to control the numerical
#'   procedures for calculating the first and second derivatives. Some values are
#'   suggested by default. Refer to functions \code{grad} and \code{hessian} under
#'   the package \code{numDeriv} for more information.
#' @param progress (optional) A logical value indicating whether the
#'   fitting progress should be displayed; TRUE by default.
#'
#' @details The function can fit mixtures via the contaminated normal distribution,
#' generalized hyperbolic distribution, and special and limiting cases of the generalized
#' hyperbolic distribution. Available distributions include
#' \itemize{
#'   \item CN - Contaminated Normal
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
#' @return A list with
#'   \item{best_mod}{An object of class \code{MixtureMissing} corresponding to the best model.}
#'   \item{all_mod}{A list of objects of class \code{MixtureMissing} corresponding to all models of consideration.
#'     The list is in the order of \code{model}.}
#'   \item{criterion}{A numeric vector containing the chosen information criterion values of all models of consideration.
#'     The vector is in the order of best-to-worst models.}
#' Each object of class \code{MixtureMissing} have slots depending on the fitted model. See
#' the returned value of \link{MCNM} and \link{MGHM}.
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
#' mod <- select_mixture(X, G = 2, model = c('CN', 'GH', 'St'), criterion = 'BIC', max_iter = 10)
#'
#' #++++ With missing values ++++#
#'
#' set.seed(1234)
#'
#' X <- hide_values(bankruptcy[, 2:3], prop_cases = 0.1)
#' mod <- select_mixture(X, G = 2, model = c('CN', 'GH', 'St'), criterion = 'BIC', max_iter = 10)
#'
#' @import numDeriv Bessel
#' @importFrom stats complete.cases cov cutree dist dnorm hclust kmeans
#'   mahalanobis pchisq rmultinom runif var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
select_mixture <- function(
    X,
    G,
    model          = c('CN', 'GH', 'NIG', 'SNIG', 'SC', 'C', 'St', 't', 'N', 'SGH', 'HUM', 'H', 'SH'),
    criterion      = c('BIC', 'AIC', 'KIC', 'KICc', 'AIC3', 'CAIC', 'AICc', 'ICL', 'AWE', 'CLC'),
    max_iter       = 20,
    epsilon        = 0.01,
    init_method    = c('kmedoids', 'kmeans', 'hierarchical', 'manual'),
    clusters       = NULL,
    eta_min        = 1.001,
    outlier_cutoff = 0.95,
    deriv_ctrl     = list(eps = 1e-8, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-7),
                          r = 6, v = 2, show.details = FALSE),
    progress       = TRUE
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

  #--------------------------------------------#
  #    Model list and information criterion    #
  #--------------------------------------------#

  criterion <- match.arg(criterion)
  best_info <- Inf
  best_mod  <- NULL

  m_codes <- c('CN', 'GH',
               'NIG', 'SNIG',
               'SC', 'C',
               'St', 't', 'N',
               'SGH', 'HUM',
               'H', 'SH')

  m_names <- c('Contaminated Normal', 'Generalized Hyperbolic',
               'Normal-Inverse Gaussian', 'Symmetric Normal-Inverse Gaussian',
               'Skew-Cauchy', 'Cauchy',
               'Skew-t', 't', 'Normal',
               'Symmetric Generalized Hyperbolic', 'Hyperbolic Univariate Marginals',
               'Hyperbolic', 'Symmetric Hyperbolic')

  model <- unique(model)

  if ( !all(model %in% m_codes) ) {
    stop('Unavailable model(s) input. Please check again.')
  }

  m_names <- m_names[match(model, m_codes)]
  m_codes <- model

  n_models     <- length(model)
  infos        <- rep(NA_real_, n_models)
  names(infos) <- model

  #--------------------------------#
  #    Parameter Initialization    #
  #--------------------------------#

  init_method <- match.arg(init_method)

  if (progress) {
    cat('\nData Set:', ifelse(any(is.na(X)), 'Incomplete', 'Complete'), '\n')
    cat('\nInformation Criterion:', criterion, '\n')
    cat('\nInitialization:', init_method, '\n')
  }

  X_imp <- X

  if (G == 1) {

    max_iter <- 1

  } else {

    X_imp <- mean_impute(X)
    init <- initialize_clusters(
      X           = X_imp,
      G           = G,
      init_method = init_method,
      clusters    = clusters
    )

    clusters <- init$clusters

  }

  #---------------------#
  #    Model Fitting    #
  #---------------------#

  if (progress) {
    cat('\nModel Fitting:\n')
  }

  allmod <- vector('list', n_models)

  for (j in 1:n_models) {

    #++++ Fit each model ++++#

    if (model[j] == 'CN') {

      mod <- tryCatch({
        MCNM(
          X           = X,
          G           = G,
          max_iter    = max_iter,
          epsilon     = epsilon,
          init_method = 'manual',
          clusters    = clusters,
          eta_min     = eta_min,
          progress    = FALSE
        )
      }, error = function(err) { return(NULL) })

    } else {

      mod <- tryCatch({
        MGHM(
          X              = X,
          G              = G,
          model          = model[j],
          max_iter       = max_iter,
          epsilon        = epsilon,
          init_method    = 'manual',
          clusters       = clusters,
          outlier_cutoff = outlier_cutoff,
          deriv_ctrl     = deriv_ctrl,
          progress       = FALSE
        )
      }, error = function(err) { return(NULL) })

    }

    allmod[[j]] <- mod

    #++++ Compare to the current best model ++++#

    if (!is.null(mod)) {
      infos[j] <- mod[[criterion]]

      if (best_info > infos[j]) {
        best_info <- infos[j]
        best_mod  <- mod
      }
    }

    # if (progress) {
    #   if (is.null(mod)) {
    #     cat('  ', j, '. ', m_names[match(model[j], m_codes)], ': Failed\n', sep = '')
    #   } else {
    #     cat('  ', j, '. ', m_names[match(model[j], m_codes)], ': ', infos[j], '\n', sep = '')
    #   }
    # }

    #++++ Update progress ++++#

    if (progress) {
      cat('  ')
      if (is.null(mod)) {
        cat(model[j], ': Failed  ', sep = '')
        #cat( m_names[match(model[j], m_codes)], ': Failed', sep = ' ')
      } else {
        # cat( m_names[match(model[j], m_codes)],  sep = ' ')
        cat(model[j], ' ', sep = ' ')
      }
    }

  }    # end for (j in 1:n_models)

  #--------------------------------------------#
  #    Summarize Results and Prepare Output    #
  #--------------------------------------------#

  # if (progress) {
  #   if (sum(is.na(infos)) == n_models) {
  #     cat('\nThe best mixture model cannot be determined\n')
  #   } else {
  #     cat('\nAccording to ', criterion, ', the best mixture model is based on the ',
  #         m_names[which.min(infos)], ' distribution\n', sep = '')
  #   }
  #   cat('\n')
  # }

  if (progress) {
    if (sum(is.na(infos)) == n_models) {
      cat('\n\nThe best mixture model cannot be determined\n')
    } else {
      cat('\n\nAccording to ', criterion, ', the best mixture model is based on the ',
          m_names[which.min(infos)], ' distribution\n', sep = '')
    }

    cat('\nModel rank according to ', criterion, ':', sep = '')
    infos <- sort(infos, na.last = TRUE)
    for(j in 1:length(infos)){
      cat('\n')
      if (!is.na(infos[j])) {
        cat('  ', j, '. ', m_names[match(names(infos[j]), m_codes)], ': ', round(infos[j], digits = 4), sep = '')
      } else {
        cat('  ', j, '. ', m_names[match(names(infos[j]), m_codes)], ': Failed', sep = '')
      }
    }
    cat('\n\n')
  }

  # output <- list(
  #   best_model = best_mod,
  #   infos      = infos
  # )
  #
  # return(output)

  names(allmod) <- model

  output <- list(
    best_mod  = best_mod,
    all_mod   = allmod,
    criterion = infos
  )

  return(output)

}
