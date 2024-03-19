#########################################################
###                                                   ###
###       Extractor Function for MixtureMissing       ###
###                                                   ###
#########################################################

#' Extractor function for MixtureMissing
#'
#' Extract values from \code{MixtureMissing} objects or from outputs of
#' \link{select_mixture}.
#'
#' @param object A \code{MixtureMissing} object or an output of \link{select_mixture}.
#' @param what The specific value to be extracted. See the return section for possible
#'   values.
#' @param criterion If \code{what = "information"}, \code{criterion} is a vector of
#'   desired information criteria. All criteria will be extracted by default. Duplicate values
#'   in the vector will not be shown again. See the details section for a list of available
#'   information criteria.
#' @param m_code Only used in the case when \code{object} is an output of \link{select_mixture}.
#'   If \code{m_code = NULL}, extracting will be based on the best model. If \code{m_code} is one of 'CN', 'GH',
#'   'NIG', 'SNIG', 'SC', 'C', 'St', 't', 'N', 'SGH', 'HUM', 'H', and 'SH', the function will look
#'   for this specific model and extract accordingly.
#'
#' @return One of the following depending on \code{what}
#' \itemize{
#'   \item If \code{what = "model"} - A data frame showing the component distribution
#'     and its abbreviation, number of clusters, and whether the data set is complete
#'     or incomplete.
#'   \item If \code{what = "parameters"} - A list containing the relevant parameters.
#'   \item If \code{what = "cluster"} - A numeric vector of length \eqn{n} indicating cluster
#'     memberships determined by the model.
#'   \item If \code{what = "posterior"} - An \eqn{n} by \eqn{G} matrix where each
#'     row indicates the expected probabilities that the corresponding observation
#'     belongs to each cluster.
#'   \item If \code{what = "outlier"} - A logical vector of length \eqn{n} indicating observations that are outliers.
#'     Only available if \code{model} is CN or t; NULL otherwise with a warning.
#'   \item If \code{what = "missing"} - A data frame showing how many observations (cases)
#'     have missing values and the number of missing values per variables.
#'   \item If \code{what = "imputed"} - The original data set if it is complete; otherwise, this is
#'     the data set with missing values imputed by appropriate expectations.
#'   \item If \code{what = "complete"} - An \eqn{n} by \eqn{d} logical matrix indicating which cells have no missing values.
#'   \item If \code{what = "information"} - A data frame showing the number of clusters, final observed
#'     log-likelihood value, number of parameters, and desired information criteria.
#' }
#'
#' @details Available information criteria include
#' \itemize{
#'   \item AIC - Akaike information criterion
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
#' @examples
#'
#' #++++ With no missing values ++++#
#'
#' X <- iris[, 1:4]
#' mod <- MGHM(X, G = 2, model = 'GH', max_iter = 10)
#' extract(mod, what = "model")
#' extract(mod, what = "parameters")
#' extract(mod, what = "cluster")
#'
#' #++++ With missing values ++++#
#'
#' set.seed(123)
#' X <- hide_values(iris[, 1:4], n_cases = 20)
#' mod <- MGHM(X, G = 2, model = 'GH', max_iter = 10)
#' extract(mod, what = "outlier")
#' extract(mod, what = "missing")
#' extract(mod, what = "imputed")
#'
#' @export
extract <- function(
    object,
    what      = c("model", "parameters", "cluster", "posterior", "outlier", "missing", "imputed", "complete", "information"),
    criterion = c("AIC", "BIC", "KIC", "KICc", "AIC3", "CAIC", "AICc", "ICL", "AWE", "CLC"),
    m_code    = NULL
) {



  what <- match.arg(what, several.ok = FALSE)

  #++++ Model ++++#

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

  if(!inherits(object, "MixtureMissing")) {
    object2 <- try(object$best_mod, silent = TRUE)

    if(!inherits(object2, "MixtureMissing")){
      stop("object not of class 'MixtureMissing'")
    } else {

      if (is.null(m_code)) {
        object <- object$best_mod
      } else {
        models <- names(object$all_mod)
        if (m_code %in% models) {
          object <- object$all_mod[[match.arg(m_code, models)]]

          if (is.null(object)) {
            stop(paste0(m_code, ' was failed to fit so nothing can be extracted'))
          }
        } else {
          if (m_code %in% m_codes) {
            stop(paste0(m_code, ' was not fitted so nothing can be extracted'))
          } else {
            stop(paste0(m_code, ' is not an available distribution'))
          }
        }
      }

    }
  }

  incomplete_data <- grepl('incomplete', object$model)

  if (incomplete_data) {
    model_code <- gsub('_incomplete_data', '', object$model)
  } else {
    model_code <- gsub('_complete_data', '', object$model)
  }

  res <- data.frame(
    Distribution = m_names[match(model_code, m_codes)],
    Code         = model_code,
    G            = length(object$pi),
    Data         = ifelse(incomplete_data, 'Incomplete', 'Complete')
  )

  #++++ Parameters ++++#

  if (what == 'parameters') {

    ### Contaminated Normal

    if (model_code == 'CN') {
      res <- list(
        pi    = object$pi,
        mu    = object$mu,
        Sigma = object$Sigma,
        alpha = object$alpha,
        eta   = object$eta
      )
    }

    ### Generalized Hyperbolic
    ### Normal-Inverse Gaussian
    ### Skew Normal-Inverse Gaussian
    ### Symmetric Generalized Hyperbolic
    ### Hyperbolic Univarate Marginals
    ### Hyperbolic
    ### Symmetric Hyperbolic

    if (model_code %in% c('GH', 'NIG', 'SNIG', 'SGH', 'HUM', 'H', 'SH')) {
      res <- list(
        pi     = object$pi,
        mu     = object$mu,
        Sigma  = object$Sigma,
        beta   = object$beta,
        lambda = object$lambda,
        omega  = object$omega
      )
    }

    ### Skew-t
    ### Skew-Cauchy

    if (model_code %in% c('St', 'SC')) {
      res <- list(
        pi    = object$pi,
        mu    = object$mu,
        Sigma = object$Sigma,
        beta  = object$beta
      )

      if (model_code == 'St') {
        res$df <- object$df
      }
    }

    ### t
    ### Cauchy
    ### Normal

    if (model_code %in% c('t', 'C', 'N')) {
      res <- list(
        pi    = object$pi,
        mu    = object$mu,
        Sigma = object$Sigma
      )

      if (model_code == 't') {
        res$df <- object$df
      }
    }

  }    # end if (what == 'parameter')

  #++++ Cluster Membership ++++#

  if (what == 'cluster') {
    res <- object$clusters
  }

  #++++ Posterior Probabilities ++++#

  if (what == 'posterior') {
    res <- object$z_tilde
  }

  #++++ Outlier Detection ++++#

  if (what == 'outlier') {
    res <- NULL
    if (model_code %in% c('CN', 't')) {
      res <- object$outliers
    } else {
      warning("NULL returned. Outlier detection is only available for contaminated normal and t distributions.")
    }

  }

  #++++ Missing Values Information ++++#

  if (what == 'missing') {
    complete_obs   <- apply(object$complete, 1, all)
    incomplete_obs <- !complete_obs

    complete_matr <- object$complete
    incomplete_matr <- !complete_matr

    if (is.null(colnames(incomplete_matr))) {
      colnames(incomplete_matr) <- paste0('V', 1:ncol(incomplete_matr))
    }

    var_names <- colnames(incomplete_matr)

    res <- data.frame(
      Cases = sum(incomplete_obs)
    )

    for (j in 1:ncol(incomplete_matr)) {
      res[[var_names[j]]] <- sum(incomplete_matr[, j])
    }

  }

  #++++ Original Data with Values Imputed ++++#

  if (what == 'imputed') {
    res <- object$data
  }

  #++++ Complete Data Indicator Matrix ++++#

  if (what == 'complete') {
    res <- object$complete
  }

  #++++ Information Criteria ++++#

  if (what == 'information') {
    criterion <- unique(criterion)
    criterion <- match.arg(criterion, several.ok = TRUE)

    res <- data.frame(
      G          = length(object$pi),
      Loglik     = object$final_loglik,
      Parameters = object$npar$total
    )

    for (crit in criterion) {
      res[[crit]] <- object[[crit]]
    }

  }

  return(res)

}
