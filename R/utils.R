###########################################################################
###                                                                     ###
###                    Cluster Memberships to Matrix                    ###
###                                                                     ###
###########################################################################

clusters_to_matrix <- function(clusters, G = max(clusters)) {

  if (G %% 1 != 0) {
    stop('Number of clusters G must be an integer')
  }

  if (G < max(clusters)) {
    stop('G must be at least max(clusters)')
  }

  n <- length(clusters)

  if (G == 1) {

    matr <- matrix(1, nrow = n, ncol = 1)

  } else {

    matr <- matrix(0, nrow = n, ncol = G)

    for (g in 1:G) {
      matr[clusters == g, g] <- 1
    }

  }

  return(matr)

}

###########################################################################
###                                                                     ###
###                    Sample Statistics of Clusters                    ###
###                                                                     ###
###########################################################################

cluster_pars <- function(
    X,
    clusters
) {

  #---------------------#
  #    Input Checking   #
  #---------------------#

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X)) {
    X <- matrix(X, nrow = length(X), ncol = 1)
  }

  if (!is.numeric(X)) {
    stop('X must be a numeric matrix, data frame or vector')
  }

  if (length(clusters) != nrow(X)) {
    stop('clusters must match the number of observations')
  }

  if ( !is.vector(clusters) | !is.numeric(clusters) ) {
    stop('clusters must be a numeric vector')
  }

  # if (!all(clusters %in% 1:G)) {
  #   stop('All cluster memberships 1:G must be present')
  # }

  #--------------------------#
  #    Extract Parameters    #
  #--------------------------#

  n <- nrow(X)
  d <- ncol(X)
  G <- max(clusters)

  N <- table( factor(clusters, levels = 1:G) )

  if (any(N == 0)) {
    stop('At least one cluster has no observations')
  }

  py <- N / n

  incomplete <- any(is.na(X))

  if (incomplete) {
    cc       <- complete.cases(X)
    X        <- X[cc, ]
    clusters <- clusters[cc]
  }

  mu    <- matrix(NA, nrow = G, ncol = d)
  Sigma <- array(NA, dim = c(d, d, G))

  for (g in 1:G) {
    mu[g, ] <- colMeans(X[clusters == g, , drop = FALSE])

    if (N[g] == 1) {
      Sigma[, , g] <- matrix(0, nrow = d, ncol = d)
    } else {
      Sigma[, , g] <- var(X[clusters == g, , drop = FALSE])
    }

  }

  #---------------------------------#
  #    Prepare and return output    #
  #---------------------------------#

  output <- list(
    size  = N,
    pi    = py,
    mu    = mu,
    Sigma = Sigma
  )

  return(output)

}

###########################################################################
###                                                                     ###
###           Cluster Initialization using a Heuristic Method           ###
###                                                                     ###
###########################################################################

#' Cluster Initialization using a Heuristic Method
#'
#' Initialize cluster memberships and component parameters to start the EM algorithm
#' using a heuristic clustering method or user-defined labels.
#'
#' @param X An \eqn{n} x \eqn{d} matrix or data frame where \eqn{n} is the number of
#'   observations and \eqn{d} is the number of columns or variables. Alternately,
#'   \code{X} can be a vector of \eqn{n} observations.
#' @param G The number of clusters, which must be at least 1. If \code{G = 1}, then
#'   user-defined \code{clusters} is ignored.
#' @param init_method (optional) A string specifying the method to initialize
#'   the EM algorithm. "kmedoids" clustering is used by default. Alternative
#'   methods include "kmeans", "hierarchical", "manual". When
#'   "manual" is chosen, a vector \code{clusters} of length \eqn{n} must
#'   be specified. When \code{G = 1} and "kmedoids" clustering is used, the medoid
#'   will be returned, not the sample mean.
#' @param clusters A numeric vector of length \eqn{n} that specifies the initial
#'   cluster memberships of the user when \code{init_method} is set to "manual".
#'   This argument is NULL by default, so that it is ignored whenever other given
#'   initialization methods are chosen.
#'
#' @details Available heuristic methods include k-medoids clustering, k-means clustering,
#'   and hierarchical clustering. Alternately, the user can also enter pre-specified
#'   cluster memberships, making other initialization methods possible. If the given
#'   data set contains missing values, only observations with complete records will
#'   be used to initialize clusters. However, in this case, except when \code{G = 1}, the resulting cluster
#'   memberships will be set to \code{NULL} since they represent those complete records
#'   rather than the original data set as a whole.
#'
#' @return A list with the following slots:
#'   \item{pi}{Component mixing proportions.}
#'   \item{mu}{A \eqn{G} by \eqn{d} matrix where each row is the component mean vector.}
#'   \item{Sigma}{A \eqn{G}-dimensional array where each \eqn{d} by \eqn{d} matrix
#'     is the component covariance matrix.}
#'   \item{clusters}{An numeric vector with values from 1 to \eqn{G} indicating
#'     initial cluster memberships if \code{X} is a complete data set; NULL otherwise.}
#'
#' @references
#' Everitt, B., Landau, S., Leese, M., and Stahl, D. (2011). \emph{Cluster Analysis}. John Wiley & Sons. \cr \cr
#' Kaufman, L. and Rousseeuw, P. J. (2009). \emph{Finding  groups  in  data:  an
#'   introduction  to  cluster analysis}, volume 344. John Wiley & Sons. \cr \cr
#' Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A K-means clustering
#'  algorithm. \emph{Applied Statistics}, \strong{28}, 100-108. doi: 10.2307/2346830.
#'
#' @examples
#'
#' #++++ Initialization using a heuristic method ++++#
#'
# set.seed(1234)
#
# init <- initialize_clusters(iris[1:4], G = 3)
# init <- initialize_clusters(iris[1:4], G = 3, init_method = 'kmeans')
# init <- initialize_clusters(iris[1:4], G = 3, init_method = 'hierarchical')
#'
#' #++++ Initialization using user-defined labels ++++#
#'
#' init <- initialize_clusters(iris[1:4], G = 3, init_method = 'manual',
#'                             clusters = as.numeric(iris$Species))
#'
#' #++++ Initial parameters and pairwise scatterplot showing the mapping ++++#
#'
#' init$pi
#' init$mu
#' init$Sigma
#' init$clusters
#'
#' pairs(iris[1:4], col = init$clusters, pch = 16)
#'
#' @import cluster
#' @export
initialize_clusters <- function(
    X,
    G,
    init_method = c("kmedoids", "kmeans", "hierarchical", "manual"),
    clusters    = NULL
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
    stop('X must be a numeric matrix')
  }

  init_method <- match.arg(init_method)

  #------------------------------#
  #    Cluster Initialization    #
  #------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  incomplete <- any(is.na(X))

  if (incomplete) {
    cc <- complete.cases(X)
    X  <- X[cc, ]

    if (init_method == 'manual') {
      clusters <- clusters[cc]
    }
  }

  if (G == 1) {

    clusters <- rep(1, n)

    if (init_method == 'kmedoids') {
      kmed <- pam(X, G)
    }

  } else {

    #-----------------#
    #    K-Medoids    #
    #-----------------#

    if (init_method == 'kmedoids') {
      kmed     <- pam(X, G)
      clusters <- kmed$clustering
    }

    #---------------#
    #    K-Means    #
    #---------------#

    if (init_method == 'kmeans') {
      km       <- kmeans(X, G)
      clusters <- km$cluster
    }

    #-------------------------------#
    #    Hierarchical Clustering    #
    #-------------------------------#

    if (init_method == 'hierarchical') {
      hc       <- hclust(dist(X), method = "ward.D")
      clusters <- cutree(hc, k = G)
    }

    #-------------------------------#
    #    User-Defined Clustering    #
    #-------------------------------#

    if (init_method == 'manual') {

      if (is.null(clusters)) {
        stop('clusters must be specified if manual initialization is chosen')
      }

      if ( !is.vector(clusters) | !is.numeric(clusters) ) {
        stop('clusters must be a numeric vector')
      }

      if (length(clusters) != nrow(X)) {
        stop('clusters must match the number of observations')
      }

      if (!all(clusters %in% 1:G)) {
        stop('All cluster memberships 1:G must be present')
      }

    }

  }

  pars <- cluster_pars(X, clusters)

  pi <- pars$pi

  if (init_method == 'kmedoids') {
    mu <- kmed$medoids
  } else {
    mu <- pars$mu
  }

  Sigma <- pars$Sigma

  #---------------------#
  #    Prepare Output   #
  #---------------------#

  if (incomplete & G > 1) {
    clusters <- NULL
  }

  output <- list(
    pi       = pi,
    mu       = mu,
    Sigma    = Sigma,
    clusters = clusters
  )

  return(output)

}

###########################################################################
###                                                                     ###
###                  Mean Imputation Based on Clusters                  ###
###                                                                     ###
###########################################################################

#' Imputation using Cluster Means
#'
#' Replace missing values within each cluster with the corresponding cluster mean
#' obtained by other observed values. In other words, a separate mean imputation
#' is applied for every cluster.
#'
#' @param X An \eqn{n} x \eqn{d} matrix or data frame where \eqn{n} is the number of
#'   observations and \eqn{d} is the number of columns or variables. Alternately,
#'   \code{X} can be a vector of \eqn{n} observations.
#' @param clusters A numeric vector containing cluster memberships. Every integer
#'   from 1 to \code{G} must be present.
#'
#' @return A complete data matrix with missing values imputed accordingly.
#'
#' @examples
#'
#' X <- matrix(nrow = 6, ncol = 3, byrow = TRUE, c(
#'   NA,  2,  2,
#'    3, NA,  5,
#'    4,  3,  2,
#'   NA, NA,  3,
#'    7,  2, NA,
#'   NA,  4,  2
#' ))
#'
#' cluster_impute(X, clusters = c(1, 1, 1, 2, 2, 2))
#'
#' @export
cluster_impute <- function(
    X,
    clusters
) {

  #---------------------#
  #    Input Checking   #
  #---------------------#

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X)) {
    X <- matrix(X, nrow = length(X), ncol = 1)
  }

  if (!is.numeric(X)) {
    stop('X must be a numeric matrix, data frame or vector')
  }

  if (length(clusters) != nrow(X)) {
    stop('clusters must match the number of observations')
  }

  if ( !is.vector(clusters) | !is.numeric(clusters) ) {
    stop('clusters must be a numeric vector')
  }

  G <- max(clusters)

  if (!all(clusters %in% 1:G)) {
    stop('All cluster memberships 1:G must be present')
  }

  #-------------------------------#
  #    Cluster Mean Imputation    #
  #-------------------------------#

  n <- nrow(X)
  d <- ncol(X)

  for (g in 1:G) {

    Ig       <- clusters == g
    centroid <- colMeans(X[Ig, , drop = FALSE], na.rm = TRUE)

    for (h in 1:d) {
      X[Ig, h][!complete.cases(X[Ig, h])] <- centroid[h]
    }

  }

  return(X)

}

###########################################################################
###                                                                     ###
###                           Mean Imputation                           ###
###                                                                     ###
###########################################################################

#' Mean Imputation
#'
#' Replace missing values of data set by the mean of other observed values.
#'
#' @param X An \eqn{n} x \eqn{d} matrix or data frame where \eqn{n} is the number of
#'   observations and \eqn{d} is the number of columns or variables. Alternately,
#'   \code{X} can be a vector of \eqn{n} observations.
#'
#' @return A complete data matrix with missing values imputed accordingly.
#'
#' @references
#' Schafer, J. L. and Graham, J. W. (2002). Missing data: our view of the state of the art.
#'   \emph{Psychological Methods}, 7(2):147–177.  \cr \cr
#' Little, R. J. A. and Rubin, D. B. (2020). \emph{Statistical analysis with missing data}.
#'   Wiley Series in Probability and Statistics. Wiley, Hoboken, NJ, 3rd edition
#'
#' @examples
#'
#' X <- matrix(nrow = 6, ncol = 3, byrow = TRUE, c(
#'   NA,  2,  2,
#'    3, NA,  5,
#'    4,  3,  2,
#'   NA, NA,  3,
#'    7,  2, NA,
#'   NA,  4,  2
#' ))
#'
#' mean_impute(X)
#'
#' @export
mean_impute <- function(X) {

  #---------------------#
  #    Input Checking   #
  #---------------------#

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X)) {
    X <- matrix(X, nrow = length(X), ncol = 1)
  }

  if (!is.numeric(X)) {
    stop('X must be a numeric matrix, data frame or vector')
  }

  #-----------------------#
  #    Mean Imputation    #
  #-----------------------#

  n <- nrow(X)
  d <- ncol(X)

  centroid <- colMeans(X, na.rm = TRUE)

  for (h in 1:d) {
    X[!complete.cases(X[, h]), h] <- centroid[h]
  }

  return(X)

}


###################################################################
###                                                             ###
###       Randomly introduce missing values to a data set       ###
###                                                             ###
###################################################################

#' Missing Values Generation
#'
#' A convenient function that randomly introduces missing values to an at-least-bivariate
#' data set. The user can specify either the proportion of observations that contain some
#' missing values or the exact number of observations that contain some missing values.
#' Note that the function does not guarantee that underlying missing-data
#' mechanism to be missing at random (MAR).
#'
#' @param X An \eqn{n} by \eqn{d} matrix or data frame where \eqn{n} is the number of
#'   observations and \eqn{d} is the number of columns or variables. \code{X} must
#'   have at least 2 rows and 2 columns.
#' @param prop_cases (optional) Proportion of observations that contain some missing values.
#'   \code{prop_cases} must be a number in \eqn{(0, 1)}. \code{prop_cases = 0.1}
#'   by default, but will be ignored if \code{n_cases} is specified.
#' @param n_cases (optional) Number of observations that contain some missing values.
#'   \code{n_cases} must be an integer ranging from 1 to \code{nrow(X) - 1}.
#'
#' @return The orginal \eqn{n} by \eqn{d} matrix or data frame with missing values.
#'
#' @details If subject to missingness, an observation can have at least 1 and at
#'   most \code{ncol(X) - 1} missing values. Depending on the data
#'   set, it is not guaranteed that the resulting matrix will have the number of
#'   rows with missing values matches the specified proportion.
#'
#' @examples
#' set.seed(1234)
#'
#' hide_values(iris[1:4])
#' hide_values(iris[1:4], prop_cases = 0.5)
#' hide_values(iris[1:4], n_cases = 80)
#'
#' @export
hide_values <- function(
    X,                # numeric data matrix or data frame
    prop_cases = 0.1, # proportion of observations subject to missingness
    n_cases = NULL    # number of observations with missing values; if specified, prop_cases is ignored
) {

  if (!is.matrix(X) & !is.data.frame(X)) {
    stop('X must be a matrix or data frame')
  }

  n <- nrow(X)
  d <- ncol(X)

  if (n < 2) {
    stop('X must have at least 2 rows')
  }

  if (d < 2) {
    stop('X must have at least 2 columns')
  }

  if (is.null(n_cases)) {
    if (prop_cases <= 0 | prop_cases >= 1) {
      stop('prop_cases must be in (0, 1)')
    }

    n_cases <- floor(prop_cases * n)

    if (n_cases < 1) {
      stop('prop_cases is too small')
    }
  } else {
    if (n_cases < 1 | n_cases > n - 1 | n_cases %% 1 != 0) {
      stop('n_cases must be an integer in (1, n - 1)')
    }
  }

  indices_to_hide <- sample(1:n, n_cases)

  for (i in indices_to_hide) {
    n_cols_to_hide     <- sample(1:(d - 1), 1)
    cols_to_hide       <- sample(1:d, n_cols_to_hide)
    X[i, cols_to_hide] <- NA
  }

  return(X)
}

##########################################################
###                                                    ###
###       Generate general missing-data patterns       ###
###                                                    ###
##########################################################

#' Missing-Data Pattern Generation
#'
#' Generate all possible missing patterns in a multivariate data set. The function
#' can be used to complement the function \code{ampute()} from package \code{mice}
#' in which a matrix of patterns is needed to allow for general missing-data
#' patterns with missing-data mechanism missing at random (MAR). Using this
#' function, each observation can have more than one missing value.
#'
#' @param d The number of variables or columns of the data set. \code{d} must be
#'   an integer greater than 1.
#'
#' @return A matrix where 0 indicates that a variable should have missing values
#'   and 1 indicates that a variable should remain complete. This matrix has \code{d}
#'   columns and \eqn{2^d - 2} rows.
#'
#' @details An observation cannot have all values missing values. A complete observation
#'   is not qualified for missing-data pattern. Note that a large value of \code{d} may
#'   result in memory allocation error.
#'
#' @examples
#' generate_patterns(4)
#'
#' #++++ To use with the function ampute() from package mice ++++#
#' library(mice)
#'
#' patterns_matr <- generate_patterns(4)
#' data_missing <- ampute(iris[1:4], prop = 0.5, patterns = patterns_matr)$amp
#'
#' @export
generate_patterns <- function(d) {

  if (d < 2) {
    stop('Number of variables must be at least 2')
  }

  if (d %% 1 != 0) {
    stop('Number of variables d must be an integer')
  }

  matr <- expand.grid(replicate(d, 0:1, simplify = FALSE))[c(-1, -2^d), ]
  matr <- as.matrix(matr)
  matr <- matr[order(rowSums(matr), decreasing = F), ]
  matr <- matr[nrow(matr):1, ]

  matr           <- as.matrix(matr)
  rownames(matr) <- NULL
  colnames(matr) <- NULL

  return(matr)

}

############################################################################
###                                                                      ###
###                   Binary Classification Evaluation                   ###
###                                                                      ###
############################################################################

#' Binary Classification Evaluation
#'
#' Evaluate the performance of a classification model by comparing its predicted
#' labels to the true labels. Various metrics are returned to give an insight on
#' how well the model classifies the observations. This function is added to aid
#' outlier detection evaluation of MCNM and MtM in case that true outliers are
#' known in advance.
#'
#' @param true_labels An 0-1 or logical vector denoting the true labels. The
#'   meaning of 0 and 1 (or TRUE and FALSE) is up to the user.
#' @param pred_labels An 0-1 or logical vector denoting the true labels. The
#'   meaning of 0 and 1 (or TRUE and FALSE) is up to the user.
#'
#' @return A list with the following slots:
#'   \item{matr}{The confusion matrix built upon true labels and predicted labels.}
#'   \item{TN}{True negative.}
#'   \item{FP}{False positive (type I error).}
#'   \item{FN}{False negative (type II error).}
#'   \item{TP}{True positive.}
#'   \item{TPR}{True positive rate (sensitivy).}
#'   \item{FPR}{False positive rate.}
#'   \item{TNR}{True negative rate (specificity).}
#'   \item{FNR}{False negative rate.}
#'   \item{precision}{Precision or positive predictive value (PPV).}
#'   \item{accuracy}{Accuracy.}
#'   \item{error_rate}{Error rate.}
#'   \item{FDR}{False discovery rate.}
#'
#' @examples
#'
#' #++++ Inputs are 0-1 vectors ++++#
#'
#' evaluation_metrics(
#'   true_labels = c(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1),
#'   pred_labels = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1)
#' )
#'
#' #++++ Inputs are logical vectors ++++#
#'
#' evaluation_metrics(
#'   true_labels = c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
#'   pred_labels = c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE)
#' )
#'
#' @export
evaluation_metrics <- function(
    true_labels,   # a 0-1 or logical vector
    pred_labels    # a 0-1 or logical vector
) {

  #----------------------#
  #    Input checking    #
  #----------------------#

  if (!is.vector(true_labels) | !is.vector(pred_labels)) {
    stop('true_labels and pred_labels must be vectors')
  }

  if (any(is.na(true_labels)) | any(is.na(pred_labels))) {
    stop('true_labels and pred_labels must not contain NA')
  }

  if (is.logical(true_labels)) {
    true_labels <- as.numeric(true_labels)
  }

  if (!is.numeric(true_labels)) {
    stop('true_labels must be numeric or logical')
  }

  if (is.logical(pred_labels)) {
    pred_labels <- as.numeric(pred_labels)
  }

  if (!is.numeric(pred_labels)) {
    stop('pred_labels must be numeric or logical')
  }

  if (length(true_labels) != length(pred_labels)) {
    stop('true_labels and pred_labels must have the same length')
  }

  #------------------------#
  #    Confusion matrix    #
  #------------------------#

  # cat('
  #   Structure of the confusion matrix:
  #
  #                     pred_labels
  #                         0     1
  #   true_labels        0 TN    FP
  #                      1 FN    TP\n\n')

  matr <- matrix(0, nrow = 2, ncol = 2,
                 dimnames = list(c('true_0', 'true_1'),
                                 c('pred_0', 'pred_1')))

  tn <- sum(true_labels == 0 & pred_labels == 0)    # true negative
  fp <- sum(true_labels == 0 & pred_labels == 1)    # false positive
  fn <- sum(true_labels == 1 & pred_labels == 0)    # false negative
  tp <- sum(true_labels == 1 & pred_labels == 1)    # true positive

  matr[1, 1] <- tn
  matr[1, 2] <- fp
  matr[2, 1] <- fn
  matr[2, 2] <- tp

  p   <- fn + tp    # real positive cases
  n   <- tn + fp    # real negative cases

  tpr <- tp / (tp + fn)    # true positive rate (sensitivy)
  fpr <- fp / (fp + tn)    # false positive rate (specificity)
  tnr <- tn / (tn + fp)    # true negative rate
  fnr <- fn / (fn + tp)    # false negative rate

  preci <- tp / (tp + fp)         # precision
  acc   <- (tn + tp) / (p + n)    # accuracy
  err   <- 1 - acc                # error rate
  fdr   <- fp / (fp + tp)         # false discovery rate

  output <- list(
    matr       = matr,
    TN         = tn,
    FP         = fp,
    FN         = fn,
    TP         = tp,
    TPR        = tpr,
    FPR        = fpr,
    TNR        = tnr,
    FNR        = fnr,
    precision  = preci,
    accuracy   = acc,
    error_rate = err,
    FDR        = fdr
  )

  return(output)
}


############################################################################
###                                                                      ###
###       Approximate the Asymptotic Maximum of the Log-Likelihood       ###
###                                                                      ###
############################################################################

getall <- function(loglik) {
  if (length(loglik) < 3) {
    return(Inf)
  }
  n       <- length(loglik)
  lm1     <- loglik[n]
  lm      <- loglik[(n-1)]
  lm_1    <- loglik[(n-2)]
  am      <- (lm1 - lm)/(lm - lm_1)
  lm1.Inf <- lm + (lm1 - lm)/(1-am)
  val     <- lm1.Inf - lm

  if (is.nan(val)) {
    val = 0
  }

  if (val < 0) {
    val = 1
  }

  return( val )
}

############################################################################
###                                                                      ###
###                       Trace of a Square Matrix                       ###
###                                                                      ###
############################################################################

Tr <- function(X) {

  if (!is.matrix(X)) {
    stop('X must be a matrix')
  }

  if (nrow(X) != ncol(X)) {
    stop('X must be a square matrix')
  }

  return(
    sum(diag(X))
  )

}
