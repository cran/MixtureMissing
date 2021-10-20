###########################################################################
###                                                                     ###
###           Cluster initialization using a heuristic method           ###
###                                                                     ###
###########################################################################

#' Cluster Initialization
#'
#' Initialize cluster memberships and component parameters to start the EM algorithm
#' using a heuristic clustering method or user-defined labels.
#'
#' @param X An \eqn{n} by \eqn{d} matrix or data frame where \eqn{n} is the number of
#'   observations and \eqn{d} is the number of columns or variables. Alternately,
#'   \code{X} can be a vector of \eqn{n} observations.
#' @param G The number of clusters.
#' @param init_method (optional) A string specifying the method to initialize
#'   the EM algorithm. "kmedoids" clustering is used by default. Alternative
#'   methods include "kmeans", "hierarchical", "manual", "soft", "hard". When
#'   "manual" is chosen, a vector \code{manual_clusters} of length \eqn{n} must
#'   be specified.
#' @param manual_clusters A vector of length \eqn{n} that specifies the initial
#'   cluster memberships of the user when \code{init_method} is set to "manual".
#'   Both numeric and character vectors are acceptable. This argument is NULL by
#'   default, so that it is ignored whenever other given initialization methods
#'   are chosen.
#'
#' @details Available heuristic methods include k-medoids clustering, k-means clustering,
#'   hierarchical clustering, soft and hard clustering. Alternately, the user can also
#'   enter pre-specified cluster memberships, making other initialization methods possible.
#'
#' @return A list with the following slots:
#'   \item{z}{Mapping probabilities in the form of an \eqn{n} by \eqn{G} matrix.}
#'   \item{pi}{Component mixing proportions.}
#'   \item{mu}{If \code{X} is a matrix or data frame, \code{mu} is an \eqn{G} by \eqn{d}
#'     matrix where each row is the component mean vector. If \code{X} is a vector, \code{mu}
#'     is a vector of \eqn{G} component means.}
#'   \item{sigma}{If \code{X} is a matrix or data frame, \code{sigma} is a \eqn{G}-dimensional
#'     array where each \eqn{d} by \eqn{d} matrix is the component covariance matrix. If
#'     \code{X} is a vector, \code{sigma} is a vector of \eqn{G} component variances.}
#'
#' @details
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
#' set.seed(1234)
#'
#' init <- initialize_clusters(iris[1:4], G = 3)
#' init <- initialize_clusters(iris[1:4], G = 3, init_method = 'kmeans')
#' init <- initialize_clusters(iris[1:4], G = 3, init_method = 'hierarchical')
#' init <- initialize_clusters(iris[1:4], G = 3, init_method = 'soft')
#' init <- initialize_clusters(iris[1:4], G = 3, init_method = 'hard')
#'
#' #++++ Initialization using user-defined labels ++++#
#'
#' init <- initialize_clusters(iris[1:4], G = 3, init_method = 'manual',
#'                             manual_clusters = iris$Species)
#'
#' #++++ Initial parameters and pairwise scatterplot showing the mapping ++++#
#'
#' init$z
#' init$pi
#' init$mu
#' init$sigma
#'
#' pairs(iris[1:4], col = init$clusters, pch = 16)
#'
#' @import cluster
#' @export
initialize_clusters <- function(
  X,  # numeric data matrix; data frame will be converted to matrix
  G,  # number of clusters
  init_method = c("kmedoids", "kmeans", "hierarchical", "manual", "soft", "hard"),
  manual_clusters = NULL
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

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X)) {
    X <- matrix(X, nrow = length(X), ncol = 1)
  }

  if (!is.numeric(X)) {
    stop('X must be a numeric matrix')
  }

  if (any(is.na(X))) {
    X <- X[complete.cases(X), ]
  }

  init_method <- match.arg(init_method)

  n <- nrow(X) # number of observations
  d <- ncol(X) # number of variables/dimensions

  z <- matrix(0, nrow = n, ncol = G)

  #-----------------------------------------------#
  #  Initialization of z and cluster memberships  #
  #-----------------------------------------------#

  if (init_method == 'kmedoids') {
    kmed <- pam(X, G)
    clusters <- kmed$clustering

    for (g in 1:G) {
      z[clusters == g, g] <- 1
    }
  }

  if (init_method == 'kmeans') {
    km <- kmeans(X, G)
    clusters <- km$cluster

    for (g in 1:G) {
      z[clusters == g, g] <- 1
    }
  }

  if (init_method == 'hierarchical') {
    hc <- hclust(dist(X), method = "ward.D")
    clusters <- cutree(hc, k = G)
    for (g in 1:G) {
      z[clusters == g, g] <- 1
    }
  }

  if (init_method == 'manual') {
    if (is.null(manual_clusters)) {
      stop('Vector manual_clusters must be specified if manual initialization is chosen')
    }

    if (is.factor(manual_clusters)) {
      manual_clusters <- as.character(manual_clusters)
    }

    if (!is.vector(manual_clusters)) {
      stop('manual_clusters must be a vector')
    }

    if (length(manual_clusters) != n) {
      stop('The length of manual_clusters must be the same as the number of rows')
    }

    clusters <- match(manual_clusters, sort(unique(manual_clusters)))
    for (g in 1:G) {
      z[clusters == g, g] <- 1
    }
  }

  if (init_method == 'soft') {
    z  <- matrix(runif(n * G), nrow = n, ncol = G)
    z  <- z / rowSums(z)             #
    clusters <- apply(z, 1, which.max)
  }

  if (init_method == 'hard') {
    z <- t(rmultinom(n, size = 1, prob = rep(1/G, G)))
    clusters <- apply(z, 1, which.max)
  }

  #--------------------------------------------------------------#
  #  Parameters according the initial z and cluster memberships  #
  #--------------------------------------------------------------#

  py <- colSums(z) / n

  if (d > 1) {
    mu <- t(sapply(1:G, function(g) {
      X_g <- X[clusters == g, ]
      if (is.vector(X_g)) {
        return(X_g)
      }
      return(colMeans(X_g))
    }))

    sigma <- sapply(1:G, function(g) {
      X_g <- X[clusters == g, ]
      if (is.vector(X_g)) {
        return(matrix(0, nrow = d, ncol = d))
      }
      return(cov(X[clusters == g, ]))
    }, simplify = 'array')

    # rownames(mu) <- paste("cluster_", 1:G, sep = "")
  } else {
    mu <- sapply(1:G, function(g) mean(X[clusters == g, ]))
    sigma <- sapply(1:G, function(g) var(X[clusters == g, ]))
  }

  # names(sigma) <- paste("cluster_", 1:G, sep = "")

  return(list(z = z, clusters = clusters, pi = py, mu = mu, sigma = sigma))
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

  n <- nrow(X) # number of observations
  d <- ncol(X) # number of variables/dimensions

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

  matr <- as.matrix(matr)
  rownames(matr) <- NULL
  colnames(matr) <- NULL

  return(matr)
}

############################################################################
###                                                                      ###
###       Approximate the asymptotic maximum of the log-likelihood       ###
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
