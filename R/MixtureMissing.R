#######################################
###                                 ###
###       Plot MixtureMissing       ###
###                                 ###
#######################################

#' Mixture Missing Plotting
#'
#' Provide a parallel plot of up to the first 10 variables of a multivariate
#' data sets, and a line plot showing log-likelihood values at every iteration
#' during the EM algorithm. When applicable, pairwise scatter plots highlighting
#' outliers denoted by triangles and/or observations whose values are missing but are replaced by
#' expectations obtained in the EM algorithm will be included.
#'
#' @param x A \code{MixtureMissing} object.
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#'
#' @return No return value, called to visualize the fitted model's results
#'
#' @examples
#'
#' data('auto')
#'
#' #++++ With no missing values ++++#
#'
#' X <- auto[, c('engine_size', 'city_mpg', 'highway_mpg')]
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#' plot(mod)
#'
#' #++++ With missing values ++++#
#'
#' X <- auto[, c('normalized_losses', 'horsepower', 'highway_mpg', 'price')]
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#' plot(mod)
#'
#' @importFrom graphics axis barplot legend matplot plot
#' @importFrom graphics pairs par
#' @importFrom MASS parcoord
#' @export plot.MixtureMissing
#' @export
plot.MixtureMissing <- function(x, ...) {

  pch_good       <- 1
  pch_bad        <- 2
  col_complete   <- 'black'
  col_incomplete <- '#e41a1c'

  model    <- x$model
  clusters <- x$clusters
  G        <- length(x$pi)
  dat      <- x$data
  d        <- ncol(dat)
  d        <- ifelse(is.null(d), 1, d)

  outliers <- x$outliers
  if (is.null(outliers)) {
    outliers <- rep(FALSE, length(clusters))
  }

  if (G <= 8) {
    cols <- c('#377eb8', '#e41a1c', '#4daf4a', '#984ea3',
              '#ff7f00', '#a65628', '#f781bf', '#ffff33')
    col_clusters <- rep(NA, G)
    for (g in 1:G) {
      col_clusters[clusters == g] <- cols[g]
    }
  } else {
    col_clusters <- clusters
  }

  if (d > 1) {

    if (d < 10) {
      # print(
      #   ggparcoord(data = cbind(dat, clusters), groupColumn = d + 1, columns = 1:d) +
      #   theme_bw() + ggtitle('Parallel Coordinate Plot')
      # )

      parcoord(dat, col= col_clusters, main = 'Parallel Coordinate Plot')

      if (d > 2) {
        pairs(dat, col = col_clusters, pch = ifelse(outliers, pch_bad, pch_good),
              main = 'Component Memberships')
      } else {
        plot(dat, col = col_clusters, pch = ifelse(outliers, pch_bad, pch_good),
             main = 'Component Memberships')
      }

      if ( grepl('incomplete', model) ) {
        if (d > 2) {
          pairs(dat, col = ifelse(x$complete, col_complete, col_incomplete),
                pch = ifelse(outliers, pch_bad, pch_good),
                main = 'Complete v.s. Incomplete Observations')
        } else {
          plot(dat, col = ifelse(x$complete, col_complete, col_incomplete),
               pch = ifelse(outliers, pch_bad, pch_good),
               main = 'Complete v.s. Incomplete Observations')
        }

      }

    } else {

      # print(
      #   ggparcoord(data = cbind(clusters, dat), groupColumn = 1, columns = 2:11) +
      #   theme_bw() + ggtitle('Parallel Coordinate Plot - First 10 varaibles')
      # )

      parcoord(dat[, 1:10], col = col_clusters, main = 'Parallel Coordinate Plot')

    }
  }

  #++++ Log-likelihood over iterations ++++#

  plot(x$loglik, type = 'b', pch = 16, xlab = 'Iteration', ylab = 'Log-Likelihood',
       main = 'Log-Likelihood Values over EM Iterations')

}

############################################
###                                      ###
###       Summarize MixtureMissing       ###
###                                      ###
############################################

#' Summary for Mixture Missing
#'
#' Summarizes main information regarding a \code{MixtureMissing} object.
#'
#' @param object A \code{MixtureMissing} object.
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#'
#' @return No return value, called to summarize the fitted model's results
#'
#' @details Information includes the model used to fit the data set, initialization
#'   method, clustering table, total outliers, outliers per cluster, mixing proportions,
#'   component means and variances, final log-likelihood value, information criteria.
#'
#' @examples
#'
#' #++++ With no missing values ++++#
#'
#' X <- auto[, c('horsepower', 'highway_mpg', 'price')]
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#' # summary(mod)
#'
#' #++++ With missing values ++++#
#'
#' X <- auto[, c('normalized_losses', 'horsepower', 'highway_mpg', 'price')]
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#' summary(mod)
#'
#' @export summary.MixtureMissing
#' @export
summary.MixtureMissing <- function(object, ...) {

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

  incomplete_data <- grepl('incomplete', object$model)

  G       <- length(object$pi)
  G_comps <- paste(G, '-Component', sep = '')

  if (incomplete_data) {
    model_code <- gsub('_incomplete_data', '', object$model)
    cat('\nModel:', G_comps, m_names[match(model_code, m_codes)], 'Mixture with Incomplete Data\n')
  } else {
    model_code <- gsub('_complete_data', '', object$model)
    cat('\nModel:', G_comps, m_names[match(model_code, m_codes)], 'Mixture with Complete Data\n')
  }

  if (incomplete_data) {
    cat('\nObservations with missing values:', sum(!object$complete), '/', nrow(object$data), '\n')
  }

  cat('\nIterations:', object$iter_stop, '/', object$max_iter)



  if (G == 1) {
    cat('\n\nInitialization: None')
  } else {
    cat('\n\nInitialization:', object$init_method)
  }

  cat("\n\nComponent frequency table:\n")
  tab           <- matrix(NA, nrow = 1, ncol = G)
  rownames(tab) <- ''
  colnames(tab) <- paste('comp', 1:G, sep = '')
  for (g in 1:G) {
    tab[, g] <- sum(object$clusters == g)
  }
  print(tab)

  if (model_code %in% c('CN', 't')) {
    cat('\nTotal outliers:', sum(object$outliers), '\n')
    cat('\nOutliers per component:\n')

    outliers <- object$outliers
    if (is.null(outliers)) {
      outliers <- rep(FALSE, length(object$clusters))
    }

    for (g in 1:G) {
      tab[, g] <- sum(object$clusters == g & outliers)
    }

    print(tab)
  }

  cat('\nMixing proportions:\n')
  print(object$pi)

  cat('\nComponent means:\n')
  print(object$mu)

  cat('\nComponent variances:\n')
  print(object$Sigma)

  cat('\nFinal log-likelihood:', object$final_loglik, '\n')

  cat('\nTotal parameters:', object$npar$total, '\n')

  cat('\nInformation Criteria:\n')
  print(
    data.frame(
      AIC  = object$AIC,
      BIC  = object$BIC,
      KIC  = object$KIC,
      KICc = object$KICc,
      AIC3 = object$AIC3,
      CAIC = object$CAIC,
      AICc = object$AICc,
      ICL  = object$ICL,
      AWE  = object$AWE,
      CLC  = object$CLC,
      row.names = ''
    )
  )

  cat('\n')

}
