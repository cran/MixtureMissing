#' Mixture Missing Plotting
#'
#' @param x A \code{MixtureMissing} object.
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#'
#' @return No return value, called to visualize the fitted model's results
#'
#' @examples
#'
#' data('nm_5_noise_close_100')
#'
#' #++++ With no missing values ++++#
#'
#' X <- nm_5_noise_close_100[, 1:2]
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#' plot(mod)
#'
#' #++++ With missing values ++++#
#'
#' set.seed(1234)
#'
#' X <- hide_values(nm_5_noise_close_100[, 1:2], prop_cases = 0.1)
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#' plot(mod)
#'
#' @importFrom graphics axis barplot legend matplot plot
#' @importFrom GGally ggparcoord
#' @importFrom ggplot2 ggplot ggtitle theme_bw
#' @importFrom graphics pairs par
#' @export plot.MixtureMissing
#' @export
plot.MixtureMissing <- function(x, ...) {

  pch_good       <- 16
  pch_bad        <- 8
  col_complete   <- 'black'
  col_incomplete <- 'tomato'

  model    <- x$model
  clusters <- as.character(x$clusters)
  outliers <- x$outliers
  dat      <- x$data
  d        <- ncol(dat)
  d        <- ifelse(is.null(d), 1, d)

  if (d > 1) {
    if (d < 10) {
      ggparcoord(data = cbind(dat, clusters), groupColumn = d + 1, columns = 1:d) +
        theme_bw() + ggtitle('Parallel Coordinate Plot')

      if (d > 2) {
        pairs(dat, col = clusters, pch = ifelse(outliers, pch_bad, pch_good),
              main = 'Cluster Memberships')
      } else {
        plot(dat, col = clusters, pch = ifelse(outliers, pch_bad, pch_good),
             main = 'Cluster Memberships')
      }

      if (model %in% c('MCNM_incomplete_data', 'MtM_incomplete_data', 'MNM_incomplete_data')) {
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
      ggparcoord(data = cbind(clusters, dat), groupColumn = 1, columns = 2:11) +
        theme_bw() + ggtitle('Parallel Coordinate Plot - First 10 varaibles')
    }
  }

  #++++ Log-likelihood over iterations ++++#
  plot(x$loglik, type = 'b', pch = 16, xlab = 'Iteration', ylab = 'Log-Likelihood')
}

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
#'   component means and variances.
#'
#' @examples
#'
#' data('nm_5_noise_close_100')
#'
#' #++++ With no missing values ++++#
#'
#' # X <- nm_5_noise_close_100[, 1:2]
#' # mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#' # summary(mod)
#'
#' #++++ With missing values ++++#
#'
#' set.seed(1234)
#'
#' X <- hide_values(nm_5_noise_close_100[, 1:2], prop_cases = 0.1)
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#' summary(mod)
#'
#' @export summary.MixtureMissing
#' @export
summary.MixtureMissing <- function(object, ...) {

  cat('\nModel:', object$model)

  if (object$model %in% c('MCNM_incomplete_data', 'MtM_incomplete_data', 'MNM_incomplete_data')) {
    cat('\nNumber of observations with missing values:', sum(!object$complete), '/', nrow(object$data))
  }

  cat('\nIterations:', object$iter_stop, '/', object$max_iter)

  cat('\nInitialization method:', object$init_method)

  cat("\n\nClustering table:")
  print(table(object$clusters))

  if (object$model %in% c('MCNM_incomplete_data', 'MtM_incomplete_data', 'MCNM_complete_data', 'MtM_complete_data')) {
    cat('\nTotal outliers', sum(object$outliers))
    cat('\nOutliers per cluster:\n')
    print(table(object$clusters, object$outliers)[, 2])
  }

  cat('\nMixing proportions:\n')
  print(object$pi)

  cat('\nComponent means:\n')
  print(object$mu)

  cat('\nComponent variances:\n')
  print(object$sigma)

  cat('\nInformation Criteria:\n')
  print(data.frame(
    AIC  = object$AIC,
    BIC  = object$BIC,
    KIC  = object$KIC,
    KICc = object$KICc,
    AIC3 = object$AIC3,
    CAIC = object$CAIC,
    AICc = object$AICc,
    ICL  = object$ICL,
    AWE  = object$AWE,
    CLC  = object$CLC
  ))
}