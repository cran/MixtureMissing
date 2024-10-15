############################################
###                                      ###
###       Summarize MixtureMissing       ###
###                                      ###
############################################

#' Summary for MixtureMissing
#'
#' Summarizes main information regarding a \code{MixtureMissing} object.
#'
#' @param object A \code{MixtureMissing} object or an output of \link[MixtureMissing]{select_mixture}.
#'   In the latter, only the best model will be considered.
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
#' summary(mod)
#'
#' #++++ With missing values ++++#
#'
#' X <- auto[, c('normalized_losses', 'horsepower', 'highway_mpg', 'price')]
#' mod <- MCNM(X, G = 2, init_method = 'kmedoids', max_iter = 10)
#' summary(mod)
#'
#' @export
summary.MixtureMissing <- function(object, ...) {

  if (!inherits(object, "MixtureMissing")) {
    object <- try(object$best_mod, silent = TRUE)

    if (!inherits(object, "MixtureMissing")) {
      stop("object not of class 'MixtureMissing'")
    }
  }

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
    mis_cases <- apply( object$complete, 1, function(ro) any(!ro) )
    cat('\nObservations with missing values:', sum(mis_cases), '/', nrow(object$data), '\n')
    cat('\nMissing values per variable:\n')
    complete_matr <- object$complete
    if (is.null(colnames(complete_matr))) {
      colnames(complete_matr) <- paste0('V', 1:ncol(complete_matr))
    }
    print(colSums(!complete_matr))
  }

  cat('\nIterations:', object$iter_stop, '/', object$max_iter)

  if (G == 1) {
    cat('\n\nInitialization: None')
  } else {
    cat('\n\nInitialization:', object$init_method)
  }

  cat("\n\nComponent frequency table:\n")
  tab           <- matrix(NA_real_, nrow = 1, ncol = G)
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

  cat('\nComponent location vectors:\n')
  print(object$mu)

  cat('\nComponent dispersion matrices:\n')
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

########################################
###                                  ###
###       Print MixtureMissing       ###
###                                  ###
########################################

#' Print for MixtureMissing
#'
#' Print \code{MixtureMissing} object.
#'
#' @param x A \code{MixtureMissing} object or an output of \link[MixtureMissing]{select_mixture}.
#'   In the latter, only the best model will be considered.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return No return value, called to print the fitted model's description.
#'
#' @details The description includes information on the complete or incomplete data,
#'   number of clusters, and component distribution.
#'
#' @examples
#'
#' #++++ With no missing values ++++#
#'
#' X <- iris[, 1:4]
#' mod <- MGHM(X, G = 2, model = 'GH', max_iter = 10)
#' print(mod)
#'
#' #++++ With missing values ++++#
#'
#' set.seed(123)
#' X <- hide_values(iris[, 1:4], n_cases = 20)
#' mod <- MGHM(X, G = 2, model = 'GH', max_iter = 10)
#' print(mod)
#'
#' @export
print.MixtureMissing <- function(x, ...) {

  object <- x

  if(!inherits(object, "MixtureMissing")) {
    object <- try(x$best_mod, silent = TRUE)

    if (!inherits(object, "MixtureMissing")) {
      stop("object not of class 'MixtureMissing'")
    }
  }

  G <- length(object$pi)

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

  if (incomplete_data) {
    model_code <- gsub('_incomplete_data', '', object$model)
  } else {
    model_code <- gsub('_complete_data', '', object$model)
  }

  cat('\nThe fitted mixture model is based on ', ifelse(incomplete_data, 'incomplete', 'complete'), ' data, G = ', G,
      ', and the ', m_names[match(model_code, m_codes)], ' distribution.\n\n', sep = '')

}

#######################################
###                                 ###
###       Plot MixtureMissing       ###
###                                 ###
#######################################

#' MixtureMissing Plotting
#'
#' Provide four model-based clustering plots for a \code{MixtureMissing} object. The options
#' include (1) pairwise scatter plots showing cluster memberships and highlighting outliers denoted by triangles;
#' (2) pairwise scatter plots highlighting in red observations whose values are missing but are replaced by
#' expectations obtained in the EM algorithm; (3) parallel plot of up to the first 10 variables of a multivariate
#' data set; and (4) plots of estimated density in the form of contours. A single or multiple options
#' can be specified. In the latter case, interactive mode will be triggered for the user to choose.
#'
#' @param x A \code{MixtureMissing} object or an output of \link[MixtureMissing]{select_mixture}.
#'   In the latter, only the best model will be considered.
#' @param what A string or a character vector specifying the desired plots. See the details section for
#'   a list of available plots.
#' @param nlevels Number of contour levels desired; 15 by default.
#' @param drawlabels Contour levels are labelled if \code{TRUE}.
#' @param addpoints Colored points showing cluster memberships are added if \code{TRUE}.
#' @param cex.point A numerical value giving the amount by which data points should be magnified relative to the default.
#' @param cex.axis The magnification to be used for axis annotation.
#' @param cex.labels A numerical value to control the character size of variable labels.
#' @param lwd The contour line width, a positive number, defaulting to 1.
#' @param col_line The color of contour; "gray" by default.
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#'
#' @return No return value, called to visualize the fitted model's results
#'
#' @details The plots that can be retrieved include
#' \itemize{
#'   \item If \code{what = "classification"} - Pairwise scatter plots showing cluster memberships
#'      and highlighting outliers denoted by triangles.
#'   \item If \code{what = "missing"} - Pairwise scatter plots highlighting in red observations
#'     whose values are missing but are replaced by expectations obtained in the EM algorithm.
#'   \item If \code{what = "parallel"} - Parallel plot of up to the first 10 variables of a multivariate
#'     data set.
#'   \item If \code{what = "density"} - Plots of estimated density in the form of contours.
#' }
#'
#' @examples
#'
#' set.seed(123)
#' X <- hide_values(iris[, 1:4], n_cases = 20)
#' mod <- MCNM(X, G = 2, max_iter = 10)
#' plot(mod, what = 'classification')
#'
#' @importFrom graphics axis barplot legend matplot plot layout points text contour
#' @importFrom graphics pairs par
#' @importFrom MASS parcoord
#' @importFrom utils menu
#' @export
plot.MixtureMissing <- function(
    x,
    what       = c("classification", "missing", "parallel", "density"),
    nlevels    = 15,
    drawlabels = TRUE,
    addpoints  = TRUE,
    cex.point  = 1,
    cex.axis   = 1,
    cex.labels = 2,
    lwd        = 1,
    col_line   = "gray",
    ...
) {

  #++++ Save current graphical parameters to reset ++++#

  oldpar <- par(no.readonly = TRUE)

  #++++ Extract object information ++++#

  object <- x

  if (!inherits(object, "MixtureMissing")) {
    object <- try(x$best_mod, silent = TRUE)

    if (!inherits(object, "MixtureMissing")) {
      stop("object not of class 'MixtureMissing'")
    }
  }

  pch_good       <- 16
  pch_bad        <- 17
  col_complete   <- 'black'
  col_incomplete <- '#e41a1c'

  X_imp <- object$data
  d     <- ncol(X_imp)
  d     <- ifelse(is.null(d), 1, d)

  if (d == 1) {
    stop('The data set must be at least bivariate')
  }

  G            <- length(object$pi)
  clusters     <- object$clusters
  model        <- object$model
  complete_obs <- apply(object$complete, 1, all)

  outliers <- object$outliers
  if (is.null(outliers)) {
    outliers <- rep(FALSE, length(clusters))
  }

  #++++ Prepare color palette ++++#

  if (G <= 8) {
    cols <- c('#377eb8', '#e41a1c', '#4daf4a', '#984ea3',
              '#ff7f00', '#a65628', '#f781bf', '#ffff33')
    col_clusters <- rep(NA_real_, G)
    for (g in 1:G) {
      col_clusters[clusters == g] <- cols[g]
    }
  } else {
    col_clusters <- clusters
  }

  #++++ Plot according to user choice ++++#

  what <- unique(what)
  what <- match.arg(what, several.ok = TRUE)

  if (interactive()) {

    if (length(what) == 1) {
      choice <- 1
    } else {
      choice <- menu(what, graphics = FALSE, title = "\nModel-based clustering plots:")
    }

    while (choice != 0) {

      ### Classification

      if (what[choice] == 'classification') {

        if (d > 2) {
          pairs(X_imp, col = col_clusters, pch = ifelse(outliers, pch_bad, pch_good),
                main = '', cex = cex.point, cex.axis = cex.axis, cex.labels = cex.labels)
        } else {
          plot(X_imp, col = col_clusters, pch = ifelse(outliers, pch_bad, pch_good),
               main = '', cex = cex.point, cex.axis = cex.axis)
        }

      }

      ### Missing values imputed

      if (what[choice] == 'missing') {

        if ( grepl('incomplete', model) ) {

          if (d > 2) {
            pairs(X_imp, col = ifelse(complete_obs, col_complete, col_incomplete),
                  pch = ifelse(outliers, pch_bad, pch_good),
                  main = '', cex = cex.point, cex.axis = cex.axis, cex.labels = cex.labels)
          } else {
            plot(X_imp, col = ifelse(complete_obs, col_complete, col_incomplete),
                 pch = ifelse(outliers, pch_bad, pch_good),
                 main = '', cex = cex.point, cex.axis = cex.axis)
          }

        } else {

          if (d > 2) {
            pairs(X_imp, col = col_complete, pch = pch_good,
                  main = '', cex = cex.point, cex.axis = cex.axis, cex.labels = cex.labels)
          } else {
            plot(X_imp, col = col_complete, pch = pch_good,
                 main = '', cex = cex.point, cex.axis = cex.axis)
          }

        }

      }

      ### Parallel plot

      if (what[choice] == 'parallel') {

        if (d < 10) {
          parcoord(X_imp, col= col_clusters, main = '', cex.axis = cex.axis, lwd = lwd)
        } else{
          parcoord(X_imp[, 1:10], col = col_clusters, main = '', cex.axis = cex.axis, lwd = lwd)
        }

      }

      ### Density contour plot

      if (what[choice] == 'density') {

        # if (d > 10) {
        #   X_imp <- X_imp[, 1:10]
        #   d     <- 10
        # }

        ### Obtain model code

        # m_codes <- c('CN', 'GH',
        #              'NIG', 'SNIG',
        #              'SC', 'C',
        #              'St', 't', 'N',
        #              'SGH', 'HUM',
        #              'H', 'SH')

        incomplete_data <- grepl('incomplete', object$model)

        if (incomplete_data) {
          model_code <- gsub('_incomplete_data', '', object$model)
        } else {
          model_code <- gsub('_complete_data', '', object$model)
        }

        #++++ Fill the upper triangle with numbers starting from 1 ++++#

        matr  <- matrix(0, nrow = d, ncol = d)
        m_upp <- matrix(0, nrow = d, ncol = d)
        m_low <- matrix(0, nrow = d, ncol = d)

        len <- length(matr[upper.tri(matr, diag = FALSE)])

        # m_upp[upper.tri(m_upp)] <- seq(from = 1, by = 2, length.out = len)
        # m_low[lower.tri(m_low)] <- seq(from = 2, by = 2, length.out = len)

        m_upp[lower.tri(m_upp)] <- seq(from = 1, by = 2, length.out = len)
        m_upp <- t(m_upp)

        m_low[lower.tri(m_low)] <- seq(from = 2, by = 2, length.out = len)

        matr <- m_upp + m_low

        diag(matr) <- (max(matr) + 1):(max(matr) + d)

        below <- matr[d, seq(from = 1, by = 2, length.out = d %/% 2)]
        left  <- matr[seq(from = 2, by = 2, length.out = d %/% 2), 1]
        above <- matr[1, seq(from = 2, by = 2, length.out = d %/% 2)]
        right <- matr[seq(from = 1, by = 2, length.out = d %/% 2), d]

        layout(matr)
        par(
          oma = c(4, 4, 4, 4),
          mar = c(0.5, 0.5, 0.5, 0.5)
        )

        #++++ Plot according to the model ++++#

        zero_matr <- matrix(0, nrow = 60, ncol = 60)

        iter <- 0

        ### Contaminated Normal

        if (model_code == 'CN') {

          prior <- object$pi
          mu    <- object$mu
          Sigma <- object$Sigma
          alpha <- object$alpha
          eta   <- object$eta

          for (ii in 1:(d - 1)) {
            for (jj in (ii + 1):d) {

              in1 <- jj
              in2 <- ii
              x.v <- seq(min(X_imp[, in1]), max(X_imp[, in1]), length.out = 60)
              y.v <- seq(min(X_imp[, in2]), max(X_imp[, in2]), length.out = 60)

              XY <- cbind(x.v, y.v)

              xyS <- zero_matr

              for (g in 1:G) {
                dens_g <- outer(
                  X   = x.v,
                  Y   = y.v,
                  FUN = bivariate_dCN,
                  mu = mu[g, c(in1, in2)], Sigma = Sigma[c(in1, in2), c(in1, in2), g],
                  alpha = alpha[g], eta = eta[g]
                )

                xyS <- xyS + prior[g] * dens_g
              }

              iter <- iter + 1

              contour(x = x.v, y = y.v, z = xyS, col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in1, in2)], pch = 16, col = col_clusters, cex = cex.point)
              }

              iter <- iter + 1

              contour(x = y.v, y = x.v, z = t(xyS), col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in2, in1)], pch = 16, col = col_clusters, cex = cex.point)
              }

            }
          }    # end for loop

        }    # end if (model_code == 'CN')

        ### Generalized Hyperbolic
        ### Normal-Inverse Gaussian
        ### Skew Normal-Inverse Gaussian
        ### Symmetric Generalized Hyperbolic
        ### Hyperbolic Univarate Marginals
        ### Hyperbolic
        ### Symmetric Hyperbolic

        if (model_code %in% c('GH', 'NIG', 'SNIG', 'SGH', 'HUM', 'H', 'SH')) {

          prior  <- object$pi
          mu     <- object$mu
          Sigma  <- object$Sigma
          beta   <- object$beta
          lambda <- object$lambda
          omega  <- object$omega

          for (ii in 1:(d - 1)) {
            for (jj in (ii + 1):d) {

              in1 <- jj
              in2 <- ii
              x.v <- seq(min(X_imp[, in1]), max(X_imp[, in1]), length.out = 60)
              y.v <- seq(min(X_imp[, in2]), max(X_imp[, in2]), length.out = 60)

              XY <- cbind(x.v, y.v)

              xyS <- zero_matr

              for (g in 1:G) {
                dens_g <- outer(
                  X   = x.v,
                  Y   = y.v,
                  FUN = bivariate_dGH,
                  mu = mu[g, c(in1, in2)], Sigma = Sigma[c(in1, in2), c(in1, in2), g],
                  beta = beta[g, c(in1, in2)], lambda = lambda[g], omega = omega[g]
                )

                xyS <- xyS + prior[g] * dens_g
              }

              iter <- iter + 1

              contour(x = x.v, y = y.v, z = xyS, col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in1, in2)], pch = 16, col = col_clusters, cex = cex.point)
              }

              iter <- iter + 1

              contour(x = y.v, y = x.v, z = t(xyS), col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in2, in1)], pch = 16, col = col_clusters, cex = cex.point)
              }

            }
          }    # end for loop

        }    # end if (model_code %in% c('GH', 'NIG', 'SNIG', 'SGH', 'HUM', 'H', 'SH'))

        ### Skew-t
        ### Skew-Cauchy

        if (model_code %in% c('St', 'SC')) {

          prior  <- object$pi
          mu     <- object$mu
          Sigma  <- object$Sigma
          beta   <- object$beta

          if (model_code == 'St') {
            df <- object$df
          } else {
            df <- rep(1, G)
          }

          for (ii in 1:(d - 1)) {
            for (jj in (ii + 1):d) {

              in1 <- jj
              in2 <- ii
              x.v <- seq(min(X_imp[, in1]), max(X_imp[, in1]), length.out = 60)
              y.v <- seq(min(X_imp[, in2]), max(X_imp[, in2]), length.out = 60)

              XY <- cbind(x.v, y.v)

              xyS <- zero_matr

              for (g in 1:G) {
                dens_g <- outer(
                  X   = x.v,
                  Y   = y.v,
                  FUN = bivariate_dSt,
                  mu = mu[g, c(in1, in2)], Sigma = Sigma[c(in1, in2), c(in1, in2), g],
                  beta = beta[g, c(in1, in2)], df = df[g]
                )

                xyS <- xyS + prior[g] * dens_g
              }

              iter <- iter + 1

              contour(x = x.v, y = y.v, z = xyS, col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in1, in2)], pch = 16, col = col_clusters, cex = cex.point)
              }

              iter <- iter + 1

              contour(x = y.v, y = x.v, z = t(xyS), col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in2, in1)], pch = 16, col = col_clusters, cex = cex.point)
              }

            }
          }    # end for loop

        }    # end if (model_code %in% c('St', 'SC'))

        ### t
        ### Cauchy

        if (model_code %in% c('t', 'C')) {

          prior  <- object$pi
          mu     <- object$mu
          Sigma  <- object$Sigma

          if (model_code == 't') {
            df <- object$df
          } else {
            df <- rep(1, G)
          }

          for (ii in 1:(d - 1)) {
            for (jj in (ii + 1):d) {

              in1 <- jj
              in2 <- ii
              x.v <- seq(min(X_imp[, in1]), max(X_imp[, in1]), length.out = 60)
              y.v <- seq(min(X_imp[, in2]), max(X_imp[, in2]), length.out = 60)

              XY <- cbind(x.v, y.v)

              xyS <- zero_matr

              for (g in 1:G) {
                dens_g <- outer(
                  X   = x.v,
                  Y   = y.v,
                  FUN = bivariate_dmvt,
                  delta = mu[g, c(in1, in2)], sigma = Sigma[c(in1, in2), c(in1, in2), g], df = df[g]
                )

                xyS <- xyS + prior[g] * dens_g
              }

              iter <- iter + 1

              contour(x = x.v, y = y.v, z = xyS, col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in1, in2)], pch = 16, col = col_clusters, cex = cex.point)
              }

              iter <- iter + 1

              contour(x = y.v, y = x.v, z = t(xyS), col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in2, in1)], pch = 16, col = col_clusters, cex = cex.point)
              }

            }
          }    # end for loop

        }    # end if (model_code %in% c('t', 'C'))

        ### Normal

        if (model_code == 'N') {

          prior  <- object$pi
          mu     <- object$mu
          Sigma  <- object$Sigma

          for (ii in 1:(d - 1)) {
            for (jj in (ii + 1):d) {

              in1 <- jj
              in2 <- ii
              x.v <- seq(min(X_imp[, in1]), max(X_imp[, in1]), length.out = 60)
              y.v <- seq(min(X_imp[, in2]), max(X_imp[, in2]), length.out = 60)

              XY <- cbind(x.v, y.v)

              xyS <- zero_matr

              for (g in 1:G) {
                dens_g <- outer(
                  X   = x.v,
                  Y   = y.v,
                  FUN = bivariate_dmvnorm,
                  mean = mu[g, c(in1, in2)], sigma = Sigma[c(in1, in2), c(in1, in2), g]
                )

                xyS <- xyS + prior[g] * dens_g
              }

              iter <- iter + 1

              contour(x = x.v, y = y.v, z = xyS, col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in1, in2)], pch = 16, col = col_clusters, cex = cex.point)
              }

              iter <- iter + 1

              contour(x = y.v, y = x.v, z = t(xyS), col = col_line, lwd = lwd,
                      nlevels = nlevels, drawlabels = drawlabels,
                      xaxt = 'n', yaxt = 'n', ann = FALSE)

              if (iter %in% below) axis(1, cex.axis = cex.axis)
              if (iter %in% left)  axis(2, cex.axis = cex.axis)
              if (iter %in% above) axis(3, cex.axis = cex.axis)
              if (iter %in% right) axis(4, cex.axis = cex.axis)

              if (addpoints) {
                points(X_imp[, c(in2, in1)], pch = 16, col = col_clusters, cex = cex.point)
              }

            }
          }    # end for loop

        }    # end if (model_code == 'N')


        for(j in 1:d){
          plot(1, xaxt = 'n', yaxt = 'n', ann = FALSE,
               xlim = c(0, 1), ylim = c(0, 1), col = 'white')
          text(0.5, 0.5, colnames(X_imp)[j], cex = cex.labels)
        }

      }

      par(oldpar)

      if (length(what) == 1) {
        choice <- 0
      } else {
        choice <- menu(what, graphics = FALSE, title = "\nModel-based clustering plots:")
      }

    }

  }    # while (choice != 0)

  on.exit(par(oldpar))

}

###########################################################################
###                                                                     ###
###                Bivariate Contaminated Normal Density                ###
###                                                                     ###
###########################################################################

bivariate_dCN <- function(
    x,
    y,
    mu    = rep(0, 2),    # location
    Sigma = diag(2),      # dispersion
    alpha = 0.99,         # proportion of good observations
    eta   = 1.01,         # degree of contamination
    log   = FALSE
) {

  #----------------------#
  #    Input Checking    #
  #----------------------#

  if (length(mu) != 2) {
    stop('mu must be a vector of length 2')
  }

  if (nrow(Sigma) != 2 | ncol(Sigma) != 2) {
    stop('Sigma must be a 2 x 2 matrix')
  }

  if (alpha < 0 | alpha > 1) {
    stop('alpha must be between 0 and 1')
  }

  if (eta <= 0) {
    stop('eta must be greater than 0')
  }

  xy <- cbind(x, y)

  good_norm <- exp( mvtnorm::dmvnorm(xy, mu, Sigma, log = TRUE) )
  bad_norm  <- exp( mvtnorm::dmvnorm(xy, mu, eta * Sigma, log = TRUE) )

  den <- alpha * good_norm + (1 - alpha) * bad_norm
  den[den <= 10^(-323)] <- 10^(-323)

  return(den)

}

############################################################################
###                                                                      ###
###               Bivariate Generalized Hyperbolic Density               ###
###                                                                      ###
############################################################################


bivariate_dGH <- function(
    x,
    y,
    mu     = rep(0, 2),    # location
    Sigma  = diag(2),      # dispersion
    beta   = rep(0, 2),    # skewness
    lambda = 0.5,          # index
    omega  = 1,            # concentration
    log    = FALSE
) {

  #----------------------#
  #    Input Checking    #
  #----------------------#

  if (length(mu) != 2) {
    stop('mu must be a vector of length 2')
  }

  if (nrow(Sigma) != 2 | ncol(Sigma) != 2) {
    stop('Sigma must be a 2 x 2 matrix')
  }

  if (length(beta) != 2) {
    stop('beta must be a vector of length 2')
  }

  if (omega <= 0) {
    stop('omega must be positive')
  }

  if (is.nan(log(det(Sigma)))) {
    Sigma <- diag(2)
  }

  X <- cbind(x, y)
  d  <- 2

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

}

############################################################################
###                                                                      ###
###                       Bivariate Skew t Density                       ###
###                                                                      ###
############################################################################

bivariate_dSt <- function(
    x,
    y,
    mu    = rep(0, 2),    # location
    Sigma = diag(2),      # dispersion
    beta  = rep(0, 2),    # skewness
    df    = 10,           # degree of freedom
    log   = FALSE
) {

  #----------------------#
  #    Input Checking    #
  #----------------------#

  if (length(mu) != 2) {
    stop('mu must be a vector of length 2')
  }

  if (nrow(Sigma) != 2 | ncol(Sigma) != 2) {
    stop('Sigma must be a 2 x 2 matrix')
  }

  if (length(beta) != 2) {
    stop('beta must be a vector of length 2')
  }

  X <- cbind(x, y)
  d  <- 2

  X_centrd  <- sweep(X, 2, mu, '-')
  Sigma_inv <- solve(Sigma)

  psi <- c( t(beta) %*% Sigma_inv %*% beta )
  chi <- df + mahalanobis(X, center = mu, cov = Sigma_inv, inverted = TRUE)

  s1 <- sqrt(psi * chi)

  res <- (-df/2 - d/2) * log(s1) + (df/2 + d/2) * log(psi)
  res <- res + (df / 2) * log(df) + log(besselK(s1, nu = -df/2 - d/2, expon.scaled = TRUE)) - s1

  if (is.nan(log(det(Sigma)))) {
    Sigma <- diag(d)
  }

  res <- res - d/2 * (log(2) + log(pi)) - 1/2 * log(det(Sigma))
  res <- res - lgamma(df / 2) - (df/2 - 1) * log(2)
  res <- res + X_centrd %*% Sigma_inv %*% beta

  if (!log) {
    res <- c(exp(res))
    res[res <= 10^(-323)] <- 10^(-323)
  }

  return(res)

}

############################################################################
###                                                                      ###
###                       Bivariate Normal Density                       ###
###                                                                      ###
############################################################################

bivariate_dmvnorm <- function(
    x,
    y,
    mean          = rep(0, 2),
    sigma         = diag(2),
    log           = FALSE,
    checkSymmetry = TRUE
) {

  res <- mvtnorm::dmvnorm(cbind(x, y), mean = mean, sigma = sigma, log = log, checkSymmetry = checkSymmetry)

  return(res)

}

###########################################################################
###                                                                     ###
###                         Bivariate t Density                         ###
###                                                                     ###
###########################################################################

bivariate_dmvt <- function(
    x,
    y,
    delta         = rep(0, 2),
    sigma         = diag(2),
    df            = 1,
    log           = FALSE,
    type          = "shifted",
    checkSymmetry = TRUE
) {

  res <- mvtnorm::dmvt(cbind(x, y), delta = delta, sigma = sigma, df = df, log = log, type = type, checkSymmetry = checkSymmetry)

  return(res)

}
