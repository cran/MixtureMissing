#' Automobile Data Set
#'
#' This data set consists of three types of entities: (a) the specification of
#' an auto in terms of various characteristics, (b) its assigned insurance risk
#' rating, (c) its normalized losses in use as compared to other cars. The second
#' rating corresponds to the degree to which the auto is more risky than its
#' price indicates. Cars are initially assigned a risk factor symbol associated
#' with its price. Then, if it is more risky (or less), this symbol is adjusted
#' by moving it up (or down) the scale. Actuarians call this process "symboling".
#' A value of +3 indicates that the auto is risky, -3 that it is probably pretty safe.
#'
#' @format A data frame with 205 rows and 26 variables. The first 15 variables are
#'   continuous, while the last 11 variables are categorical. There are 45 rows
#'   with missing values.
#' \describe{
#'   \item{normalized_losses}{continuous from 65 to 256.}
#'   \item{wheel_base}{continuous from 86.6 120.9.}
#'   \item{length}{continuous from 141.1 to 208.1.}
#'   \item{width}{continuous from 60.3 to 72.3.}
#'   \item{height}{continuous from 47.8 to 59.8.}
#'   \item{curb_weight}{continuous from 1488 to 4066.}
#'   \item{engine_size}{continuous from 61 to 326.}
#'   \item{bore}{continuous from 2.54 to 3.94.}
#'   \item{stroke}{continuous from 2.07 to 4.17.}
#'   \item{compression_ratio}{continuous from 7 to 23.}
#'   \item{horsepower}{continuous from 48 to 288.}
#'   \item{peak_rpm}{continuous from 4150 to 6600.}
#'   \item{city_mpg}{continuous from 13 to 49.}
#'   \item{highway_mpg}{continuous from 16 to 54.}
#'   \item{price}{continuous from 5118 to 45400.}
#'   \item{symboling}{-3, -2, -1, 0, 1, 2, 3.}
#'   \item{make}{alfa-romero, audi, bmw, chevrolet, dodge, honda, isuzu, jaguar,
#'     mazda, mercedes-benz, mercury, mitsubishi, nissan, peugot, plymouth, porsche,
#'     renault, saab, subaru, toyota, volkswagen, volvo}
#'   \item{fuel_type}{diesel, gas.}
#'   \item{aspiration}{std, turbo.}
#'   \item{num_doors}{four, two.}
#'   \item{body_style}{hardtop, wagon, sedan, hatchback, convertible.}
#'   \item{drive_wheels}{4wd, fwd, rwd.}
#'   \item{engine_location}{front, rear.}
#'   \item{engine_type}{dohc, dohcv, l, ohc, ohcf, ohcv, rotor.}
#'   \item{num_cylinders}{eight, five, four, six, three, twelve, two.}
#'   \item{fuel_system}{1bbl, 2bbl, 4bbl, idi, mfi, mpfi, spdi, spfi.}
#' }
#' @source Kibler, D., Aha, D.W., & Albert,M. (1989). Instance-based prediction of real-valued attributes. Computational Intelligence, Vol 5, 51--57.
#' \url{https://archive.ics.uci.edu/ml/datasets/automobile}
"auto"

#' Bankruptcy Data Set
#'
#' The data set contain the ratio of retained earnings (RE) to total assets, and
#' the ratio of earnings before interests and taxes (EBIT) to total assets of 66
#' American firms recorded in the form of ratios. Half of the selected firms had
#' filed for bankruptcy.
#'
#' @format A data frame with 66 rows and 3 variables:
#' \describe{
#'   \item{Y}{Status of the firm: 0 for bankruptcy and 1 for financially sound.}
#'   \item{RE}{Ratio of retained earnings.}
#'   \item{EBIT}{ratio of earnings before interests and taxes.}
#' }
#' @source Altman E.I. (1968). Financial ratios, discriminant analysis and the prediction of corporate bankruptcy. \emph{J Finance} 23(4): 589-609
#' \url{https://www.jstor.org/stable/2978933}
"bankruptcy"

#' A Mixture of Two Far Student's \eqn{t} Distributions - 500 Observations
#'
#' A simulated mixture of two far Student's \eqn{t} distributions. Refer to Punzo
#' and McNicholas (2016) for more information about the underlying distribution
#' that generates this data set.
#'
#' @format A matrix with 500 rows and 3 variables. The first two variables make
#'   the bivariate data, while the last variable refers to cluster memberships.
#'   The first 150 rows belong to cluster 1, and the last 350 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"tm_far_500"

#' A Mixture of Two Close Student's \eqn{t} Distributions - 500 Observations
#'
#' A simulated mixture of two close Student's \eqn{t} distributions. Refer to Punzo
#' and McNicholas (2016) for more information about the underlying distribution
#' that generates this data set.
#'
#' @format A matrix with 500 rows and 3 variables. The first two variables make
#'   the bivariate data, while the last variable refers to cluster memberships.
#'   The first 150 rows belong to cluster 1, and the last 350 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"tm_close_500"

#' A Mixture of Two Far Contaminated Normal Distributions - 500 Observations
#'
#' A simulated mixture of two far contaminated normal distributions. Refer to Punzo
#' and McNicholas (2016) for more information about the underlying distribution
#' that generates this data set.
#'
#' @format A matrix with 500 rows and 3 variables. The first two variables make
#'   the bivariate data, while the last variable refers to cluster memberships.
#'   The first 150 rows belong to cluster 1, and the last 350 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"cnm_far_500"

#' A Mixture of Two Close Contaminated Normal Distributions - 500 Observations
#'
#' A simulated mixture of two close contaminated normal distributions. Refer to Punzo
#' and McNicholas (2016) for more information about the underlying distribution
#' that generates this data set.
#'
#' @format A matrix with 500 rows and 3 variables. The first two variables make
#'   the bivariate data, while the last variable refers to cluster memberships.
#'   The first 150 rows belong to cluster 1, and the last 350 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"cnm_close_500"

#' A Mixture of Two Far Normal Distributions with 1% of Points Randomly Substituted
#' by High Atypical Points - 500 Observations
#'
#' A simulated mixture of two far normal distributions with 1% of points randomly
#' substituted by high atypical points. Refer to Punzo and McNicholas (2016) for
#' more information about the underlying distribution that generates this data set.
#'
#' @format A matrix with 500 rows and 3 variables. The first two variables make
#'   the bivariate data, while the last variable refers to cluster memberships.
#'   The first 150 rows belong to cluster 1, and the last 350 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"nm_1_noise_far_500"

#' A Mixture of Two Close Normal Distributions with 1% of Points Randomly Substituted
#' by High Atypical Points - 500 Observations
#'
#' A simulated mixture of two close normal distributions with 1% of points randomly
#' substituted by high atypical points. Refer to Punzo and McNicholas (2016) for
#' more information about the underlying distribution that generates this data set.
#'
#' @format A matrix with 500 rows and 4 variables. The first two variables make
#'   the bivariate data; the last variable refers to cluster memberships; and the
#'   last variable refers to outlier indication (1 means outlier, 0 otherwise).
#'   The first 150 rows belong to cluster 1, and the last 350 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#'   \item{outlier}{outlier indication}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"nm_1_noise_close_500"

#' A Mixture of Two Far Normal Distributions with 5% of Points Randomly Substituted
#' by Noise - 500 Observations
#'
#' A simulated mixture of two far normal distributions with 1% of points randomly
#' substituted by noise. Refer to Punzo and McNicholas (2016) for more information
#' about the underlying distribution that generates this data set.
#'
#' @format A matrix with 500 rows and 4 variables. The first two variables make
#'   the bivariate data; the last variable refers to cluster memberships; and the
#'   last variable refers to outlier indication (1 means outlier, 0 otherwise).
#'   The first 150 rows belong to cluster 1, and the last 350 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#'   \item{outlier}{outlier indication}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"nm_5_noise_far_500"

#' A Mixture of Two Close Normal Distributions with 5% of Points Randomly Substituted
#' by Noise - 500 Observations
#'
#' A simulated mixture of two close normal distributions with 1% of points randomly
#' substituted by noise. Refer to Punzo and McNicholas (2016) for more information
#' about the underlying distribution that generates this data set.
#'
#' @format A matrix with 500 rows and 4 variables. The first two variables make
#'   the bivariate data; the last variable refers to cluster memberships; and the
#'   last variable refers to outlier indication (1 means outlier, 0 otherwise).
#'   The first 150 rows belong to cluster 1, and the last 350 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#'   \item{outlier}{outlier indication}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"nm_5_noise_close_500"

#' A Mixture of Two Far Student's \eqn{t} Distributions - 100 Observations
#'
#' A simulated mixture of two far Student's \eqn{t} distributions. Refer to Punzo
#' and McNicholas (2016) for more information about the underlying distribution
#' that generates this data set.
#'
#' @format A matrix with 500 rows and 4 variables. The first two variables make
#'   the bivariate data; the last variable refers to cluster memberships; and the
#'   last variable refers to outlier indication (1 means outlier, 0 otherwise).
#'   The first 30 rows belong to cluster 1, and the last 70 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#'   \item{outlier}{outlier indication}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"tm_far_100"

#' A Mixture of Two Close Student's \eqn{t} Distributions - 100 Observations
#'
#' A simulated mixture of two close Student's \eqn{t} distributions. Refer to Punzo
#' and McNicholas (2016) for more information about the underlying distribution
#' that generates this data set.
#'
#' @format A matrix with 100 rows and 3 variables. The first two variables make
#'   the bivariate data, while the last variable refers to cluster memberships.
#'   The first 30 rows belong to cluster 1, and the last 70 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"tm_close_100"

#' A Mixture of Two Far Contaminated Normal Distributions - 100 Observations
#'
#' A simulated mixture of two far contaminated normal distributions. Refer to Punzo
#' and McNicholas (2016) for more information about the underlying distribution
#' that generates this data set.
#'
#' @format A matrix with 100 rows and 3 variables. The first two variables make
#'   the bivariate data, while the last variable refers to cluster memberships.
#'   The first 30 rows belong to cluster 1, and the last 70 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"cnm_far_100"

#' A Mixture of Two Close Contaminated Normal Distributions - 100 Observations
#'
#' A simulated mixture of two close contaminated normal distributions. Refer to Punzo
#' and McNicholas (2016) for more information about the underlying distribution
#' that generates this data set.
#'
#' @format A matrix with 100 rows and 3 variables. The first two variables make
#'   the bivariate data, while the last variable refers to cluster memberships.
#'   The first 30 rows belong to cluster 1, and the last 70 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"cnm_close_100"

#' A Mixture of Two Far Normal Distributions with 1% of Points Randomly Substituted
#' by High Atypical Points - 100 Observations
#'
#' A simulated mixture of two far normal distributions with 1% of points randomly
#' substituted by high atypical points. Refer to Punzo and McNicholas (2016) for
#' more information about the underlying distribution that generates this data set.
#'
#' @format A matrix with 100 rows and 3 variables. The first two variables make
#'   the bivariate data, while the last variable refers to cluster memberships.
#'   The first 30 rows belong to cluster 1, and the last 70 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"nm_1_noise_far_100"

#' A Mixture of Two Close Normal Distributions with 1% of Points Randomly Substituted
#' by High Atypical Points - 100 Observations
#'
#' A simulated mixture of two close normal distributions with 1% of points randomly
#' substituted by high atypical points. Refer to Punzo and McNicholas (2016) for
#' more information about the underlying distribution that generates this data set.
#'
#' @format A matrix with 100 rows and 4 variables. The first two variables make
#'   the bivariate data; the last variable refers to cluster memberships; and the
#'   last variable refers to outlier indication (1 means outlier, 0 otherwise).
#'   The first 30 rows belong to cluster 1, and the last 70 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#'   \item{outlier}{outlier indication}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"nm_1_noise_close_100"

#' A Mixture of Two Far Normal Distributions with 5% of Points Randomly Substituted
#' by Noise - 100 Observations
#'
#' A simulated mixture of two far normal distributions with 1% of points randomly
#' substituted by noise. Refer to Punzo and McNicholas (2016) for more information
#' about the underlying distribution that generates this data set.
#'
#' @format A matrix with 100 rows and 4 variables. The first two variables make
#'   the bivariate data; the last variable refers to cluster memberships; and the
#'   last variable refers to outlier indication (1 means outlier, 0 otherwise).
#'   The first 30 rows belong to cluster 1, and the last 70 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#'   \item{outlier}{outlier indication}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"nm_5_noise_far_100"

#' A Mixture of Two Close Normal Distributions with 5% of Points Randomly Substituted
#' by Noise - 100 Observations
#'
#' A simulated mixture of two close normal distributions with 1% of points randomly
#' substituted by noise. Refer to Punzo and McNicholas (2016) for more information
#' about the underlying distribution that generates this data set.
#'
#' @format A matrix with 100 rows and 4 variables. The first two variables make
#'   the bivariate data; the last variable refers to cluster memberships; and the
#'   last variable refers to outlier indication (1 means outlier, 0 otherwise).
#'   The first 30 rows belong to cluster 1, and the last 70 rows belong to cluster
#'   2
#' \describe{
#'   \item{d1}{variable 1.}
#'   \item{d2}{variable 2.}
#'   \item{cluster}{cluster memberships}
#'   \item{outlier}{outlier indication}
#' }
#' @source Punzo, A. and McNicholas, P.D., 2016. Parsimonious mixtures of multivariate
#'   contaminated normal distributions. \emph{Biometrical Journal, 58}(6), pp.1506-1537.
"nm_5_noise_close_100"

