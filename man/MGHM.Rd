% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MGHM.R
\name{MGHM}
\alias{MGHM}
\title{Multivariate Generalized Hyperbolic Mixture (MCNM)}
\usage{
MGHM(
  X,
  G,
  max_iter = 20,
  epsilon = 0.01,
  init_method = c("kmedoids", "kmeans", "hierarchical", "manual", "emEM", "RndEM"),
  clusters = NULL,
  impute = FALSE,
  equal_prop = FALSE,
  identity_cov = FALSE,
  deriv_ctrl = list(eps = 1e-08, d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07), r
    = 6, v = 2, show.details = FALSE),
  progress = TRUE,
  n_run = 100,
  n_short = NULL,
  short_eps = 0.1
)
}
\arguments{
\item{X}{An \eqn{n} x \eqn{d} matrix or data frame where \eqn{n} is the number of
observations and \eqn{d} is the number of variables.}

\item{G}{The number of clusters, which must be at least 1. If \code{G = 1}, then
both \code{init_method} and \code{clusters} are ignored.}

\item{max_iter}{(optional) A numeric value giving the maximum number of
iterations each EM algorithm is allowed to use; 20 by default.}

\item{epsilon}{(optional) A number specifying the epsilon value for the
Aitken-based stopping criterion used in the EM algorithm: 0.01 by default.}

\item{init_method}{(optional) A string specifying the method to initialize
the EM algorithm. "kmedoids" clustering is used by default. Alternative
methods include "kmeans", "hierarchical", "manual", "emEM",
and "RndEM". When "manual" is chosen, a vector \code{clusters} of
length \eqn{n} must be specified.}

\item{clusters}{(optional) A numeric vector of length \eqn{n} that specifies the initial
cluster memberships of the user when \code{init_method} is set to "manual".
This argument is NULL by default, so that it is ignored whenever other given
initialization methods are chosen.}

\item{impute}{(optional) A logical value indicating whether missing values should
be imputed for initialization. It is FALSE by default, in which only complete
observations are used for obtaining initial parameters. When it is TRUE, imputation
varies depending on the initialization method selected. For "emEM" and "RndEM",
after observations are randomly assigned cluster memberships, missing values
are replaced by the corresponding cluster means. For other heuristic methods,
mean imputation is applied on the whole data set as a pre-processing step.}

\item{equal_prop}{(optional) A logical value indicating whether mixing
proportions should be equal when initialized with emEM or RndEM; FALSE by
default.}

\item{identity_cov}{(optional) A logical value indicating whether covariance
matrices should be set to identity matrices when initialized with emEM or RndEM;
FALSE by default.}

\item{deriv_ctrl}{(optional) A list containing arguments to control the numerical
procedures for calculating the first and second derivatives. Some values are
suggested by default. Refer to functions \code{grad} and \code{hessian} under
the package \code{numDeriv} for more information.}

\item{progress}{(optional) A logical value indicating whether the
fitting progress should be displayed; TRUE by default.}

\item{n_run}{(optional) Number of random sets to consider for initialization
if \code{init_method = "emEM"} or \code{init_method = "RndEM"}; 100 by default.}

\item{n_short}{(optional) Number of iterations in each run of the short EM
phase if \code{init_method = "emEM"}. It is ignored when another initialization
method is used. When \code{init_method = "emEM"}, emEM reduces to RndEM. It is
NULL by default.}

\item{short_eps}{(optional) The epsilon value for the Aitken-based stopping criterion
used the short EM phase. The value is ignored if \code{n_short} is specified (not NULL).
By default, it is 0.1.}
}
\value{
An object of class \code{MixtureMissing} with:
  \item{model}{The model used to fit the data set}
  \item{pi}{Mixing proportions.}
  \item{mu}{Component mean vectors (location).}
  \item{Sigma}{Component covariance matrices (dispersion).}
  \item{alpha}{Component skewness vectors.}
  \item{lambda}{Component index parameters.}
  \item{omega}{Component concentration parameters.}
  \item{z_tilde}{An \eqn{n} by \eqn{G} matrix where each row indicates the expected
    probabilities that the corresponding observation belongs to each cluster.}
  \item{clusters}{A numeric vector of length \eqn{n} indicating cluster
    memberships determined by the model.}
  \item{data}{The original data set if it is complete; otherwise, this is
    the data set with missing values imputed by appropriate expectations.}
  \item{complete}{A logical vector of length \eqn{n} indicating which observation(s)
    have no missing values.}
  \item{npar}{The breakdown of the number of parameters to estimate.}
  \item{max_iter}{Maximum number of iterations allowed in the EM algorithm.}
  \item{iter_stop}{The actual number of iterations needed when fitting the
    data set.}
  \item{final_loglik}{The final value of log-likelihood.}
  \item{loglik}{All the values of log-likelihood.}
  \item{AIC}{Akaike information criterion.}
  \item{BIC}{Bayesian information criterion.}
  \item{KIC}{Kullback information criterion.}
  \item{KICc}{Corrected Kullback information criterion.}
  \item{AIC3}{Modified AIC.}
  \item{CAIC}{Bozdogan's consistent AIC.}
  \item{AICc}{Small-sample version of AIC.}
  \item{ent}{Entropy}
  \item{ICL}{Integrated Completed Likelihood criterion.}
  \item{AWE}{Approximate weight of evidence.}
  \item{CLC}{Classification likelihood criterion.}
  \item{init_method}{The initialization method used in model fitting.}
  \item{n_run}{Number of random sets considered for initialization if emEM or RndEM is used.}
  \item{n_short}{Number of iterations used in each run of the short EM phase.}
  \item{short_eps}{The epsilon value for the Aitken-based stopping criterion used the short EM phase.}
}
\description{
Carries out model-based clustering using a multivariate generalized hyperbolic
mixture (MGHM). The function will determine itself if the data set is
complete or incomplete and fit the appropriate model accordingly. In the incomplete
case, the data set must be at least bivariate, and missing values are assumed to
be missing at random (MAR).
}
\examples{

data('bankruptcy')

#++++ With no missing values ++++#

X <- bankruptcy[, 2:3]
mod <- MGHM(X, G = 2, init_method = 'kmedoids', max_iter = 10)

summary(mod)
plot(mod)

#++++ With missing values ++++#

set.seed(1234)

X <- hide_values(bankruptcy[, 2:3], prop_cases = 0.1)
mod <- MGHM(X, G = 2, init_method = 'kmedoids', max_iter = 10)

summary(mod)
plot(mod)

}
\references{
Browne, R. P. and McNicholas, P. D. (2015). A mixture of generalized hyperbolic distributions.
  \emph{Canadian Journal of Statistics}, 43(2):176–198. \cr \cr
Wei, Y., Tang, Y., and McNicholas, P. D. (2019). Mixtures of generalized hyperbolic
  distributions and mixtures of skew-\emph{t} distributions for model-based clustering
   with incomplete data. \emph{Computational Statistics & Data Analysis}, 130:18–41.
}