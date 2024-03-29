% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bankruptcy}
\alias{bankruptcy}
\title{Bankruptcy Data Set}
\format{
A data frame with 66 rows and 3 variables:
\describe{
  \item{Y}{Status of the firm: 0 for bankruptcy and 1 for financially sound.}
  \item{RE}{Ratio of retained earnings.}
  \item{EBIT}{Ratio of earnings before interests and taxes.}
}
}
\source{
Altman E.I. (1968). Financial ratios, discriminant analysis and the prediction of corporate bankruptcy. \emph{J Finance} 23(4): 589-609
\url{https://www.jstor.org/stable/2978933}
}
\usage{
bankruptcy
}
\description{
The data set contains the ratio of retained earnings (RE) to total assets, and
the ratio of earnings before interests and taxes (EBIT) to total assets of 66
American firms recorded in the form of ratios. Half of the selected firms had
filed for bankruptcy.
}
\keyword{datasets}
