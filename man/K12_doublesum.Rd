% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/K12_doublesum.R
\name{K12_doublesum}
\alias{K12_doublesum}
\title{Doublesum in K12 estimation
The double sum of indicators in K12 estimator. No edge correction or anything fancy.}
\usage{
K12_doublesum(xy, mark, r)
}
\arguments{
\item{xy}{coordinates}

\item{mark}{1-base integer vector giving tpyes}

\item{r}{ranges}
}
\description{
Doublesum in K12 estimation
The double sum of indicators in K12 estimator. No edge correction or anything fancy.
}
\details{
Use double_integral_ball to compute the global edge correction factor.
}
