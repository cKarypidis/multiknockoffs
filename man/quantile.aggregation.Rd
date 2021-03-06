% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile.aggregation.R
\name{quantile.aggregation}
\alias{quantile.aggregation}
\title{P-values quantile aggregation.}
\usage{
\method{quantile}{aggregation}(pvals, gamma = NULL)
}
\arguments{
\item{pvals}{B x p matrix of p-values for all variables. B is the number of p-values
per variable.}

\item{gamma}{value between (0,1) which defines the quantile value used for the aggregation.
If \code{gamma = NULL}, the adaptive search by Meinshausen et al. (2009) is used.
Default: \code{NULL}.}
}
\value{
Aggregated p-values.
}
\description{
This function aggregates p-values by the quantile aggregation proposed by
Meinshausen et al. (2009).
}
\details{
The function is a helper function that it used in the p-value knockoff filter.
}
\references{
Meinshausen, Meier and Buehlmann (2009). \emph{p-Values for High-Dimensional Regression}.
Journal of the American Statistical Association 104(488), 1671-1681.
}
